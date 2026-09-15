module ML_Surrogate_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use iso_c_binding
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! The network shape (n_in -> n_h0 -> n_h1 -> n_h2 -> n_out) is read from the
  ! weight file itself (dataset extents), so a retrain with a different width
  ! or a different number of inputs/outputs needs no recompile.  Two output
  ! stages are supported, selected by the root attribute 'output_stage_id'
  ! (absent -> 0, which keeps every pre-existing weight file valid):
  !
  !   0  ML_STAGE_PHI_SINH  : rate[1/y] = asinh_scale * sinh(o(1)*y_scale + y_mean)
  !                           (B-huber / B-thuber / NODE-phi files; 7 inputs, 1 output)
  !   1  ML_STAGE_NODE_MIN  : 'node_min_v1' physics-bounded multi-state stage
  !                           (NODE-min files; 10 inputs, 5 outputs).  The host
  !                           integrates four mineral-volume states and derives
  !                           the porosity from them; see MLPForwardNodeMin.
  !
  ! Weight-file datasets (unchanged contract): /scalers/X_mean, X_scale
  ! [/scalers/y_mean, y_scale for stage 0], /transforms/log10_mask, clip_lo,
  ! clip_hi, /layers/layer_00..03/{W,b} with W stored PyTorch-style (n_out,n_in)
  ! in C order, read into Fortran (n_in,n_out) so the forward is matmul(h, W).
  PetscInt, parameter, public :: ML_STAGE_PHI_SINH = 0
  PetscInt, parameter, public :: ML_STAGE_NODE_MIN = 1
  PetscInt, parameter, public :: ML_RAMP_SMOOTHSTEP2 = 0
  PetscInt, parameter, public :: ML_RAMP_LINEAR = 1

  ! Output rate transform of stage 0: rate[1/y] = ASINH_SCALE * sinh(z).
  PetscReal, parameter, public :: ML_ASINH_SCALE = 1.d-5

  type, public :: ml_mlp_type
    PetscInt :: n_in
    PetscInt :: n_h0, n_h1, n_h2
    PetscInt :: n_out
    PetscInt :: stage                  ! ML_STAGE_*
    PetscReal, pointer :: X_mean(:)
    PetscReal, pointer :: X_scale(:)
    PetscReal, pointer :: clip_lo(:)
    PetscReal, pointer :: clip_hi(:)
    PetscReal, pointer :: log10_mask(:)  ! >0.5 -> apply log10 after clip
    PetscReal :: y_mean
    PetscReal :: y_scale
    PetscReal :: asinh_scale
    ! node_min_v1 output-stage constants (root attributes of the weight file)
    PetscReal :: rate_scale
    PetscReal :: phi_res
    PetscReal :: ramp_width
    PetscInt  :: ramp_form             ! ML_RAMP_*
    PetscReal :: k_exch
    PetscReal :: cs_min
    PetscReal :: cs_max
    PetscReal :: phi0
    PetscReal :: rate_max
    PetscReal :: csh_vm_a              ! Vm_CSH(cs) = csh_vm_a + csh_vm_b*cs  [cm^3/mol]
    PetscReal :: csh_vm_b
    PetscReal :: m4ac_vm               ! [cm^3/mol]
    PetscReal, pointer :: W0(:,:), b0(:)
    PetscReal, pointer :: W1(:,:), b1(:)
    PetscReal, pointer :: W2(:,:), b2(:)
    PetscReal, pointer :: W3(:,:), b3(:)
  end type ml_mlp_type

  ! Everything the node_min_v1 stage returns for one cell and one call.
  type, public :: ml_node_min_out_type
    PetscReal :: d_vf_csh      ! d vf_CSH_total / dt        [1/y]  (>= 0)
    PetscReal :: d_s_csh       ! d S_CSH / dt, S = vf*meanCS [1/y]  (signed, bounded)
    PetscReal :: d_vf_m4ac     ! d vf_M4Ac / dt              [1/y]  (>= 0)
    PetscReal :: d_vf_other    ! d vf_other / dt             [1/y]  (>= 0)
    PetscReal :: dphi          ! = -(d_vf_csh + d_vf_m4ac + d_vf_other)  (<= 0)
    PetscReal :: cs_new        ! Ca/Si of the CSH being formed, in [cs_min, cs_max]
    PetscReal :: sink_Si       ! aqueous sinks [mol / m^3 bulk / y], positive = consumed
    PetscReal :: sink_Ca
    PetscReal :: sink_Mg
    PetscReal :: sink_Al
    PetscReal :: sink_HCO3
    PetscReal :: sink_TotalH   ! negative: H+ is released on precipitation
  end type ml_node_min_out_type

  public :: MLPCreate, &
            MLPLoad, &
            MLPLoadFromFileId, &
            MLPTrunk, &
            MLPForward, &
            MLPForwardNodeMin, &
            MLPDestroy, &
            MLPRelu

contains

! ************************************************************************** !

subroutine MLPCreate(mlp)
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none
  type(ml_mlp_type) :: mlp

  mlp%n_in = 0; mlp%n_h0 = 0; mlp%n_h1 = 0; mlp%n_h2 = 0; mlp%n_out = 0
  mlp%stage = ML_STAGE_PHI_SINH
  mlp%y_mean = 0.d0
  mlp%y_scale = 1.d0
  mlp%asinh_scale = ML_ASINH_SCALE
  mlp%rate_scale = 0.d0; mlp%phi_res = 0.d0; mlp%ramp_width = 1.d0
  mlp%ramp_form = ML_RAMP_SMOOTHSTEP2; mlp%k_exch = 0.d0
  mlp%cs_min = 0.d0; mlp%cs_max = 0.d0; mlp%phi0 = 0.d0; mlp%rate_max = 0.d0
  mlp%csh_vm_a = 0.d0; mlp%csh_vm_b = 0.d0; mlp%m4ac_vm = 0.d0
  nullify(mlp%X_mean); nullify(mlp%X_scale)
  nullify(mlp%clip_lo); nullify(mlp%clip_hi); nullify(mlp%log10_mask)
  nullify(mlp%W0); nullify(mlp%b0)
  nullify(mlp%W1); nullify(mlp%b1)
  nullify(mlp%W2); nullify(mlp%b2)
  nullify(mlp%W3); nullify(mlp%b3)

end subroutine MLPCreate

! ************************************************************************** !

subroutine MLPLoad(mlp, filename, option)
  !
  ! Opens the weight file through PFLOTRAN's HDF5 helper and loads it.  Any
  ! contract violation is fatal (PrintErrMsg) with the offending item named.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24 (generalized to file-sized layers and two output stages,
  !       Claude Fable 5.1 for Heeho Park, 2026-09-14)
  !
  use HDF5_Aux_module
  use Option_module

  implicit none

  type(ml_mlp_type) :: mlp
  character(len=*) :: filename
  type(option_type) :: option

  integer(HID_T) :: file_id
  integer :: hdf5_err, ierr
  character(len=MAXSTRINGLENGTH) :: errmsg

  call HDF5FileOpenReadOnly(filename, file_id, PETSC_FALSE, &
                            'ML surrogate weights', option)

  call MLPLoadFromFileId(mlp, file_id, ierr, errmsg)
  if (ierr /= 0) then
    option%io_buffer = 'ML surrogate weights file "' // trim(filename) // &
                       '": ' // trim(errmsg)
    call PrintErrMsg(option)
  endif

  call h5fclose_f(file_id, hdf5_err)
  if (hdf5_err /= 0) then
    option%io_buffer = 'Error closing ML surrogate weights file "' // &
                       trim(filename) // '".'
    call PrintErrMsg(option)
  endif

  write(option%io_buffer,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
    'ML surrogate loaded: n_in=', mlp%n_in, ' hidden=', mlp%n_h0, '/', &
    mlp%n_h1, '/', mlp%n_h2, ' n_out=', mlp%n_out, ' output_stage=', mlp%stage
  call PrintMsg(option)

end subroutine MLPLoad

! ************************************************************************** !

subroutine MLPLoadFromFileId(mlp, file_id, ierr, errmsg)
  !
  ! Option-free loader shared by MLPLoad and the standalone validator: reads
  ! scalers, transforms, the four layers (allocating them from the dataset
  ! extents) and the output-stage attributes, and checks the shapes agree.
  ! ierr /= 0 with errmsg set on any violation.
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  type(ml_mlp_type) :: mlp
  integer(HID_T) :: file_id
  integer :: ierr
  character(len=*) :: errmsg

  PetscInt :: itmp
  PetscBool :: found

  ierr = 0; errmsg = ''

  call ReadVectorAlloc(file_id, '/scalers/X_mean',  mlp%X_mean,  ierr, errmsg); if (ierr /= 0) return
  call ReadVectorAlloc(file_id, '/scalers/X_scale', mlp%X_scale, ierr, errmsg); if (ierr /= 0) return
  call ReadVectorAlloc(file_id, '/transforms/clip_lo',    mlp%clip_lo,    ierr, errmsg); if (ierr /= 0) return
  call ReadVectorAlloc(file_id, '/transforms/clip_hi',    mlp%clip_hi,    ierr, errmsg); if (ierr /= 0) return
  call ReadVectorAlloc(file_id, '/transforms/log10_mask', mlp%log10_mask, ierr, errmsg); if (ierr /= 0) return
  mlp%n_in = size(mlp%X_mean)
  if (size(mlp%X_scale) /= mlp%n_in .or. size(mlp%clip_lo) /= mlp%n_in .or. &
      size(mlp%clip_hi) /= mlp%n_in .or. size(mlp%log10_mask) /= mlp%n_in) then
    ierr = 1; errmsg = '/scalers and /transforms vectors do not all have length n_in.'; return
  endif

  call ReadLayerAlloc(file_id, '/layers/layer_00', mlp%W0, mlp%b0, ierr, errmsg); if (ierr /= 0) return
  call ReadLayerAlloc(file_id, '/layers/layer_01', mlp%W1, mlp%b1, ierr, errmsg); if (ierr /= 0) return
  call ReadLayerAlloc(file_id, '/layers/layer_02', mlp%W2, mlp%b2, ierr, errmsg); if (ierr /= 0) return
  call ReadLayerAlloc(file_id, '/layers/layer_03', mlp%W3, mlp%b3, ierr, errmsg); if (ierr /= 0) return
  mlp%n_h0 = size(mlp%b0); mlp%n_h1 = size(mlp%b1); mlp%n_h2 = size(mlp%b2); mlp%n_out = size(mlp%b3)
  if (size(mlp%W0,1) /= mlp%n_in  .or. size(mlp%W0,2) /= mlp%n_h0 .or. &
      size(mlp%W1,1) /= mlp%n_h0  .or. size(mlp%W1,2) /= mlp%n_h1 .or. &
      size(mlp%W2,1) /= mlp%n_h1  .or. size(mlp%W2,2) /= mlp%n_h2 .or. &
      size(mlp%W3,1) /= mlp%n_h2  .or. size(mlp%W3,2) /= mlp%n_out) then
    ierr = 1
    errmsg = 'layer shapes are inconsistent (expect W_k as (n_out,n_in) C-order chained n_in->h0->h1->h2->n_out).'
    return
  endif

  ! ---- output stage ------------------------------------------------------
  call ReadAttrInt(file_id, 'output_stage_id', itmp, found)
  mlp%stage = ML_STAGE_PHI_SINH
  if (found) mlp%stage = itmp

  select case(mlp%stage)
    case(ML_STAGE_PHI_SINH)
      if (mlp%n_out /= 1) then
        ierr = 1; errmsg = 'output_stage 0 (sinh rate) requires n_out = 1.'; return
      endif
      call ReadScalarDataset(file_id, '/scalers/y_mean',  mlp%y_mean,  ierr, errmsg); if (ierr /= 0) return
      call ReadScalarDataset(file_id, '/scalers/y_scale', mlp%y_scale, ierr, errmsg); if (ierr /= 0) return
      mlp%asinh_scale = ML_ASINH_SCALE
    case(ML_STAGE_NODE_MIN)
      if (mlp%n_out /= 5 .or. mlp%n_in /= 10) then
        ierr = 1; errmsg = 'output_stage 1 (node_min_v1) requires n_in = 10 and n_out = 5.'; return
      endif
      call ReadAttrReal(file_id, 'rate_scale', mlp%rate_scale, found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'phi_res',    mlp%phi_res,    found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'ramp_width', mlp%ramp_width, found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'k_exch',     mlp%k_exch,     found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'cs_min',     mlp%cs_min,     found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'cs_max',     mlp%cs_max,     found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'phi0',       mlp%phi0,       found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'rate_max',   mlp%rate_max,   found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'csh_vm_a',   mlp%csh_vm_a,   found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'csh_vm_b',   mlp%csh_vm_b,   found); if (.not.found) ierr = 1
      call ReadAttrReal(file_id, 'm4ac_vm',    mlp%m4ac_vm,    found); if (.not.found) ierr = 1
      call ReadAttrInt (file_id, 'ramp_form_id', mlp%ramp_form, found); if (.not.found) ierr = 1
      if (ierr /= 0) then
        errmsg = 'node_min_v1 requires the root attributes rate_scale, phi_res, ramp_width, k_exch, cs_min, cs_max, ' // &
                 'phi0, rate_max, csh_vm_a, csh_vm_b, m4ac_vm, ramp_form_id.'
        return
      endif
      if (mlp%ramp_form /= ML_RAMP_SMOOTHSTEP2 .and. mlp%ramp_form /= ML_RAMP_LINEAR) then
        ierr = 1; errmsg = 'ramp_form_id must be 0 (smoothstep2) or 1 (linear).'; return
      endif
    case default
      ierr = 1; errmsg = 'unknown output_stage_id (supported: 0 = sinh rate, 1 = node_min_v1).'; return
  end select

end subroutine MLPLoadFromFileId

! ************************************************************************** !

subroutine ReadLayerAlloc(file_id, group_path, W, b, ierr, errmsg)
  !
  ! Reads layer_NN/W and layer_NN/b, allocating from the on-disk extents.
  ! W is stored on disk in PyTorch (n_out, n_in) C-order; it is held in
  ! Fortran as (n_in, n_out) so the bytes line up without a transpose and the
  ! forward pass uses MATMUL(h, W).
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: group_path
  PetscReal, pointer :: W(:,:)
  PetscReal, pointer :: b(:)
  integer :: ierr
  character(len=*) :: errmsg

  call ReadMatrixAlloc(file_id, trim(group_path) // '/W', W, ierr, errmsg)
  if (ierr /= 0) return
  call ReadVectorAlloc(file_id, trim(group_path) // '/b', b, ierr, errmsg)

end subroutine ReadLayerAlloc

! ************************************************************************** !

subroutine ReadVectorAlloc(file_id, dataset_path, vec, ierr, errmsg)
  !
  ! Reads a 1-D float64 dataset into a freshly allocated vector.
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscReal, pointer :: vec(:)
  integer :: ierr
  character(len=*) :: errmsg

  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T) :: dims(1), max_dims(1)
  integer :: hdf5_err, rank

  ierr = 0
  call h5dopen_f(file_id, dataset_path, dset_id, hdf5_err)
  if (hdf5_err < 0) then
    ierr = 1; errmsg = 'dataset "' // trim(dataset_path) // '" not found.'; return
  endif
  call h5dget_space_f(dset_id, space_id, hdf5_err)
  call h5sget_simple_extent_ndims_f(space_id, rank, hdf5_err)
  if (rank /= 1) then
    ierr = 1; errmsg = 'dataset "' // trim(dataset_path) // '" must be 1-D.'
    call h5sclose_f(space_id, hdf5_err); call h5dclose_f(dset_id, hdf5_err); return
  endif
  call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdf5_err)
  allocate(vec(int(dims(1))))
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vec, dims, hdf5_err)
  if (hdf5_err /= 0) then
    ierr = 1; errmsg = 'error reading dataset "' // trim(dataset_path) // '".'
  endif
  call h5sclose_f(space_id, hdf5_err)
  call h5dclose_f(dset_id, hdf5_err)

end subroutine ReadVectorAlloc

! ************************************************************************** !

subroutine ReadMatrixAlloc(file_id, dataset_path, mat, ierr, errmsg)
  !
  ! Reads a 2-D float64 dataset whose on-disk (C-order) shape is (n_out, n_in)
  ! into a freshly allocated Fortran (n_in, n_out) array (no transpose copy).
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscReal, pointer :: mat(:,:)
  integer :: ierr
  character(len=*) :: errmsg

  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T) :: dims(2), max_dims(2)
  integer :: hdf5_err, rank

  ierr = 0
  call h5dopen_f(file_id, dataset_path, dset_id, hdf5_err)
  if (hdf5_err < 0) then
    ierr = 1; errmsg = 'dataset "' // trim(dataset_path) // '" not found.'; return
  endif
  call h5dget_space_f(dset_id, space_id, hdf5_err)
  call h5sget_simple_extent_ndims_f(space_id, rank, hdf5_err)
  if (rank /= 2) then
    ierr = 1; errmsg = 'dataset "' // trim(dataset_path) // '" must be 2-D.'
    call h5sclose_f(space_id, hdf5_err); call h5dclose_f(dset_id, hdf5_err); return
  endif
  ! The Fortran HDF5 API reports dims in Fortran (fastest-first) order, i.e.
  ! dims(1) = n_in, dims(2) = n_out for a C-order (n_out, n_in) dataset.
  call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdf5_err)
  allocate(mat(int(dims(1)), int(dims(2))))
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mat, dims, hdf5_err)
  if (hdf5_err /= 0) then
    ierr = 1; errmsg = 'error reading dataset "' // trim(dataset_path) // '".'
  endif
  call h5sclose_f(space_id, hdf5_err)
  call h5dclose_f(dset_id, hdf5_err)

end subroutine ReadMatrixAlloc

! ************************************************************************** !

subroutine ReadScalarDataset(file_id, dataset_path, scalar, ierr, errmsg)
  !
  ! Reads a length-1 dataset into a scalar.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscReal :: scalar
  integer :: ierr
  character(len=*) :: errmsg

  PetscReal, pointer :: tmp(:)

  nullify(tmp)
  call ReadVectorAlloc(file_id, dataset_path, tmp, ierr, errmsg)
  if (ierr /= 0) return
  if (size(tmp) /= 1) then
    ierr = 1; errmsg = 'dataset "' // trim(dataset_path) // '" must have length 1.'
  else
    scalar = tmp(1)
  endif
  deallocate(tmp)

end subroutine ReadScalarDataset

! ************************************************************************** !

subroutine ReadAttrReal(file_id, name, val, found)
  !
  ! Reads a scalar float64 root attribute; found = .false. if absent.
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: name
  PetscReal, target :: val
  PetscBool :: found

  integer(HID_T) :: attr_id
  integer :: hdf5_err
  logical :: exists
  type(c_ptr) :: ptr

  found = PETSC_FALSE
  call h5aexists_f(file_id, name, exists, hdf5_err)
  if (hdf5_err /= 0 .or. .not.exists) return
  call h5aopen_f(file_id, name, attr_id, hdf5_err)
  if (hdf5_err /= 0) return
  ptr = c_loc(val)
  call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, ptr, hdf5_err)
  found = (hdf5_err == 0)
  call h5aclose_f(attr_id, hdf5_err)

end subroutine ReadAttrReal

! ************************************************************************** !

subroutine ReadAttrInt(file_id, name, val, found)
  !
  ! Reads a scalar integer root attribute (any on-disk integer width; HDF5
  ! converts to the native 4-byte integer); found = .false. if absent.
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: name
  PetscInt :: val
  PetscBool :: found

  integer(HID_T) :: attr_id
  integer :: hdf5_err
  integer(kind=4), target :: ival
  logical :: exists
  type(c_ptr) :: ptr

  found = PETSC_FALSE
  call h5aexists_f(file_id, name, exists, hdf5_err)
  if (hdf5_err /= 0 .or. .not.exists) return
  call h5aopen_f(file_id, name, attr_id, hdf5_err)
  if (hdf5_err /= 0) return
  ptr = c_loc(ival)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, ptr, hdf5_err)
  found = (hdf5_err == 0)
  if (found) val = int(ival, kind(val))
  call h5aclose_f(attr_id, hdf5_err)

end subroutine ReadAttrInt

! ************************************************************************** !

subroutine MLPTrunk(mlp, x_in, o)
  !
  ! Feature pre-processing and the network up to (and including) the linear
  ! output layer: per-feature clip, log10 where masked, standardize, three
  ! ReLU hidden layers, linear output.  o(1:n_out) is the RAW output vector;
  ! the output stage is applied by the caller.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24 (split out of MLPForward, 2026-09-14)
  !
  implicit none

  type(ml_mlp_type) :: mlp
  PetscReal :: x_in(:)
  PetscReal :: o(:)

  PetscReal :: x(mlp%n_in)
  PetscReal :: h_in(mlp%n_in)
  PetscReal :: h0(mlp%n_h0)
  PetscReal :: h1(mlp%n_h1)
  PetscReal :: h2(mlp%n_h2)
  PetscInt :: i

  do i = 1, mlp%n_in
    x(i) = max(mlp%clip_lo(i), min(mlp%clip_hi(i), x_in(i)))
    if (mlp%log10_mask(i) > 0.5d0) then
      x(i) = log10(max(x(i), 1.d-300))
    endif
  enddo

  h_in = (x - mlp%X_mean) / mlp%X_scale

  ! W stored as (n_in, n_out), so "Wh + b" is matmul(h, W) + b in Fortran.
  h0 = MLPRelu(MATMUL(h_in, mlp%W0) + mlp%b0)
  h1 = MLPRelu(MATMUL(h0,   mlp%W1) + mlp%b1)
  h2 = MLPRelu(MATMUL(h1,   mlp%W2) + mlp%b2)
  o(1:mlp%n_out) = MATMUL(h2, mlp%W3) + mlp%b3

end subroutine MLPTrunk

! ************************************************************************** !

subroutine MLPForward(mlp, x_in, rate)
  !
  ! Output stage 0: rate[1/y] = asinh_scale * sinh(o(1)*y_scale + y_mean).
  ! The signed rate is returned as-is; the caller applies the clogging clamp,
  ! envelope clamp and Euler update.  Matches chuber_reference.py exactly.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  type(ml_mlp_type) :: mlp
  PetscReal :: x_in(:)
  PetscReal :: rate

  PetscReal :: o(mlp%n_out)
  PetscReal :: z

  call MLPTrunk(mlp, x_in, o)
  z = o(1) * mlp%y_scale + mlp%y_mean
  rate = mlp%asinh_scale * sinh(z)

end subroutine MLPForward

! ************************************************************************** !

subroutine MLPForwardNodeMin(mlp, x_in, out)
  !
  ! Output stage 1, 'node_min_v1' (executable specification:
  ! reference_node_min() in node_development/export_node.py).  Inputs
  ! x_in(1:10) = porosity, aperture, TotalH, Ca, Si, Al, Mg, vf_CSH_total,
  ! CSH_meanCS (S/vf, 0 if vf <= 1e-6), vf_M4Ac.  Returns the four state rates,
  ! dphi/dt, the Ca/Si of newly formed CSH and the stoichiometric aqueous sinks.
  !
  !   s      = clip((phi - phi_res)/ramp_width, 0, 1)
  !   ramp   = s (linear) | (3 s^2 - 2 s^3)^2 (smoothstep2)
  !   r_csh  = rate_scale*softplus(o1)*ramp ;  cs_new = cs_min + (cs_max-cs_min)*sigmoid(o2)
  !   r_m4ac = rate_scale*softplus(o3)*ramp ;  r_oth  = rate_scale*softplus(o4)*ramp
  !   tot > rate_max -> scale the three rates by rate_max/tot
  !   u = tanh(o5) ; r_x = k_exch*(max(u,0)*(cs_max*vf - S) - max(-u,0)*(S - cs_min*vf))
  !   d vf_csh = r_csh ; d S = cs_new*r_csh + r_x ; d vf_m4ac = r_m4ac ; d vf_other = r_oth
  !   dphi = -(r_csh + r_m4ac + r_oth)
  !   n_csh = r_csh/Vm(cs_new) ; n_m4 = r_m4ac/Vm_m4ac  (Vm in m^3/mol)
  !   sinks: Si = n_csh ; Ca = cs_new*n_csh + r_x/Vm(meanCS) ; Mg = 4 n_m4 ; Al = 2 n_m4 ;
  !          HCO3 = n_m4 ; TotalH = -(2 cs_new n_csh + 13 n_m4)
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  type(ml_mlp_type) :: mlp
  PetscReal :: x_in(:)
  type(ml_node_min_out_type) :: out

  PetscReal :: o(mlp%n_out)
  PetscReal :: phi, vf, meancs, s_state, s, ramp, r_csh, r_m4, r_oth, tot, g
  PetscReal :: cs_new, u, r_x, vm_new, vm_cur, n_csh, n_m4

  call MLPTrunk(mlp, x_in, o)

  phi = x_in(1); vf = x_in(8); meancs = x_in(9)
  s_state = vf * meancs

  s = min(max((phi - mlp%phi_res) / mlp%ramp_width, 0.d0), 1.d0)
  if (mlp%ramp_form == ML_RAMP_LINEAR) then
    ramp = s
  else
    ramp = (3.d0*s*s - 2.d0*s*s*s)**2
  endif

  r_csh = mlp%rate_scale * Softplus(o(1)) * ramp
  cs_new = mlp%cs_min + (mlp%cs_max - mlp%cs_min) * Sigmoid(o(2))
  r_m4  = mlp%rate_scale * Softplus(o(3)) * ramp
  r_oth = mlp%rate_scale * Softplus(o(4)) * ramp
  tot = r_csh + r_m4 + r_oth
  if (tot > mlp%rate_max) then
    g = mlp%rate_max / tot
    r_csh = r_csh * g; r_m4 = r_m4 * g; r_oth = r_oth * g; tot = mlp%rate_max
  endif
  u = tanh(o(5))
  r_x = mlp%k_exch * (max(u, 0.d0) * (mlp%cs_max*vf - s_state) - &
                      max(-u, 0.d0) * (s_state - mlp%cs_min*vf))

  out%d_vf_csh   = r_csh
  out%d_s_csh    = cs_new * r_csh + r_x
  out%d_vf_m4ac  = r_m4
  out%d_vf_other = r_oth
  out%dphi       = -tot
  out%cs_new     = cs_new

  vm_new = (mlp%csh_vm_a + mlp%csh_vm_b * cs_new) * 1.d-6
  if (vf > 1.d-6) then
    vm_cur = (mlp%csh_vm_a + mlp%csh_vm_b * meancs) * 1.d-6
  else
    vm_cur = vm_new
  endif
  n_csh = r_csh / vm_new
  n_m4  = r_m4 / (mlp%m4ac_vm * 1.d-6)
  out%sink_Si     = n_csh
  out%sink_Ca     = cs_new * n_csh + r_x / vm_cur
  out%sink_Mg     = 4.d0 * n_m4
  out%sink_Al     = 2.d0 * n_m4
  out%sink_HCO3   = n_m4
  out%sink_TotalH = -(2.d0 * cs_new * n_csh) - 13.d0 * n_m4

end subroutine MLPForwardNodeMin

! ************************************************************************** !

elemental function Softplus(x) result(y)
  !
  ! log(1 + exp(x)), overflow-safe and accurate for very negative x
  ! (matches numpy log1p(exp(-|x|)) + max(x,0) to round-off).
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  PetscReal, intent(in) :: x
  PetscReal :: y

  y = Log1p(exp(-abs(x))) + max(x, 0.d0)

end function Softplus

! ************************************************************************** !

elemental function Log1p(x) result(y)
  !
  ! log(1 + x) without cancellation for small x (series below 1e-4).
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  PetscReal, intent(in) :: x
  PetscReal :: y

  if (abs(x) < 1.d-4) then
    y = x * (1.d0 - x * (0.5d0 - x * (1.d0/3.d0 - x * 0.25d0)))
  else
    y = log(1.d0 + x)
  endif

end function Log1p

! ************************************************************************** !

elemental function Sigmoid(x) result(y)
  !
  ! 1 / (1 + exp(-x)), overflow-safe.
  !
  ! Author: Claude Fable 5.1 (for Heeho Park)
  ! Date: 2026-09-14
  !
  implicit none

  PetscReal, intent(in) :: x
  PetscReal :: y

  if (x >= 0.d0) then
    y = 1.d0 / (1.d0 + exp(-x))
  else
    y = exp(x) / (1.d0 + exp(x))
  endif

end function Sigmoid

! ************************************************************************** !

elemental function MLPRelu(x) result(y)
  !
  ! Plain ReLU, matching torch.nn.ReLU.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  PetscReal, intent(in) :: x
  PetscReal :: y

  y = max(x, 0.d0)

end function MLPRelu

! ************************************************************************** !

subroutine MLPDestroy(mlp)
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  type(ml_mlp_type) :: mlp

  if (associated(mlp%X_mean)) deallocate(mlp%X_mean)
  if (associated(mlp%X_scale)) deallocate(mlp%X_scale)
  if (associated(mlp%clip_lo)) deallocate(mlp%clip_lo)
  if (associated(mlp%clip_hi)) deallocate(mlp%clip_hi)
  if (associated(mlp%log10_mask)) deallocate(mlp%log10_mask)
  if (associated(mlp%W0)) deallocate(mlp%W0)
  if (associated(mlp%b0)) deallocate(mlp%b0)
  if (associated(mlp%W1)) deallocate(mlp%W1)
  if (associated(mlp%b1)) deallocate(mlp%b1)
  if (associated(mlp%W2)) deallocate(mlp%W2)
  if (associated(mlp%b2)) deallocate(mlp%b2)
  if (associated(mlp%W3)) deallocate(mlp%W3)
  if (associated(mlp%b3)) deallocate(mlp%b3)
  nullify(mlp%X_mean); nullify(mlp%X_scale)
  nullify(mlp%clip_lo); nullify(mlp%clip_hi); nullify(mlp%log10_mask)
  nullify(mlp%W0); nullify(mlp%b0)
  nullify(mlp%W1); nullify(mlp%b1)
  nullify(mlp%W2); nullify(mlp%b2)
  nullify(mlp%W3); nullify(mlp%b3)

end subroutine MLPDestroy

! ************************************************************************** !

end module ML_Surrogate_Aux_module
