module ML_Surrogate_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! B-huber architecture (fixed by trained model): 7 -> 256 -> 128 -> 64 -> 1
  ! ReLU on the three hidden layers, linear output.  Features are porosity,
  ! aperture, TotalH, Ca, Si, Al, Mg -- no mineral state (deadlock-free).
  ! See the ML side's handoff docs for the full contract.
  PetscInt, parameter, public :: ML_N_IN  = 7
  PetscInt, parameter, public :: ML_N_H0  = 256
  PetscInt, parameter, public :: ML_N_H1  = 128
  PetscInt, parameter, public :: ML_N_H2  = 64
  PetscInt, parameter, public :: ML_N_OUT = 1

  ! Output rate transform: rate[1/y] = ASINH_SCALE * sinh(z).  Matches the
  ! 'asinh_scale' root attribute in the weights file; the golden-vector
  ! validator flags any divergence.
  PetscReal, parameter, public :: ML_ASINH_SCALE = 1.d-5

  type, public :: ml_mlp_type
    PetscReal :: X_mean(ML_N_IN)
    PetscReal :: X_scale(ML_N_IN)
    PetscReal :: clip_lo(ML_N_IN)
    PetscReal :: clip_hi(ML_N_IN)
    PetscReal :: log10_mask(ML_N_IN)   ! >0.5 -> apply log10 after clip
    PetscReal :: y_mean
    PetscReal :: y_scale
    PetscReal :: asinh_scale
    PetscReal, pointer :: W0(:,:), b0(:)
    PetscReal, pointer :: W1(:,:), b1(:)
    PetscReal, pointer :: W2(:,:), b2(:)
    PetscReal, pointer :: W3(:,:), b3(:)
  end type ml_mlp_type

  public :: MLPCreate, &
            MLPLoad, &
            MLPForward, &
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

  mlp%X_mean = 0.d0
  mlp%X_scale = 1.d0
  mlp%clip_lo = -1.d30
  mlp%clip_hi = 1.d30
  mlp%log10_mask = 0.d0
  mlp%y_mean = 0.d0
  mlp%y_scale = 1.d0
  mlp%asinh_scale = ML_ASINH_SCALE
  nullify(mlp%W0); nullify(mlp%b0)
  nullify(mlp%W1); nullify(mlp%b1)
  nullify(mlp%W2); nullify(mlp%b2)
  nullify(mlp%W3); nullify(mlp%b3)

end subroutine MLPCreate

! ************************************************************************** !

subroutine MLPLoad(mlp, filename, option)
  !
  ! Reads C-huber weights, scalers and per-feature transforms from the HDF5
  ! file written by the NUMO chemistry surrogate training pipeline.  Layer
  ! matrices are stored in PyTorch (n_out, n_in) convention, so Fortran
  ! MATMUL(h, W) + b needs no transpose.  All ranks read the small file.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use HDF5_Aux_module
  use Option_module

  implicit none

  type(ml_mlp_type) :: mlp
  character(len=*) :: filename
  type(option_type) :: option

  integer(HID_T) :: file_id
  integer :: hdf5_err

  call HDF5FileOpenReadOnly(filename, file_id, PETSC_FALSE, &
                            'ML surrogate weights', option)

  call ReadFixedVector(file_id, '/scalers/X_mean',  ML_N_IN, mlp%X_mean,  option)
  call ReadFixedVector(file_id, '/scalers/X_scale', ML_N_IN, mlp%X_scale, option)
  call ReadScaler(file_id, '/scalers/y_mean',  mlp%y_mean,  option)
  call ReadScaler(file_id, '/scalers/y_scale', mlp%y_scale, option)

  ! Per-feature pre-transforms (int32 mask read as double via HDF5 conversion).
  call ReadFixedVector(file_id, '/transforms/clip_lo',    ML_N_IN, mlp%clip_lo,    option)
  call ReadFixedVector(file_id, '/transforms/clip_hi',    ML_N_IN, mlp%clip_hi,    option)
  call ReadFixedVector(file_id, '/transforms/log10_mask', ML_N_IN, mlp%log10_mask, option)

  mlp%asinh_scale = ML_ASINH_SCALE

  call ReadLayer(file_id, '/layers/layer_00', ML_N_H0, ML_N_IN, mlp%W0, mlp%b0, option)
  call ReadLayer(file_id, '/layers/layer_01', ML_N_H1, ML_N_H0, mlp%W1, mlp%b1, option)
  call ReadLayer(file_id, '/layers/layer_02', ML_N_H2, ML_N_H1, mlp%W2, mlp%b2, option)
  call ReadLayer(file_id, '/layers/layer_03', ML_N_OUT, ML_N_H2, mlp%W3, mlp%b3, option)

  call h5fclose_f(file_id, hdf5_err)
  if (hdf5_err /= 0) then
    option%io_buffer = 'Error closing ML surrogate weights file "' // &
                       trim(filename) // '".'
    call PrintErrMsg(option)
  endif

end subroutine MLPLoad

! ************************************************************************** !

subroutine ReadLayer(file_id, group_path, n_out, n_in, W, b, option)
  !
  ! Reads layer_NN/W and layer_NN/b.  W is stored on disk in PyTorch
  ! (n_out, n_in) C-order; held in Fortran as (n_in, n_out) so the bytes
  ! line up without a transpose, and the forward pass uses MATMUL(h, W).
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Option_module

  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: group_path
  PetscInt :: n_out, n_in
  PetscReal, pointer :: W(:,:)
  PetscReal, pointer :: b(:)
  type(option_type) :: option

  allocate(W(n_in, n_out))
  allocate(b(n_out))
  call ReadFixedMatrix(file_id, trim(group_path) // '/W', n_out, n_in, W, option)
  call ReadFixedVector(file_id, trim(group_path) // '/b', n_out,       b, option)

end subroutine ReadLayer

! ************************************************************************** !

subroutine ReadFixedVector(file_id, dataset_path, n, vec, option)
  !
  ! Reads a 1D dataset of fixed expected length n into vec(1:n).  Verifies
  ! the on-disk shape.  Integer datasets are converted to double on read.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Option_module

  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscInt :: n
  PetscReal :: vec(n)
  type(option_type) :: option

  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T) :: dims(1), max_dims(1)
  integer :: hdf5_err, rank

  call OpenDataset(file_id, dataset_path, dset_id, option)

  call h5dget_space_f(dset_id, space_id, hdf5_err)
  call h5sget_simple_extent_ndims_f(space_id, rank, hdf5_err)
  if (rank /= 1) then
    option%io_buffer = 'ML surrogate dataset "' // trim(dataset_path) // &
                       '" must be 1-D.'
    call PrintErrMsg(option)
  endif
  call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdf5_err)
  if (dims(1) /= int(n, HSIZE_T)) then
    option%io_buffer = 'ML surrogate dataset "' // trim(dataset_path) // &
                       '" has unexpected length.'
    call PrintErrMsg(option)
  endif

  dims(1) = int(n, HSIZE_T)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vec, dims, hdf5_err)
  if (hdf5_err /= 0) then
    option%io_buffer = 'Error reading ML surrogate dataset "' // &
                       trim(dataset_path) // '".'
    call PrintErrMsg(option)
  endif

  call h5sclose_f(space_id, hdf5_err)
  call h5dclose_f(dset_id, hdf5_err)

end subroutine ReadFixedVector

! ************************************************************************** !

subroutine ReadFixedMatrix(file_id, dataset_path, n_out, n_in, mat, option)
  !
  ! Reads a 2D dataset whose on-disk (C-order) shape is (n_out, n_in) into a
  ! Fortran array declared (n_in, n_out).  C row-major bytes read column-major
  ! give the element-for-element (n_in, n_out) array with no transpose copy.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Option_module

  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscInt :: n_out, n_in
  PetscReal :: mat(n_in, n_out)
  type(option_type) :: option

  integer(HID_T) :: dset_id, space_id
  integer(HSIZE_T) :: dims(2), max_dims(2)
  integer :: hdf5_err, rank

  call OpenDataset(file_id, dataset_path, dset_id, option)

  call h5dget_space_f(dset_id, space_id, hdf5_err)
  call h5sget_simple_extent_ndims_f(space_id, rank, hdf5_err)
  if (rank /= 2) then
    option%io_buffer = 'ML surrogate dataset "' // trim(dataset_path) // &
                       '" must be 2-D.'
    call PrintErrMsg(option)
  endif
  call h5sget_simple_extent_dims_f(space_id, dims, max_dims, hdf5_err)
  if (dims(1) /= int(n_in, HSIZE_T) .or. dims(2) /= int(n_out, HSIZE_T)) then
    option%io_buffer = 'ML surrogate dataset "' // trim(dataset_path) // &
                       '" has unexpected shape.'
    call PrintErrMsg(option)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mat, dims, hdf5_err)
  if (hdf5_err /= 0) then
    option%io_buffer = 'Error reading ML surrogate dataset "' // &
                       trim(dataset_path) // '".'
    call PrintErrMsg(option)
  endif

  call h5sclose_f(space_id, hdf5_err)
  call h5dclose_f(dset_id, hdf5_err)

end subroutine ReadFixedMatrix

! ************************************************************************** !

subroutine ReadScaler(file_id, dataset_path, scalar, option)
  !
  ! Reads a length-1 dataset into a scalar.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Option_module

  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscReal :: scalar
  type(option_type) :: option

  PetscReal :: tmp(1)

  call ReadFixedVector(file_id, dataset_path, 1, tmp, option)
  scalar = tmp(1)

end subroutine ReadScaler

! ************************************************************************** !

subroutine OpenDataset(file_id, dataset_path, dset_id, option)
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Option_module

  implicit none

  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  integer(HID_T) :: dset_id
  type(option_type) :: option

  integer :: hdf5_err

  call h5dopen_f(file_id, dataset_path, dset_id, hdf5_err)
  if (hdf5_err < 0) then
    option%io_buffer = 'ML surrogate dataset "' // trim(dataset_path) // &
                       '" not found.'
    call PrintErrMsg(option)
  endif

end subroutine OpenDataset

! ************************************************************************** !

subroutine MLPForward(mlp, x_in, rate)
  !
  ! 10-input, 1-output MLP forward pass.  Applies per-feature clip then
  ! log10 (where masked), standardizes, runs three ReLU hidden layers with
  ! a linear output, then the asinh rate transform:
  !   rate[1/y] = asinh_scale * sinh(y_norm*y_scale + y_mean)
  ! The signed rate is returned as-is; the caller applies the clogging clamp,
  ! envelope clamp and Euler update.  Matches chuber_reference.py exactly.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  type(ml_mlp_type) :: mlp
  PetscReal :: x_in(ML_N_IN)
  PetscReal :: rate

  PetscReal :: x(ML_N_IN)
  PetscReal :: h_in(ML_N_IN)
  PetscReal :: h0(ML_N_H0)
  PetscReal :: h1(ML_N_H1)
  PetscReal :: h2(ML_N_H2)
  PetscReal :: z
  PetscInt :: i

  do i = 1, ML_N_IN
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

  z = (dot_product(h2, mlp%W3(:,1)) + mlp%b3(1)) * mlp%y_scale + mlp%y_mean

  rate = mlp%asinh_scale * sinh(z)

end subroutine MLPForward

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

  if (associated(mlp%W0)) deallocate(mlp%W0)
  if (associated(mlp%b0)) deallocate(mlp%b0)
  if (associated(mlp%W1)) deallocate(mlp%W1)
  if (associated(mlp%b1)) deallocate(mlp%b1)
  if (associated(mlp%W2)) deallocate(mlp%W2)
  if (associated(mlp%b2)) deallocate(mlp%b2)
  if (associated(mlp%W3)) deallocate(mlp%W3)
  if (associated(mlp%b3)) deallocate(mlp%b3)
  nullify(mlp%W0); nullify(mlp%b0)
  nullify(mlp%W1); nullify(mlp%b1)
  nullify(mlp%W2); nullify(mlp%b2)
  nullify(mlp%W3); nullify(mlp%b3)

end subroutine MLPDestroy

! ************************************************************************** !

end module ML_Surrogate_Aux_module
