program validate_ml_surrogate
  !
  ! Standalone bit-verification of the C-huber MLP forward pass.  Loads the
  ! trained weights directly via the HDF5 Fortran API (no PFLOTRAN init
  ! plumbing), then runs MLPForward on every row of the golden test-vector
  ! CSV and compares the predicted rate against the expected column.  The
  ! forward pass must reproduce chuber_reference.py to ~1e-9 relative.
  !
  ! Run:
  !   make validate_ml_surrogate
  !   ./validate_ml_surrogate <weights.h5> <chuber_test_vectors.csv>
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use ML_Surrogate_Aux_module

  implicit none

  type(ml_mlp_type) :: mlp
  character(len=1024) :: filename, csvname
  character(len=4096) :: line, hdr
  PetscReal :: x_in(ML_N_IN), rate, expected, vals(64), rel, worst
  PetscErrorCode :: ierr
  integer :: nargs, hdf5_err, ios, i, nrows, ncol, erate_col
  integer(HID_T) :: file_id

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRA(ierr)

  nargs = command_argument_count()
  if (nargs < 2) then
    write(*,'(a)') 'Usage: validate_ml_surrogate <weights.h5> <test_vectors.csv>'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 1
  endif
  call get_command_argument(1, filename)
  call get_command_argument(2, csvname)

  call h5open_f(hdf5_err)
  call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdf5_err)
  if (hdf5_err < 0) then
    write(*,'(a,a,a)') 'Cannot open HDF5 file "', trim(filename), '".'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 2
  endif

  call MLPCreate(mlp)
  call LoadWeights(file_id, mlp)

  call h5fclose_f(file_id, hdf5_err)
  call h5close_f(hdf5_err)

  open(unit=10, file=trim(csvname), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,'(a,a,a)') 'Cannot open CSV "', trim(csvname), '".'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 3
  endif

  ! Header-driven column layout: features are cols 4..3+ML_N_IN; the target
  ! column is located by name so C-huber and B-huber CSVs both work.
  read(10,'(a)',iostat=ios) hdr    ! header (comma-separated)
  ncol = 1
  do i = 1, len_trim(hdr)
    if (hdr(i:i) == ',') ncol = ncol + 1
  enddo
  erate_col = 1
  do i = 1, index(hdr,'expected_rate_per_yr') - 1
    if (hdr(i:i) == ',') erate_col = erate_col + 1
  enddo

  worst = 0.d0
  nrows = 0
  do
    read(10,'(a)',iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    do i = 1, len_trim(line)
      if (line(i:i) == ',') line(i:i) = ' '
    enddo
    read(line,*,iostat=ios) vals(1:ncol)
    if (ios /= 0) cycle

    x_in = vals(4:3+ML_N_IN)         ! features, idx order 1..ML_N_IN
    expected = vals(erate_col)       ! expected_rate_per_yr
    call MLPForward(mlp, x_in, rate)

    rel = abs(rate - expected) / max(abs(expected), 1.d-30)
    worst = max(worst, rel)
    nrows = nrows + 1
  enddo
  close(10)

  write(*,'(i0,a,es12.4)') nrows, ' vectors checked; worst relative error = ', worst
  if (worst < 1.d-9) then
    write(*,'(a)') 'PASS -- Fortran forward pass matches the golden vectors.'
  else
    write(*,'(a)') 'FAIL -- transform chain broken (check log10 mask, clip '// &
                   'bounds, feature order, TotalH sign, asinh_scale).'
  endif

  call MLPDestroy(mlp)
  call PetscFinalize(ierr); CHKERRA(ierr)
  if (worst >= 1.d-9) stop 4

contains

! ************************************************************************** !

subroutine LoadWeights(file_id, mlp)
  !
  ! Validator-private raw-HDF5 loader so MLPForward can be exercised without
  ! OptionCreate/Comm/Driver plumbing.  Mirrors MLPLoad's dataset contract.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  integer(HID_T) :: file_id
  type(ml_mlp_type) :: mlp

  call ReadVec(file_id, '/scalers/X_mean',  ML_N_IN, mlp%X_mean)
  call ReadVec(file_id, '/scalers/X_scale', ML_N_IN, mlp%X_scale)
  call ReadScalar(file_id, '/scalers/y_mean',  mlp%y_mean)
  call ReadScalar(file_id, '/scalers/y_scale', mlp%y_scale)

  call ReadVec(file_id, '/transforms/clip_lo',    ML_N_IN, mlp%clip_lo)
  call ReadVec(file_id, '/transforms/clip_hi',    ML_N_IN, mlp%clip_hi)
  call ReadVec(file_id, '/transforms/log10_mask', ML_N_IN, mlp%log10_mask)
  mlp%asinh_scale = ML_ASINH_SCALE

  call ReadLayer(file_id, '/layers/layer_00', ML_N_H0, ML_N_IN, mlp%W0, mlp%b0)
  call ReadLayer(file_id, '/layers/layer_01', ML_N_H1, ML_N_H0, mlp%W1, mlp%b1)
  call ReadLayer(file_id, '/layers/layer_02', ML_N_H2, ML_N_H1, mlp%W2, mlp%b2)
  call ReadLayer(file_id, '/layers/layer_03', ML_N_OUT, ML_N_H2, mlp%W3, mlp%b3)

end subroutine LoadWeights

! ************************************************************************** !

subroutine ReadLayer(file_id, group_path, n_out, n_in, W, b)
  implicit none
  integer(HID_T) :: file_id
  character(len=*) :: group_path
  PetscInt :: n_out, n_in
  PetscReal, pointer :: W(:,:)
  PetscReal, pointer :: b(:)

  allocate(W(n_in, n_out))
  allocate(b(n_out))
  call ReadMat(file_id, trim(group_path) // '/W', n_out, n_in, W)
  call ReadVec(file_id, trim(group_path) // '/b', n_out,       b)
end subroutine ReadLayer

! ************************************************************************** !

subroutine ReadVec(file_id, dataset_path, n, vec)
  implicit none
  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscInt :: n
  PetscReal :: vec(n)

  integer(HID_T) :: dset_id
  integer(HSIZE_T) :: dims(1)
  integer :: hdf5_err

  dims(1) = int(n, HSIZE_T)
  call h5dopen_f(file_id, dataset_path, dset_id, hdf5_err)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vec, dims, hdf5_err)
  call h5dclose_f(dset_id, hdf5_err)
end subroutine ReadVec

! ************************************************************************** !

subroutine ReadMat(file_id, dataset_path, n_out, n_in, mat)
  implicit none
  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscInt :: n_out, n_in
  ! On-disk shape is (n_out, n_in) C-order; Fortran (col-major) reads
  ! those bytes into a (n_in, n_out) array element-for-element.
  PetscReal :: mat(n_in, n_out)

  integer(HID_T) :: dset_id
  integer(HSIZE_T) :: dims(2)
  integer :: hdf5_err

  dims(1) = int(n_in,  HSIZE_T)
  dims(2) = int(n_out, HSIZE_T)
  call h5dopen_f(file_id, dataset_path, dset_id, hdf5_err)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mat, dims, hdf5_err)
  call h5dclose_f(dset_id, hdf5_err)
end subroutine ReadMat

! ************************************************************************** !

subroutine ReadScalar(file_id, dataset_path, scalar)
  implicit none
  integer(HID_T) :: file_id
  character(len=*) :: dataset_path
  PetscReal :: scalar

  PetscReal :: tmp(1)
  call ReadVec(file_id, dataset_path, 1, tmp)
  scalar = tmp(1)
end subroutine ReadScalar

end program validate_ml_surrogate
