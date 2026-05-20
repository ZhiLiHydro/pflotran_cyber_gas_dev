program validate_ml_surrogate
  !
  ! Standalone validator for the ML surrogate forward pass.  Reads the
  ! trained weights directly via the HDF5 Fortran API to keep the driver
  ! independent of PFLOTRAN init plumbing, populates a ml_mlp_type, then
  ! calls MLPForward on the three canned inputs from the design doc.
  ! A companion Python script reproduces the forward pass with h5py +
  ! numpy; |Delta phi| should be below 1e-10.
  !
  ! Run:
  !   make validate_ml_surrogate
  !   ./validate_ml_surrogate <path/to/model_mlp_weights.h5>
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use ML_Surrogate_Aux_module

  implicit none

  type(ml_mlp_type) :: mlp
  character(len=1024) :: filename
  PetscReal :: x_in(ML_N_IN), phi
  PetscErrorCode :: ierr
  integer :: nargs, hdf5_err
  integer(HID_T) :: file_id

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRA(ierr)

  nargs = command_argument_count()
  if (nargs < 1) then
    write(*,'(a)') 'Usage: validate_ml_surrogate <path/to/model_mlp_weights.h5>'
    call PetscFinalize(ierr); CHKERRA(ierr)
    stop 1
  endif
  call get_command_argument(1, filename)

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

  ! Canned inputs from design doc section 8.
  ! Order: log10(v_darcy), log10(aperture), time_years, r_norm
  x_in = [-9.0d0, -3.0d0, 500.0d0, 0.3d0]
  call MLPForward(mlp, x_in, phi)
  write(*,'(a,4(f12.4,1x),a,es25.17)') 'input ', x_in, '-> phi = ', phi

  x_in = [-6.0d0, -2.0d0, 100.0d0, 0.1d0]
  call MLPForward(mlp, x_in, phi)
  write(*,'(a,4(f12.4,1x),a,es25.17)') 'input ', x_in, '-> phi = ', phi

  x_in = [-12.0d0, -5.0d0, 900.0d0, 0.9d0]
  call MLPForward(mlp, x_in, phi)
  write(*,'(a,4(f12.4,1x),a,es25.17)') 'input ', x_in, '-> phi = ', phi

  call MLPDestroy(mlp)
  call PetscFinalize(ierr); CHKERRA(ierr)

contains

! ************************************************************************** !

subroutine LoadWeights(file_id, mlp)
  !
  ! Validator-private raw-HDF5 weight loader so we can exercise MLPForward
  ! without dragging in OptionCreate/Comm/Driver plumbing.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
  !
  implicit none

  integer(HID_T) :: file_id
  type(ml_mlp_type) :: mlp

  call ReadVec(file_id, '/scalers/X_mean',  ML_N_IN, mlp%X_mean)
  call ReadVec(file_id, '/scalers/X_scale', ML_N_IN, mlp%X_scale)
  call ReadScalar(file_id, '/scalers/y_mean',  mlp%y_mean)
  call ReadScalar(file_id, '/scalers/y_scale', mlp%y_scale)

  call ReadLayer(file_id, '/layers/layer_00', ML_N_H0, ML_N_IN, mlp%W0, mlp%b0)
  call ReadLayer(file_id, '/layers/layer_01', ML_N_H1, ML_N_H0, mlp%W1, mlp%b1)
  call ReadLayer(file_id, '/layers/layer_02', ML_N_H2, ML_N_H1, mlp%W2, mlp%b2)
  call ReadLayer(file_id, '/layers/layer_03', ML_N_H3, ML_N_H2, mlp%W3, mlp%b3)
  call ReadLayer(file_id, '/layers/layer_04', ML_N_OUT, ML_N_H3, mlp%W4, mlp%b4)

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
