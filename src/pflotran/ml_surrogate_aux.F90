module ML_Surrogate_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use hdf5
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! Architecture is fixed by the trained model: 4 -> 64 -> 128 -> 128 -> 64 -> 1
  PetscInt, parameter, public :: ML_N_IN  = 4
  PetscInt, parameter, public :: ML_N_H0  = 64
  PetscInt, parameter, public :: ML_N_H1  = 128
  PetscInt, parameter, public :: ML_N_H2  = 128
  PetscInt, parameter, public :: ML_N_H3  = 64
  PetscInt, parameter, public :: ML_N_OUT = 1

  PetscReal, parameter :: ML_PHI_MIN = 0.02d0
  PetscReal, parameter :: ML_PHI_MAX = 0.99d0

  type, public :: ml_mlp_type
    PetscReal :: X_mean(ML_N_IN)
    PetscReal :: X_scale(ML_N_IN)
    PetscReal :: y_mean
    PetscReal :: y_scale
    PetscReal, pointer :: W0(:,:), b0(:)
    PetscReal, pointer :: W1(:,:), b1(:)
    PetscReal, pointer :: W2(:,:), b2(:)
    PetscReal, pointer :: W3(:,:), b3(:)
    PetscReal, pointer :: W4(:,:), b4(:)
  end type ml_mlp_type

  public :: MLPCreate, &
            MLPLoad, &
            MLPForward, &
            MLPDestroy, &
            MLPGelu

contains

! ************************************************************************** !

subroutine MLPCreate(mlp)
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
  !
  implicit none
  type(ml_mlp_type) :: mlp

  mlp%X_mean = 0.d0
  mlp%X_scale = 1.d0
  mlp%y_mean = 0.d0
  mlp%y_scale = 1.d0
  nullify(mlp%W0); nullify(mlp%b0)
  nullify(mlp%W1); nullify(mlp%b1)
  nullify(mlp%W2); nullify(mlp%b2)
  nullify(mlp%W3); nullify(mlp%b3)
  nullify(mlp%W4); nullify(mlp%b4)

end subroutine MLPCreate

! ************************************************************************** !

subroutine MLPLoad(mlp, filename, option)
  !
  ! Reads MLP weights and scalers from an HDF5 file written by the NUMO
  ! fracture-clogging surrogate training pipeline.  Layer matrices are stored
  ! in PyTorch (n_out, n_in) convention, so Fortran MATMUL(W, h) + b needs no
  ! transpose.  All ranks read the (small) file independently.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
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

  call ReadLayer(file_id, '/layers/layer_00', ML_N_H0, ML_N_IN, mlp%W0, mlp%b0, option)
  call ReadLayer(file_id, '/layers/layer_01', ML_N_H1, ML_N_H0, mlp%W1, mlp%b1, option)
  call ReadLayer(file_id, '/layers/layer_02', ML_N_H2, ML_N_H1, mlp%W2, mlp%b2, option)
  call ReadLayer(file_id, '/layers/layer_03', ML_N_H3, ML_N_H2, mlp%W3, mlp%b3, option)
  call ReadLayer(file_id, '/layers/layer_04', ML_N_OUT, ML_N_H3, mlp%W4, mlp%b4, option)

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
  ! Reads layer_NN/W and layer_NN/b from an open H5 file.  Allocates W
  ! and b.  W is stored on disk in PyTorch (n_out, n_in) C-order layout;
  ! we hold it in Fortran as (n_in, n_out) so the byte interpretation
  ! lines up without an explicit transpose.  Math-wise this means
  ! W_fortran(j,i) == W_pytorch(i,j); the forward pass uses MATMUL with
  ! h on the left to compensate.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
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
  ! Reads a 1D float64 dataset of fixed expected length n into vec(1:n).
  ! Verifies the on-disk shape; errors out via PrintErrMsg on mismatch.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
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
  ! Reads a 2D float64 dataset whose on-disk (C-order) shape is
  ! (n_out, n_in) into a Fortran array declared (n_in, n_out).  HDF5
  ! writes bytes in C row-major order and Fortran reads them in
  ! column-major order, so a (n_out, n_in) C-order dataset becomes a
  ! (n_in, n_out) Fortran array element-for-element.  No transpose copy.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
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
  ! Reads a length-1 float64 dataset (y_mean / y_scale) into a scalar.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
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
  ! Date: 2026-05-18
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

subroutine MLPForward(mlp, x_in, phi)
  !
  ! 4-input, 1-output MLP forward pass with GELU on the four hidden layers
  ! and no activation on the output.  Output is denormalized via
  ! (phi = phi_norm * y_scale + y_mean) and clipped to [0.02, 0.99].
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
  !
  implicit none

  type(ml_mlp_type) :: mlp
  PetscReal :: x_in(ML_N_IN)
  PetscReal :: phi

  PetscReal :: h_in(ML_N_IN)
  PetscReal :: h0(ML_N_H0)
  PetscReal :: h1(ML_N_H1)
  PetscReal :: h2(ML_N_H2)
  PetscReal :: h3(ML_N_H3)
  PetscReal :: phi_norm

  h_in = (x_in - mlp%X_mean) / mlp%X_scale

  ! W is stored as (n_in, n_out) (see ReadFixedMatrix), so the math
  ! "Wh + b" becomes matmul(h, W) + b in Fortran.
  h0 = MLPGelu(MATMUL(h_in, mlp%W0) + mlp%b0)
  h1 = MLPGelu(MATMUL(h0,   mlp%W1) + mlp%b1)
  h2 = MLPGelu(MATMUL(h1,   mlp%W2) + mlp%b2)
  h3 = MLPGelu(MATMUL(h2,   mlp%W3) + mlp%b3)

  phi_norm = dot_product(h3, mlp%W4(:,1)) + mlp%b4(1)

  phi = phi_norm * mlp%y_scale + mlp%y_mean
  phi = max(ML_PHI_MIN, min(ML_PHI_MAX, phi))

end subroutine MLPForward

! ************************************************************************** !

elemental function MLPGelu(x) result(y)
  !
  ! Tanh-approximation GELU, matching PyTorch's nn.GELU(approximate='tanh').
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
  !
  implicit none

  PetscReal, intent(in) :: x
  PetscReal :: y

  PetscReal, parameter :: c1 = 0.79788456080286535587d0  ! sqrt(2/pi)
  PetscReal, parameter :: c2 = 0.044715d0

  y = 0.5d0 * x * (1.d0 + tanh(c1 * (x + c2 * x**3)))

end function MLPGelu

! ************************************************************************** !

subroutine MLPDestroy(mlp)
  !
  ! Author: Heeho Park
  ! Date: 2026-05-18
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
  if (associated(mlp%W4)) deallocate(mlp%W4)
  if (associated(mlp%b4)) deallocate(mlp%b4)
  nullify(mlp%W0); nullify(mlp%b0)
  nullify(mlp%W1); nullify(mlp%b1)
  nullify(mlp%W2); nullify(mlp%b2)
  nullify(mlp%W3); nullify(mlp%b3)
  nullify(mlp%W4); nullify(mlp%b4)

end subroutine MLPDestroy

! ************************************************************************** !

end module ML_Surrogate_Aux_module
