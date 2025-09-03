module MPI_Utility_module

 ! Wrapper routines for simplifying the interface between PFLOTRAN and
 ! MPI

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  interface MUAllreduce
    module procedure :: MUAllreduceLogical
    module procedure :: MUAllreduceLogicalArray
    module procedure :: MUAllreduceLogicalArray2D
  end interface

  interface MUBcast
    module procedure :: MUBcastLogical
  end interface

  public :: MUAllReduce, &
            MUBcast

contains

! ************************************************************************** !

subroutine MUAllreduceLogical(lflag,operation,communicator)
  !
  ! Calculates a global reduction for a logical
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/25
  !
  implicit none

  PetscBool :: lflag
  PetscMPIInt :: operation
  PetscMPIInt :: communicator

  logical :: logical_flag
  PetscMPIInt, parameter :: one_int = 1
  PetscErrorCode :: ierr

  logical_flag = lflag
  call MPI_Allreduce(MPI_IN_PLACE,logical_flag,one_int,MPI_LOGICAL, &
                     operation,communicator,ierr);CHKERRQ(ierr)
  lflag = logical_flag

end subroutine MUAllreduceLogical

! ************************************************************************** !

subroutine MUAllreduceLogicalArray(lflag,length,operation,communicator)
  !
  ! Calculates a global reduction for an array of logicals
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/25
  !
  implicit none

  PetscBool :: lflag(:)
  PetscInt :: length
  PetscMPIInt :: operation
  PetscMPIInt :: communicator

  logical, allocatable :: logical_flag(:)
  PetscMPIInt :: mpi_int
  PetscErrorCode :: ierr

  mpi_int = length
  allocate(logical_flag(length))
  logical_flag(1:length) = lflag(1:length)
  call MPI_Allreduce(MPI_IN_PLACE,logical_flag,mpi_int,MPI_LOGICAL, &
                     operation,communicator,ierr);CHKERRQ(ierr)
  lflag(1:length) = logical_flag(1:length)
  deallocate(logical_flag)

end subroutine MUAllreduceLogicalArray

! ************************************************************************** !

subroutine MUAllreduceLogicalArray2D(lflag,length,operation,communicator)
  !
  ! Calculates a global reduction for an array of logicals
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/25
  !
  implicit none

  PetscBool :: lflag(:,:)
  PetscInt :: length
  PetscMPIInt :: operation
  PetscMPIInt :: communicator

  logical, allocatable :: logical_flag(:)
  PetscInt :: len1, len2, i, j
  PetscMPIInt :: mpi_int
  PetscErrorCode :: ierr

  len1 = size(lflag,1)
  len2 = size(lflag,2)
  if (len1*len2 /= length) then
    stop 'error in MUAllreduceLogicalArray2D'
  endif

  mpi_int = length
  allocate(logical_flag(length))
  do j = 1, len1
    do i = 1, len2
      logical_flag((j-1)*len2+i) = lflag(i,j)
    enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,logical_flag,mpi_int,MPI_LOGICAL, &
                     operation,communicator,ierr);CHKERRQ(ierr)
  do j = 1, len1
    do i = 1, len2
      lflag(i,j) = logical_flag((j-1)*len2+i)
    enddo
  enddo
  deallocate(logical_flag)

end subroutine MUAllreduceLogicalArray2D

! ************************************************************************** !

subroutine MUBcastLogical(lflag,root,communicator)
  !
  ! Broadcasts a logical to all processes in communicator
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/25
  !
  implicit none

  PetscBool :: lflag
  PetscMPIInt :: root
  PetscMPIInt :: communicator

  logical :: logical_flag
  PetscMPIInt, parameter :: one_int = 1
  PetscErrorCode :: ierr

  logical_flag = lflag
  call MPI_Bcast(logical_flag,one_int, &
                 MPI_LOGICAL,root,communicator,ierr);CHKERRQ(ierr)
  lflag = logical_flag

end subroutine MUBcastLogical

end module MPI_Utility_module
