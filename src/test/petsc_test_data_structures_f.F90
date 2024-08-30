program test
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  PetscMPIInt :: size
  PetscMPIInt :: rank
  
  AO :: test_ao
  Vec :: test_vec

  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test program'

  test_vec = PETSC_NULL_VEC 
  if (test_vec == PETSC_NULL_VEC) print *, 'null vec'

  test_ao = PETSC_NULL_AO
  if (test_ao == PETSC_NULL_AO) print *, 'null ao'

  if (test_vec == 0) print *, 'null vec'

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr);CHKERRQ(ierr)
 
end program test
