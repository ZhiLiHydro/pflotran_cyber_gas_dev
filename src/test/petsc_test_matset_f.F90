program test

  use petscmat

  implicit none

#include "petsc/finclude/petscmat.h"

  PetscMPIInt :: size
  PetscMPIInt :: rank
  Mat :: A
  PetscInt, parameter :: nx = 4
  PetscInt, parameter :: ny = 4
  PetscInt, parameter :: n_local = 8
  PetscInt :: id
  PetscReal :: value1
  PetscReal :: values1(1)
  PetscReal :: values2(1,1)
  PetscInt :: rows(1), cols(1)
  PetscErrorCode :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr);CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr);CHKERRQ(ierr)

  if (rank == 0) print *, 'Beginning of Fortran90 test mat program'

  call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATBAIJ,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,n_local,n_local,PETSC_DETERMINE,PETSC_DETERMINE, &
                   ierr);CHKERRQ(ierr)

  id = 1
  value1 = 1.d0
  values1 = 1.d0
  values2 = 1.d0
  rows = 0
  cols = 0
  call MatSetValue(A,id-1,id-1,value1,INSERT_VALUES,ierr);CHKERRQ(ierr)

  call MatSetValues(A,1,id-1,1,id-1,value1,INSERT_VALUES,ierr);CHKERRQ(ierr)
  ! doesn't like
  call MatSetValues(A,1,id-1,1,id-1,values1,INSERT_VALUES,ierr);CHKERRQ(ierr)
  ! doesn't like
  call MatSetValues(A,1,id-1,1,id-1,values2,INSERT_VALUES,ierr);CHKERRQ(ierr)

  call MatSetValues(A,1,rows,1,cols,values1,INSERT_VALUES,ierr);CHKERRQ(ierr)
  call MatSetValues(A,1,rows,1,cols,values2,INSERT_VALUES,ierr);CHKERRQ(ierr)

  call MatSetValuesBlocked(A,1,rows,1,cols,values1,INSERT_VALUES,ierr);CHKERRQ(ierr)
  ! doesn't like
  call MatSetValuesBlocked(A,1,rows,1,cols,values2,INSERT_VALUES,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)


  call MatDestroy(A,ierr);CHKERRA(ierr)

  if (rank == 0) print *, 'End of Fortran90 test program'
 
  call PetscFinalize(ierr)
 
end program test
