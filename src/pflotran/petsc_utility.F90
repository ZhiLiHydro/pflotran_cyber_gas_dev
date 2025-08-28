module Petsc_Utility_module

 ! Wrapper routines for simplifying the interface between PFLOTRAN and
 ! PETSc

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module

  implicit none

  private

  interface PetUtilLoadVec
    module procedure PetUtilLoadVecInt
    module procedure PetUtilLoadVecReal
  end interface

  interface PetUtilUnloadVec
    module procedure PetUtilUnloadVecInt
    module procedure PetUtilUnloadVecReal
  end interface

  interface PUMSetValues
    module procedure PUMSetValues1
    module procedure PUMSetValues2
    module procedure PUMSetValues3
  end interface

  interface PUMSetValuesLocal
    module procedure PUMSetValuesLocal1
    module procedure PUMSetValuesLocal2
    module procedure PUMSetValuesLocal3
  end interface

  public :: PetUtilMatSVBL, &
            PetUtilVecSVBL, &
            PetUtilLoadVec, &
            PetUtilUnloadVec, &
            PetscUtilCompareMatrices

  !TODO(geh324): clean up the below, potentially eliminating many
  public :: PUMSetValue, &
            PUMSetValues, &
            PUMSetValuesLocal, &
            PUMSetValuesBlocked, &
            PUMSetValuesBlockedLocal

  public :: PUVGetArray, &
            PUVRestoreArray, &
            PUVGetArrayRead, &
            PUVRestoreArrayRead, &
            PUMSeqAIJGetArray, &
            PUMSeqAIJRestoreArray, &
            PUMGetRowIJ, &
            PUMRestoreRowIJ, &
            PUMGetValuesLocal, &
            PUIGetIndices, &
            PUIRestoreIndices, &
            PUVDuplicateVecs, &
            PUVDestroyVecs, &
            PetscTestFile!, &
!            PetscTime

contains

! ************************************************************************** !

subroutine PetUtilMatSVBL(A,irow,icol,matrix_block,ndof)
  !
  ! Maps values from ndof*ndof sized array into ncell*ndof blocked matrix
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/22
  !
  implicit none

  Mat :: A
  PetscInt :: irow
  PetscInt :: icol
  PetscReal :: matrix_block(:,:)
  PetscInt :: ndof

  PetscReal :: ndof_mat(ndof,ndof)
  PetscErrorCode :: ierr

  ndof_mat = matrix_block(1:ndof,1:ndof)

  call MatSetValuesBlockedLocal(A,1,[irow-1],1,[icol-1], &
                                reshape(ndof_mat,(/size(ndof_mat)/)), &
                                ADD_VALUES,ierr);CHKERRQ(ierr)

end subroutine PetUtilMatSVBL

! ************************************************************************** !

subroutine PetUtilVecSVBL(array,icell,array_block,ndof,replace)
  !
  ! Maps values from ndof sized array into ncell*ndof blocked array
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/22
  !
  implicit none

  PetscInt :: ndof
  PetscReal :: array(:)
  PetscInt :: icell
  PetscReal :: array_block(:)
  PetscBool :: replace

  PetscInt :: istart, iend

  iend = icell*ndof
  istart = iend-ndof+1
  if (replace) then
    array(istart:iend) = array_block(1:ndof)
  else
    array(istart:iend) = array(istart:iend) + array_block(1:ndof)
  endif

end subroutine PetUtilVecSVBL

! ************************************************************************** !

subroutine PetUtilLoadVecInt(vec,iarray)
  !
  ! Loads values from array into a vec (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscInt :: iarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArray(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = iarray(:)
  call VecRestoreArray(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilLoadVecInt

! ************************************************************************** !

subroutine PetUtilLoadVecReal(vec,rarray)
  !
  ! Loads values from array into a vec (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscReal :: rarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArray(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = rarray(:)
  call VecRestoreArray(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilLoadVecReal

! ************************************************************************** !

subroutine PetUtilUnloadVecInt(vec,iarray)
  !
  ! Unloads values from a vec into array (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscInt :: iarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)
  iarray(:) = nint(vec_ptr(:))
  call VecRestoreArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilUnloadVecInt

! ************************************************************************** !

subroutine PetUtilUnloadVecReal(vec,rarray)
  !
  ! Unloads values from a vec into array (must be sized identically)
  !
  ! Author: Glenn Hammond
  ! Date: 06/03/24
  !
  implicit none

  Vec :: vec
  PetscReal :: rarray(:)

  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  call VecGetArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)
  rarray(:) = vec_ptr(:)
  call VecRestoreArrayRead(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PetUtilUnloadVecReal

! ************************************************************************** !

subroutine PetscUtilCompareMatrices(A,B,nL2G,nG2A,rtol,row_rtol,option)
  !
  ! Compares values in matrices pinpointing differences.
  !
  ! Author: Glenn Hammond
  ! Date: 01/02/25
  !
  use Option_module
  use String_module

  implicit none

  Mat :: A, B
  PetscInt :: nL2G(:) ! from grid%nL2G
  PetscInt :: nG2A(:)
  PetscReal :: rtol
  PetscReal :: row_rtol
  type(option_type) :: option

  PetscInt :: nrow, ncol
  PetscInt, pointer :: irow_ptr(:), icol_ptr(:)
  PetscReal :: values_a(1000), values_b(1000)
  PetscReal :: value_a, value_b, scale, row_scale
  PetscBool :: success
  PetscInt :: i, j
  PetscInt :: irow(1)
  PetscInt :: irow_local0, icol_local0
  PetscInt :: irow_cell_local, icol_cell_local
  PetscInt :: row_cell_natural, col_cell_natural
  PetscInt :: irow_natural, icol_natural
  PetscInt :: block_size
  PetscErrorCode :: ierr

  call MatGetRowIJ(A,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                   nrow,irow_ptr,icol_ptr,success,ierr);CHKERRQ(ierr)
  if (.not.success) then
    option%io_buffer = 'Error returned from MatGetRowIJF90() in &
      &PetscUtilCompareMatrices. I believe that the type of matrix A is not &
      &supported by that routine.'
    call PrintErrMsg(option)
  endif
  call MatGetBlockSize(A,block_size,ierr);CHKERRQ(ierr)

  do i = 1, nrow
    irow(1) = i-1
    irow_local0 = i-1
    irow_cell_local = irow_local0/block_size+1
    row_cell_natural = nG2A(nL2G(irow_cell_local))
    irow_natural = (row_cell_natural-1)*block_size+mod(i-1,block_size)+1
    ncol = irow_ptr(i+1)-irow_ptr(i)
    call MatGetValuesLocal(A,ONE_INTEGER,[irow],ncol, &
                           icol_ptr(irow_ptr(i)+1:irow_ptr(i+1)), &
                           values_a,ierr);CHKERRQ(ierr)
    call MatGetValuesLocal(B,ONE_INTEGER,[irow],ncol, &
                           icol_ptr(irow_ptr(i)+1:irow_ptr(i+1)), &
                           values_b,ierr);CHKERRQ(ierr)
    do j = 1, ncol
      value_a = values_a(j)
      value_b = values_b(j)
      row_scale = max(maxval(dabs(values_a(1:ncol))), &
                      maxval(dabs(values_b(1:ncol))))
      scale = max(dabs(value_a),dabs(value_b))
      if (scale > 0.d0) then
        if (dabs(value_a - value_b)/scale > rtol .and. &
            dabs(value_a - value_b)/row_scale > row_rtol) then
          icol_local0 = (icol_ptr(irow_ptr(i)+j)+1)-1
          icol_cell_local = icol_local0/block_size+1
          col_cell_natural = nG2A(nL2G(icol_cell_local))
          icol_natural = (col_cell_natural-1)*block_size+mod(j-1,block_size)+1
          option%io_buffer = 'Large difference - (' // &
            StringWrite(irow_natural) // ',' // &
            StringWrite(icol_natural) // ') : ' // &
            StringWrite(value_a) // ' vs ' // StringWrite(value_b)
          call PrintMsgByRank(option)
        endif
      endif
    enddo
  enddo
  call MatRestoreRowIJ(A,ZERO_INTEGER,PETSC_FALSE,PETSC_FALSE, &
                       nrow,irow_ptr,icol_ptr,success,ierr);CHKERRQ(ierr)

end subroutine PetscUtilCompareMatrices

! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !

#if (PETSC_VERSION_MINOR > 21)
#define REFACTOR_323
#endif

subroutine PUVGetArray(vector,vec_ptr,ierr)
  !
  ! Maps to VecGetArray()/VecGetArrayF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  Vec :: vector
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call VecGetArray(vector,vec_ptr,ierr);CHKERRQ(ierr)
#else
  call VecGetArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUVGetArray

! ************************************************************************** !

subroutine PUVRestoreArray(vector,vec_ptr,ierr)
  !
  ! Maps to VecGetArray()/VecGetArrayF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  Vec :: vector
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call VecRestoreArray(vector,vec_ptr,ierr);CHKERRQ(ierr)
#else
  call VecRestoreArrayF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUVRestoreArray

! ************************************************************************** !

subroutine PUVGetArrayRead(vector,vec_ptr,ierr)
  !
  ! Maps to VecGetArray()/VecGetArrayF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  Vec :: vector
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call VecGetArrayRead(vector,vec_ptr,ierr);CHKERRQ(ierr)
#else
  call VecGetArrayReadF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUVGetArrayRead

! ************************************************************************** !

subroutine PUVRestoreArrayRead(vector,vec_ptr,ierr)
  !
  ! Maps to VecGetArray()/VecGetArrayF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  Vec :: vector
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call VecRestoreArrayRead(vector,vec_ptr,ierr);CHKERRQ(ierr)
#else
  call VecRestoreArrayReadF90(vector,vec_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUVRestoreArrayRead

! ************************************************************************** !

subroutine PUMSeqAIJGetArray(matrix,r_ptr,ierr)
  !
  ! Maps to PUMatSeqAIJGetArray()/PUMatSeqAIJGetArrayF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: matrix
  PetscReal, pointer :: r_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSeqAIJGetArray(matrix,r_ptr,ierr);CHKERRQ(ierr)
#else
  call MatSeqAIJGetArrayF90(matrix,r_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSeqAIJGetArray

! ************************************************************************** !

subroutine PUMSeqAIJRestoreArray(matrix,r_ptr,ierr)
  !
  ! Maps to PUMatSeqAIJRestoreArray()/PUMatSeqAIJRestoreArrayF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: matrix
  PetscReal, pointer :: r_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSeqAIJGetArray(matrix,r_ptr,ierr);CHKERRQ(ierr)
#else
  call MatSeqAIJGetArrayF90(matrix,r_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSeqAIJRestoreArray

! ************************************************************************** !

subroutine PUMGetRowIJ(matrix,shift,is_symmetric,decompressed, &
                       n,ia_p,ja_p,done,ierr)
  !
  ! Maps to PUMatGetRowIJ()/PUMatGetRowIJF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: matrix
  PetscInt :: shift
  PetscBool :: is_symmetric
  PetscBool :: decompressed
  PetscInt :: n
  PetscInt, pointer :: ia_p(:)
  PetscInt, pointer :: ja_p(:)
  PetscBool :: done
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatGetRowIJ(matrix,shift,is_symmetric,decompressed, &
                      n,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
#else
  call MatGetRowIJF90(matrix,shift,is_symmetric,decompressed, &
                      n,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMGetRowIJ

! ************************************************************************** !

subroutine PUMRestoreRowIJ(matrix,shift,is_symmetric,decompressed, &
                           n,ia_p,ja_p,done,ierr)
  !
  ! Maps to PUMatGetRowIJ()/PUMatGetRowIJF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: matrix
  PetscInt :: shift
  PetscBool :: is_symmetric
  PetscBool :: decompressed
  PetscInt :: n
  PetscInt, pointer :: ia_p(:)
  PetscInt, pointer :: ja_p(:)
  PetscBool :: done
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatRestoreRowIJ(matrix,shift,is_symmetric,decompressed, &
                          n,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
#else
  call MatRestoreRowIJF90(matrix,shift,is_symmetric,decompressed, &
                          n,ia_p,ja_p,done,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMRestoreRowIJ

! ************************************************************************** !

subroutine PUVDuplicateVecs(template_vec,num_vecs,vec_array_ptr,ierr)
  !
  ! Maps to PUMatGetRowIJ()/PUMatGetRowIJF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  Vec :: template_vec
  PetscInt :: num_vecs
  Vec, pointer :: vec_array_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call VecDuplicateVecs(template_vec,num_vecs,vec_array_ptr, &
                           ierr);CHKERRQ(ierr)
#else
  call VecDuplicateVecsF90(template_vec,num_vecs,vec_array_ptr, &
                           ierr);CHKERRQ(ierr)
#endif

end subroutine PUVDuplicateVecs

! ************************************************************************** !

subroutine PUVDestroyVecs(num_vecs,vec_array_ptr,ierr)
  !
  ! Maps to PUMatGetRowIJ()/PUMatGetRowIJF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscvec.h"
  use petscvec
  implicit none

  PetscInt :: num_vecs
  Vec, pointer :: vec_array_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call VecDestroyVecs(num_vecs,vec_array_ptr,ierr);CHKERRQ(ierr)
#else
  call VecDestroyVecsF90(num_vecs,vec_array_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUVDestroyVecs

! ************************************************************************** !

subroutine PUIGetIndices(is,is_ptr,ierr)
  !
  ! Maps to PUMatGetRowIJ()/PUMatGetRowIJF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscis.h"
  use petscis
  implicit none

  IS :: is
  PetscInt, pointer :: is_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call ISGetIndices(is,is_ptr,ierr);CHKERRQ(ierr)
#else
  call ISGetIndicesF90(is,is_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUIGetIndices

! ************************************************************************** !

subroutine PUIRestoreIndices(is,is_ptr,ierr)
  !
  ! Maps to PUMatGetRowIJ()/PUMatGetRowIJF90
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscis.h"
  use petscis
  implicit none

  IS :: is
  PetscInt, pointer :: is_ptr(:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call ISRestoreIndices(is,is_ptr,ierr);CHKERRQ(ierr)
#else
  call ISRestoreIndicesF90(is,is_ptr,ierr);CHKERRQ(ierr)
#endif

end subroutine PUIRestoreIndices

! ************************************************************************** !

subroutine PUMSetValue(A,irow,icol,scalar, &
                       insert_mode,ierr)
  !
  ! Maps to MatSetValue()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25
#include "petsc/finclude/petscmat.h"
  use petscmat

  implicit none

  Mat :: A
  PetscInt :: irow
  PetscInt :: icol
  PetscReal :: scalar
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSetValue(A,irow,icol,scalar,insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValue(A,irow,icol,scalar,insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValue

! ************************************************************************** !

subroutine PUMSetValues1(A,nrow,irow,ncol,icol,scalar, &
                         insert_mode,ierr)
  !
  ! Maps to MatSetValues()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: scalar
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

  PetscInt :: irow1d(1)
  PetscInt :: icol1d(1)
  PetscReal :: array1d(1)

  irow1d(:) = irow
  icol1d(:) = icol
  array1d(1) = scalar

#if defined(REFACTOR_323)
stop 'PUMSetValues1'
#else
  call MatSetValues(A,nrow,irow1d,ncol,icol1d,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValues1

! ************************************************************************** !

subroutine PUMSetValues2(A,nrow,irow,ncol,icol,array1d, &
                         insert_mode,ierr)
  !
  ! Maps to MatSetValues()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: array1d(:)
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSetValues(A,nrow,[irow],ncol,icol,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValues(A,nrow,irow1d,ncol,icol,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValues2

! ************************************************************************** !

subroutine PUMSetValues3(A,nrow,irow,ncol,icol,array1d, &
                         insert_mode,ierr)
  !
  ! Maps to MatSetValues()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow(:)
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: array1d(:)
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSetValues(A,nrow,irow,ncol,icol,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValues(A,nrow,irow,ncol,icol,array1d, &
                    insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValues3

! ************************************************************************** !

subroutine PUMSetValuesLocal1(A,nrow,irow,ncol,icol,scalar, &
                              insert_mode,ierr)
  !
  ! Maps to MatSetValuesLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: scalar
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

  PetscReal :: array2d(1,1)

  array2d(1,1) = scalar

#if defined(REFACTOR_323)
  call MatSetValuesLocal(A,nrow,[irow],ncol,[icol],[array2d(1,1)], &
                         insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValuesLocal(A,nrow,irow,ncol,icol,array2d, &
                         insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValuesLocal1

! ************************************************************************** !

subroutine PUMSetValuesLocal2(A,nrow,irow,ncol,icol,array2d, &
                              insert_mode,ierr)
  !
  ! Maps to MatSetValuesLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: array2d(:,:)
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSetValuesLocal(A,nrow,[irow],ncol,[icol],[array2d(1,1)], &
                         insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValuesLocal(A,nrow,irow,ncol,icol,array2d, &
                         insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValuesLocal2

! ************************************************************************** !

subroutine PUMSetValuesLocal3(A,nrow,irow,ncol,icol,array1d, &
                              insert_mode,ierr)
  !
  ! Maps to MatSetValuesLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: array1d(:)
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

  PetscReal :: irow1d(1)
  PetscReal :: array2d(1,ncol)

  irow1d = irow
  array2d(1,:) = array1d(:)

  stop 'Error in PUMSetValuesLocal3'

#if defined(REFACTOR_323)
  stop 'PUMSetValuesLocal3'
#else
!  call MatSetValuesLocal(A,nrow,irow1d,ncol,icol,array2d, &
!                         insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValuesLocal3

! ************************************************************************** !

subroutine PUMSetValuesBlocked(A,nrow,irow,ncol,icol,matrix_block, &
                                    insert_mode,ierr)
  !
  ! Maps to MatSetValuesBlocked()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: matrix_block(:,:)
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSetValuesBlocked(A,nrow,[irow],ncol,[icol], &
                                reshape(matrix_block,(/size(matrix_block)/)), &
                                insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValuesBlocked(A,nrow,irow,ncol,icol,matrix_block, &
                           insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValuesBlocked

! ************************************************************************** !

subroutine PUMSetValuesBlockedLocal(A,nrow,irow,ncol,icol,matrix_block, &
                                    insert_mode,ierr)
  !
  ! Maps to MatSetValuesBlockedLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow
  PetscInt :: ncol
  PetscInt :: icol
  PetscReal :: matrix_block(:,:)
#if defined(REFACTOR_323)
  type(einsertmode) :: insert_mode
#else
  PetscInt :: insert_mode
#endif
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatSetValuesBlockedLocal(A,nrow,[irow],ncol,[icol], &
                                reshape(matrix_block,(/size(matrix_block)/)), &
                                insert_mode,ierr);CHKERRQ(ierr)
#else
  call MatSetValuesBlockedLocal(A,nrow,irow,ncol,icol,matrix_block, &
                                insert_mode,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMSetValuesBlockedLocal

! ************************************************************************** !

subroutine PUMGetValuesLocal(A,nrow,irow,ncol,icol,values,ierr)
  !
  ! Maps to MatSetValuesBlockedLocal()
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/25

#include "petsc/finclude/petscmat.h"
  use petscmat
  implicit none

  Mat :: A
  PetscInt :: nrow
  PetscInt :: irow(:)
  PetscInt :: ncol
  PetscInt :: icol(:)
  PetscReal :: values(:,:)
  PetscErrorCode :: ierr

#if defined(REFACTOR_323)
  call MatGetValuesLocal(A,nrow,irow,ncol,icol,[values(1,1)],ierr);CHKERRQ(ierr)
#else
  call MatGetValuesLocal(A,nrow,irow,ncol,icol,values,ierr);CHKERRQ(ierr)
#endif

end subroutine PUMGetValuesLocal

! ************************************************************************** !

subroutine PetscTestFile(filename,char,flag,ierr)

  implicit none

  character(len=*) :: filename
  character(len=*) :: char
  PetscBool :: flag
  PetscErrorCode :: ierr

  inquire(file=trim(filename), exist=flag)
  ierr = 0

end subroutine PetscTestFile

#if 0
! ************************************************************************** !

subroutine PetscTime(time,ierr)

  implicit none

  PetscLogDouble :: time
  PetscErrorCode :: ierr

  time = 0.d0
  ierr = 0

end subroutine PetscTime
#endif

end module Petsc_Utility_module
