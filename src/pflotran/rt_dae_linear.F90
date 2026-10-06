module RT_DAE_Linear_module
  !
  ! Condensed linear algebra of the GIRT_DAE model.
  !
  ! Per cell the TS state is [c; N; K] (nc, nc and nk entries), and the
  ! shifted Jacobian has the rows
  !   F_c:  -dA_dc x_c + x_N - dA_dk x_K
  !   F_N:  (G x_c)    + a x_N + dR_dk x_K
  !   F_K:  -djk_dc x_c        + D x_K,        D = a I - djk_dk
  ! where a is the shift and G (block size nc) holds the flux couplings
  ! between cells plus the cell's own reaction block dR_dc. Rows F_c, F_K
  ! and columns N, K are cell-local, so N and K are eliminated cell by
  ! cell, leaving a system in c alone,
  !   S = G + a dA_dc + E D^-1 djk_dc,    E = a dA_dk + dR_dk,
  ! the classic GIRT Jacobian with the kinetic inventories folded in, of
  ! the classic size whatever nk.
  !
  ! The operator is a shell holding G and the per-cell blocks; it applies
  ! J and J^T exactly at any shift. The preconditioner is this elimination
  ! with a solve in S, exact when that solve is, and has a transpose for
  ! adjoint solves; D is factored and S formed at preconditioner setup.
  ! Inactive cells and cells with a non-finite residual carry a diagonal
  ! Jacobian.
  !
  ! Transposed S solves use their own KSP on an explicit S^T: some PETSc
  ! preconditioners lack a transpose (e.g. hypre), and PETSc 3.24's block
  ! ILU/LU of block size 15 in natural ordering applies the forward solve.
  !
  ! Each operator owns its linearization's data, freed with the matrix.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: rt_dae_operator_type
    PetscInt :: nc = 0
    PetscInt :: nk = 0
    PetscInt :: nlocal = 0
    PetscInt :: ndof = 0
    PetscReal :: shift = 0.d0
    ! flux and reaction couplings of rows N to columns c; block size nc
    Mat :: G
    ! the condensed matrix, same pattern as G
    Mat :: S
    PetscBool :: S_created = PETSC_FALSE
    KSP :: ksp_S
    PetscBool :: ksp_S_ready = PETSC_FALSE
    ! the explicit transpose of S and its solver (adjoint solves)
    Mat :: ST
    PetscBool :: ST_created = PETSC_FALSE
    PetscBool :: ST_current = PETSC_FALSE
    KSP :: ksp_ST
    PetscBool :: ksp_ST_ready = PETSC_FALSE
    PetscBool :: shift_inblocks = PETSC_FALSE
    Vec :: xc
    Vec :: yc
    PetscReal, allocatable :: dA_dc(:,:,:)   ! (nc,nc,nlocal)
    PetscReal, allocatable :: dA_dk(:,:,:)   ! (nc,nk,nlocal)
    PetscReal, allocatable :: dR_dk(:,:,:)   ! (nc,nk,nlocal)
    PetscReal, allocatable :: djk_dc(:,:,:)  ! (nk,nc,nlocal)
    PetscReal, allocatable :: D(:,:,:)       ! (nk,nk,nlocal)
    ! preconditioner data: LU factors of D and D^-1 djk_dc
    PetscReal, allocatable :: Dlu(:,:,:)
    PetscInt, allocatable :: ipiv(:,:)
    PetscReal, allocatable :: Dinv_djk_dc(:,:,:)
    ! cells whose Jacobian is diagonal, and its entries
    PetscBool, allocatable :: diag_cell(:)
    PetscReal, allocatable :: diag(:,:)      ! (ndof,nlocal)
    PetscBool :: factored = PETSC_FALSE
    PetscBool :: singular = PETSC_FALSE
    PetscMPIInt :: comm
  end type rt_dae_operator_type

  public :: RTDAEOperatorCreate, &
            RTDAEOperatorGet, &
            RTDAEOperatorIsCondensed, &
            RTDAEOperatorBeginLinearization, &
            RTDAEOperatorSetCell, &
            RTDAEOperatorSetDiagonalCell, &
            RTDAEOperatorEndLinearization, &
            RTDAEOperatorSetupKSP

contains

! ************************************************************************** !

subroutine RTDAEOperatorCreate(G,nk,nlocal,comm,A)
  !
  ! Creates the shell operator A; it takes ownership of G (block size nc,
  ! flux stencil).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  Mat :: G
  PetscInt :: nk
  PetscInt :: nlocal
  PetscMPIInt :: comm
  Mat :: A

  type(rt_dae_operator_type), pointer :: op
  PetscInt :: nc, nloc
  PetscErrorCode :: ierr

  allocate(op)
  call MatGetBlockSize(G,nc,ierr);CHKERRQ(ierr)
  op%nc = nc
  op%nk = nk
  op%nlocal = nlocal
  op%ndof = 2*nc + nk
  op%comm = comm
  op%G = G
  call MatCreateVecs(G,op%xc,op%yc,ierr);CHKERRQ(ierr)
  allocate(op%dA_dc(nc,nc,nlocal))
  allocate(op%dA_dk(nc,nk,nlocal))
  allocate(op%dR_dk(nc,nk,nlocal))
  allocate(op%djk_dc(nk,nc,nlocal))
  allocate(op%D(nk,nk,nlocal))
  allocate(op%Dlu(nk,nk,nlocal))
  allocate(op%ipiv(nk,nlocal))
  allocate(op%Dinv_djk_dc(nk,nc,nlocal))
  allocate(op%diag_cell(nlocal))
  allocate(op%diag(op%ndof,nlocal))
  op%dA_dc = 0.d0
  op%dA_dk = 0.d0
  op%dR_dk = 0.d0
  op%djk_dc = 0.d0
  op%D = 0.d0
  op%Dlu = 0.d0
  op%ipiv = 0
  op%Dinv_djk_dc = 0.d0
  op%diag_cell = PETSC_FALSE
  op%diag = 0.d0

  call KSPCreate(comm,op%ksp_S,ierr);CHKERRQ(ierr)
  call KSPSetOptionsPrefix(op%ksp_S,'tran_schur_',ierr);CHKERRQ(ierr)
  call KSPCreate(comm,op%ksp_ST,ierr);CHKERRQ(ierr)
  call KSPSetOptionsPrefix(op%ksp_ST,'tran_schur_',ierr);CHKERRQ(ierr)

  nloc = nlocal*op%ndof
  call MatCreateShell(comm,nloc,nloc,PETSC_DETERMINE,PETSC_DETERMINE,op,A, &
                      ierr);CHKERRQ(ierr)
  call MatShellSetOperation(A,MATOP_MULT,RTDAEOperatorMult,ierr);CHKERRQ(ierr)
  call MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,RTDAEOperatorMultTranspose, &
                            ierr);CHKERRQ(ierr)
  call MatShellSetOperation(A,MATOP_DESTROY,RTDAEOperatorDestroy, &
                            ierr);CHKERRQ(ierr)

end subroutine RTDAEOperatorCreate

! ************************************************************************** !

subroutine RTDAEOperatorGet(A,op)
  !
  ! The operator data of shell A.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  implicit none

  Mat :: A
  type(rt_dae_operator_type), pointer :: op

  type(c_ptr) :: ctx_c
  PetscErrorCode :: ierr

  call MatShellGetContext(A,ctx_c,ierr);CHKERRQ(ierr)
  call c_f_pointer(ctx_c,op)

end subroutine RTDAEOperatorGet

! ************************************************************************** !

function RTDAEOperatorIsCondensed(A)
  !
  ! Whether A is the condensed shell rather than an assembled matrix.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  Mat :: A
  PetscBool :: RTDAEOperatorIsCondensed

  MatType :: mat_type
  PetscErrorCode :: ierr

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  RTDAEOperatorIsCondensed = (mat_type == MATSHELL)

end function RTDAEOperatorIsCondensed

! ************************************************************************** !

subroutine RTDAEOperatorBeginLinearization(op,shift)
  !
  ! The caller refills G and every cell; the factors are stale.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(rt_dae_operator_type) :: op
  PetscReal :: shift

  PetscErrorCode :: ierr

  call MatZeroEntries(op%G,ierr);CHKERRQ(ierr)
  op%shift = shift
  op%factored = PETSC_FALSE

end subroutine RTDAEOperatorBeginLinearization

! ************************************************************************** !

subroutine RTDAEOperatorSetCell(op,icell,dA_dc,dA_dk,dR_dk,djk_dc,djk_dk)
  !
  ! The local blocks of cell icell (see the module header); dR_dc and the
  ! fluxes go into G.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(rt_dae_operator_type) :: op
  PetscInt :: icell
  PetscReal :: dA_dc(:,:), dA_dk(:,:), dR_dk(:,:), djk_dc(:,:), djk_dk(:,:)

  PetscInt :: i

  op%diag_cell(icell) = PETSC_FALSE
  op%dA_dc(:,:,icell) = dA_dc
  if (op%nk > 0) then
    op%dA_dk(:,:,icell) = dA_dk
    op%dR_dk(:,:,icell) = dR_dk
    op%djk_dc(:,:,icell) = djk_dc
    op%D(:,:,icell) = -djk_dk
    do i = 1, op%nk
      op%D(i,i,icell) = op%D(i,i,icell) + op%shift
    enddo
  endif

end subroutine RTDAEOperatorSetCell

! ************************************************************************** !

subroutine RTDAEOperatorSetDiagonalCell(op,icell,diag)
  !
  ! Cell icell has the diagonal Jacobian diag (all entries nonzero); the
  ! caller adds a zero diagonal block to G.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(rt_dae_operator_type) :: op
  PetscInt :: icell
  PetscReal :: diag(:)

  op%diag_cell(icell) = PETSC_TRUE
  op%diag(:,icell) = diag

end subroutine RTDAEOperatorSetDiagonalCell

! ************************************************************************** !

subroutine RTDAEOperatorEndLinearization(op)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  type(rt_dae_operator_type) :: op

  PetscErrorCode :: ierr

  call MatAssemblyBegin(op%G,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(op%G,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine RTDAEOperatorEndLinearization

! ************************************************************************** !

subroutine RTDAEOperatorSetupKSP(ksp,A,pc_type,shift_inblocks)
  !
  ! Installs the elimination preconditioner on ksp (S solves: see
  ! RTDAEConfigureInnerKSP). An iterative S solve makes it nonlinear and
  ! needs a flexible outer KSP.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  KSP :: ksp
  Mat :: A
  PCType :: pc_type
  PetscBool :: shift_inblocks

  type(rt_dae_operator_type), pointer :: op
  PC :: pc
  PetscErrorCode :: ierr

  call RTDAEOperatorGet(A,op)
  call RTDAEConfigureInnerKSP(op%ksp_S,pc_type,op%comm)
  call RTDAEConfigureInnerKSP(op%ksp_ST,pc_type,op%comm)
  op%shift_inblocks = shift_inblocks
  op%ksp_S_ready = PETSC_FALSE
  op%ksp_ST_ready = PETSC_FALSE

  call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
  call PCSetType(pc,PCSHELL,ierr);CHKERRQ(ierr)
  call PCShellSetSetUp(pc,RTDAEPCSetUp,ierr);CHKERRQ(ierr)
  call PCShellSetApply(pc,RTDAEPCApply,ierr);CHKERRQ(ierr)
  call PCShellSetApplyTranspose(pc,RTDAEPCApplyTranspose,ierr);CHKERRQ(ierr)
  call PCShellSetName(pc,'GIRT_DAE elimination of N and K',ierr)
  CHKERRQ(ierr)

end subroutine RTDAEOperatorSetupKSP

! ************************************************************************** !

subroutine RTDAEConfigureInnerKSP(ksp,pc_type,comm)
  !
  ! One application of pc_type, or with a blank pc_type PETSc's default
  ! (ILU on one process, block Jacobi on several); -tran_schur_ksp_type
  ! and -tran_schur_pc_type override.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  KSP :: ksp
  PCType :: pc_type
  PetscMPIInt :: comm

  PC :: pc
  PetscMPIInt :: comm_size
  PetscErrorCode :: ierr

  call KSPSetType(ksp,KSPPREONLY,ierr);CHKERRQ(ierr)
  call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
  if (len_trim(pc_type) > 0) then
    call PCSetType(pc,pc_type,ierr);CHKERRQ(ierr)
  else
    call MPI_Comm_size(comm,comm_size,ierr);CHKERRQ(ierr)
    if (comm_size == 1) then
      call PCSetType(pc,PCILU,ierr);CHKERRQ(ierr)
    else
      call PCSetType(pc,PCBJACOBI,ierr);CHKERRQ(ierr)
    endif
  endif
  call KSPSetFromOptions(ksp,ierr);CHKERRQ(ierr)

end subroutine RTDAEConfigureInnerKSP

! ************************************************************************** !

subroutine RTDAEInnerShiftInBlocks(ksp)
  !
  ! PETSc does not shift zero pivots by default; shift within blocks as
  ! SolverSetupPCShiftAndPivoting does. Must precede the setup that factors
  ! ILU/LU; block Jacobi needs a setup to create its blocks, which factor
  ! at the first solve. Needs the operators set.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  KSP :: ksp

  KSP, pointer :: sub_ksps(:)
  PC :: pc, pc_sub
  PCType :: pc_type
  PetscInt :: nsub, first_sub, i
  PetscErrorCode :: ierr

  call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
  call PCFactorSetShiftType(pc,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
  call PCGetType(pc,pc_type,ierr);CHKERRQ(ierr)
  if (pc_type == PCBJACOBI) then
    call KSPSetUp(ksp,ierr);CHKERRQ(ierr)
    nullify(sub_ksps)
    call PCBJacobiGetSubKSP(pc,nsub,first_sub,sub_ksps,ierr);CHKERRQ(ierr)
    do i = 1, nsub
      call KSPGetPC(sub_ksps(i),pc_sub,ierr);CHKERRQ(ierr)
      call PCFactorSetShiftType(pc_sub,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
    enddo
    nullify(sub_ksps)
  endif

end subroutine RTDAEInnerShiftInBlocks

! ************************************************************************** !

subroutine RTDAEPCOperator(pc,op)
  !
  ! The operator data of the preconditioner's Amat.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  type(rt_dae_operator_type), pointer :: op

  Mat :: A, B
  PetscErrorCode :: ierr

  call PCGetOperators(pc,A,B,ierr);CHKERRQ(ierr)
  call RTDAEOperatorGet(A,op)

end subroutine RTDAEPCOperator

! ************************************************************************** !

subroutine RTDAEPCSetUp(pc,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  PC :: pc
  PetscErrorCode :: ierr

  type(rt_dae_operator_type), pointer :: op

  ! PETSc never clears a PCSHELL's failure reason; reset it so a past
  ! failure does not fail later solves
  call PCSetFailedReason(pc,PC_NOERROR,ierr);CHKERRQ(ierr)
  call RTDAEPCOperator(pc,op)
  call RTDAEOperatorFactor(op)
  ierr = 0

end subroutine RTDAEPCSetUp

! ************************************************************************** !

subroutine RTDAEOperatorFactor(op)
  !
  ! Factors D per cell and forms S = G + a dA_dc + E D^-1 djk_dc.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(rt_dae_operator_type) :: op

  PetscReal :: corr(op%nc,op%nc), E(op%nc,op%nk)
  PetscInt :: icell, info, nc, nk, i, rstart, rend, row
  PetscBool :: singular
  PetscErrorCode :: ierr

  if (op%factored) return
  nc = op%nc
  nk = op%nk
  singular = PETSC_FALSE
  if (op%S_created) then
    call MatCopy(op%G,op%S,SAME_NONZERO_PATTERN,ierr);CHKERRQ(ierr)
  else
    call MatDuplicate(op%G,MAT_COPY_VALUES,op%S,ierr);CHKERRQ(ierr)
    call MatSetOption(op%S,MAT_ROW_ORIENTED,PETSC_FALSE,ierr);CHKERRQ(ierr)
    call MatSetOption(op%S,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE, &
                      ierr);CHKERRQ(ierr)
    op%S_created = PETSC_TRUE
  endif
  call MatGetOwnershipRange(op%S,rstart,rend,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    row = rstart/nc + icell - 1
    corr = 0.d0
    if (op%diag_cell(icell)) then
      ! an identity row: the right-hand side carries the solution
      do i = 1, nc
        corr(i,i) = 1.d0
      enddo
    else
      corr = op%shift*op%dA_dc(:,:,icell)
      if (nk > 0) then
        op%Dlu(:,:,icell) = op%D(:,:,icell)
        call DGETRF(nk,nk,op%Dlu(:,:,icell),nk,op%ipiv(:,icell),info)
        if (info /= 0) then
          singular = PETSC_TRUE
        else
          op%Dinv_djk_dc(:,:,icell) = op%djk_dc(:,:,icell)
          call DGETRS('N',nk,nc,op%Dlu(:,:,icell),nk,op%ipiv(:,icell), &
                      op%Dinv_djk_dc(:,:,icell),nk,info)
          E = op%shift*op%dA_dk(:,:,icell) + op%dR_dk(:,:,icell)
          corr = corr + matmul(E,op%Dinv_djk_dc(:,:,icell))
        endif
      endif
    endif
    call MatSetValuesBlocked(op%S,1,[row],1,[row],reshape(corr,[nc*nc]), &
                             ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo
  call MatAssemblyBegin(op%S,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(op%S,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,singular,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,op%comm,ierr);CHKERRQ(ierr)
  op%singular = singular
  call KSPSetOperators(op%ksp_S,op%S,op%S,ierr);CHKERRQ(ierr)
  if (.not.op%ksp_S_ready) then
    if (op%shift_inblocks) call RTDAEInnerShiftInBlocks(op%ksp_S)
    op%ksp_S_ready = PETSC_TRUE
  endif
  op%ST_current = PETSC_FALSE
  op%factored = PETSC_TRUE

end subroutine RTDAEOperatorFactor

! ************************************************************************** !

subroutine RTDAEOperatorMult(A,x,y,ierr)
  !
  ! y = J x
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  Mat :: A
  Vec :: x, y
  PetscErrorCode :: ierr

  type(rt_dae_operator_type), pointer :: op
  PetscReal, pointer :: x_p(:), y_p(:), xc_p(:), yc_p(:)
  PetscInt :: icell, o, oc, nc, nk

  call RTDAEOperatorGet(A,op)
  nc = op%nc
  nk = op%nk
  call VecGetArrayRead(x,x_p,ierr);CHKERRQ(ierr)
  call VecGetArray(op%xc,xc_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    xc_p((icell-1)*nc+1:icell*nc) = x_p(o+1:o+nc)
  enddo
  call VecRestoreArray(op%xc,xc_p,ierr);CHKERRQ(ierr)
  call MatMult(op%G,op%xc,op%yc,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(y,y_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    oc = (icell-1)*nc
    if (op%diag_cell(icell)) then
      y_p(o+1:o+op%ndof) = op%diag(:,icell)*x_p(o+1:o+op%ndof)
      cycle
    endif
    y_p(o+1:o+nc) = -matmul(op%dA_dc(:,:,icell),x_p(o+1:o+nc)) + &
                    x_p(o+nc+1:o+2*nc)
    y_p(o+nc+1:o+2*nc) = yc_p(oc+1:oc+nc) + op%shift*x_p(o+nc+1:o+2*nc)
    if (nk > 0) then
      y_p(o+1:o+nc) = y_p(o+1:o+nc) - &
        matmul(op%dA_dk(:,:,icell),x_p(o+2*nc+1:o+op%ndof))
      y_p(o+nc+1:o+2*nc) = y_p(o+nc+1:o+2*nc) + &
        matmul(op%dR_dk(:,:,icell),x_p(o+2*nc+1:o+op%ndof))
      y_p(o+2*nc+1:o+op%ndof) = &
        -matmul(op%djk_dc(:,:,icell),x_p(o+1:o+nc)) + &
        matmul(op%D(:,:,icell),x_p(o+2*nc+1:o+op%ndof))
    endif
  enddo
  call VecRestoreArray(y,y_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(x,x_p,ierr);CHKERRQ(ierr)
  ierr = 0

end subroutine RTDAEOperatorMult

! ************************************************************************** !

subroutine RTDAEOperatorMultTranspose(A,x,y,ierr)
  !
  ! y = J^T x
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  Mat :: A
  Vec :: x, y
  PetscErrorCode :: ierr

  type(rt_dae_operator_type), pointer :: op
  PetscReal, pointer :: x_p(:), y_p(:), xc_p(:), yc_p(:)
  PetscInt :: icell, o, oc, nc, nk

  call RTDAEOperatorGet(A,op)
  nc = op%nc
  nk = op%nk
  call VecGetArrayRead(x,x_p,ierr);CHKERRQ(ierr)
  call VecGetArray(op%xc,xc_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    xc_p((icell-1)*nc+1:icell*nc) = x_p(o+nc+1:o+2*nc)
  enddo
  call VecRestoreArray(op%xc,xc_p,ierr);CHKERRQ(ierr)
  call MatMultTranspose(op%G,op%xc,op%yc,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(y,y_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    oc = (icell-1)*nc
    if (op%diag_cell(icell)) then
      y_p(o+1:o+op%ndof) = op%diag(:,icell)*x_p(o+1:o+op%ndof)
      cycle
    endif
    y_p(o+1:o+nc) = yc_p(oc+1:oc+nc) - &
      matmul(transpose(op%dA_dc(:,:,icell)),x_p(o+1:o+nc))
    y_p(o+nc+1:o+2*nc) = x_p(o+1:o+nc) + op%shift*x_p(o+nc+1:o+2*nc)
    if (nk > 0) then
      y_p(o+1:o+nc) = y_p(o+1:o+nc) - &
        matmul(transpose(op%djk_dc(:,:,icell)),x_p(o+2*nc+1:o+op%ndof))
      y_p(o+2*nc+1:o+op%ndof) = &
        -matmul(transpose(op%dA_dk(:,:,icell)),x_p(o+1:o+nc)) + &
        matmul(transpose(op%dR_dk(:,:,icell)),x_p(o+nc+1:o+2*nc)) + &
        matmul(transpose(op%D(:,:,icell)),x_p(o+2*nc+1:o+op%ndof))
    endif
  enddo
  call VecRestoreArray(y,y_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(x,x_p,ierr);CHKERRQ(ierr)
  ierr = 0

end subroutine RTDAEOperatorMultTranspose

! ************************************************************************** !

subroutine RTDAEPCApply(pc,r,z,ierr)
  !
  ! z = J^-1 r by elimination of N and K:
  !   S z_c = r_N - a r_c - E D^-1 r_K
  !   z_K = D^-1 r_K + D^-1 djk_dc z_c
  !   z_N = r_c + dA_dc z_c + dA_dk z_K
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  Vec :: r, z
  PetscErrorCode :: ierr

  type(rt_dae_operator_type), pointer :: op
  PetscReal, pointer :: r_p(:), z_p(:), xc_p(:), yc_p(:)
  PetscReal, allocatable :: tK(:), E(:,:)
  PetscInt :: icell, o, oc, nc, nk, info

  call PCSetFailedReason(pc,PC_NOERROR,ierr);CHKERRQ(ierr)
  call RTDAEPCOperator(pc,op)
  call RTDAEOperatorFactor(op)
  if (RTDAEPCFailed(pc,op,z)) then
    ierr = 0
    return
  endif
  nc = op%nc
  nk = op%nk
  allocate(tK(nk),E(nc,nk))

  call VecGetArrayRead(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArray(op%xc,xc_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    oc = (icell-1)*nc
    if (op%diag_cell(icell)) then
      xc_p(oc+1:oc+nc) = r_p(o+1:o+nc)/op%diag(1:nc,icell)
      cycle
    endif
    xc_p(oc+1:oc+nc) = r_p(o+nc+1:o+2*nc) - op%shift*r_p(o+1:o+nc)
    if (nk > 0) then
      tK = r_p(o+2*nc+1:o+op%ndof)
      call DGETRS('N',nk,1,op%Dlu(:,:,icell),nk,op%ipiv(:,icell),tK,nk,info)
      E = op%shift*op%dA_dk(:,:,icell) + op%dR_dk(:,:,icell)
      xc_p(oc+1:oc+nc) = xc_p(oc+1:oc+nc) - matmul(E,tK)
    endif
  enddo
  call VecRestoreArray(op%xc,xc_p,ierr);CHKERRQ(ierr)

  call RTDAESolveS(pc,op,PETSC_FALSE)

  call VecGetArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(z,z_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    oc = (icell-1)*nc
    if (op%diag_cell(icell)) then
      z_p(o+1:o+op%ndof) = r_p(o+1:o+op%ndof)/op%diag(:,icell)
      cycle
    endif
    z_p(o+1:o+nc) = yc_p(oc+1:oc+nc)
    z_p(o+nc+1:o+2*nc) = r_p(o+1:o+nc) + &
      matmul(op%dA_dc(:,:,icell),yc_p(oc+1:oc+nc))
    if (nk > 0) then
      tK = r_p(o+2*nc+1:o+op%ndof)
      call DGETRS('N',nk,1,op%Dlu(:,:,icell),nk,op%ipiv(:,icell),tK,nk,info)
      tK = tK + matmul(op%Dinv_djk_dc(:,:,icell),yc_p(oc+1:oc+nc))
      z_p(o+2*nc+1:o+op%ndof) = tK
      z_p(o+nc+1:o+2*nc) = z_p(o+nc+1:o+2*nc) + &
        matmul(op%dA_dk(:,:,icell),tK)
    endif
  enddo
  call VecRestoreArray(z,z_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(r,r_p,ierr);CHKERRQ(ierr)
  deallocate(tK,E)
  call RTDAEPCCheckFinite(pc,z)
  ierr = 0

end subroutine RTDAEPCApply

! ************************************************************************** !

subroutine RTDAEPCApplyTranspose(pc,r,z,ierr)
  !
  ! z = J^-T r (adjoint solves):
  !   S^T z_N = r_c + dA_dc^T r_N + (D^-1 djk_dc)^T (r_K + dA_dk^T r_N)
  !   z_c = r_N - a z_N
  !   z_K = D^-T (r_K + dA_dk^T z_c - dR_dk^T z_N)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  Vec :: r, z
  PetscErrorCode :: ierr

  type(rt_dae_operator_type), pointer :: op
  PetscReal, pointer :: r_p(:), z_p(:), xc_p(:), yc_p(:)
  PetscReal, allocatable :: tK(:), zc(:)
  PetscInt :: icell, o, oc, nc, nk, info

  call PCSetFailedReason(pc,PC_NOERROR,ierr);CHKERRQ(ierr)
  call RTDAEPCOperator(pc,op)
  call RTDAEOperatorFactor(op)
  if (RTDAEPCFailed(pc,op,z)) then
    ierr = 0
    return
  endif
  nc = op%nc
  nk = op%nk
  allocate(tK(nk),zc(nc))

  call VecGetArrayRead(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArray(op%xc,xc_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    oc = (icell-1)*nc
    if (op%diag_cell(icell)) then
      xc_p(oc+1:oc+nc) = r_p(o+1:o+nc)/op%diag(1:nc,icell)
      cycle
    endif
    xc_p(oc+1:oc+nc) = r_p(o+1:o+nc) + &
      matmul(transpose(op%dA_dc(:,:,icell)),r_p(o+nc+1:o+2*nc))
    if (nk > 0) then
      tK = r_p(o+2*nc+1:o+op%ndof) + &
        matmul(transpose(op%dA_dk(:,:,icell)),r_p(o+nc+1:o+2*nc))
      xc_p(oc+1:oc+nc) = xc_p(oc+1:oc+nc) + &
        matmul(transpose(op%Dinv_djk_dc(:,:,icell)),tK)
    endif
  enddo
  call VecRestoreArray(op%xc,xc_p,ierr);CHKERRQ(ierr)

  call RTDAESolveS(pc,op,PETSC_TRUE)

  call VecGetArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(z,z_p,ierr);CHKERRQ(ierr)
  do icell = 1, op%nlocal
    o = (icell-1)*op%ndof
    oc = (icell-1)*nc
    if (op%diag_cell(icell)) then
      z_p(o+1:o+op%ndof) = r_p(o+1:o+op%ndof)/op%diag(:,icell)
      cycle
    endif
    ! yc holds z_N
    zc = r_p(o+nc+1:o+2*nc) - op%shift*yc_p(oc+1:oc+nc)
    z_p(o+1:o+nc) = zc
    z_p(o+nc+1:o+2*nc) = yc_p(oc+1:oc+nc)
    if (nk > 0) then
      tK = r_p(o+2*nc+1:o+op%ndof) + &
        matmul(transpose(op%dA_dk(:,:,icell)),zc) - &
        matmul(transpose(op%dR_dk(:,:,icell)),yc_p(oc+1:oc+nc))
      call DGETRS('T',nk,1,op%Dlu(:,:,icell),nk,op%ipiv(:,icell),tK,nk,info)
      z_p(o+2*nc+1:o+op%ndof) = tK
    endif
  enddo
  call VecRestoreArray(z,z_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(op%yc,yc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(r,r_p,ierr);CHKERRQ(ierr)
  deallocate(tK,zc)
  call RTDAEPCCheckFinite(pc,z)
  ierr = 0

end subroutine RTDAEPCApplyTranspose

! ************************************************************************** !

function RTDAEPCFailed(pc,op,z)
  !
  ! A singular D admits no elimination: fail the preconditioner so the
  ! step is retried smaller (larger shift). KSP reads the failed reason
  ! only when a norm is not finite, so z is flagged with infinity; every
  ! process flags z (singular is reduced in RTDAEOperatorFactor).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  type(rt_dae_operator_type) :: op
  Vec :: z
  PetscBool :: RTDAEPCFailed

  PetscErrorCode :: ierr

  RTDAEPCFailed = op%singular
  if (.not.op%singular) return
  call VecFlag(z,ONE_INTEGER,ierr);CHKERRQ(ierr)
  call PCSetFailedReason(pc,PC_SUBPC_ERROR,ierr);CHKERRQ(ierr)

end function RTDAEPCFailed

! ************************************************************************** !

subroutine RTDAEPCCheckFinite(pc,z)
  !
  ! An elimination that overflows fails the preconditioner as a singular D
  ! does (RTDAEPCFailed): with LOG_FORMULATION, a vanishing species whose
  ! inventory no concentration matches has a step of order (N - A)/A.
  !
  ! Author: Hansol Suh
  ! Date: 10/06/26
  !
  implicit none

  PC :: pc
  Vec :: z

  PetscReal, pointer :: z_p(:)
  PetscBool :: finite
  PetscErrorCode :: ierr

  call VecGetArrayRead(z,z_p,ierr);CHKERRQ(ierr)
  finite = all(abs(z_p) <= huge(1.d0))
  call VecRestoreArrayRead(z,z_p,ierr);CHKERRQ(ierr)
  if (finite) return
  call VecFlag(z,ONE_INTEGER,ierr);CHKERRQ(ierr)
  call PCSetFailedReason(pc,PC_SUBPC_ERROR,ierr);CHKERRQ(ierr)

end subroutine RTDAEPCCheckFinite

! ************************************************************************** !

subroutine RTDAESolveS(pc,op,transposed)
  !
  ! yc = S^-1 xc (S^-T xc when transposed).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  type(rt_dae_operator_type) :: op
  PetscBool :: transposed

  KSPConvergedReason :: reason
  PetscErrorCode :: ierr

  if (transposed) then
    if (.not.op%ST_current) then
      if (op%ST_created) then
        call MatTranspose(op%S,MAT_REUSE_MATRIX,op%ST,ierr);CHKERRQ(ierr)
      else
        call MatTranspose(op%S,MAT_INITIAL_MATRIX,op%ST,ierr);CHKERRQ(ierr)
        op%ST_created = PETSC_TRUE
      endif
      call KSPSetOperators(op%ksp_ST,op%ST,op%ST,ierr);CHKERRQ(ierr)
      if (.not.op%ksp_ST_ready) then
        if (op%shift_inblocks) call RTDAEInnerShiftInBlocks(op%ksp_ST)
        op%ksp_ST_ready = PETSC_TRUE
      endif
      op%ST_current = PETSC_TRUE
    endif
    call KSPSolve(op%ksp_ST,op%xc,op%yc,ierr);CHKERRQ(ierr)
    call KSPGetConvergedReason(op%ksp_ST,reason,ierr);CHKERRQ(ierr)
  else
    call KSPSolve(op%ksp_S,op%xc,op%yc,ierr);CHKERRQ(ierr)
    call KSPGetConvergedReason(op%ksp_S,reason,ierr);CHKERRQ(ierr)
  endif
  if (reason%v < 0) then
    call PCSetFailedReason(pc,PC_SUBPC_ERROR,ierr);CHKERRQ(ierr)
  endif

end subroutine RTDAESolveS

! ************************************************************************** !

subroutine RTDAEOperatorDestroy(A,ierr)
  !
  ! Frees the operator data with the matrix.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  Mat :: A
  PetscErrorCode :: ierr

  type(rt_dae_operator_type), pointer :: op

  call RTDAEOperatorGet(A,op)
  call MatDestroy(op%G,ierr);CHKERRQ(ierr)
  if (op%S_created) then
    call MatDestroy(op%S,ierr);CHKERRQ(ierr)
  endif
  call VecDestroy(op%xc,ierr);CHKERRQ(ierr)
  call VecDestroy(op%yc,ierr);CHKERRQ(ierr)
  call KSPDestroy(op%ksp_S,ierr);CHKERRQ(ierr)
  if (op%ST_created) then
    call MatDestroy(op%ST,ierr);CHKERRQ(ierr)
  endif
  call KSPDestroy(op%ksp_ST,ierr);CHKERRQ(ierr)
  deallocate(op)
  ierr = 0

end subroutine RTDAEOperatorDestroy

end module RT_DAE_Linear_module
