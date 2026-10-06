module TS_Model_Checks_module
  !
  ! Development checks of an equation model (TS_Model_class) at a state,
  ! through the model interface only; they leave the run unchanged:
  !   TSModelContractTest  the model contract, -ts_model_contract_test
  !   TSModelColumnCheck   J(1/dt) entry by entry against central
  !                        differences, -ts_model_column_check
  !   TSModelAdjointTest   a discrete adjoint Taylor test on backward
  !                        Euler, -ts_model_adjoint_test
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use PFLOTRAN_Constants_module
  use TS_Model_class
  use Option_module

  implicit none

  private

  ! the run the checks mirror
  type, public :: ts_model_checks_type
    class(ts_model_type), pointer :: model => null()
    type(option_type), pointer :: option => null()
    ! state the checks start from; restored where a check moves it
    Vec :: Y
    ! the run's absolute local error tolerances (weights and FD steps)
    Vec :: vatol
    ! the run's Jacobian, duplicated unless the model creates its own
    Mat :: jacobian
    ! the model supplies the Jacobian operator and its preconditioner
    PetscBool :: model_jacobian = PETSC_FALSE
    ! as in TSModelSetupLinearSolver
    PCType :: pc_type = ''
    PetscBool :: shift_inblocks = PETSC_FALSE
  end type ts_model_checks_type

  ! what PETSc hands the callbacks of the checks' own TS
  type :: ts_model_checks_ctx_type
    class(ts_model_type), pointer :: model => null()
  end type ts_model_checks_ctx_type

  public :: TSModelContractTest, &
            TSModelColumnCheck, &
            TSModelAdjointTest

contains

! ************************************************************************** !

subroutine TSModelChecksIFunction(ts,time,Y,Ydot,F,ctx,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: Y, Ydot, F
  type(ts_model_checks_ctx_type) :: ctx
  PetscErrorCode :: ierr

  SNES :: snes
  PetscInt :: status

  call ctx%model%Evaluate(time,Y,Ydot,F,status)
  if (status /= TS_MODEL_OK) then
    call TSGetSNES(ts,snes,ierr);CHKERRQ(ierr)
    call SNESSetFunctionDomainError(snes,ierr);CHKERRQ(ierr)
  endif
  ierr = 0

end subroutine TSModelChecksIFunction

! ************************************************************************** !

subroutine TSModelChecksIJacobian(ts,time,Y,Ydot,shift,A,B,ctx,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: Y, Ydot
  PetscReal :: shift
  Mat :: A, B
  type(ts_model_checks_ctx_type) :: ctx
  PetscErrorCode :: ierr

  SNES :: snes
  PetscInt :: status

  call ctx%model%Linearize(time,Y,Ydot,shift,A,B,status)
  if (status /= TS_MODEL_OK) then
    call TSGetSNES(ts,snes,ierr);CHKERRQ(ierr)
    call SNESSetJacobianDomainError(snes,ierr);CHKERRQ(ierr)
  endif
  ierr = 0

end subroutine TSModelChecksIJacobian

! ************************************************************************** !

subroutine TSModelChecksNewJacobian(checks,J)
  !
  ! A new Jacobian operator of the run's kind, independent of the run's.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_model_checks_type) :: checks
  Mat :: J

  PetscBool :: created
  PetscErrorCode :: ierr

  call checks%model%CreateJacobian(J,created)
  if (created) return
  call MatDuplicate(checks%jacobian,MAT_DO_NOT_COPY_VALUES,J, &
                    ierr);CHKERRQ(ierr)
  call MatSetOption(J,MAT_ROW_ORIENTED,PETSC_FALSE,ierr);CHKERRQ(ierr)

end subroutine TSModelChecksNewJacobian

! ************************************************************************** !

subroutine TSModelChecksForward(checks,ctx,time,h,nsteps,Y,save_trajectory, &
                                ts_out)
  !
  ! nsteps fixed backward Euler steps of size h from time, on a TS of its
  ! own; Y becomes the final state. With save_trajectory the TS is returned
  ! for an adjoint solve (ctx must outlive it), otherwise destroyed.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_model_checks_type) :: checks
  type(ts_model_checks_ctx_type) :: ctx
  PetscReal :: time
  PetscReal :: h
  PetscInt :: nsteps
  Vec :: Y
  PetscBool :: save_trajectory
  TS :: ts_out

  TS :: ts
  TSAdapt :: adapt
  SNES :: snes
  SNESLineSearch :: linesearch
  KSP :: ksp
  Vec :: r
  Mat :: J
  TSConvergedReason :: reason
  PetscErrorCode :: ierr

  call VecDuplicate(Y,r,ierr);CHKERRQ(ierr)
  call TSModelChecksNewJacobian(checks,J)
  call TSCreate(checks%option%mycomm,ts,ierr);CHKERRQ(ierr)
  call TSSetOptionsPrefix(ts,'adjoint_test_',ierr);CHKERRQ(ierr)
  call TSSetType(ts,TSBEULER,ierr);CHKERRQ(ierr)
  call TSSetEquationType(ts,TS_EQ_DAE_IMPLICIT_INDEX1,ierr);CHKERRQ(ierr)
  call TSSetIFunction(ts,r,TSModelChecksIFunction,ctx,ierr);CHKERRQ(ierr)
  call TSSetIJacobian(ts,J,J,TSModelChecksIJacobian,ctx,ierr);CHKERRQ(ierr)
  call TSSetTime(ts,time,ierr);CHKERRQ(ierr)
  call TSSetTimeStep(ts,h,ierr);CHKERRQ(ierr)
  call TSSetMaxSteps(ts,nsteps,ierr);CHKERRQ(ierr)
  call TSSetMaxTime(ts,time+nsteps*h*(1.d0+1.d-12),ierr);CHKERRQ(ierr)
  call TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER,ierr);CHKERRQ(ierr)
  call TSGetAdapt(ts,adapt,ierr);CHKERRQ(ierr)
  call TSAdaptSetType(adapt,TSADAPTNONE,ierr);CHKERRQ(ierr)
  call TSSetErrorIfStepFails(ts,PETSC_FALSE,ierr);CHKERRQ(ierr)
  call TSGetSNES(ts,snes,ierr);CHKERRQ(ierr)
  ! solve each stage to roundoff: the Taylor test needs exact forward maps
  call SNESSetTolerances(snes,1.d-40,1.d-14,1.d-16,50,10000, &
                         ierr);CHKERRQ(ierr)
  ! no divergence test (PETSC_UNLIMITED), see PMCTSDriverSetupSolvers
  call SNESSetDivergenceTolerance(snes,-3.d0,ierr);CHKERRQ(ierr)
  call SNESGetLineSearch(snes,linesearch,ierr);CHKERRQ(ierr)
  call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC,ierr);CHKERRQ(ierr)
  if (checks%model_jacobian) then
    call SNESGetKSP(snes,ksp,ierr);CHKERRQ(ierr)
    call checks%model%SetupLinearSolver(ksp,J,checks%pc_type, &
                                        checks%shift_inblocks)
  endif
  if (save_trajectory) then
    call TSSetSaveTrajectory(ts,ierr);CHKERRQ(ierr)
  endif
  call TSSetFromOptions(ts,ierr);CHKERRQ(ierr)
  call TSSolve(ts,Y,ierr);CHKERRQ(ierr)
  call TSGetConvergedReason(ts,reason,ierr);CHKERRQ(ierr)
  if (reason%v < 0) then
    call PrintMsg(checks%option,' TS model adjoint test: forward solve &
                  &failed')
  endif
  call VecDestroy(r,ierr);CHKERRQ(ierr)
  call MatDestroy(J,ierr);CHKERRQ(ierr)
  if (save_trajectory) then
    ts_out = ts
  else
    call TSDestroy(ts,ierr);CHKERRQ(ierr)
  endif

end subroutine TSModelChecksForward

! ************************************************************************** !

subroutine TSModelAdjointTest(checks,time,dt,t_end)
  !
  ! Discrete adjoint check on backward Euler: J(Y0) = w . Y(T) over a few
  ! fixed steps inside the first forcing segment, dJ/dY0 by
  ! TSAdjointSolve, and the Taylor test
  !   |J(Y0 + e v) - J(Y0) - e dJ/dY0 . v| = O(e^2).
  ! Exercises the zero-shift Jacobian, transposed solves and the absence
  ! of hidden state. Prints the remainders and their ratios.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_model_checks_type) :: checks
  PetscReal :: time
  PetscReal :: dt
  PetscReal :: t_end

  type(ts_model_checks_ctx_type) :: ctx
  TS :: ts
  Vec :: Y0, Y, w, lambda(1), v
  PetscRandom :: rand
  PetscReal :: h, J0, Jp, dJv, eps, rem(6)
  PetscReal, pointer :: v_p(:), y_p(:), w_p(:), a_p(:)
  PetscInt :: nsteps, i, n, ieps
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  ctx%model => checks%model
  nsteps = 4
  h = min(dt,(t_end-time)/dble(nsteps))
  call VecDuplicate(checks%Y,Y0,ierr);CHKERRQ(ierr)
  call VecDuplicate(checks%Y,Y,ierr);CHKERRQ(ierr)
  call VecDuplicate(checks%Y,w,ierr);CHKERRQ(ierr)
  call VecDuplicate(checks%Y,v,ierr);CHKERRQ(ierr)
  call VecDuplicate(checks%Y,lambda(1),ierr);CHKERRQ(ierr)
  call VecCopy(checks%Y,Y0,ierr);CHKERRQ(ierr)

  call PetscRandomCreate(checks%option%mycomm,rand,ierr);CHKERRQ(ierr)
  call PetscRandomSetInterval(rand,0.5d0,1.5d0,ierr);CHKERRQ(ierr)
  call PetscRandomSetSeed(rand,7654321_8,ierr);CHKERRQ(ierr)
  call PetscRandomSeed(rand,ierr);CHKERRQ(ierr)
  call VecSetRandom(w,rand,ierr);CHKERRQ(ierr)
  call PetscRandomSetInterval(rand,-1.d0,1.d0,ierr);CHKERRQ(ierr)
  call VecSetRandom(v,rand,ierr);CHKERRQ(ierr)
  call PetscRandomDestroy(rand,ierr);CHKERRQ(ierr)
  ! w weighs differential entries (atol > 0) by 1/max(|Y0|,atol); v is
  ! relative to |Y0|
  call VecGetLocalSize(v,n,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Y0,y_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(checks%vatol,a_p,ierr);CHKERRQ(ierr)
  call VecGetArray(v,v_p,ierr);CHKERRQ(ierr)
  call VecGetArray(w,w_p,ierr);CHKERRQ(ierr)
  do i = 1, n
    v_p(i) = v_p(i)*max(abs(y_p(i)),1.d-30)
    if (a_p(i) > 0.d0) then
      w_p(i) = w_p(i)/max(abs(y_p(i)),a_p(i))
    else
      w_p(i) = 0.d0
    endif
  enddo
  call VecRestoreArray(w,w_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(v,v_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(checks%vatol,a_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y0,y_p,ierr);CHKERRQ(ierr)

  ! forward with trajectory, then the adjoint
  call VecCopy(Y0,Y,ierr);CHKERRQ(ierr)
  call TSModelChecksForward(checks,ctx,time,h,nsteps,Y,PETSC_TRUE,ts)
  call VecDot(w,Y,J0,ierr);CHKERRQ(ierr)
  call VecCopy(w,lambda(1),ierr);CHKERRQ(ierr)
  call TSSetCostGradients(ts,1,lambda,PETSC_NULL_VEC_ARRAY,ierr);CHKERRQ(ierr)
  call TSAdjointSolve(ts,ierr);CHKERRQ(ierr)
  call VecDot(lambda(1),v,dJv,ierr);CHKERRQ(ierr)
  call TSDestroy(ts,ierr);CHKERRQ(ierr)

  do ieps = 1, 6
    eps = 10.d0**(-ieps)
    call VecWAXPY(Y,eps,v,Y0,ierr);CHKERRQ(ierr)
    call TSModelChecksForward(checks,ctx,time,h,nsteps,Y,PETSC_FALSE,ts)
    call VecDot(w,Y,Jp,ierr);CHKERRQ(ierr)
    rem(ieps) = abs(Jp-J0-eps*dJv)
    write(string,'(a,es9.2,a,es12.4,a,es12.4)') &
      ' TS model adjoint test: eps = ',eps,'  |dJ| = ',abs(Jp-J0), &
      '  Taylor remainder = ',rem(ieps)
    if (ieps > 1) then
      if (rem(ieps) > 0.d0) &
        write(string,'(a,a,f8.2)') trim(string),'  ratio = ', &
          rem(ieps-1)/rem(ieps)
    endif
    call PrintMsg(checks%option,trim(string))
  enddo
  write(string,'(a,es12.4,a,es12.4)') ' TS model adjoint test: J = ',J0, &
    ', dJ/dY0 . v = ',dJv
  call PrintMsg(checks%option,trim(string))

  ! the state of the run proper is untouched
  call VecCopy(Y0,checks%Y,ierr);CHKERRQ(ierr)
  call VecDestroy(Y0,ierr);CHKERRQ(ierr)
  call VecDestroy(Y,ierr);CHKERRQ(ierr)
  call VecDestroy(w,ierr);CHKERRQ(ierr)
  call VecDestroy(v,ierr);CHKERRQ(ierr)
  call VecDestroy(lambda(1),ierr);CHKERRQ(ierr)

end subroutine TSModelAdjointTest

! ************************************************************************** !

subroutine TSModelContractTest(checks,time,dt)
  !
  ! Checks of the model contract at checks%Y, one summary line each:
  !   1. evaluation is repeatable, with no stale reuse by a cache keyed on
  !      time alone or on the vector object;
  !   2. J at A is independent of earlier evaluations/linearizations;
  !   3. J(a) v against central differences, a = 0 and 1/dt;
  !   4. (model operators, serial) the model's preconditioner with an exact
  !      condensed solve inverts J and J^T;
  !   5. the run's preconditioner and its transposed application (adjoint
  !      solves) are transposes of each other.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_model_checks_type) :: checks
  PetscReal :: time
  PetscReal :: dt

  type(option_type), pointer :: option
  Vec :: Y, Ydot, F0, F1, F2, V, Yp, Ydotp, Jv, Jv2, Fd
  Mat :: J
  KSP :: ksp
  PC :: pc
  PCType :: pc_type
  PetscRandom :: rand
  PetscReal :: shifts(2), eps, err, best_err, norm_jv, norm_diff
  PetscReal, pointer :: y_p(:), v_p(:)
  PetscInt :: ishift, ieps, status, i, n
  PetscBool :: same, same2
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  option => checks%option
  Y = checks%Y
  call VecDuplicate(Y,Ydot,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,F0,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,F1,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,F2,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,V,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,Yp,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,Ydotp,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,Jv,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,Jv2,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,Fd,ierr);CHKERRQ(ierr)
  call TSModelChecksNewJacobian(checks,J)

  ! a deterministic perturbation direction, relative to |Y|
  call PetscRandomCreate(option%mycomm,rand,ierr);CHKERRQ(ierr)
  call PetscRandomSetInterval(rand,-1.d0,1.d0,ierr);CHKERRQ(ierr)
  call PetscRandomSetSeed(rand,1234567_8,ierr);CHKERRQ(ierr)
  call PetscRandomSeed(rand,ierr);CHKERRQ(ierr)
  call VecSetRandom(V,rand,ierr);CHKERRQ(ierr)
  call PetscRandomDestroy(rand,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(V,n,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArray(V,v_p,ierr);CHKERRQ(ierr)
  do i = 1, n
    v_p(i) = v_p(i)*max(abs(y_p(i)),1.d-30)
  enddo
  call VecRestoreArray(V,v_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  ! a nonzero rate so that the Ydot dependence is exercised
  call VecCopy(V,Ydot,ierr);CHKERRQ(ierr)
  call VecScale(Ydot,1.d0/dt,ierr);CHKERRQ(ierr)

  ! 1. an evaluation right after one at another time is a reference (no
  !    reuse possible). Then F(A) after F(B) at the same time, and F of a
  !    vector (Jv2) changed in place from A to B.
  call VecWAXPY(Yp,1.d-3,V,Y,ierr);CHKERRQ(ierr)
  call checks%model%Evaluate(time,Y,Ydot,F0,status)
  call checks%model%Evaluate(time+0.5d0*dt,Yp,Ydot,F1,status)
  call checks%model%Evaluate(time,Y,Ydot,F2,status)
  call VecEqual(F0,F2,same,ierr);CHKERRQ(ierr)
  call checks%model%Evaluate(time+0.5d0*dt,Y,Ydot,F1,status)
  call checks%model%Evaluate(time,Yp,Ydot,F1,status)
  call VecCopy(Y,Jv2,ierr);CHKERRQ(ierr)
  call checks%model%Evaluate(time,Jv2,Ydot,F2,status)
  call VecEqual(F0,F2,same2,ierr);CHKERRQ(ierr)
  same = same .and. same2
  call VecCopy(Yp,Jv2,ierr);CHKERRQ(ierr)
  call checks%model%Evaluate(time,Jv2,Ydot,F2,status)
  call VecEqual(F1,F2,same2,ierr);CHKERRQ(ierr)
  if (same .and. same2) then
    string = 'PASS'
  else
    string = 'FAIL'
  endif
  call PrintMsg(option,' TS model contract: repeatable evaluation, &
                &state changed in place: ' // trim(string))

  ! 2. reference J at A (after an evaluation at another time), then J at A
  !    after: F and J at B; F(B) at the same time; F(A), which a model may
  !    reuse
  call checks%model%Evaluate(time+0.5d0*dt,Yp,Ydot,F1,status)
  call checks%model%Linearize(time,Y,Ydot,1.d0/dt,J,J,status)
  call MatMult(J,V,Jv,ierr);CHKERRQ(ierr)
  call VecNorm(Jv,NORM_INFINITY,norm_jv,ierr);CHKERRQ(ierr)
  norm_diff = 0.d0
  do i = 1, 3
    select case(i)
      case(1)
        call checks%model%Evaluate(time+0.5d0*dt,Yp,Ydot,F1,status)
        call checks%model%Linearize(time+0.5d0*dt,Yp,Ydot,0.d0,J,J,status)
      case(2)
        ! B at time is computed (the time changed before it)
        call checks%model%Evaluate(time+0.5d0*dt,Y,Ydot,F1,status)
        call checks%model%Evaluate(time,Yp,Ydot,F1,status)
      case(3)
        call checks%model%Evaluate(time+0.5d0*dt,Yp,Ydot,F1,status)
        call checks%model%Evaluate(time,Y,Ydot,F1,status)
    end select
    call checks%model%Linearize(time,Y,Ydot,1.d0/dt,J,J,status)
    call MatMult(J,V,Jv2,ierr);CHKERRQ(ierr)
    call VecAXPY(Jv2,-1.d0,Jv,ierr);CHKERRQ(ierr)
    call VecNorm(Jv2,NORM_INFINITY,err,ierr);CHKERRQ(ierr)
    norm_diff = max(norm_diff,err)
  enddo
  if (norm_diff == 0.d0) then
    string = 'PASS'
  else
    string = 'FAIL'
  endif
  write(option%io_buffer,'(a,es10.3,2a)') ' TS model contract: &
    &order-independent linearization: max |J1 v - J2 v| = ', &
    norm_diff/max(norm_jv,1.d-300),': ',trim(string)
  call PrintMsg(option)

  ! 3. directional finite differences at two shifts
  shifts(1) = 0.d0
  shifts(2) = 1.d0/dt
  do ishift = 1, 2
    call checks%model%Linearize(time,Y,Ydot,shifts(ishift),J,J,status)
    call MatMult(J,V,Jv,ierr);CHKERRQ(ierr)
    call VecNorm(Jv,NORM_2,norm_jv,ierr);CHKERRQ(ierr)
    best_err = MAX_DOUBLE
    ! central differences: O(eps^2) truncation
    do ieps = 1, 5
      eps = 10.d0**(-1-ieps)
      call VecWAXPY(Yp,eps,V,Y,ierr);CHKERRQ(ierr)
      call VecWAXPY(Ydotp,eps*shifts(ishift),V,Ydot,ierr);CHKERRQ(ierr)
      call checks%model%Evaluate(time,Yp,Ydotp,F1,status)
      call VecWAXPY(Yp,-eps,V,Y,ierr);CHKERRQ(ierr)
      call VecWAXPY(Ydotp,-eps*shifts(ishift),V,Ydot,ierr);CHKERRQ(ierr)
      call checks%model%Evaluate(time,Yp,Ydotp,F2,status)
      call VecWAXPY(Fd,-1.d0,F2,F1,ierr);CHKERRQ(ierr)
      call VecScale(Fd,0.5d0/eps,ierr);CHKERRQ(ierr)
      call VecAXPY(Fd,-1.d0,Jv,ierr);CHKERRQ(ierr)
      call VecNorm(Fd,NORM_2,norm_diff,ierr);CHKERRQ(ierr)
      err = norm_diff/max(norm_jv,1.d-300)
      best_err = min(best_err,err)
    enddo
    if (best_err < 1.d-7) then
      string = 'PASS'
    else
      string = 'FAIL'
    endif
    write(option%io_buffer,'(a,es10.3,a,es10.3,2a)') &
      ' TS model contract: J(a)v vs finite differences, a = ', &
      shifts(ishift),': rel. error ',best_err,': ',trim(string)
    call PrintMsg(option)
  enddo

  ! 4. exact preconditioner
  if (checks%model_jacobian .and. option%comm%size == 1) then
    call checks%model%Linearize(time,Y,Ydot,1.d0/dt,J,J,status)
    call KSPCreate(PETSC_COMM_SELF,ksp,ierr);CHKERRQ(ierr)
    call KSPSetOperators(ksp,J,J,ierr);CHKERRQ(ierr)
    call KSPSetType(ksp,KSPPREONLY,ierr);CHKERRQ(ierr)
    pc_type = PCLU
    call checks%model%SetupLinearSolver(ksp,J,pc_type,PETSC_FALSE)
    call MatMult(J,V,Jv,ierr);CHKERRQ(ierr)
    call VecNorm(Jv,NORM_2,norm_jv,ierr);CHKERRQ(ierr)
    call KSPSolve(ksp,Jv,Jv2,ierr);CHKERRQ(ierr)
    call MatMult(J,Jv2,Fd,ierr);CHKERRQ(ierr)
    call VecAXPY(Fd,-1.d0,Jv,ierr);CHKERRQ(ierr)
    call VecNorm(Fd,NORM_2,norm_diff,ierr);CHKERRQ(ierr)
    err = norm_diff/max(norm_jv,1.d-300)
    call MatMultTranspose(J,V,Jv,ierr);CHKERRQ(ierr)
    call VecNorm(Jv,NORM_2,norm_jv,ierr);CHKERRQ(ierr)
    call KSPSolveTranspose(ksp,Jv,Jv2,ierr);CHKERRQ(ierr)
    call MatMultTranspose(J,Jv2,Fd,ierr);CHKERRQ(ierr)
    call VecAXPY(Fd,-1.d0,Jv,ierr);CHKERRQ(ierr)
    call VecNorm(Fd,NORM_2,norm_diff,ierr);CHKERRQ(ierr)
    best_err = norm_diff/max(norm_jv,1.d-300)
    if (max(err,best_err) < 1.d-10) then
      string = 'PASS'
    else
      string = 'FAIL'
    endif
    write(option%io_buffer,'(a,es10.3,a,es10.3,2a)') &
      ' TS model contract: exact preconditioner: |J z - b|/|b| = ',err, &
      ', transpose ',best_err,': ',trim(string)
    call PrintMsg(option)
    call KSPDestroy(ksp,ierr);CHKERRQ(ierr)
  endif

  ! 5. y.(P^-1 x) = (P^-T y).x for the run's preconditioner
  call MatDestroy(J,ierr);CHKERRQ(ierr)
  call TSModelChecksNewJacobian(checks,J)
  call checks%model%Linearize(time,Y,Ydot,1.d0/dt,J,J,status)
  call KSPCreate(option%mycomm,ksp,ierr);CHKERRQ(ierr)
  call KSPSetOperators(ksp,J,J,ierr);CHKERRQ(ierr)
  call KSPSetType(ksp,KSPPREONLY,ierr);CHKERRQ(ierr)
  if (checks%model_jacobian) then
    call checks%model%SetupLinearSolver(ksp,J,checks%pc_type, &
                                        checks%shift_inblocks)
  else if (len_trim(checks%pc_type) > 0) then
    call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
    call PCSetType(pc,checks%pc_type,ierr);CHKERRQ(ierr)
  endif
  ! x = V, y = F0 (the initial residual): independent
  call KSPSolve(ksp,V,Jv,ierr);CHKERRQ(ierr)
  call KSPSolveTranspose(ksp,F0,Jv2,ierr);CHKERRQ(ierr)
  call VecDot(F0,Jv,err,ierr);CHKERRQ(ierr)
  call VecDot(Jv2,V,best_err,ierr);CHKERRQ(ierr)
  call VecNorm(F0,NORM_2,norm_jv,ierr);CHKERRQ(ierr)
  call VecNorm(Jv,NORM_2,norm_diff,ierr);CHKERRQ(ierr)
  err = abs(err-best_err)/max(norm_jv*norm_diff,1.d-300)
  if (err < 1.d-10) then
    string = 'PASS'
  else
    string = 'FAIL'
  endif
  write(option%io_buffer,'(a,es10.3,2a)') &
    ' TS model contract: transposed preconditioner: |y.P^-1 x - &
    &P^-T y.x| / |y||P^-1 x| = ',err,': ',trim(string)
  call PrintMsg(option)
  call KSPDestroy(ksp,ierr);CHKERRQ(ierr)

  call VecDestroy(Ydot,ierr);CHKERRQ(ierr)
  call VecDestroy(F0,ierr);CHKERRQ(ierr)
  call VecDestroy(F1,ierr);CHKERRQ(ierr)
  call VecDestroy(F2,ierr);CHKERRQ(ierr)
  call VecDestroy(V,ierr);CHKERRQ(ierr)
  call VecDestroy(Yp,ierr);CHKERRQ(ierr)
  call VecDestroy(Ydotp,ierr);CHKERRQ(ierr)
  call VecDestroy(Jv,ierr);CHKERRQ(ierr)
  call VecDestroy(Jv2,ierr);CHKERRQ(ierr)
  call VecDestroy(Fd,ierr);CHKERRQ(ierr)
  call MatDestroy(J,ierr);CHKERRQ(ierr)

end subroutine TSModelContractTest

! ************************************************************************** !

subroutine TSModelColumnCheck(checks,time,dt)
  !
  ! Entry-wise check of J(1/dt) against central differences, one column at
  ! a time (serial); a norm-wise check hides wrong small entries. Counts,
  ! and prints the first 60 of, the entries with relative error > 1e-4
  ! that are not negligible in their column and not FD noise (the
  ! differences at two step sizes agree), and that larger steps do not
  ! resolve: at a zero entry the step (1e-6 of its atol) can be below the
  ! roundoff of a large row, where the differences at both steps can round
  ! to the same wrong value.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_model_checks_type) :: checks
  PetscReal :: time
  PetscReal :: dt

  type(option_type), pointer :: option
  Vec :: Y, Ydot, Je, E
  Vec :: Fdiff(4), work(3)
  Mat :: J
  PetscReal, pointer :: y_p(:), fd1_p(:), fd2_p(:), fd3_p(:), fd4_p(:)
  PetscReal, pointer :: je_p(:), e_p(:), a_p(:)
  PetscReal :: a, h, h0, fd, fd2, an, colmax, worst
  PetscInt :: n, i, k, nbad, ndof, istep, nunresolved, status
  PetscBool :: larger
  PetscErrorCode :: ierr

  option => checks%option
  if (option%comm%size > 1) return
  ndof = checks%model%ndof
  Y = checks%Y
  a = 1.d0/dt
  call VecDuplicate(Y,Ydot,ierr);CHKERRQ(ierr)
  do i = 1, 4
    call VecDuplicate(Y,Fdiff(i),ierr);CHKERRQ(ierr)
  enddo
  call VecDuplicate(Y,Je,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,E,ierr);CHKERRQ(ierr)
  do i = 1, 3
    call VecDuplicate(Y,work(i),ierr);CHKERRQ(ierr)
  enddo
  call VecZeroEntries(Ydot,ierr);CHKERRQ(ierr)
  call TSModelChecksNewJacobian(checks,J)
  call checks%model%Linearize(time,Y,Ydot,a,J,J,status)
  call VecGetLocalSize(Y,n,ierr);CHKERRQ(ierr)
  nbad = 0
  nunresolved = 0
  worst = 0.d0
  do k = 1, n
    call VecZeroEntries(E,ierr);CHKERRQ(ierr)
    call VecGetArray(E,e_p,ierr);CHKERRQ(ierr)
    e_p(k) = 1.d0
    call VecRestoreArray(E,e_p,ierr);CHKERRQ(ierr)
    call MatMult(J,E,Je,ierr);CHKERRQ(ierr)
    ! step relative to the entry, or to its atol when (nearly) zero: a step
    ! below roundoff differences to exactly zero
    call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
    call VecGetArrayRead(checks%vatol,a_p,ierr);CHKERRQ(ierr)
    h = 1.d-6*max(abs(y_p(k)),a_p(k),1.d-20)
    h0 = h
    call VecRestoreArrayRead(checks%vatol,a_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
    do istep = 1, 2
      call TSModelChecksCentralDifference(checks,time,Y,Ydot,a,E,h, &
                                          Fdiff(istep),work)
      h = 10.d0*h
    enddo
    call VecGetArrayRead(Fdiff(1),fd1_p,ierr);CHKERRQ(ierr)
    call VecGetArrayRead(Fdiff(2),fd2_p,ierr);CHKERRQ(ierr)
    call VecGetArrayRead(Je,je_p,ierr);CHKERRQ(ierr)
    colmax = 0.d0
    do i = 1, n
      colmax = max(colmax,abs(je_p(i)),abs(fd1_p(i)))
    enddo
    larger = PETSC_FALSE
    do i = 1, n
      fd = fd1_p(i)
      fd2 = fd2_p(i)
      an = je_p(i)
      if (max(abs(fd),abs(an)) <= 1.d-8*colmax) cycle
      if (abs(fd-an) <= 1.d-4*max(abs(fd),abs(an))) cycle
      ! unresolved: the perturbation is below the roundoff of the row
      if (fd == 0.d0 .and. fd2 == 0.d0) then
        nunresolved = nunresolved + 1
        cycle
      endif
      ! noise: the two differences disagree with each other
      if (abs(fd-fd2) > 1.d-2*abs(fd-an)) cycle
      ! roundoff: differences at 1e4 and 1e5 times the step agree
      if (.not.larger) then
        do istep = 3, 4
          call TSModelChecksCentralDifference(checks,time,Y,Ydot,a,E, &
                                              10.d0**(istep+1)*h0, &
                                              Fdiff(istep),work)
        enddo
        call VecGetArrayRead(Fdiff(3),fd3_p,ierr);CHKERRQ(ierr)
        call VecGetArrayRead(Fdiff(4),fd4_p,ierr);CHKERRQ(ierr)
        larger = PETSC_TRUE
      endif
      if (abs(fd3_p(i)-an) <= 1.d-4*max(abs(fd3_p(i)),abs(an)) .and. &
          abs(fd4_p(i)-an) <= 1.d-4*max(abs(fd4_p(i)),abs(an))) cycle
      nbad = nbad + 1
      worst = max(worst,abs(fd-an)/max(abs(fd),abs(an)))
      if (nbad <= 60) then
        if (ndof > 0) then
          write(option%io_buffer,'(a,4i6,3es13.5)') &
            ' TS model column check: cell,row,cell,col,J,FD,rel = ', &
            (i-1)/ndof+1,mod(i-1,ndof)+1,(k-1)/ndof+1,mod(k-1,ndof)+1, &
            an,fd,abs(fd-an)/max(abs(fd),abs(an))
        else
          ! not laid out cell by cell: local entry indices
          write(option%io_buffer,'(a,2i8,3es13.5)') &
            ' TS model column check: row,col,J,FD,rel = ',i,k, &
            an,fd,abs(fd-an)/max(abs(fd),abs(an))
        endif
        call PrintMsg(option)
      endif
    enddo
    if (larger) then
      call VecRestoreArrayRead(Fdiff(4),fd4_p,ierr);CHKERRQ(ierr)
      call VecRestoreArrayRead(Fdiff(3),fd3_p,ierr);CHKERRQ(ierr)
    endif
    call VecRestoreArrayRead(Je,je_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayRead(Fdiff(2),fd2_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayRead(Fdiff(1),fd1_p,ierr);CHKERRQ(ierr)
  enddo
  write(option%io_buffer,'(a,i8,a,es10.3,a,i8,a)') &
    ' TS model column check: ',nbad,' entries off by more than 1e-4; &
    &worst ',worst,'; ',nunresolved,' entries unresolved (the difference &
    &is below roundoff)'
  call PrintMsg(option)
  call VecDestroy(Ydot,ierr);CHKERRQ(ierr)
  do i = 1, 4
    call VecDestroy(Fdiff(i),ierr);CHKERRQ(ierr)
  enddo
  call VecDestroy(Je,ierr);CHKERRQ(ierr)
  call VecDestroy(E,ierr);CHKERRQ(ierr)
  do i = 1, 3
    call VecDestroy(work(i),ierr);CHKERRQ(ierr)
  enddo
  call MatDestroy(J,ierr);CHKERRQ(ierr)

end subroutine TSModelColumnCheck

! ************************************************************************** !

subroutine TSModelChecksCentralDifference(checks,time,Y,Ydot,a,E,h,Fd,work)
  !
  ! Fd = (F(Y+hE, Ydot+ahE) - F(Y-hE, Ydot-ahE)) / 2h.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_model_checks_type) :: checks
  PetscReal :: time
  Vec :: Y, Ydot, E, Fd
  PetscReal :: a, h
  Vec :: work(3)

  PetscInt :: status
  PetscErrorCode :: ierr

  call VecWAXPY(work(1),h,E,Y,ierr);CHKERRQ(ierr)
  call VecWAXPY(work(2),a*h,E,Ydot,ierr);CHKERRQ(ierr)
  call checks%model%Evaluate(time,work(1),work(2),Fd,status)
  call VecWAXPY(work(1),-h,E,Y,ierr);CHKERRQ(ierr)
  call VecWAXPY(work(2),-a*h,E,Ydot,ierr);CHKERRQ(ierr)
  call checks%model%Evaluate(time,work(1),work(2),work(3),status)
  call VecAXPY(Fd,-1.d0,work(3),ierr);CHKERRQ(ierr)
  call VecScale(Fd,0.5d0/h,ierr);CHKERRQ(ierr)

end subroutine TSModelChecksCentralDifference

end module TS_Model_Checks_module
