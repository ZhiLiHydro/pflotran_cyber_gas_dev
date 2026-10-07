module TS_Adjoint_class
  !
  ! Discrete adjoint (TSAdjointSolve) of a TS_Model_class model integrated
  ! by PMC_TS_Driver_class. Configure prepares the forward TS; after the
  ! final time, Solve runs the backward sweep while the model replays its
  ! recorded forcing, so the sweep linearizes the forward steps.
  !
  ! A concrete adjoint binds costs and parameters to one kind of model
  ! (Setup, TerminalGradients, ParameterJacobian, Report, Get/SetParameters).
  ! Report keeps the values it prints (KeepResult) for the regression file.
  ! Solve returns a status; quiet (a rerun of the forward run, see
  ! PMC_TS_Driver_class) only returns it, with the costs and gradients
  ! (GetCostValues, GetCostGradients), and prints or keeps nothing.
  !
  ! Only backward Euler (PETSc has no BDF adjoint) and the basic (disk) or
  ! memory trajectory are supported. Runs with regime transitions (TSEvent)
  ! skip the sweep: their jump maps and event-time sensitivities are not
  ! included.
  !
  ! The forward run is one TSSolve per waypoint interval. For the memory
  ! trajectory EndInterval stores each interval's last step, which PETSc
  ! takes for the run's end, by setting the converged reason and calling
  ! TSTrajectorySet itself: recheck it after a PETSc upgrade (the memory
  ! trajectory check of make verify compares it with the disk one).
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

  ! Solve: a gradient, or why there is none
  PetscInt, parameter, public :: TS_ADJ_OK = 0
  PetscInt, parameter, public :: TS_ADJ_SKIPPED_EVENTS = 1
  PetscInt, parameter, public :: TS_ADJ_SKIPPED_TERMINAL_TIME = 2
  PetscInt, parameter, public :: TS_ADJ_INCOMPLETE = 3
  ! the forward run did not complete (set by the driver)
  PetscInt, parameter, public :: TS_ADJ_FORWARD_FAILED = 4

  ! what PETSc hands the parameter Jacobian callback
  type, public :: ts_adjoint_ctx_type
    class(ts_adjoint_type), pointer :: adjoint => null()
    class(ts_model_type), pointer :: model => null()
  end type ts_adjoint_ctx_type

  type, public, abstract :: ts_adjoint_type
    type(option_type), pointer :: option => null()
    PetscInt :: ncosts = 0
    PetscInt :: nparams = 0
    ! J_i at the final time of the last sweep (TerminalGradients)
    PetscReal, allocatable :: cost_value(:)
    ! a rerun: Solve prints, reports and keeps nothing
    PetscBool :: quiet = PETSC_FALSE
    ! dJ_i/dY, laid out as the TS solution, and dJ_i/dp
    Vec, allocatable :: lambda(:)
    Vec, allocatable :: mu(:)
    ! dF/dp: a dense matrix with one column per parameter
    Mat :: Jacp
    ! the forward run's initial state and time, and its final time
    Vec :: Y0
    PetscReal :: t0 = 0.d0
    PetscReal :: final_time = 0.d0
    ! -ts_trajectory_type memory, and the last step EndInterval stored
    PetscBool :: memory_trajectory = PETSC_FALSE
    PetscInt :: interval_end_step = -1
    type(ts_adjoint_ctx_type), pointer :: ctx => null()
    ! the reported costs and gradients, for the regression file
    character(len=MAXWORDLENGTH), allocatable :: result_name(:)
    PetscReal, allocatable :: result_value(:)
  contains
    procedure(TSAdjSetupIface), deferred :: Setup
    procedure(TSAdjTerminalTimeIface), deferred :: CheckTerminalTime
    procedure(TSAdjTerminalGradientsIface), deferred :: TerminalGradients
    procedure(TSAdjParameterJacobianIface), deferred :: ParameterJacobian
    procedure(TSAdjReportIface), deferred :: Report
    procedure(TSAdjGetParametersIface), deferred :: GetParameters
    procedure(TSAdjSetParametersIface), deferred :: SetParameters
    procedure, public :: Configure => TSAdjConfigure
    procedure, public :: CheckMethod => TSAdjCheckMethod
    procedure, public :: RecordInitialState => TSAdjRecordInitialState
    procedure, public :: EndInterval => TSAdjEndInterval
    procedure, public :: Solve => TSAdjSolveBackward
    procedure, public :: BeginRerun => TSAdjBeginRerun
    procedure, public :: GetCostValues => TSAdjGetCostValues
    procedure, public :: GetCostGradients => TSAdjGetCostGradients
    procedure, public, non_overridable :: ResetResults => TSAdjResetResults
    procedure, public, non_overridable :: KeepResult => TSAdjKeepResult
    procedure, public :: RegressionOutput => TSAdjRegressionOutput
    procedure, public :: Destroy => TSAdjDestroy
  end type ts_adjoint_type

  abstract interface

    subroutine TSAdjSetupIface(this,model,final_time)
      ! Resolves the costs and parameters for model and sets ncosts and
      ! nparams; rejects configurations the adjoint does not support.
      use petscts
      import :: ts_adjoint_type, ts_model_type
      implicit none
      class(ts_adjoint_type) :: this
      class(ts_model_type), target :: model
      PetscReal :: final_time
    end subroutine TSAdjSetupIface

    subroutine TSAdjTerminalTimeIface(this,time,valid)
      ! Whether the costs are defined at the final time of the run.
      use petscts
      import :: ts_adjoint_type
      implicit none
      class(ts_adjoint_type) :: this
      PetscReal :: time
      PetscBool :: valid
    end subroutine TSAdjTerminalTimeIface

    subroutine TSAdjTerminalGradientsIface(this,model,time,Y)
      ! lambda(i) = dJ_i/dY at the final state Y.
      use petscts
      import :: ts_adjoint_type, ts_model_type
      implicit none
      class(ts_adjoint_type) :: this
      class(ts_model_type), target :: model
      PetscReal :: time
      Vec :: Y
    end subroutine TSAdjTerminalGradientsIface

    subroutine TSAdjParameterJacobianIface(this,model,time,Y,Ydot,Jacp,status)
      ! Jacp = dF/dp at (time,Y,Ydot), assembled.
      use petscts
      import :: ts_adjoint_type, ts_model_type
      implicit none
      class(ts_adjoint_type) :: this
      class(ts_model_type), target :: model
      PetscReal :: time
      Vec :: Y
      Vec :: Ydot
      Mat :: Jacp
      PetscInt :: status
    end subroutine TSAdjParameterJacobianIface

    subroutine TSAdjReportIface(this,model,time0,Y0)
      ! Prints the costs and their gradients; lambda(i) holds dJ_i/dY(t0)
      ! at the initial state Y0.
      use petscts
      import :: ts_adjoint_type, ts_model_type
      implicit none
      class(ts_adjoint_type) :: this
      class(ts_model_type), target :: model
      PetscReal :: time0
      Vec :: Y0
    end subroutine TSAdjReportIface

    subroutine TSAdjGetParametersIface(this,p)
      ! p(j): the current value of parameter j (deck order), on every
      ! process.
      use petscts
      import :: ts_adjoint_type
      implicit none
      class(ts_adjoint_type) :: this
      PetscReal :: p(:)
    end subroutine TSAdjGetParametersIface

    subroutine TSAdjSetParametersIface(this,model,p)
      ! Sets parameter j to p(j) in the data model reads, and expires the
      ! model's cache.
      use petscts
      import :: ts_adjoint_type, ts_model_type
      implicit none
      class(ts_adjoint_type) :: this
      class(ts_model_type) :: model
      PetscReal :: p(:)
    end subroutine TSAdjSetParametersIface

  end interface

  public :: TSAdjInit, &
            TSAdjDestroy

contains

! ************************************************************************** !

subroutine TSAdjInit(this,option)
  !
  ! Initializes the base; call first from a concrete adjoint's Create.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  type(option_type), pointer :: option

  this%option => option
  PetscObjectNullify(this%Jacp)
  PetscObjectNullify(this%Y0)

end subroutine TSAdjInit

! ************************************************************************** !

subroutine TSAdjConfigure(this,ts,model,Y,final_time)
  !
  ! Backward Euler, a saved trajectory, dF/dp and forcing recording. Call
  ! after TSSetIJacobian and before TSSetFromOptions; CheckMethod after it.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_adjoint_type), target :: this
  TS :: ts
  class(ts_model_type), pointer :: model
  Vec :: Y
  PetscReal :: final_time

  PetscInt :: nlocal, nglobal, i
  PetscErrorCode :: ierr

  call this%Setup(model,final_time)
  this%final_time = final_time
  allocate(this%cost_value(this%ncosts))
  this%cost_value = 0.d0

  allocate(this%ctx)
  this%ctx%adjoint => this
  this%ctx%model => model

  call VecGetLocalSize(Y,nlocal,ierr);CHKERRQ(ierr)
  call VecGetSize(Y,nglobal,ierr);CHKERRQ(ierr)
  allocate(this%lambda(this%ncosts))
  allocate(this%mu(this%ncosts))
  do i = 1, this%ncosts
    call VecDuplicate(Y,this%lambda(i),ierr);CHKERRQ(ierr)
  enddo
  call MatCreate(this%option%mycomm,this%Jacp,ierr);CHKERRQ(ierr)
  call MatSetSizes(this%Jacp,nlocal,PETSC_DECIDE,nglobal,this%nparams, &
                   ierr);CHKERRQ(ierr)
  call MatSetType(this%Jacp,MATDENSE,ierr);CHKERRQ(ierr)
  call MatSetUp(this%Jacp,ierr);CHKERRQ(ierr)
  do i = 1, this%ncosts
    call MatCreateVecs(this%Jacp,this%mu(i),PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  enddo
  call TSSetIJacobianP(ts,this%Jacp,TSAdjIJacobianP,this%ctx, &
                       ierr);CHKERRQ(ierr)

  call TSSetType(ts,TSBEULER,ierr);CHKERRQ(ierr)
  call TSSetSaveTrajectory(ts,ierr);CHKERRQ(ierr)
  call model%RecordForcing()

  call PrintMsg(this%option,'  TSAdjoint: backward Euler, saved trajectory, &
    &parameter Jacobian and recorded forcing')

end subroutine TSAdjConfigure

! ************************************************************************** !

subroutine TSAdjCheckMethod(this,ts)
  !
  ! Rejects command-line choices outside the supported scope (see module
  ! header), including a memory trajectory that drops steps or stages.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  TS :: ts

  character(len=*), parameter :: limits(4) = &
    [character(len=32) :: '-ts_trajectory_max_cps_ram', &
     '-ts_trajectory_max_cps_disk','-ts_trajectory_max_units_ram', &
     '-ts_trajectory_max_units_disk']
  TSType :: ts_type
  TSTrajectory :: trajectory
  PetscBool :: basic, solution_only, flag
  PetscInt :: i
  PetscErrorCode :: ierr

  call TSGetType(ts,ts_type,ierr);CHKERRQ(ierr)
  if (trim(ts_type) /= TSBEULER) then
    this%option%io_buffer = 'ADJOINT requires backward Euler; the TS type "' &
      // trim(ts_type) // '" set on the command line is not supported (BDF &
      &has no discrete adjoint in PETSc).'
    call PrintErrMsg(this%option)
  endif
  call TSGetTrajectory(ts,trajectory,ierr);CHKERRQ(ierr)
  call PetscObjectTypeCompare(trajectory,TSTRAJECTORYBASIC,basic, &
                              ierr);CHKERRQ(ierr)
  call PetscObjectTypeCompare(trajectory,TSTRAJECTORYMEMORY, &
                              this%memory_trajectory,ierr);CHKERRQ(ierr)
  if (.not.(basic .or. this%memory_trajectory)) then
    call PrintErrMsg(this%option,'ADJOINT supports the TS trajectories basic &
      &(the default) and memory (-ts_trajectory_type).')
  endif
  ! solution-only mode: PETSc rejects it for the basic trajectory, but only
  ! after the forward run; the memory trajectory drops the stored step that
  ! ends each TSSolve (one per waypoint interval)
  call TSTrajectoryGetSolutionOnly(trajectory,solution_only, &
                                   ierr);CHKERRQ(ierr)
  if (solution_only) then
    call PrintErrMsg(this%option,'ADJOINT does not support &
      &-ts_trajectory_solution_only.')
  endif
  if (this%memory_trajectory) then
    ! checkpointing recomputes steps: not verified
    do i = 1, size(limits)
      call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                               trim(limits(i)),flag,ierr);CHKERRQ(ierr)
      if (flag) then
        call PrintErrMsg(this%option,'ADJOINT keeps every step in memory: &
          &' // trim(limits(i)) // ' (checkpointing) is not supported.')
      endif
    enddo
  endif
  ! TSHistory rejects a step set twice (EndInterval sets a step again, a
  ! rerun every step; the basic trajectory rewrites their files); the
  ! sweep reads stages only
  call TSTrajectorySetUseHistory(trajectory,PETSC_FALSE,ierr);CHKERRQ(ierr)

end subroutine TSAdjCheckMethod

! ************************************************************************** !

subroutine TSAdjRecordInitialState(this,time,Y)
  !
  ! Keeps the forward run's initial state for the report.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  PetscReal :: time
  Vec :: Y

  PetscErrorCode :: ierr

  if (PetscObjectIsNull(this%Y0)) then
    call VecDuplicate(Y,this%Y0,ierr);CHKERRQ(ierr)
  endif
  call VecCopy(Y,this%Y0,ierr);CHKERRQ(ierr)
  this%t0 = time

end subroutine TSAdjRecordInitialState

! ************************************************************************** !

subroutine TSAdjEndInterval(this,ts,Y)
  !
  ! After a successful TSSolve of a waypoint interval. The memory
  ! trajectory skips the step that ends a TSSolve, taking it for the end of
  ! the run; store it here for intervals that end before the final time.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  TS :: ts
  Vec :: Y

  TSTrajectory :: trajectory
  TSConvergedReason :: reason
  PetscReal :: time
  PetscInt :: step
  PetscErrorCode :: ierr

  if (.not.this%memory_trajectory) return
  call TSGetTime(ts,time,ierr);CHKERRQ(ierr)
  if (time >= this%final_time - &
      64.d0*epsilon(1.d0)*max(1.d0,abs(this%final_time))) return
  call TSGetStepNumber(ts,step,ierr);CHKERRQ(ierr)
  ! an interval that took no step
  if (step == this%interval_end_step .or. step == 0) return
  this%interval_end_step = step
  call TSGetTrajectory(ts,trajectory,ierr);CHKERRQ(ierr)
  call TSGetConvergedReason(ts,reason,ierr);CHKERRQ(ierr)
  call TSSetConvergedReason(ts,TS_CONVERGED_ITERATING,ierr);CHKERRQ(ierr)
  call TSTrajectorySet(trajectory,ts,step,time,Y,ierr);CHKERRQ(ierr)
  call TSSetConvergedReason(ts,reason,ierr);CHKERRQ(ierr)

end subroutine TSAdjEndInterval

! ************************************************************************** !

subroutine TSAdjSolveBackward(this,ts,model,Y,num_events,status)
  !
  ! Backward sweep and report after a complete forward run with final
  ! state Y; leaves the TS and model at that state. num_events: this
  ! process's regime transitions; any on any process skips the sweep.
  ! status: TS_ADJ_OK, or why there is no gradient. Unless quiet, a skip
  ! sets EXIT_ADJOINT_SKIPPED and an incomplete sweep stops the run.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  TS :: ts
  class(ts_model_type), target :: model
  Vec :: Y
  PetscInt :: num_events
  PetscInt :: status

  type(option_type), pointer :: option
  TSConvergedReason :: reason
  Vec :: Y_final
  PetscReal :: final_time, final_dt
  PetscInt :: final_step, step, i, events
  PetscBool :: valid
  PetscErrorCode :: ierr

  option => this%option
  status = TS_ADJ_OK
  ! collective: the backward sweep is
  call MPI_Allreduce(num_events,events,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  if (events > 0) then
    status = TS_ADJ_SKIPPED_EVENTS
    if (this%quiet) return
    write(option%io_buffer,'(a,i0,a)') 'ADJOINT cannot run after ', &
      events,' regime transition(s) (TSEvent): the adjoint does not &
      &include their jump maps and event-time sensitivities.'
    call PrintMsg(option)
    call PrintMsg(option,'ADJOINT: forward solve completed; adjoint skipped.')
    option%driver%exit_code = EXIT_ADJOINT_SKIPPED
    return
  endif
  call TSGetTime(ts,final_time,ierr);CHKERRQ(ierr)
  call this%CheckTerminalTime(final_time,valid)
  if (.not.valid) then
    status = TS_ADJ_SKIPPED_TERMINAL_TIME
    if (this%quiet) return
    call PrintMsg(option,'ADJOINT: a cost is not defined at the final time &
      &of the run; forward solve completed; adjoint skipped.')
    option%driver%exit_code = EXIT_ADJOINT_SKIPPED
    return
  endif

  ! the sweep rewinds the TS to t0
  call VecDuplicate(Y,Y_final,ierr);CHKERRQ(ierr)
  call VecCopy(Y,Y_final,ierr);CHKERRQ(ierr)
  call TSGetStepNumber(ts,final_step,ierr);CHKERRQ(ierr)
  call TSGetTimeStep(ts,final_dt,ierr);CHKERRQ(ierr)

  call this%TerminalGradients(model,final_time,Y)
  ! TSAdjointSolve accumulates into mu; dJ/dp(T) = 0 for a terminal cost
  do i = 1, this%ncosts
    call VecZeroEntries(this%mu(i),ierr);CHKERRQ(ierr)
  enddo
  call TSSetCostGradients(ts,this%ncosts,this%lambda,this%mu, &
                          ierr);CHKERRQ(ierr)

  call model%BeginReplay()
  if (.not.this%quiet) &
    call PrintMsg(option,'=== TSAdjointSolve: running backward ===')
  ! CHKERRA: CHKERRQ would return silently and the run would report success
  call TSAdjointSolve(ts,ierr);CHKERRA(ierr)
  call TSGetConvergedReason(ts,reason,ierr);CHKERRQ(ierr)
  call TSGetStepNumber(ts,step,ierr);CHKERRQ(ierr)
  if (reason%v <= 0 .or. step /= 0) then
    status = TS_ADJ_INCOMPLETE
    if (.not.this%quiet) then
      write(option%io_buffer,'(a,i0,a,i0,a)') &
        'ADJOINT incomplete: backward solve stopped at step ',step, &
        ' with TSConvergedReason ',reason%v,'; no gradient will be &
        &reported.'
      call PrintErrMsg(option)
    endif
  endif

  if (.not.this%quiet) then
    call this%ResetResults()
    call this%Report(model,this%t0,this%Y0)
  endif
  call model%EndReplay()

  ! the final forward state, for output and the regression file
  call VecCopy(Y_final,Y,ierr);CHKERRQ(ierr)
  call VecDestroy(Y_final,ierr);CHKERRQ(ierr)
  call TSSetTime(ts,final_time,ierr);CHKERRQ(ierr)
  call TSSetStepNumber(ts,final_step,ierr);CHKERRQ(ierr)
  call TSSetTimeStep(ts,final_dt,ierr);CHKERRQ(ierr)

end subroutine TSAdjSolveBackward

! ************************************************************************** !

subroutine TSAdjBeginRerun(this,ts)
  !
  ! Before a rerun of the forward run: forgets the last sweep's costs and
  ! the trajectory's steps. The memory trajectory is created again: set up
  ! a second time (TSTrajectoryReset), it keeps the container its stack
  ! grew to but sizes it for the first waypoint interval, and loses
  ! elements (PETSc 3.24.5). The basic one (on disk) rewrites the files of
  ! the steps the rerun sets; it cannot be set up again in its directory.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  TS :: ts

  TSTrajectory :: trajectory
  PetscErrorCode :: ierr

  this%interval_end_step = -1
  if (allocated(this%cost_value)) this%cost_value = 0.d0
  if (this%memory_trajectory) then
    call TSResetTrajectory(ts,ierr);CHKERRQ(ierr)
    call TSGetTrajectory(ts,trajectory,ierr);CHKERRQ(ierr)
    call TSTrajectorySetType(trajectory,ts,TSTRAJECTORYMEMORY, &
                             ierr);CHKERRQ(ierr)
    call TSTrajectorySetFromOptions(trajectory,ts,ierr);CHKERRQ(ierr)
    call TSTrajectorySetUseHistory(trajectory,PETSC_FALSE,ierr);CHKERRQ(ierr)
  endif

end subroutine TSAdjBeginRerun

! ************************************************************************** !

subroutine TSAdjGetCostValues(this,J)
  !
  ! J(i) = J_i of the last sweep, on every process.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  PetscReal :: J(:)

  J = this%cost_value

end subroutine TSAdjGetCostValues

! ************************************************************************** !

subroutine TSAdjGetCostGradients(this,dJ_dp)
  !
  ! dJ_dp(j,i) = dJ_i/dp_j of the last sweep, on every process.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  PetscReal :: dJ_dp(:,:)

  PetscReal, pointer :: mu_p(:)
  PetscInt :: i, j, col_start, col_end
  PetscErrorCode :: ierr

  dJ_dp = 0.d0
  do i = 1, this%ncosts
    call VecGetOwnershipRange(this%mu(i),col_start,col_end, &
                              ierr);CHKERRQ(ierr)
    call VecGetArrayRead(this%mu(i),mu_p,ierr);CHKERRQ(ierr)
    do j = col_start+1, col_end
      dJ_dp(j,i) = mu_p(j-col_start)
    enddo
    call VecRestoreArrayRead(this%mu(i),mu_p,ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,dJ_dp,size(dJ_dp),MPI_DOUBLE_PRECISION, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)

end subroutine TSAdjGetCostGradients

! ************************************************************************** !

subroutine TSAdjIJacobianP(ts,time,Y,Ydot,shift,Jacp,ctx,ierr)
  !
  ! dF/dp for TSAdjointSolve; dF/dp does not depend on the shift.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: Y, Ydot
  PetscReal :: shift
  Mat :: Jacp
  type(ts_adjoint_ctx_type) :: ctx
  PetscErrorCode :: ierr

  PetscInt :: status

  call ctx%adjoint%ParameterJacobian(ctx%model,time,Y,Ydot,Jacp,status)
  if (status /= TS_MODEL_OK) then
    call PrintErrMsg(ctx%adjoint%option,'ADJOINT: the parameter Jacobian &
      &could not be evaluated at a state of the forward run.')
  endif
  ierr = 0

end subroutine TSAdjIJacobianP

! ************************************************************************** !

subroutine TSAdjKeepResult(this,name,value)
  !
  ! Keeps a reported value for the regression file.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  character(len=*) :: name
  PetscReal :: value

  if (.not.allocated(this%result_name)) then
    allocate(this%result_name(0))
    allocate(this%result_value(0))
  endif
  this%result_name = [character(len=MAXWORDLENGTH) :: this%result_name, name]
  this%result_value = [this%result_value, value]

end subroutine TSAdjKeepResult

! ************************************************************************** !

subroutine TSAdjResetResults(this)
  !
  ! Forgets the kept values.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_adjoint_type) :: this

  if (allocated(this%result_name)) deallocate(this%result_name)
  if (allocated(this%result_value)) deallocate(this%result_value)

end subroutine TSAdjResetResults

! ************************************************************************** !

subroutine TSAdjRegressionOutput(this,fid_out)
  !
  ! Writes the kept values as a section of the regression file.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
  implicit none

  class(ts_adjoint_type) :: this
  PetscInt :: fid_out

  PetscInt :: i

  if (.not.allocated(this%result_name)) return
  if (.not.OptionIsIORank(this%option)) return
  write(fid_out,'(''-- GENERIC: Adjoint --'')')
  do i = 1, size(this%result_name)
    write(fid_out,'(3x,a,'': '',es21.13)') trim(this%result_name(i)), &
      this%result_value(i)
  enddo

end subroutine TSAdjRegressionOutput

! ************************************************************************** !

subroutine TSAdjDestroy(this)
  !
  ! Frees the gradients and the parameter Jacobian. A concrete adjoint's
  ! Destroy calls this last.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Petsc_Utility_module, only : PUVecDestroy

  implicit none

  class(ts_adjoint_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  if (allocated(this%lambda)) then
    do i = 1, size(this%lambda)
      call PUVecDestroy(this%lambda(i))
    enddo
    deallocate(this%lambda)
  endif
  if (allocated(this%mu)) then
    do i = 1, size(this%mu)
      call PUVecDestroy(this%mu(i))
    enddo
    deallocate(this%mu)
  endif
  if (.not.PetscObjectIsNull(this%Jacp)) then
    call MatDestroy(this%Jacp,ierr);CHKERRQ(ierr)
  endif
  call PUVecDestroy(this%Y0)
  if (associated(this%ctx)) then
    nullify(this%ctx%adjoint)
    nullify(this%ctx%model)
    deallocate(this%ctx)
    nullify(this%ctx)
  endif
  call this%ResetResults()
  if (allocated(this%cost_value)) deallocate(this%cost_value)

end subroutine TSAdjDestroy

end module TS_Adjoint_class
