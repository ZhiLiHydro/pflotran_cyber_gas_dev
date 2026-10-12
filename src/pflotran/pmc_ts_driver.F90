module PMC_TS_Driver_class
  !
  ! Process model coupler that integrates an equation model (TS_Model_class)
  ! with PETSc TS: one TSSolve per waypoint interval, each with its own
  ! forcing segment; the TS history restarts where the forcing jumps. The
  ! driver never writes the TS solution between steps. With an adjoint
  ! (TS_Adjoint_class) the backward sweep runs at the final time.
  ! TSDriverRerun repeats the run from its initial state with new
  ! parameters, an optimizer's evaluation (-ts_model_rerun_check checks it).
  ! With ADJOINT OPTIMIZATION (TS_Optimizer_class) TAO's evaluations are
  ! reruns, and a last one at the optimum writes the output.
  !
  ! RunToTime replaces PMCBaseRunToTime: changes to that loop do not reach
  ! this one. TSAdapt chooses the steps; of TIMESTEPPER only MAX_STEPS,
  ! MAX_TS_CUTS, the step sizes and their factors apply, the rest are
  ! rejected, as are child and peer couplers.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use PFLOTRAN_Constants_module
  use PMC_Subsurface_class
  use PM_Base_class
  use TS_Model_class
  use TS_Adjoint_class
  use TS_Optimizer_class
  use Timestepper_Base_class, only : TS_CONTINUE, TS_STOP_FAILURE, &
                                     TS_STOP_WALLCLOCK_EXCEEDED
  use Option_module

  implicit none

  private

  ! failed stages in an interval after which BDF continues at order 1
  PetscInt, parameter :: STAGE_FAILURES_TO_LOWER = 3

  abstract interface
    subroutine TSDriverGetModelInterface(pm,model)
      import :: pm_base_type, ts_model_type
      implicit none
      class(pm_base_type) :: pm
      class(ts_model_type), pointer :: model
    end subroutine TSDriverGetModelInterface
  end interface

  ! what PETSc hands the callbacks
  type, public :: ts_driver_ctx_type
    class(ts_model_type), pointer :: model => null()
    type(option_type), pointer :: option => null()
    ! the coupler, for the output of accepted steps (TSDriverStepOutput)
    class(pmc_ts_driver_type), pointer :: pmc => null()
    ! time of the last residual or Jacobian evaluation: the current stage's
    PetscReal :: stage_time = 0.d0
    ! output of an interval's last step, which RunToTime writes with the
    ! waypoint's: any due, and the snap/obs/msbl/cons flags
    PetscBool :: output_due = PETSC_FALSE
    PetscBool :: output_flags(4) = PETSC_FALSE
    ! step size ceiling of the current interval
    PetscReal :: dt_max = MAX_DOUBLE
    PetscReal :: dt_min = 0.d0
    ! step size the next interval resumes from, see TSDriverPostStep
    PetscReal :: unmatched_dt = 0.d0
    ! size of the last accepted step
    PetscReal :: dt_accepted = 0.d0
    PetscInt :: domain_errors = 0
    PetscInt :: events = 0
    ! intervals redone at BDF order 1 after a failure at a higher order
    PetscInt :: order_fallbacks = 0
    ! BDF order to restore at the interval's end after TSDriverPostStep
    ! lowered it to 1; 0 if not lowered
    PetscInt :: order_lowered = 0
    PetscInt :: order_lowerings = 0
    ! failed stages before the interval: SNES failures and domain errors
    PetscInt :: interval_start_snes_failures = 0
    PetscInt :: interval_start_domain_errors = 0
    TS :: ts
    ! fraction of the local error tolerance left to the Newton iteration,
    ! see TSDriverSNESConverged
    PetscReal :: newton_kappa = 0.1d0
    ! no output or progress messages (a rerun)
    PetscBool :: quiet = PETSC_FALSE
    ! -adj_opt_fail_eval (a test): the step that fails, 0 for none
    PetscInt :: fail_at_step = 0
  end type ts_driver_ctx_type

  ! the run's statistics, which reruns must leave as they were
  type :: ts_driver_stats_type
    PetscInt :: steps = 0
    PetscInt :: cuts = 0
    PetscInt :: newton = 0
    PetscInt :: linear = 0
    PetscReal :: solver_time = 0.d0
    PetscLogDouble :: time = 0.d0
  end type ts_driver_stats_type

  type, public, extends(pmc_subsurface_type) :: pmc_ts_driver_type
    procedure(TSDriverGetModelInterface), nopass, pointer :: &
      GetModel => null()
    class(ts_model_type), pointer :: model => null()
    type(ts_driver_ctx_type), pointer :: ctx => null()
    ! the discrete adjoint of the run, if the deck asks for one, and the
    ! calibration of its parameters
    class(ts_adjoint_type), pointer :: adjoint => null()
    type(ts_optimizer_type), pointer :: optimizer => null()
    Vec :: solution
    Vec :: residual
    PetscBool :: ts_initialized = PETSC_FALSE
    PetscBool :: contract_test = PETSC_FALSE
    PetscBool :: adjoint_test = PETSC_FALSE
    PetscBool :: column_check = PETSC_FALSE
    ! the model supplies its Jacobian operator and preconditioner
    PetscBool :: model_jacobian = PETSC_FALSE
    ! the preconditioner the input deck asks for
    PCType :: deck_pc_type = ''
    ! the run's initial time and step size, and its waypoint intervals:
    ! their end times and whether the forcing jumped at their start
    PetscReal :: t0 = 0.d0
    PetscReal :: dt0 = 0.d0
    PetscInt :: num_intervals = 0
    PetscReal, allocatable :: interval_end(:)
    PetscBool, allocatable :: interval_discontinuous(:)
    ! a rerun is in progress, and the interval it is in
    PetscBool :: evaluating = PETSC_FALSE
    PetscInt :: rerun_interval = 0
    ! status of the last adjoint sweep (TS_Adjoint_class)
    PetscInt :: adjoint_status = TS_ADJ_OK
    PetscBool :: rerun_check = PETSC_FALSE
    ! tests of the optimizer's failure paths: the evaluation whose rerun
    ! fails (-adj_opt_fail_eval) or reaches the wall-clock limit
    ! (-adj_opt_wallclock_eval)
    PetscInt :: fail_eval = 0
    PetscInt :: wallclock_eval = 0
    ! the run's stop flag if the optimization or its rerun at the solution
    ! did not complete
    PetscInt :: optimization_stop_flag = TS_CONTINUE
  contains
    procedure, public :: ReadNumericalMethods => PMCTSDriverReadNumMethods
    procedure, public :: SetupSolvers => PMCTSDriverSetupSolvers
    procedure, public :: RunToTime => PMCTSDriverRunToTime
    procedure, public :: Destroy => PMCTSDriverDestroy
  end type pmc_ts_driver_type

  ! an optimizer's evaluation: a rerun at given parameters
  type, extends(ts_evaluator_type) :: ts_driver_evaluator_type
    class(pmc_ts_driver_type), pointer :: pmc => null()
  contains
    procedure, public :: Evaluate => TSDriverEvaluate
  end type ts_driver_evaluator_type

  ! the C wrappers of petsc_fortran_interop.c: PETSc has no Fortran binding
  ! for TSSetEventHandler or SNESSetComputeInitialGuess
  interface
    subroutine pflotran_ts_set_event_handler(ts,nevents,direction, &
                                             terminate,indicator, &
                                             postevent,ctx,ierr)
      use petscts
      import :: ts_driver_ctx_type
      implicit none
      TS :: ts
      PetscInt :: nevents
      PetscInt :: direction(*)
      PetscBool :: terminate(*)
      external :: indicator, postevent
      type(ts_driver_ctx_type) :: ctx
      PetscErrorCode :: ierr
    end subroutine pflotran_ts_set_event_handler
    subroutine pflotran_snes_set_compute_initial_guess(snes,guess,ctx,ierr)
      use petscsnes
      import :: ts_driver_ctx_type
      implicit none
      SNES :: snes
      external :: guess
      type(ts_driver_ctx_type) :: ctx
      PetscErrorCode :: ierr
    end subroutine pflotran_snes_set_compute_initial_guess
  end interface

  public :: PMCTSDriverCreate

contains

! ************************************************************************** !

function PMCTSDriverCreate()
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  class(pmc_ts_driver_type), pointer :: PMCTSDriverCreate

  class(pmc_ts_driver_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  PetscObjectNullify(pmc%solution)
  PetscObjectNullify(pmc%residual)
  PMCTSDriverCreate => pmc

end function PMCTSDriverCreate

! ************************************************************************** !

subroutine PMCTSDriverReadNumMethods(this,input,pm_label)
  !
  ! NUMERICAL_METHODS configure pm_ptr, which need not be first in pm_list:
  ! with flow and transport as one model, pm_list is flow, then transport,
  ! and the block (NUMERICAL_METHODS TRANSPORT) is transport's.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Input_Aux_module

  implicit none

  class(pmc_ts_driver_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: pm_label

  class(pm_base_type), pointer :: pm_list

  pm_list => this%pm_list
  this%pm_list => this%pm_ptr%pm
  call this%pmc_subsurface_type%ReadNumericalMethods(input,pm_label)
  this%pm_list => pm_list

end subroutine PMCTSDriverReadNumMethods

! ************************************************************************** !

subroutine PMCTSDriverSetupSolvers(this)
  !
  ! Creates TS, SNES, KSP and the Jacobian for the coupler's model.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Discretization_module
  use Solver_module
  use Timestepper_TS_class
  use Waypoint_module, only : WaypointListGetFinalTime

  implicit none

  class(pmc_ts_driver_type) :: this

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver
  class(timestepper_TS_type), pointer :: timestepper
  type(discretization_type), pointer :: disc
  TSAdapt :: adapt
  SNESLineSearch :: linesearch
  Vec :: vatol, vrtol
  PetscReal :: atol, rtol
  PetscReal :: zero_pivot_tol
  PetscBool :: flag, linear_shift
  PetscErrorCode :: ierr

  option => this%option
  select type(ts => this%timestepper)
    class is(timestepper_TS_type)
      timestepper => ts
    class default
      call PrintErrMsg(option,'The TS driver requires a TS timestepper.')
  end select
  solver => timestepper%solver
  disc => this%realization%discretization

  if (.not.associated(this%GetModel)) then
    call PrintErrMsg(option,'The TS driver has no equation model.')
  endif
  nullify(this%model)
  call this%GetModel(this%pm_ptr%pm,this%model)
  if (.not.associated(this%model)) then
    call PrintErrMsg(option,'The TS driver found no equation model.')
  endif
  call TSDriverCheckTimestepper(timestepper,this%model%name,option)
  allocate(this%ctx)
  this%ctx%model => this%model
  this%ctx%option => option

  ! RunToTime never runs a child or peer coupler (e.g. waste form, UFD
  ! decay): reject rather than skip it silently
  if (associated(this%child) .or. associated(this%peer)) then
    option%io_buffer = ''
    if (associated(this%child)) option%io_buffer = trim(this%child%name)
    if (associated(this%peer)) option%io_buffer = &
      trim(option%io_buffer) // ' ' // trim(this%peer%name)
    option%io_buffer = trim(this%model%name) // ' does not run process &
      &models coupled to it: ' // trim(adjustl(option%io_buffer)) // '.'
    call PrintErrMsg(option)
  endif

  call PrintMsg(option,'  Beginning setup of ' // trim(this%model%name) // &
                ' TS driver')

  call SolverCreateTS(solver,option%mycomm)
  this%ctx%ts = solver%ts
  call TSSetOptionsPrefix(solver%ts,trim(this%model%prefix),ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  if (Uninitialized(solver%Mpre_mat_type)) solver%Mpre_mat_type = MATBAIJ
  if (Uninitialized(solver%M_mat_type)) &
    solver%M_mat_type = solver%Mpre_mat_type
  if (solver%M_mat_type == MATMFFD .or. solver%Mpre_mat_type == MATMFFD) then
    call PrintErrMsg(option,'The TS driver does not support a matrix-free &
                     &Jacobian (MFFD).')
  endif
  call this%model%CreateJacobian(solver%Mpre,this%model_jacobian)
  if (.not.this%model_jacobian) then
    call DiscretizationCreateMatrix(disc,NDAEDOF,solver%Mpre_mat_type, &
                                    solver%Mpre,PETSC_FALSE,option)
  endif
  call MatSetOptionsPrefix(solver%Mpre,trim(this%model%prefix),ierr)
  CHKERRQ(ierr)
  solver%M = solver%Mpre

  call this%model%CreateVector(this%solution,flag)
  if (.not.flag) then
    call DiscretizationCreateVector(disc,NDAEDOF,this%solution,GLOBAL,option)
  endif
  call VecDuplicate(this%solution,this%residual,ierr);CHKERRQ(ierr)
  this%pm_ptr%pm%solution_vec = this%solution
  this%pm_ptr%pm%residual_vec = this%residual

  call TSSetType(solver%ts,TSBDF,ierr);CHKERRQ(ierr)
  call TSSetEquationType(solver%ts,TS_EQ_DAE_IMPLICIT_INDEX1, &
                         ierr);CHKERRQ(ierr)
  call TSSetIFunction(solver%ts,this%residual,TSDriverIFunction,this%ctx, &
                      ierr);CHKERRQ(ierr)
  call TSSetIJacobian(solver%ts,solver%M,solver%Mpre,TSDriverIJacobian, &
                      this%ctx,ierr);CHKERRQ(ierr)
  call TSSetApplicationContext(solver%ts,this%ctx,ierr);CHKERRQ(ierr)
  call TSSetPostStep(solver%ts,TSDriverPostStep,ierr);CHKERRQ(ierr)
  call TSGetSNES(solver%ts,solver%snes,ierr);CHKERRQ(ierr)
  call SNESGetKSP(solver%snes,solver%ksp,ierr);CHKERRQ(ierr)
  call KSPGetPC(solver%ksp,solver%pc,ierr);CHKERRQ(ierr)

  ! PETSc owns the steps and their retries
  call TSGetAdapt(solver%ts,adapt,ierr);CHKERRQ(ierr)
  call TSAdaptSetType(adapt,TSADAPTBASIC,ierr);CHKERRQ(ierr)
  call TSAdaptSetClip(adapt,timestepper%time_step_reduction_factor, &
                      timestepper%time_step_max_growth_factor, &
                      ierr);CHKERRQ(ierr)
  call TSSetMaxSNESFailures(solver%ts,NEG_ONE_INTEGER,ierr);CHKERRQ(ierr)
  call TSSetMaxStepRejections(solver%ts,timestepper%max_time_step_cuts, &
                              ierr);CHKERRQ(ierr)
  call TSSetMaxSteps(solver%ts,timestepper%max_time_step,ierr);CHKERRQ(ierr)
  call TSSetErrorIfStepFails(solver%ts,PETSC_FALSE,ierr);CHKERRQ(ierr)
  call TSSetExactFinalTime(solver%ts,TS_EXACTFINALTIME_MATCHSTEP, &
                           ierr);CHKERRQ(ierr)

  call TSDriverAddForcingWaypoints(this)

  ! backward Euler, a saved trajectory, dF/dp and recorded forcing; the
  ! observation times are interval ends
  if (associated(this%adjoint)) then
    call this%adjoint%Configure(solver%ts,this%model,this%solution, &
                                WaypointListGetFinalTime(this%waypoint_list))
    call TSDriverAddObservationWaypoints(this)
  endif

  call this%model%DefaultTolerances(atol,rtol)
  call TSSetTolerances(solver%ts,atol,PETSC_NULL_VEC,rtol,PETSC_NULL_VEC, &
                       ierr);CHKERRQ(ierr)
  ! CHKERRA: CHKERRQ would return silently with a half-configured TS
  call TSSetFromOptions(solver%ts,ierr);CHKERRA(ierr)
  if (associated(this%adjoint)) call this%adjoint%CheckMethod(solver%ts)

  ! local error tolerances per entry
  call TSGetTolerances(solver%ts,atol,PETSC_NULL_VEC,rtol,PETSC_NULL_VEC, &
                       ierr);CHKERRQ(ierr)
  call VecDuplicate(this%solution,vatol,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%solution,vrtol,ierr);CHKERRQ(ierr)
  call this%model%SetTolerances(atol,rtol,vatol,vrtol)
  call TSSetTolerances(solver%ts,PETSC_CURRENT_REAL,vatol, &
                       PETSC_CURRENT_REAL,vrtol,ierr);CHKERRQ(ierr)
  call VecDestroy(vatol,ierr);CHKERRQ(ierr)
  call VecDestroy(vrtol,ierr);CHKERRQ(ierr)

  call TSDriverSetupEvents(this,solver%ts)

  ! Newton convergence test of the model
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_newton_kappa',this%ctx%newton_kappa, &
                           flag,ierr);CHKERRQ(ierr)
  call SNESSetConvergenceTest(solver%snes,TSDriverSNESConverged,this%ctx, &
                              PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)
  call pflotran_snes_set_compute_initial_guess(solver%snes, &
                                               TSDriverInitialGuess, &
                                               this%ctx,ierr);CHKERRQ(ierr)
  ! full Newton steps (PreCheck keeps them in the domain): backtracking
  ! stalls at the rate rows' roundoff floor, shift*|N|*eps
  call SNESGetLineSearch(solver%snes,linesearch,ierr);CHKERRQ(ierr)
  call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC,ierr);CHKERRQ(ierr)

  ! The deck's NEWTON_SOLVER ATOL/RTOL/STOL are 2-norm tolerances; here the
  ! rows mix algebraic and rate equations, and an absolute residual
  ! tolerance accepts stages whose algebraic equations have no solution.
  ! Convergence is TSDriverSNESConverged's update test (-<prefix>snes_*tol
  ! on the command line still apply); -ts_model_newton_kappa 0 restores
  ! the deck's tolerances. No divergence test unless the deck sets DTOL,
  ! with either test: a stage starts at a consistent state with a small
  ! residual, which the first update of the algebraic rows raises by
  ! orders of magnitude.
  ! PETSC_UNLIMITED (-3); -1 (PETSC_DETERMINE) is PETSc's default
  if (solver%newton_dtol == PETSC_DEFAULT_REAL) solver%newton_dtol = -3.d0
  if (this%ctx%newton_kappa > 0.d0) then
    solver%newton_atol = 0.d0
    solver%newton_rtol = 0.d0
    solver%newton_stol = 0.d0
    write(option%io_buffer,'(a,es9.2,a)') '  Newton convergence: the &
      &model''s update test (-ts_model_newton_kappa ',this%ctx%newton_kappa, &
      '); NEWTON_SOLVER ATOL/RTOL/STOL do not apply'
    call PrintMsg(option)
  endif
  ! A model operator brings its own (PCSHELL) preconditioner, which applies
  ! the deck's (or -<prefix>pc_type) to the model's condensed system.
  this%deck_pc_type = solver%pc_type
  if (this%model_jacobian) then
    solver%pc_type = PCSHELL
    linear_shift = solver%linear_shift
    zero_pivot_tol = solver%linear_zero_pivot_tol
    solver%linear_shift = PETSC_FALSE
    solver%linear_zero_pivot_tol = UNINITIALIZED_DOUBLE
  endif
  call SolverSetSNESOptions(solver,option)
  if (this%model_jacobian) then
    if (solver%pc_type /= PCSHELL) this%deck_pc_type = solver%pc_type
    solver%pc_type = PCSHELL
    solver%linear_shift = linear_shift
    solver%linear_zero_pivot_tol = zero_pivot_tol
    if (Initialized(zero_pivot_tol)) then
      call PrintWrnMsg(option,'The LINEAR_SOLVER zero pivot tolerance does &
        &not apply to the condensed system of ' // trim(this%model%name) // &
        '; use -tran_schur_pc_factor_zeropivot.')
    endif
    call this%model%SetupLinearSolver(solver%ksp,solver%M, &
                                      this%deck_pc_type,solver%linear_shift)
  endif
  if (solver%snes_type == SNESNEWTONLS) then
    call SNESGetLineSearch(solver%snes,linesearch,ierr);CHKERRQ(ierr)
    call SNESLineSearchSetPreCheck(linesearch,TSDriverPreCheck,this%ctx, &
                                   ierr);CHKERRQ(ierr)
  endif

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_contract_test',flag,ierr);CHKERRQ(ierr)
  this%contract_test = flag
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_adjoint_test',flag,ierr);CHKERRQ(ierr)
  this%adjoint_test = flag
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_column_check',flag,ierr);CHKERRQ(ierr)
  this%column_check = flag
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_rerun_check',flag,ierr);CHKERRQ(ierr)
  this%rerun_check = flag
  if (this%rerun_check .and. .not.associated(this%adjoint)) then
    call PrintErrMsg(option,'-ts_model_rerun_check requires an ADJOINT &
      &block: a rerun replays the forcing recorded for the adjoint.')
  endif
  if (associated(this%optimizer)) then
    if (this%rerun_check) then
      call PrintErrMsg(option,'-ts_model_rerun_check does not run with &
        &ADJOINT OPTIMIZATION.')
    endif
    call TSDriverSetupOptimizer(this)
  endif

  call PrintMsg(option,'  Finished setting up ' // trim(this%model%name) // &
                ' TS driver')

end subroutine PMCTSDriverSetupSolvers

! ************************************************************************** !

subroutine TSDriverAddForcingWaypoints(this)
  !
  ! Makes each knot of the model's forcing a waypoint: a segment holds or
  ! linearly interpolates the forcing of its ends, so it would miss an
  ! interior knot, and PFLOTRAN adds flow condition knots only for rates,
  ! fluxes and SYNC_TIMESTEP_WITH_UPDATE.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Waypoint_module

  implicit none

  class(pmc_ts_driver_type) :: this

  type(waypoint_type), pointer :: waypoint
  PetscReal, allocatable :: times(:)
  PetscReal :: final_time
  PetscInt :: i

  final_time = WaypointListGetFinalTime(this%waypoint_list)
  call this%model%ForcingTimes(final_time,times)
  do i = 1, size(times)
    if (times(i) <= 0.d0 .or. times(i) > final_time) cycle
    waypoint => WaypointCreate()
    waypoint%time = times(i)
    waypoint%update_conditions = PETSC_TRUE
    call WaypointInsertInList(waypoint,this%waypoint_list,this%option)
  enddo
  deallocate(times)

end subroutine TSDriverAddForcingWaypoints

! ************************************************************************** !

subroutine TSDriverAddObservationWaypoints(this)
  !
  ! Makes each observation time of the adjoint a waypoint, so that the run
  ! has a step state there.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  use Waypoint_module

  implicit none

  class(pmc_ts_driver_type) :: this

  type(waypoint_type), pointer :: waypoint
  PetscReal, allocatable :: times(:)
  PetscInt :: i

  call this%adjoint%ObservationTimes(times)
  do i = 1, size(times)
    waypoint => WaypointCreate()
    waypoint%time = times(i)
    call WaypointInsertInList(waypoint,this%waypoint_list,this%option)
  enddo
  deallocate(times)

end subroutine TSDriverAddObservationWaypoints

! ************************************************************************** !

recursive subroutine PMCTSDriverRunToTime(this,sync_time,stop_flag)
  !
  ! Integrates to sync_time, one TSSolve per waypoint interval, keeping
  ! PETSc's solution and history throughout.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Solver_module
  use Timestepper_Base_class
  use Timestepper_TS_class
  use Waypoint_module

  implicit none

  class(pmc_ts_driver_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag

  class(timestepper_TS_type), pointer :: timestepper
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  type(waypoint_type), pointer :: waypoint
  TSAdapt :: adapt
  TSConvergedReason :: reason
  PetscReal :: time, end_time, previous_time, next_dt, tolerance
  PetscReal :: dt_max
  PetscInt :: steps_before, steps_after, newton_before, newton_after
  PetscInt :: linear_before, linear_after
  PetscInt :: rejected_before, rejected_after
  PetscInt :: incoming_stop_flag
  PetscLogDouble :: start_clock, end_clock
  PetscBool :: discontinuous, at_waypoint, redo
  PetscBool :: print_snap_output, print_obs_output
  PetscBool :: print_msbl_output, print_cons_output
  PetscBool :: wallclock_stop
  PetscInt :: order_saved
  PetscInt :: counts(2)
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE .or. &
      stop_flag == TS_STOP_WALLCLOCK_EXCEEDED) return
  incoming_stop_flag = stop_flag
  select type(ts => this%timestepper)
    class is(timestepper_TS_type)
      timestepper => ts
  end select
  solver => timestepper%solver
  option => this%option
  tolerance = 64.d0*epsilon(1.d0)*max(1.d0,abs(sync_time))
  if (timestepper%target_time >= sync_time-tolerance) return
  ! for the output of accepted steps
  this%ctx%pmc => this

  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif
  call PetscTime(start_clock,ierr);CHKERRQ(ierr)
  call this%GetAuxData()
  call TSGetAdapt(solver%ts,adapt,ierr);CHKERRQ(ierr)

  time = timestepper%target_time
  waypoint => this%waypoint_list%first
  do while (associated(waypoint))
    if (waypoint%time > time+tolerance) exit
    waypoint => waypoint%next
  enddo

  stop_flag = TS_CONTINUE
  order_saved = 0
  redo = PETSC_FALSE
  do while (time < sync_time-tolerance)
    end_time = sync_time
    dt_max = timestepper%dt_max
    if (associated(waypoint)) then
      end_time = min(end_time,waypoint%time)
      dt_max = waypoint%dt_max
    endif

    ! forcing segment of this interval; a redo keeps the interval's, a
    ! rerun replays the recorded one of the same interval
    if (redo) then
      discontinuous = PETSC_FALSE
      redo = PETSC_FALSE
    else if (this%evaluating) then
      this%rerun_interval = this%rerun_interval + 1
      if (.not.TSDriverIntervalRecorded(this,this%rerun_interval, &
                                        end_time)) then
        stop_flag = TS_STOP_FAILURE
        exit
      endif
      discontinuous = this%interval_discontinuous(this%rerun_interval)
    else
      call this%model%BeginSegment(time,end_time,discontinuous)
      call TSDriverRecordInterval(this,end_time,discontinuous)
    endif
    if (.not.this%ts_initialized) then
      call this%model%InitialState(time,this%solution)
      this%t0 = time
      this%dt0 = timestepper%dt
      call TSSetTime(solver%ts,time,ierr);CHKERRQ(ierr)
      call TSSetTimeStep(solver%ts,timestepper%dt,ierr);CHKERRQ(ierr)
      call TSSetStepNumber(solver%ts,timestepper%steps,ierr);CHKERRQ(ierr)
      call TSSetSolution(solver%ts,this%solution,ierr);CHKERRQ(ierr)
      if (this%contract_test .or. this%column_check .or. &
          this%adjoint_test) then
        call TSDriverModelChecks(this,time,timestepper%dt,end_time)
      endif
      if (associated(this%adjoint)) &
        call this%adjoint%RecordInitialState(time,this%solution)
      this%ts_initialized = PETSC_TRUE
    else if (discontinuous) then
      call TSRestartStep(solver%ts,ierr);CHKERRQ(ierr)
    endif

    this%ctx%dt_max = dt_max
    this%ctx%dt_min = timestepper%dt_min
    call TSDriverSetStepLimits(adapt,timestepper%dt_min,this%ctx%dt_max)
    call TSGetTimeStep(solver%ts,next_dt,ierr);CHKERRQ(ierr)
    next_dt = min(next_dt,this%ctx%dt_max)
    call TSSetTimeStep(solver%ts,next_dt,ierr);CHKERRQ(ierr)
    this%ctx%unmatched_dt = next_dt
    call TSSetMaxTime(solver%ts,end_time,ierr);CHKERRQ(ierr)

    call TSGetStepNumber(solver%ts,steps_before,ierr);CHKERRQ(ierr)
    call TSGetSNESIterations(solver%ts,newton_before,ierr);CHKERRQ(ierr)
    call TSGetKSPIterations(solver%ts,linear_before,ierr);CHKERRQ(ierr)
    call TSGetStepRejections(solver%ts,rejected_before,ierr);CHKERRQ(ierr)
    call TSGetSNESFailures(solver%ts,this%ctx%interval_start_snes_failures, &
                           ierr);CHKERRQ(ierr)
    this%ctx%interval_start_domain_errors = this%ctx%domain_errors
    this%ctx%output_due = PETSC_FALSE
    this%ctx%output_flags = PETSC_FALSE

    ! failure unless TSSolve returns with a positive reason
    stop_flag = TS_STOP_FAILURE
    call TSSolve(solver%ts,this%solution,ierr);CHKERRQ(ierr)
    call TSGetConvergedReason(solver%ts,reason,ierr);CHKERRQ(ierr)
    call TSGetTime(solver%ts,time,ierr);CHKERRQ(ierr)
    call TSGetPrevTime(solver%ts,previous_time,ierr);CHKERRQ(ierr)
    call TSGetTimeStep(solver%ts,next_dt,ierr);CHKERRQ(ierr)
    call TSGetStepNumber(solver%ts,steps_after,ierr);CHKERRQ(ierr)
    call TSGetSNESIterations(solver%ts,newton_after,ierr);CHKERRQ(ierr)
    call TSGetKSPIterations(solver%ts,linear_after,ierr);CHKERRQ(ierr)
    call TSGetStepRejections(solver%ts,rejected_after,ierr);CHKERRQ(ierr)
    ! the adjoint's trajectory must keep the interval's last step
    if (associated(this%adjoint) .and. reason%v > 0) &
      call this%adjoint%EndInterval(solver%ts,this%model,this%solution)

    ! resume from the step size used before MATCHSTEP shortened the last
    ! step(s), unless TSAdapt asked for a reduction
    if (reason%v >= 0 .and. next_dt >= time-previous_time) then
      next_dt = max(next_dt,this%ctx%unmatched_dt)
      call TSSetTimeStep(solver%ts,next_dt,ierr);CHKERRQ(ierr)
    endif

    timestepper%steps = timestepper%steps + steps_after-steps_before
    timestepper%target_time = time
    timestepper%dt = next_dt
    timestepper%cumulative_time_step_cuts = &
      timestepper%cumulative_time_step_cuts + rejected_after-rejected_before
    timestepper%cumulative_newton_iterations = &
      timestepper%cumulative_newton_iterations + newton_after-newton_before
    timestepper%cumulative_linear_iterations = &
      timestepper%cumulative_linear_iterations + linear_after-linear_before
    option%time = time
    option%dt = time-previous_time
    if (this%model%integrates_flow) option%flow_dt = option%dt
    if (this%model%integrates_transport) option%tran_dt = option%dt

    if (reason%v < 0) then
      ! BDF order > 1 is not positivity preserving: a vanishing species can
      ! make every stage inadmissible. Redo the rest of the interval at
      ! order 1 (positivity preserving for decay), then restore the order.
      if (TSDriverLowerOrder(solver%ts,order_saved)) then
        this%ctx%order_fallbacks = this%ctx%order_fallbacks + 1
        write(option%io_buffer,'(a,es12.5,a)') ' ' // &
          trim(this%model%name) // ' TS: step failure at t = ',time, &
          '; redoing the interval at BDF order 1'
        call PrintMsg(option)
        call TSSetConvergedReason(solver%ts,TS_CONVERGED_ITERATING, &
                                  ierr);CHKERRQ(ierr)
        call TSRestartStep(solver%ts,ierr);CHKERRQ(ierr)
        ! the last redo: from the last accepted step size, or dt_min before
        ! any (the failed sizes, already cut by TSAdapt, may fail again)
        call TSSetTimeStep(solver%ts,max(this%ctx%dt_accepted, &
                                         timestepper%dt_min),ierr)
        CHKERRQ(ierr)
        stop_flag = TS_CONTINUE
        redo = PETSC_TRUE
        cycle
      endif
      write(option%io_buffer,'(a,i4,a,es12.5)') &
        trim(this%model%name) // ' TS failed (reason ',reason%v, &
        ') at t = ',time
      call PrintMsg(option)
      exit
    endif
    if (order_saved > 0) then
      call TSBDFSetOrder(solver%ts,order_saved,ierr);CHKERRQ(ierr)
      order_saved = 0
    endif
    if (this%ctx%order_lowered > 0) then
      call TSBDFSetOrder(solver%ts,this%ctx%order_lowered,ierr);CHKERRQ(ierr)
      this%ctx%order_lowered = 0
    endif
    stop_flag = TS_CONTINUE
    if (.not.this%ctx%quiet) then
      write(option%io_buffer, &
        '("== ",a," TS: ",i8," steps, t=",es12.5," next dt=",es12.5)') &
        trim(this%model%name),steps_after-steps_before,time,next_dt
      call PrintMsg(option)
    endif

    ! output at the waypoint, together with that of the interval's last
    ! step (see TSDriverStepOutput)
    at_waypoint = PETSC_FALSE
    if (associated(waypoint)) at_waypoint = (time >= waypoint%time-tolerance)
    if ((at_waypoint .or. this%ctx%output_due) .and. &
        .not.this%ctx%quiet) then
      call this%model%Publish(time,this%solution)
      if (associated(this%Output)) then
        ! Output() clears its flag arguments; pass copies
        print_snap_output = this%ctx%output_flags(1)
        print_obs_output = this%ctx%output_flags(2)
        print_msbl_output = this%ctx%output_flags(3)
        print_cons_output = this%ctx%output_flags(4)
        if (at_waypoint) then
          print_snap_output = print_snap_output .or. &
                              waypoint%print_snap_output
          print_obs_output = print_obs_output .or. waypoint%print_obs_output
          print_msbl_output = print_msbl_output .or. &
                              waypoint%print_msbl_output
          print_cons_output = print_cons_output .or. &
                              waypoint%print_cons_output
        endif
        call this%Output(this%pm_list%realization_base,print_snap_output, &
                         print_obs_output,print_msbl_output, &
                         print_cons_output)
      endif
    endif
    this%ctx%output_due = PETSC_FALSE
    this%ctx%output_flags = PETSC_FALSE
    if (at_waypoint) waypoint => waypoint%next
    if (this%is_master .and. option%wallclock_stop_flag) then
      wallclock_stop = timestepper%WallClockStop(option)
      call MPI_Bcast(wallclock_stop,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     ZERO_INTEGER_MPI,option%mycomm,ierr);CHKERRQ(ierr)
      if (wallclock_stop) then
        stop_flag = TS_STOP_WALLCLOCK_EXCEEDED
        exit
      endif
    endif
    if (time < end_time-tolerance) then
      ! a positive reason can mean the step limit was reached
      stop_flag = TS_STOP_MAX_TIME_STEP
      exit
    endif
  enddo

  if (stop_flag == TS_CONTINUE .and. .not.associated(waypoint)) &
    stop_flag = TS_STOP_END_SIMULATION

  ! the backward sweep of a complete run; it leaves the final state
  if (associated(this%adjoint)) then
    if (stop_flag == TS_STOP_END_SIMULATION) then
      call this%adjoint%Solve(solver%ts,this%model,this%solution, &
                              this%ctx%events,this%adjoint_status)
      if (this%rerun_check .and. .not.this%evaluating) &
        call TSDriverRerunCheck(this)
      if (associated(this%optimizer) .and. .not.this%evaluating) then
        call TSDriverOptimize(this)
        if (this%optimization_stop_flag /= TS_CONTINUE) &
          stop_flag = this%optimization_stop_flag
      endif
    else if (stop_flag /= TS_CONTINUE .and. this%evaluating) then
      this%adjoint_status = TS_ADJ_FORWARD_FAILED
    else if (stop_flag /= TS_CONTINUE) then
      write(option%io_buffer,'(a,i0,a)') 'ADJOINT: forward run stopped &
        &before completion (stop_flag = ',stop_flag,'); adjoint solve &
        &skipped.'
      call PrintMsg(option)
      ! no gradients: not a success (a failure already exits with
      ! EXIT_FAILURE, see SimSubsurfFinalizeRun)
      if (stop_flag /= TS_STOP_FAILURE) &
        option%driver%exit_code = EXIT_ADJOINT_SKIPPED
    endif
  endif

  ! publish for the parent's output at sync_time
  call this%model%Publish(time,this%solution)
  if (stop_flag /= TS_CONTINUE .and. .not.this%ctx%quiet) then
    ! events and domain errors are counted per process
    counts(1) = this%ctx%events
    counts(2) = this%ctx%domain_errors
    call MPI_Allreduce(MPI_IN_PLACE,counts,TWO_INTEGER_MPI,MPIU_INTEGER, &
                       MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'(a,4(a,i8))') ' ' // trim(this%model%name), &
      ' TS summary: events',counts(1),', BDF order fallbacks', &
      this%ctx%order_fallbacks,', order lowerings', &
      this%ctx%order_lowerings,', domain errors',counts(2)
    call PrintMsg(option)
    call this%model%PrintSummary()
  endif
  stop_flag = max(stop_flag,incoming_stop_flag)
  call this%SetAuxData()
  call PetscTime(end_clock,ierr);CHKERRQ(ierr)
  this%cumulative_time = this%cumulative_time + end_clock-start_clock
  timestepper%cumulative_solver_time = &
    timestepper%cumulative_solver_time + end_clock-start_clock
  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif

end subroutine PMCTSDriverRunToTime

! ************************************************************************** !

subroutine TSDriverRecordInterval(this,t_end,discontinuous)
  !
  ! Records a waypoint interval of the run, for a rerun.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(pmc_ts_driver_type) :: this
  PetscReal :: t_end
  PetscBool :: discontinuous

  PetscReal, allocatable :: t_grown(:)
  PetscBool, allocatable :: d_grown(:)
  PetscInt :: n

  n = this%num_intervals
  if (.not.allocated(this%interval_end)) then
    allocate(this%interval_end(16))
    allocate(this%interval_discontinuous(16))
  else if (n == size(this%interval_end)) then
    allocate(t_grown(2*n))
    allocate(d_grown(2*n))
    t_grown(1:n) = this%interval_end(1:n)
    d_grown(1:n) = this%interval_discontinuous(1:n)
    call move_alloc(t_grown,this%interval_end)
    call move_alloc(d_grown,this%interval_discontinuous)
  endif
  this%num_intervals = n + 1
  this%interval_end(n+1) = t_end
  this%interval_discontinuous(n+1) = discontinuous

end subroutine TSDriverRecordInterval

! ************************************************************************** !

function TSDriverIntervalRecorded(this,k,t_end)
  !
  ! Whether a rerun's k-th waypoint interval, ending at t_end, is the
  ! run's; reports a mismatch.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(pmc_ts_driver_type) :: this
  PetscInt :: k
  PetscReal :: t_end
  PetscBool :: TSDriverIntervalRecorded

  TSDriverIntervalRecorded = PETSC_FALSE
  if (k <= this%num_intervals) then
    TSDriverIntervalRecorded = &
      abs(this%interval_end(k)-t_end) <= &
      64.d0*epsilon(1.d0)*max(1.d0,abs(t_end))
  endif
  if (TSDriverIntervalRecorded) return
  write(this%option%io_buffer,'(a,i0,a,es12.5,a,i0,a)') ' ' // &
    trim(this%model%name) // ' TS rerun: waypoint interval ',k, &
    ' ends at t = ',t_end,', not as in the run (',this%num_intervals, &
    ' intervals)'
  call PrintMsg(this%option)

end function TSDriverIntervalRecorded

! ************************************************************************** !

subroutine TSDriverRerun(this,status,quiet)
  !
  ! Integrates the run again from its initial state with the current
  ! parameters, replaying the forcing recorded for the adjoint over the
  ! same waypoint intervals, then runs the adjoint sweep; quiet: outputs
  ! nothing and prints only a failure. InitialState and BeginSegment are
  ! not called: the
  ! first reads the output fields, which Publish has overwritten, and
  ! conditions only move forward in time. Requires a completed run with
  ! an adjoint; leaves the TS and model at the rerun's final state, and
  ! adds to the run's statistics (callers restore them, TSDriverSaveStats).
  ! status: that of the sweep (TS_Adjoint_class), TS_ADJ_FORWARD_FAILED or
  ! TS_ADJ_WALLCLOCK.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use Timestepper_Base_class
  use Timestepper_TS_class
  use Waypoint_module, only : WaypointListGetFinalTime

  implicit none

  class(pmc_ts_driver_type) :: this
  PetscInt :: status
  PetscBool :: quiet

  class(timestepper_TS_type), pointer :: timestepper
  type(ts_driver_ctx_type), pointer :: ctx
  TS :: ts
  PetscReal :: final_time
  PetscInt :: stop_flag
  PetscBool :: ctx_quiet, adjoint_quiet
  PetscErrorCode :: ierr

  select type(stepper => this%timestepper)
    class is(timestepper_TS_type)
      timestepper => stepper
  end select
  ts = timestepper%solver%ts
  ctx => this%ctx
  final_time = WaypointListGetFinalTime(this%waypoint_list)

  call VecCopy(this%adjoint%Y0,this%solution,ierr);CHKERRQ(ierr)
  call TSSetTime(ts,this%t0,ierr);CHKERRQ(ierr)
  call TSSetStepNumber(ts,ZERO_INTEGER,ierr);CHKERRQ(ierr)
  call TSSetTimeStep(ts,this%dt0,ierr);CHKERRQ(ierr)
  call TSSetConvergedReason(ts,TS_CONVERGED_ITERATING,ierr);CHKERRQ(ierr)
  call TSRestartStep(ts,ierr);CHKERRQ(ierr)
  call this%model%ResetToInitial()
  call this%model%BeginReplay()
  call this%adjoint%BeginRerun(ts)

  timestepper%target_time = this%t0
  timestepper%dt = this%dt0
  this%option%time = this%t0
  ctx%stage_time = this%t0
  ctx%output_due = PETSC_FALSE
  ctx%output_flags = PETSC_FALSE
  ctx%unmatched_dt = 0.d0
  ctx%dt_accepted = 0.d0
  ctx%domain_errors = 0
  ctx%events = 0
  ctx%order_fallbacks = 0
  ctx%order_lowered = 0
  ctx%order_lowerings = 0
  ctx%interval_start_snes_failures = 0
  ctx%interval_start_domain_errors = 0
  this%rerun_interval = 0
  this%adjoint_status = TS_ADJ_FORWARD_FAILED

  ctx_quiet = ctx%quiet
  adjoint_quiet = this%adjoint%quiet
  this%evaluating = PETSC_TRUE
  ctx%quiet = quiet
  this%adjoint%quiet = quiet
  stop_flag = TS_CONTINUE
  call this%RunToTime(final_time,stop_flag)
  this%evaluating = PETSC_FALSE
  ctx%quiet = ctx_quiet
  this%adjoint%quiet = adjoint_quiet
  call this%model%EndReplay()

  status = this%adjoint_status
  if (stop_flag == TS_STOP_WALLCLOCK_EXCEEDED) then
    status = TS_ADJ_WALLCLOCK
  else if (stop_flag /= TS_STOP_END_SIMULATION .or. &
           this%rerun_interval /= this%num_intervals) then
    status = TS_ADJ_FORWARD_FAILED
  endif

end subroutine TSDriverRerun

! ************************************************************************** !

subroutine TSDriverRerunCheck(this)
  !
  ! -ts_model_rerun_check, after the run's adjoint sweep: reruns at the
  ! deck's parameters p0, at 1.1 p0 and at p0 again. Both reruns at p0
  ! must reproduce the run's costs, gradients, adjoint initial state and
  ! final state bit for bit; the costs and gradients at 1.1 p0 are
  ! printed, for comparison with a run of a deck at those values (a regime
  ! transition there is allowed). Restores the run's final state,
  ! statistics (the times include the reruns) and kept results.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use Timestepper_TS_class

  implicit none

  class(pmc_ts_driver_type) :: this

  type(option_type), pointer :: option
  class(timestepper_TS_type), pointer :: timestepper
  class(ts_adjoint_type), pointer :: adjoint
  TS :: ts
  Vec :: Y_ref, work
  Vec, allocatable :: lambda_ref(:)
  PetscReal, allocatable :: p0(:), cost0(:), grad0(:,:), cost(:)
  PetscReal, allocatable :: grad(:,:)
  type(ts_driver_stats_type) :: stats
  PetscReal :: final_time, final_dt, target_time, dt, time, option_dt
  PetscReal :: flow_dt, tran_dt, diff(4), value
  PetscInt :: final_step, events, errors
  PetscInt :: fallbacks, lowerings, status, pass, i, j, ncosts, nparams
  PetscBool :: passed
  character(len=*), parameter :: label(3) = &
    [character(len=10) :: 'p0', '1.1 p0', 'p0 again']
  PetscErrorCode :: ierr

  option => this%option
  adjoint => this%adjoint
  select type(stepper => this%timestepper)
    class is(timestepper_TS_type)
      timestepper => stepper
  end select
  ts = timestepper%solver%ts
  if (this%adjoint_status /= TS_ADJ_OK) then
    call PrintErrMsg(option,'-ts_model_rerun_check requires a run whose &
      &adjoint sweep completes.')
  endif
  call PrintMsg(option,' TS model rerun check: rerunning at p0, 1.1 p0 and &
    &p0')

  ncosts = adjoint%ncosts
  nparams = adjoint%nparams
  allocate(p0(nparams),cost0(ncosts),cost(ncosts))
  allocate(grad0(nparams,ncosts),grad(nparams,ncosts))
  call adjoint%GetParameters(p0)
  call adjoint%GetCostValues(cost0)
  call adjoint%GetCostGradients(grad0)
  call VecDuplicate(this%solution,Y_ref,ierr);CHKERRQ(ierr)
  call VecCopy(this%solution,Y_ref,ierr);CHKERRQ(ierr)
  call VecDuplicate(this%solution,work,ierr);CHKERRQ(ierr)
  allocate(lambda_ref(ncosts))
  do i = 1, ncosts
    call VecDuplicate(adjoint%lambda(i),lambda_ref(i),ierr);CHKERRQ(ierr)
    call VecCopy(adjoint%lambda(i),lambda_ref(i),ierr);CHKERRQ(ierr)
  enddo
  ! the run's state and statistics
  call TSGetTime(ts,final_time,ierr);CHKERRQ(ierr)
  call TSGetStepNumber(ts,final_step,ierr);CHKERRQ(ierr)
  call TSGetTimeStep(ts,final_dt,ierr);CHKERRQ(ierr)
  target_time = timestepper%target_time
  dt = timestepper%dt
  call TSDriverSaveStats(this,stats)
  time = option%time
  option_dt = option%dt
  flow_dt = option%flow_dt
  tran_dt = option%tran_dt
  events = this%ctx%events
  errors = this%ctx%domain_errors
  fallbacks = this%ctx%order_fallbacks
  lowerings = this%ctx%order_lowerings

  passed = PETSC_TRUE
  do pass = 1, 3
    if (pass == 2) then
      call adjoint%SetParameters(this%model,1.1d0*p0)
    else
      call adjoint%SetParameters(this%model,p0)
    endif
    call TSDriverRerun(this,status,PETSC_TRUE)
    if (status /= TS_ADJ_OK) then
      write(option%io_buffer,'(a,a,a,i0)') ' TS model rerun check at ', &
        trim(label(pass)),': no gradient, status ',status
      call PrintMsg(option)
      ! a regime transition at 1.1 p0 skips its sweep; the rerun at p0
      ! after it checks that the run's regimes are restored
      if (pass /= 2 .or. status /= TS_ADJ_SKIPPED_EVENTS) &
        passed = PETSC_FALSE
      cycle
    endif
    call adjoint%GetCostValues(cost)
    call adjoint%GetCostGradients(grad)
    if (pass == 2) then
      do j = 1, nparams
        write(option%io_buffer,'(a,i0,a,es24.16)') &
          ' TS model rerun check at 1.1 p0: p_',j,' = ',1.1d0*p0(j)
        call PrintMsg(option)
      enddo
      do i = 1, ncosts
        write(option%io_buffer,'(a,i0,a,es24.16)') &
          ' TS model rerun check at 1.1 p0: J_',i,' = ',cost(i)
        call PrintMsg(option)
        do j = 1, nparams
          write(option%io_buffer,'(a,i0,a,i0,a,es24.16)') &
            ' TS model rerun check at 1.1 p0: dJ_',i,'/dp_',j,' = ',grad(j,i)
          call PrintMsg(option)
        enddo
      enddo
      cycle
    endif
    diff(1) = maxval(abs(cost-cost0))
    diff(2) = maxval(abs(grad-grad0))
    diff(3) = 0.d0
    do i = 1, ncosts
      call VecWAXPY(work,-1.d0,lambda_ref(i),adjoint%lambda(i), &
                    ierr);CHKERRQ(ierr)
      call VecNorm(work,NORM_INFINITY,value,ierr);CHKERRQ(ierr)
      diff(3) = max(diff(3),value)
    enddo
    call VecWAXPY(work,-1.d0,Y_ref,this%solution,ierr);CHKERRQ(ierr)
    call VecNorm(work,NORM_INFINITY,diff(4),ierr);CHKERRQ(ierr)
    write(option%io_buffer,'(a,a,a,4(a,es10.3))') ' TS model rerun check &
      &at ',trim(label(pass)),':',' max |dJ| ',diff(1),', max |d(dJ/dp)| ', &
      diff(2),', max |d lambda(t0)| ',diff(3),', max |dY(T)| ',diff(4)
    call PrintMsg(option)
    if (any(diff /= 0.d0)) passed = PETSC_FALSE
  enddo

  ! back to the run
  call adjoint%SetParameters(this%model,p0)
  call VecCopy(Y_ref,this%solution,ierr);CHKERRQ(ierr)
  call TSSetTime(ts,final_time,ierr);CHKERRQ(ierr)
  call TSSetStepNumber(ts,final_step,ierr);CHKERRQ(ierr)
  call TSSetTimeStep(ts,final_dt,ierr);CHKERRQ(ierr)
  timestepper%target_time = target_time
  timestepper%dt = dt
  call TSDriverRestoreStats(this,stats)
  option%time = time
  option%dt = option_dt
  option%flow_dt = flow_dt
  option%tran_dt = tran_dt
  this%ctx%events = events
  this%ctx%domain_errors = errors
  this%ctx%order_fallbacks = fallbacks
  this%ctx%order_lowerings = lowerings
  adjoint%cost_value = cost0
  do i = 1, ncosts
    call VecDestroy(lambda_ref(i),ierr);CHKERRQ(ierr)
  enddo
  call VecDestroy(Y_ref,ierr);CHKERRQ(ierr)
  call VecDestroy(work,ierr);CHKERRQ(ierr)

  if (.not.passed) then
    call PrintErrMsg(option,'TS model rerun check FAILED: a rerun does not &
      &reproduce the run.')
  endif
  call PrintMsg(option,' TS model rerun check: PASSED')

end subroutine TSDriverRerunCheck

! ************************************************************************** !

subroutine TSDriverSaveStats(this,stats)
  !
  ! The run's statistics, before reruns: steps, cuts, Newton and linear
  ! iterations and times (the coupler's RunToTime counts a nested rerun's
  ! time in its own).
  !
  ! Author: Hansol Suh
  ! Date: 10/11/26
  !
  use Timestepper_TS_class

  implicit none

  class(pmc_ts_driver_type) :: this
  type(ts_driver_stats_type) :: stats

  stats%steps = this%timestepper%steps
  stats%cuts = this%timestepper%cumulative_time_step_cuts
  select type(stepper => this%timestepper)
    class is(timestepper_TS_type)
      stats%newton = stepper%cumulative_newton_iterations
      stats%linear = stepper%cumulative_linear_iterations
  end select
  stats%solver_time = this%timestepper%cumulative_solver_time
  stats%time = this%cumulative_time

end subroutine TSDriverSaveStats

! ************************************************************************** !

subroutine TSDriverRestoreStats(this,stats)
  !
  ! Author: Hansol Suh
  ! Date: 10/11/26
  !
  use Timestepper_TS_class

  implicit none

  class(pmc_ts_driver_type) :: this
  type(ts_driver_stats_type) :: stats

  this%timestepper%steps = stats%steps
  this%timestepper%cumulative_time_step_cuts = stats%cuts
  select type(stepper => this%timestepper)
    class is(timestepper_TS_type)
      stepper%cumulative_newton_iterations = stats%newton
      stepper%cumulative_linear_iterations = stats%linear
  end select
  this%timestepper%cumulative_solver_time = stats%solver_time
  this%cumulative_time = stats%time

end subroutine TSDriverRestoreStats

! ************************************************************************** !

subroutine TSDriverSetupOptimizer(this)
  !
  ! The optimizer's parameters, from the adjoint; the run is quiet: its
  ! output is that of the last rerun, at the optimum.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(pmc_ts_driver_type) :: this

  character(len=MAXSTRINGLENGTH), allocatable :: name(:)
  PetscReal, allocatable :: p0(:)
  PetscBool, allocatable :: positive(:), log_scale(:)
  PetscBool :: flag
  PetscInt :: j, n
  PetscErrorCode :: ierr

  if (this%adjoint%iobjective == 0) then
    call PrintErrMsg(this%option,'ADJOINT OPTIMIZATION requires &
      &OBSERVATIONS: it minimizes their misfit.')
  endif
  n = this%adjoint%nparams
  allocate(name(n),p0(n),positive(n),log_scale(n))
  call this%adjoint%GetParameters(p0)
  do j = 1, n
    call this%adjoint%ParameterInfo(j,name(j),positive(j),log_scale(j))
  enddo
  call this%optimizer%Setup(this%option,name,positive,log_scale,p0)
  this%ctx%quiet = PETSC_TRUE
  this%adjoint%quiet = PETSC_TRUE
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                          '-adj_opt_fail_eval',this%fail_eval,flag, &
                          ierr);CHKERRQ(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                          '-adj_opt_wallclock_eval',this%wallclock_eval, &
                          flag,ierr);CHKERRQ(ierr)

end subroutine TSDriverSetupOptimizer

! ************************************************************************** !

subroutine TSDriverObjective(this,f,g)
  !
  ! The optimizer's objective, the observations' misfit of the last sweep,
  ! and its gradient.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(pmc_ts_driver_type) :: this
  PetscReal :: f
  PetscReal :: g(:)

  PetscReal :: cost(this%adjoint%ncosts)
  PetscReal :: grad(this%adjoint%nparams,this%adjoint%ncosts)

  call this%adjoint%GetCostValues(cost)
  call this%adjoint%GetCostGradients(grad)
  f = cost(this%adjoint%iobjective)
  g = grad(:,this%adjoint%iobjective)

end subroutine TSDriverObjective

! ************************************************************************** !

subroutine TSDriverEvaluate(this,p,f,g,status)
  !
  ! The objective and its gradient at parameters p: a rerun.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(ts_driver_evaluator_type) :: this
  PetscReal :: p(:)
  PetscReal :: f
  PetscReal :: g(:)
  PetscInt :: status

  class(pmc_ts_driver_type), pointer :: pmc
  PetscInt :: evaluation

  pmc => this%pmc
  ! the tests of the failure paths: a step failure in the middle of the
  ! run, or the wall-clock check after the first interval
  evaluation = pmc%optimizer%num_evaluations + 1
  if (evaluation == pmc%fail_eval) pmc%ctx%fail_at_step = 2
  if (evaluation == pmc%wallclock_eval) then
    pmc%option%wallclock_stop_flag = PETSC_TRUE
    pmc%option%wallclock_stop_time = 0.d0
  endif
  call pmc%adjoint%SetParameters(pmc%model,p)
  call TSDriverRerun(pmc,status,PETSC_TRUE)
  pmc%ctx%fail_at_step = 0
  f = 0.d0
  g = 0.d0
  select case(status)
    case(TS_ADJ_OK)
      call TSDriverObjective(pmc,f,g)
      status = TS_OPT_EVAL_OK
    case(TS_ADJ_WALLCLOCK)
      status = TS_OPT_EVAL_STOP
  end select

end subroutine TSDriverEvaluate

! ************************************************************************** !

subroutine TSDriverOptimize(this)
  !
  ! ADJOINT OPTIMIZATION, after the run's sweep: TAO from the run's
  ! objective and gradient, then a rerun at its solution with output (and
  ! the adjoint's report of it). The run must have a gradient. The
  ! statistics stay the run's; a stop at the wall-clock limit or a failed
  ! rerun at the solution sets optimization_stop_flag.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(pmc_ts_driver_type), target :: this

  type(ts_driver_evaluator_type), target :: evaluator
  type(ts_driver_stats_type) :: stats
  PetscReal :: f
  PetscReal :: g(this%adjoint%nparams)
  PetscInt :: status

  if (this%adjoint_status /= TS_ADJ_OK) then
    write(this%option%io_buffer,'(a,i0,a)') 'ADJOINT OPTIMIZATION: the &
      &run at the initial parameters has no gradient (adjoint status ', &
      this%adjoint_status,'; 1: mineral regime transitions, 2: an &
      &observation not recorded, 3: incomplete backward solve).'
    call PrintErrMsg(this%option)
  endif
  call TSDriverObjective(this,f,g)
  call this%optimizer%SetInitialEvaluation(f,g)
  call TSDriverSaveStats(this,stats)
  evaluator%pmc => this
  call this%optimizer%Solve(evaluator)
  ! at the wall-clock limit: no rerun at the solution (the fields are those
  ! of the interrupted evaluation)
  if (this%optimizer%stopped) then
    this%optimization_stop_flag = TS_STOP_WALLCLOCK_EXCEEDED
    call TSDriverRestoreStats(this,stats)
    return
  endif

  call PrintMsg(this%option,' OPTIMIZATION: rerunning at the solution &
    &with output')
  call this%adjoint%SetParameters(this%model,this%optimizer%p_opt)
  call TSDriverRerun(this,status,PETSC_FALSE)
  call TSDriverRestoreStats(this,stats)
  select case(status)
    case(TS_ADJ_OK)
    case(TS_ADJ_WALLCLOCK)
      ! the output stops where the rerun did
      call PrintWrnMsg(this%option,'ADJOINT OPTIMIZATION: the rerun at the &
        &solution stopped at the wall-clock limit; its output is &
        &incomplete.')
      this%optimization_stop_flag = TS_STOP_WALLCLOCK_EXCEEDED
    case(TS_ADJ_FORWARD_FAILED)
      ! TAO evaluated this point: its rerun failing is not deterministic
      call PrintWrnMsg(this%option,'ADJOINT OPTIMIZATION: the rerun at the &
        &solution failed; its output is incomplete.')
      this%optimization_stop_flag = TS_STOP_FAILURE
    case default
      ! the output is complete, the adjoint's report of it is not
      write(this%option%io_buffer,'(a,i0,a)') 'ADJOINT OPTIMIZATION: the &
        &rerun at the solution has no gradient (status ',status,').'
      call PrintWrnMsg(this%option)
  end select

end subroutine TSDriverOptimize

! ************************************************************************** !

function TSDriverLowerOrder(ts,order_saved)
  !
  ! Drops a TSBDF of order > 1 to order 1, returning the previous order in
  ! order_saved; false if there is nothing to lower.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  TS :: ts
  PetscInt :: order_saved
  PetscBool :: TSDriverLowerOrder

  TSType :: ts_type
  PetscInt :: order
  PetscErrorCode :: ierr

  TSDriverLowerOrder = PETSC_FALSE
  call TSGetType(ts,ts_type,ierr);CHKERRQ(ierr)
  if (trim(ts_type) /= TSBDF) return
  call TSBDFGetOrder(ts,order,ierr);CHKERRQ(ierr)
  if (order <= 1) return
  order_saved = order
  call TSBDFSetOrder(ts,1,ierr);CHKERRQ(ierr)
  TSDriverLowerOrder = PETSC_TRUE

end function TSDriverLowerOrder

! ************************************************************************** !

subroutine TSDriverCheckTimestepper(timestepper,model_name,option)
  !
  ! Rejects the TIMESTEPPER keywords of fields the driver does not use, as
  ! PETSc's adaptor chooses the steps: those set to other than their
  ! defaults. The driver uses MAX_STEPS, MAX_TS_CUTS, the initial, minimum
  ! and maximum step sizes, and the reduction and growth factors.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
  use Timestepper_TS_class

  implicit none

  class(timestepper_TS_type) :: timestepper
  character(len=*) :: model_name
  type(option_type) :: option

  class(timestepper_TS_type), pointer :: defaults
  character(len=MAXSTRINGLENGTH) :: unsupported

  defaults => TimestepperTSCreate()
  unsupported = ''
  if (timestepper%iaccel /= defaults%iaccel) &
    call TSModelListAdd(unsupported,'TS_ACCELERATION')
  if (timestepper%ntfac /= defaults%ntfac) then
    call TSModelListAdd(unsupported,'DT_FACTOR')
  else if (any(timestepper%tfac /= defaults%tfac)) then
    call TSModelListAdd(unsupported,'DT_FACTOR')
  endif
  if (timestepper%constant_time_step_threshold /= &
      defaults%constant_time_step_threshold) &
    call TSModelListAdd(unsupported,'NUM_STEPS_AFTER_TS_CUT')
  if (timestepper%max_num_contig_revert /= defaults%max_num_contig_revert) &
    call TSModelListAdd(unsupported,'MAX_NUM_CONTIGUOUS_REVERTS')
  if (timestepper%time_step_tolerance /= defaults%time_step_tolerance) &
    call TSModelListAdd(unsupported,'TIMESTEP_OVERSTEP_REL_TOLERANCE')
  if (timestepper%print_ekg) call TSModelListAdd(unsupported,'PRINT_EKG')
  deallocate(defaults%tfac)
  deallocate(defaults)
  if (len_trim(unsupported) > 0) then
    ! (no TS in the text: the preprocessor would expand it)
    option%io_buffer = trim(model_name) // ' does not support the &
      &TIMESTEPPER keywords ' // trim(unsupported) // ': the PETSc &
      &adaptor chooses the time steps.'
    call PrintErrMsg(option)
  endif

end subroutine TSDriverCheckTimestepper

! ************************************************************************** !

subroutine TSDriverSetStepLimits(adapt,hmin,hmax)
  !
  ! TSAdaptSetStepLimits requires hmax > hmin; a minimum step size wins.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  TSAdapt :: adapt
  PetscReal :: hmin
  PetscReal :: hmax

  PetscErrorCode :: ierr

  hmax = max(hmax,nearest(hmin,1.d0))
  call TSAdaptSetStepLimits(adapt,hmin,hmax,ierr);CHKERRQ(ierr)

end subroutine TSDriverSetStepLimits

! ************************************************************************** !

subroutine TSDriverIFunction(ts,time,Y,Ydot,F,ctx,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: Y, Ydot, F
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  SNES :: snes
  PetscInt :: status

  ctx%stage_time = time
  call ctx%model%Evaluate(time,Y,Ydot,F,status)
  if (status /= TS_MODEL_OK) then
    ctx%domain_errors = ctx%domain_errors + 1
    call TSGetSNES(ts,snes,ierr);CHKERRQ(ierr)
    call SNESSetFunctionDomainError(snes,ierr);CHKERRQ(ierr)
  endif
  ierr = 0

end subroutine TSDriverIFunction

! ************************************************************************** !

subroutine TSDriverIJacobian(ts,time,Y,Ydot,shift,A,B,ctx,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: Y, Ydot
  PetscReal :: shift
  Mat :: A, B
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  SNES :: snes
  PetscInt :: status

  ctx%stage_time = time
  call ctx%model%Linearize(time,Y,Ydot,shift,A,B,status)
  if (status /= TS_MODEL_OK) then
    call TSGetSNES(ts,snes,ierr);CHKERRQ(ierr)
    call SNESSetJacobianDomainError(snes,ierr);CHKERRQ(ierr)
  endif
  ierr = 0

end subroutine TSDriverIJacobian

! ************************************************************************** !

subroutine TSDriverInitialGuess(snes,X,ctx,ierr)
  !
  ! The model's ProjectInitialGuess, given the step's start state.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  SNES :: snes
  Vec :: X
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  Vec :: Xstart

  call TSGetSolution(ctx%ts,Xstart,ierr);CHKERRQ(ierr)
  call ctx%model%ProjectInitialGuess(X,Xstart)
  ierr = 0

end subroutine TSDriverInitialGuess

! ************************************************************************** !

subroutine TSDriverPreCheck(linesearch,X,dX,changed,ctx,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  SNESLineSearch :: linesearch
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  call ctx%model%PreCheck(X,dX,changed)
  ! models decide per process; SNESLineSearchPreCheck requires one value
  call MPI_Allreduce(MPI_IN_PLACE,changed,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,ctx%option%mycomm,ierr);CHKERRQ(ierr)

end subroutine TSDriverPreCheck

! ************************************************************************** !

subroutine TSDriverPostStep(ts,ierr)
  !
  ! Accepted-step bookkeeping (step size limits, BDF order, the step size
  ! the next interval resumes from) and output. Never modifies the solution.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  implicit none

  TS :: ts
  PetscErrorCode :: ierr

  type(c_ptr) :: ctx_c
  type(ts_driver_ctx_type), pointer :: ctx
  TSConvergedReason :: reason
  TSAdapt :: adapt
  Vec :: Y
  PetscReal :: time, prev_time, next_dt, dt_limit, hmin, hmax, max_time
  PetscInt :: steps, failures, domain_errors
  PetscBool :: interval_end

  call TSGetApplicationContext(ts,ctx_c,ierr);CHKERRQ(ierr)
  call c_f_pointer(ctx_c,ctx)

  ! TSSolve also calls PostStep after a failed step
  call TSGetConvergedReason(ts,reason,ierr);CHKERRQ(ierr)
  if (reason%v < 0) return

  call TSGetTime(ts,time,ierr);CHKERRQ(ierr)
  call TSGetPrevTime(ts,prev_time,ierr);CHKERRQ(ierr)
  call TSGetTimeStep(ts,next_dt,ierr);CHKERRQ(ierr)
  call TSGetSolution(ts,Y,ierr);CHKERRQ(ierr)

  ctx%dt_accepted = time-prev_time
  dt_limit = ctx%dt_max
  call TSGetAdapt(ts,adapt,ierr);CHKERRQ(ierr)
  call TSAdaptGetStepLimits(adapt,hmin,hmax,ierr);CHKERRQ(ierr)
  if (dt_limit /= hmax) then
    call TSDriverSetStepLimits(adapt,hmin,dt_limit)
  endif
  if (next_dt > dt_limit) then
    next_dt = dt_limit
    call TSSetTimeStep(ts,next_dt,ierr);CHKERRQ(ierr)
  endif
  ctx%unmatched_dt = min(ctx%unmatched_dt,dt_limit)

  ! a BDF of order > 1 whose stages keep failing (a vanishing species
  ! extrapolated below zero) continues the interval at order 1; TSAdapt's
  ! error-control rejections do not count, and PETSc counts domain errors
  ! apart from SNES failures
  if (ctx%order_lowered == 0) then
    call TSGetSNESFailures(ts,failures,ierr);CHKERRQ(ierr)
    domain_errors = ctx%domain_errors - ctx%interval_start_domain_errors
    call MPI_Allreduce(MPI_IN_PLACE,domain_errors,ONE_INTEGER_MPI, &
                       MPIU_INTEGER,MPI_MAX,ctx%option%mycomm, &
                       ierr);CHKERRQ(ierr)
    failures = failures - ctx%interval_start_snes_failures + domain_errors
    if (failures >= STAGE_FAILURES_TO_LOWER) then
      if (TSDriverLowerOrder(ts,ctx%order_lowered)) then
        ctx%order_lowerings = ctx%order_lowerings + 1
        write(ctx%option%io_buffer,'(a,i6,a,es12.5)') &
          ' ' // trim(ctx%model%name) // ' TS: ',failures, &
          ' failed stages; BDF order 1 until the next waypoint, from t = ', &
          time
        call PrintMsg(ctx%option)
      endif
    endif
  endif

  ! TSAdaptChoose shortens a proposal to hit the max time unless twice it
  ! fits, and sizes later ones from the shortened step: record the last
  ! proposal it left alone
  call TSGetMaxTime(ts,max_time,ierr);CHKERRQ(ierr)
  if (time+2.d0*next_dt < max_time) ctx%unmatched_dt = next_dt

  call TSGetStepNumber(ts,steps,ierr);CHKERRQ(ierr)
  interval_end = (time >= max_time-64.d0*epsilon(1.d0)*max(1.d0,abs(max_time)))
  call TSDriverStepOutput(ctx,time,Y,steps,interval_end)

  ! -adj_opt_fail_eval: a step failure, as from the solver
  if (ctx%fail_at_step > 0 .and. steps >= ctx%fail_at_step) then
    call TSSetConvergedReason(ts,TS_DIVERGED_STEP_REJECTED,ierr);CHKERRQ(ierr)
  endif

end subroutine TSDriverPostStep

! ************************************************************************** !

subroutine TSDriverStepOutput(ctx,time,Y,steps,interval_end)
  !
  ! Output of an accepted step, as in PMCBaseRunToTime: PERIODIC TIMESTEP
  ! files, and time-averaged variables (accumulated at every output call).
  ! An interval's last step leaves its flags to RunToTime. Output is null
  ! below the top coupler.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Output_Aux_module, only : output_option_type

  implicit none

  type(ts_driver_ctx_type) :: ctx
  PetscReal :: time
  Vec :: Y
  PetscInt :: steps
  PetscBool :: interval_end

  class(pmc_ts_driver_type), pointer :: pmc
  type(output_option_type), pointer :: output_option
  PetscBool :: flags(4), due

  pmc => ctx%pmc
  if (.not.associated(pmc)) return
  if (.not.associated(pmc%Output)) return
  if (ctx%quiet) return
  output_option => pmc%pm_list%output_option
  flags(1) = mod(steps,output_option%periodic_snap_output_ts_imod) == 0
  flags(2) = mod(steps,output_option%periodic_obs_output_ts_imod) == 0
  flags(3) = mod(steps,output_option%periodic_msbl_output_ts_imod) == 0
  flags(4) = mod(steps,output_option%periodic_cons_output_ts_imod) == 0
  due = any(flags) .or. &
        associated(output_option%aveg_output_variable_list%first)
  if (interval_end) then
    ctx%output_due = due
    ctx%output_flags = flags
    return
  endif
  if (.not.due) return
  ctx%option%time = time
  call ctx%model%Publish(time,Y)
  ! Output() clears its flag arguments: flags is a local copy
  call pmc%Output(pmc%pm_list%realization_base,flags(1),flags(2), &
                  flags(3),flags(4))

end subroutine TSDriverStepOutput

! ************************************************************************** !

subroutine TSDriverSNESConverged(snes,it,xnorm,unorm,fnorm,reason,ctx,ierr)
  !
  ! PETSc's norm tests, then from the first iteration on the model's test
  ! (-ts_model_newton_kappa, default 0.1; 0 disables it). An update the
  ! model held back in part never converges.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm, unorm, fnorm
  SNESConvergedReason :: reason
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  Vec :: X, dX, F
  PetscReal :: dt
  PetscBool :: converged

  call SNESConvergedDefault(snes,it,xnorm,unorm,fnorm,reason,0,ierr)
  CHKERRQ(ierr)
  if (reason%v > 0 .and. it > 0) then
    if (ctx%model%HoldsUpdate()) reason = SNES_CONVERGED_ITERATING
  endif
  if (reason /= SNES_CONVERGED_ITERATING) return
  if (it < 1 .or. .not.(ctx%newton_kappa > 0.d0)) return
  if (.not.(fnorm < huge(fnorm))) return

  call SNESGetSolution(snes,X,ierr);CHKERRQ(ierr)
  call SNESGetSolutionUpdate(snes,dX,ierr);CHKERRQ(ierr)
  call SNESGetFunction(snes,F,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call TSGetTimeStep(ctx%ts,dt,ierr);CHKERRQ(ierr)
  ! the residual has been evaluated: stage_time is this stage's
  call ctx%model%ConvergenceTest(ctx%stage_time,X,dX,F,dt, &
                                 ctx%newton_kappa,converged)
  call MPI_Allreduce(MPI_IN_PLACE,converged,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LAND,ctx%option%mycomm,ierr);CHKERRQ(ierr)
  if (converged) reason = SNES_CONVERGED_FNORM_RELATIVE

end subroutine TSDriverSNESConverged

! ************************************************************************** !

subroutine TSDriverSetupEvents(this,ts)
  !
  ! Registers the model's event functions with TSEvent, through the C
  ! wrapper in petsc_fortran_interop.c.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(pmc_ts_driver_type) :: this
  TS :: ts

  PetscInt :: nevents, nevents_global
  PetscInt, allocatable :: direction(:)
  PetscBool, allocatable :: terminate(:)
  PetscReal, allocatable :: vtol(:)
  PetscErrorCode :: ierr

  nevents = this%model%NumEvents()
  call MPI_Allreduce(nevents,nevents_global,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)
  if (nevents_global == 0) return
  allocate(direction(max(nevents,1)))
  allocate(terminate(max(nevents,1)))
  allocate(vtol(max(nevents,1)))
  direction = 0
  terminate = PETSC_FALSE
  vtol = 1.d-6
  if (nevents > 0) &
    call this%model%EventSetup(direction(1:nevents),vtol(1:nevents))
  call pflotran_ts_set_event_handler(ts,nevents,direction,terminate, &
                                     TSDriverEventIndicator, &
                                     TSDriverPostEvent,this%ctx,ierr)
  CHKERRQ(ierr)
  call TSSetEventTolerances(ts,PETSC_CURRENT_REAL,vtol,ierr);CHKERRQ(ierr)
  deallocate(direction,terminate,vtol)

end subroutine TSDriverSetupEvents

! ************************************************************************** !

subroutine TSDriverEventIndicator(ts,time,Y,fvalue,ctx,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: Y
  PetscReal :: fvalue(*)
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  PetscInt :: n

  n = ctx%model%NumEvents()
  if (n > 0) call ctx%model%EventIndicator(time,Y,fvalue(1:n))
  ierr = 0

end subroutine TSDriverEventIndicator

! ************************************************************************** !

subroutine TSDriverPostEvent(ts,nzero,izero,time,Y,forward,dt_post,ctx,ierr)
  !
  ! Applies the model's transitions at a located event. The first step
  ! after it is the last accepted step size, within [dt_min, dt_max], so
  ! that TSAdapt keeps control of the steps that follow.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  TS :: ts
  PetscInt :: nzero
  PetscInt :: izero(*)
  PetscReal :: time
  Vec :: Y
  PetscBool :: forward
  PetscReal :: dt_post
  type(ts_driver_ctx_type) :: ctx
  PetscErrorCode :: ierr

  PetscBool :: changed
  PetscInt, allocatable :: events(:)

  ! TSAdjointSolve calls this at recorded events (forward false); runs
  ! with events skip the adjoint, and a restored state must not be mapped
  if (.not.forward) then
    dt_post = -1.d0
    ierr = 0
    return
  endif
  allocate(events(max(nzero,0)))
  if (nzero > 0) events(1:nzero) = izero(1:nzero)
  call ctx%model%PostEvent(time,Y,events,changed)
  deallocate(events)
  call MPI_Allreduce(MPI_IN_PLACE,changed,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,ctx%option%mycomm,ierr);CHKERRQ(ierr)
  if (changed) then
    ! a new Y state makes PETSc restart the integrator and reevaluate the
    ! event signs, which may change with the regimes even where Y did not
    call PetscObjectStateIncrease(Y,ierr);CHKERRQ(ierr)
  endif
  ctx%events = ctx%events + nzero
  dt_post = min(max(ctx%dt_accepted,ctx%dt_min),ctx%dt_max)
  ierr = 0

end subroutine TSDriverPostEvent

! ************************************************************************** !

subroutine TSDriverModelChecks(this,time,dt,t_end)
  !
  ! The requested TS_Model_Checks_module checks at the initial state;
  ! t_end is the end of the first forcing segment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Solver_module
  use Timestepper_TS_class
  use TS_Model_Checks_module

  implicit none

  class(pmc_ts_driver_type) :: this
  PetscReal :: time
  PetscReal :: dt
  PetscReal :: t_end

  type(ts_model_checks_type) :: checks
  type(solver_type), pointer :: solver
  PetscReal :: atol, rtol
  Vec :: vrtol
  PetscErrorCode :: ierr

  select type(ts => this%timestepper)
    class is(timestepper_TS_type)
      solver => ts%solver
  end select
  checks%model => this%model
  checks%option => this%option
  checks%Y = this%solution
  call TSGetTolerances(solver%ts,atol,checks%vatol,rtol,vrtol, &
                       ierr);CHKERRQ(ierr)
  checks%jacobian = solver%Mpre
  checks%model_jacobian = this%model_jacobian
  checks%pc_type = this%deck_pc_type
  checks%shift_inblocks = solver%linear_shift
  if (this%contract_test) call TSModelContractTest(checks,time,dt)
  if (this%column_check) call TSModelColumnCheck(checks,time,dt)
  if (this%adjoint_test) call TSModelAdjointTest(checks,time,dt,t_end)

end subroutine TSDriverModelChecks

! ************************************************************************** !

recursive subroutine PMCTSDriverDestroy(this)
  !
  ! The coupler (and its TS) first, then the vectors and callback context
  ! the TS references. The model belongs to the process model.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Petsc_Utility_module, only : PUVecDestroy

  implicit none

  class(pmc_ts_driver_type) :: this

  call this%pmc_subsurface_type%Destroy()
  if (associated(this%optimizer)) then
    call this%optimizer%Destroy()
    deallocate(this%optimizer)
    nullify(this%optimizer)
  endif
  if (associated(this%adjoint)) then
    call this%adjoint%Destroy()
    deallocate(this%adjoint)
    nullify(this%adjoint)
  endif
  call PUVecDestroy(this%solution)
  call PUVecDestroy(this%residual)
  if (associated(this%ctx)) then
    nullify(this%ctx%model)
    nullify(this%ctx%pmc)
    deallocate(this%ctx)
    nullify(this%ctx)
  endif
  nullify(this%model)

end subroutine PMCTSDriverDestroy

end module PMC_TS_Driver_class
