module TS_Model_class
  !
  ! Equation-model interface for PETSc TS.
  !
  ! A model defines F(t,Y,Ydot) = 0 and its shifted linearization
  ! J(a) = dF/dY + a*dF/dYdot. PETSc owns time integration and its history;
  ! the driver (PMC_TS_Driver_class) owns TS/SNES/KSP, the stepping policy
  ! and accepted-step bookkeeping; the model owns its state layout,
  ! evaluation scratch and forcing segments.
  !
  ! Contract:
  !   - Evaluate, Linearize, ConvergenceTest, EventIndicator and Publish
  !     depend only on their arguments, the current forcing segment and the
  !     discrete regime: never on a global time cursor, option%tran_dt, or
  !     an earlier evaluation's scratch unless cached for the same
  !     arguments. Linearize never requires a preceding Evaluate.
  !   - Evaluation writes private scratch only; accepted physical state
  !     lives in the TS solution and nothing is integrated outside it.
  !   - Linearize is valid for any shift, including zero.
  !   - A forcing segment is immutable once BeginSegment returns; a later
  !     segment is a new one (new segment_id). BeginSegment reads the
  !     segment's forcing from the realization, and may leave its
  !     conditions, couplers and velocities at the segment's end (dataset
  !     time storage only moves forward); Publish writes the state at its
  !     time back for output. Nothing else reads them in between.
  !   - After RecordForcing every segment begun is recorded; between
  !     BeginReplay and EndReplay (an adjoint's backward sweep) an
  !     evaluation at any time of the run uses that time's recorded segment.
  !     A model holding forcing ends BeginSegment with RecordSegment and
  !     calls ReplayAt before reading its forcing.
  !
  ! ConvergenceTest and PreCheck play the roles of PM_Base's
  ! CheckConvergence and CheckUpdatePre: a model must work without a
  ! process model (Flow_RT_DAE_Model_class composes two).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use PFLOTRAN_Constants_module

  implicit none

  private

  ! evaluation status
  PetscInt, parameter, public :: TS_MODEL_OK = 0
  ! recoverable: the trial state is outside the model's domain
  PetscInt, parameter, public :: TS_MODEL_DOMAIN_ERROR = 1

  ! End times of the recorded forcing segments, in run order; the model
  ! stores segment k as record k (TSModelRecordsSelect picks one by time)
  type, public :: ts_model_records_type
    PetscInt :: n = 0
    PetscReal, allocatable :: t_end(:)
  end type ts_model_records_type

  public :: TSModelListAdd, &
            TSModelReduceStatus

  type, public, abstract :: ts_model_type
    character(len=MAXWORDLENGTH) :: name = ''
    ! dofs per grid cell of the TS solution; zero if the solution is not
    ! laid out cell by cell (see CreateVector)
    PetscInt :: ndof = 0
    ! options prefix of the TS and its solvers
    character(len=MAXWORDLENGTH) :: prefix = 'tran_'
    ! which of option%flow_dt/tran_dt the driver sets to the step size
    PetscBool :: integrates_flow = PETSC_FALSE
    PetscBool :: integrates_transport = PETSC_FALSE
    ! generation number of the current forcing segment
    PetscInt :: segment_id = 0
    ! forcing records for a replay (see the module contract); replay_record
    ! is the record in use
    PetscBool :: record = PETSC_FALSE
    PetscBool :: replay = PETSC_FALSE
    type(ts_model_records_type) :: recorded
    PetscInt :: replay_record = 0
  contains
    procedure(TSModelInitialState), deferred :: InitialState
    procedure(TSModelBeginSegment), deferred :: BeginSegment
    procedure(TSModelEvaluate), deferred :: Evaluate
    procedure(TSModelLinearize), deferred :: Linearize
    procedure, public :: CreateVector => TSModelCreateVector
    procedure, public :: CreateJacobian => TSModelCreateJacobian
    procedure, public :: SetupLinearSolver => TSModelSetupLinearSolver
    procedure(TSModelSetTolerances), deferred :: SetTolerances
    procedure, public :: ConvergenceTest => TSModelConvergenceTest
    procedure, public :: DefaultTolerances => TSModelDefaultTolerances
    procedure(TSModelPublish), deferred :: Publish
    procedure, public :: PreCheck => TSModelPreCheck
    procedure, public :: HoldsUpdate => TSModelHoldsUpdate
    procedure, public :: ProjectInitialGuess => TSModelProjectInitialGuess
    procedure, public :: NumEvents => TSModelNumEvents
    procedure, public :: EventSetup => TSModelEventSetup
    procedure, public :: EventIndicator => TSModelEventIndicator
    procedure, public :: PostEvent => TSModelPostEvent
    procedure, public :: PrintSummary => TSModelPrintSummary
    procedure, public :: ForcingTimes => TSModelForcingTimes
    procedure, public :: RecordForcing => TSModelRecordForcing
    procedure, public :: BeginReplay => TSModelBeginReplay
    procedure, public :: EndReplay => TSModelEndReplay
    procedure, public, non_overridable :: RecordSegment => TSModelRecordSegment
    procedure, public, non_overridable :: ReplayAt => TSModelReplayAt
    procedure, public :: StoreRecord => TSModelStoreRecord
    procedure, public :: RestoreRecord => TSModelRestoreRecord
    procedure, public :: InvalidateCache => TSModelInvalidateCache
    procedure, public :: Destroy => TSModelDestroy
  end type ts_model_type

  abstract interface

    subroutine TSModelInitialState(this,time,Y)
      ! Fills Y with a state consistent with the algebraic equations at time.
      use petscts
      import :: ts_model_type
      implicit none
      class(ts_model_type) :: this
      PetscReal :: time
      Vec :: Y
    end subroutine TSModelInitialState

    subroutine TSModelBeginSegment(this,t_start,t_end,discontinuous)
      ! Builds the immutable forcing of [t_start, t_end]. Called once per
      ! waypoint interval, in time order; a redo reuses the interval's
      ! segment. discontinuous: the forcing jumps at t_start, so the TS
      ! history must restart.
      use petscts
      import :: ts_model_type
      implicit none
      class(ts_model_type) :: this
      PetscReal :: t_start
      PetscReal :: t_end
      PetscBool :: discontinuous
    end subroutine TSModelBeginSegment

    subroutine TSModelEvaluate(this,time,Y,Ydot,F,status)
      use petscts
      import :: ts_model_type
      implicit none
      class(ts_model_type) :: this
      PetscReal :: time
      Vec :: Y
      Vec :: Ydot
      Vec :: F
      PetscInt :: status
    end subroutine TSModelEvaluate

    subroutine TSModelLinearize(this,time,Y,Ydot,shift,A,B,status)
      use petscts
      import :: ts_model_type
      implicit none
      class(ts_model_type) :: this
      PetscReal :: time
      Vec :: Y
      Vec :: Ydot
      PetscReal :: shift
      Mat :: A
      Mat :: B
      PetscInt :: status
    end subroutine TSModelLinearize

    subroutine TSModelSetTolerances(this,atol,rtol,vatol,vrtol)
      ! Fills the per-entry local error tolerances; algebraic entries get
      ! zero, which excludes them from the error estimate.
      use petscts
      import :: ts_model_type
      implicit none
      class(ts_model_type) :: this
      PetscReal :: atol
      PetscReal :: rtol
      Vec :: vatol
      Vec :: vrtol
    end subroutine TSModelSetTolerances

    subroutine TSModelPublish(this,time,Y)
      ! Derives the legacy output fields from an accepted state. One way:
      ! evaluation never reads what Publish writes.
      use petscts
      import :: ts_model_type
      implicit none
      class(ts_model_type) :: this
      PetscReal :: time
      Vec :: Y
    end subroutine TSModelPublish

  end interface

contains

! ************************************************************************** !

subroutine TSModelCreateVector(this,Y,created)
  !
  ! Optional: a new TS solution vector. created = PETSC_FALSE: the driver
  ! creates one with ndof entries per grid cell.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  Vec :: Y
  PetscBool :: created

  created = PETSC_FALSE

end subroutine TSModelCreateVector

! ************************************************************************** !

subroutine TSModelCreateJacobian(this,J,created)
  !
  ! Optional: a new Jacobian operator of the model's own kind (e.g. with a
  ! condensed preconditioner). created = PETSC_FALSE: the driver assembles
  ! a block matrix of ndof entries per cell.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  Mat :: J
  PetscBool :: created

  created = PETSC_FALSE

end subroutine TSModelCreateJacobian

! ************************************************************************** !

subroutine TSModelSetupLinearSolver(this,ksp,J,pc_type,shift_inblocks)
  !
  ! Configures ksp's preconditioner for a CreateJacobian operator J.
  ! pc_type: the deck's preconditioner (blank: PETSc's default);
  ! shift_inblocks: factorizations shift zero pivots within the blocks.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  KSP :: ksp
  Mat :: J
  PCType :: pc_type
  PetscBool :: shift_inblocks

end subroutine TSModelSetupLinearSolver

! ************************************************************************** !

subroutine TSModelDefaultTolerances(this,atol,rtol)
  !
  ! Scalar local error tolerances used unless -<prefix>ts_atol/_rtol are
  ! given; SetTolerances receives the effective values.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol

  atol = 1.d-4
  rtol = 1.d-4

end subroutine TSModelDefaultTolerances

! ************************************************************************** !

subroutine TSModelConvergenceTest(this,time,X,dX,F,dt,kappa,converged)
  !
  ! Optional Newton convergence test of a stage at time, after each
  ! iteration: iterate X, update dX, residual F, step size dt. The
  ! iteration error may use the fraction kappa of the local error
  ! tolerance. The default defers to the SNES norm tests.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: time
  Vec :: X
  Vec :: dX
  Vec :: F
  PetscReal :: dt
  PetscReal :: kappa
  PetscBool :: converged

  converged = PETSC_FALSE

end subroutine TSModelConvergenceTest

! ************************************************************************** !

subroutine TSModelPreCheck(this,Y,dY,changed)
  !
  ! Optional Newton update limiter (a solver strategy, not physics). The
  ! default leaves the update untouched.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  Vec :: Y
  Vec :: dY
  PetscBool :: changed

  changed = PETSC_FALSE

end subroutine TSModelPreCheck

! ************************************************************************** !

function TSModelHoldsUpdate(this)
  !
  ! Whether PreCheck held part of the last Newton update back: the
  ! iteration has not converged, whatever the norm tests say. Collective.
  ! The default holds nothing.
  !
  ! Author: Hansol Suh
  ! Date: 10/13/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscBool :: TSModelHoldsUpdate

  TSModelHoldsUpdate = PETSC_FALSE

end function TSModelHoldsUpdate

! ************************************************************************** !

subroutine TSModelProjectInitialGuess(this,X,Xstart)
  !
  ! Optional: projects a stage's Newton initial guess X (e.g. a multistep
  ! extrapolation) into the model's domain, given the step's start state
  ! Xstart. A solver strategy; the equations are unchanged.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  Vec :: X
  Vec :: Xstart

end subroutine TSModelProjectInitialGuess

! ************************************************************************** !

function TSModelNumEvents(this)
  !
  ! Number of event functions on this process, one per possible regime
  ! transition (e.g. a mineral that disappears). Default: none.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscInt :: TSModelNumEvents

  TSModelNumEvents = 0

end function TSModelNumEvents

! ************************************************************************** !

subroutine TSModelEventSetup(this,direction,vtol)
  !
  ! Zero-crossing direction (-1, 0, +1) and value tolerance of each event
  ! function.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscInt :: direction(:)
  PetscReal :: vtol(:)

  direction = 0
  vtol = 1.d-6

end subroutine TSModelEventSetup

! ************************************************************************** !

subroutine TSModelEventIndicator(this,time,Y,fvalue)
  !
  ! Values of the event functions at (time,Y); a transition happens where
  ! one crosses zero in its declared direction.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscReal :: fvalue(:)

  fvalue = 1.d0

end subroutine TSModelEventIndicator

! ************************************************************************** !

subroutine TSModelPostEvent(this,time,Y,events,changed)
  !
  ! Applies the transition map of each located event (local indices) at an
  ! accepted state. changed: Y or a regime changed, so the integrator
  ! history restarts.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscInt :: events(:)
  PetscBool :: changed

  changed = PETSC_FALSE

end subroutine TSModelPostEvent

! ************************************************************************** !

subroutine TSModelPrintSummary(this)
  !
  ! End-of-run diagnostics of the model (e.g. regime transitions).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this

end subroutine TSModelPrintSummary

! ************************************************************************** !

subroutine TSModelForcingTimes(this,final_time,times)
  !
  ! Optional: the knots of the model's forcing datasets up to final_time.
  ! The driver makes each a waypoint so that no segment spans one.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: final_time
  PetscReal, allocatable :: times(:)

  allocate(times(0))

end subroutine TSModelForcingTimes

! ************************************************************************** !

subroutine TSModelRecordForcing(this)
  !
  ! Records every forcing segment from now on, for a replay. Must precede
  ! the first BeginSegment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this

  this%record = PETSC_TRUE

end subroutine TSModelRecordForcing

! ************************************************************************** !

subroutine TSModelBeginReplay(this)
  !
  ! From now on, evaluate with the recorded segment of each time (see
  ! ReplayAt).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this

  this%replay = PETSC_TRUE
  this%replay_record = 0

end subroutine TSModelBeginReplay

! ************************************************************************** !

subroutine TSModelEndReplay(this)
  !
  ! Back to the last forcing segment of the run.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this

  if (.not.this%replay) return
  if (this%recorded%n > 0) call this%RestoreRecord(this%recorded%n)
  this%replay = PETSC_FALSE
  this%replay_record = 0

end subroutine TSModelEndReplay

! ************************************************************************** !

subroutine TSModelRecordSegment(this,t_end)
  !
  ! After RecordForcing, records the segment just begun (ending at t_end)
  ! as record recorded%n, via StoreRecord. Call at the end of BeginSegment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: t_end

  if (.not.this%record) return
  call TSModelRecordsAppend(this%recorded,t_end)
  call this%StoreRecord(this%recorded%n)

end subroutine TSModelRecordSegment

! ************************************************************************** !

subroutine TSModelReplayAt(this,time)
  !
  ! During a replay, makes the record of time current via RestoreRecord.
  ! Call before reading forcing. Without records it does nothing.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscReal :: time

  PetscInt :: k

  if (.not.this%replay .or. this%recorded%n == 0) return
  k = TSModelRecordsSelect(this%recorded,time,this%replay_record)
  if (k == this%replay_record) return
  call this%RestoreRecord(k)
  this%replay_record = k

end subroutine TSModelReplayAt

! ************************************************************************** !

subroutine TSModelStoreRecord(this,k)
  !
  ! Keeps the current forcing segment as record k (= recorded%n). The
  ! default suits a model whose forcing is the same in every segment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscInt :: k

end subroutine TSModelStoreRecord

! ************************************************************************** !

subroutine TSModelRestoreRecord(this,k)
  !
  ! Makes record k the current forcing segment, expiring anything prepared
  ! with another.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this
  PetscInt :: k

end subroutine TSModelRestoreRecord

! ************************************************************************** !

subroutine TSModelInvalidateCache(this)
  !
  ! Discards any cached evaluation: parameters the model reads (e.g. a
  ! reaction rate constant) have changed in place.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(ts_model_type) :: this

end subroutine TSModelInvalidateCache

! ************************************************************************** !

subroutine TSModelRecordsAppend(records,t_end)
  !
  ! Records a segment that ends at t_end, as record records%n.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(ts_model_records_type) :: records
  PetscReal :: t_end

  PetscReal, allocatable :: grown(:)

  if (.not.allocated(records%t_end)) allocate(records%t_end(16))
  if (records%n == size(records%t_end)) then
    allocate(grown(2*records%n))
    grown(1:records%n) = records%t_end
    call move_alloc(grown,records%t_end)
  endif
  records%n = records%n + 1
  records%t_end(records%n) = t_end

end subroutine TSModelRecordsAppend

! ************************************************************************** !

function TSModelRecordsSelect(records,time,current)
  !
  ! The record of time: the first that ends at or after it (the last one
  ! past the run's end). Backward Euler evaluates step (t_n, t_n+1] at
  ! t_n+1, so a segment boundary selects the segment ending there.
  ! current, the last selection (0: none), is tried first.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(ts_model_records_type) :: records
  PetscReal :: time
  PetscInt :: current
  PetscInt :: TSModelRecordsSelect

  PetscReal :: tol
  PetscInt :: k

  tol = 64.d0*epsilon(1.d0)*max(1.d0,abs(time))
  k = current
  if (k > 0) then
    if (time <= records%t_end(k)+tol) then
      TSModelRecordsSelect = k
      if (k == 1) return
      if (time > records%t_end(k-1)+tol) return
    endif
  endif
  do k = 1, records%n
    if (time <= records%t_end(k)+tol) exit
  enddo
  TSModelRecordsSelect = min(k,records%n)

end function TSModelRecordsSelect

! ************************************************************************** !

subroutine TSModelListAdd(list,item)
  !
  ! Appends item to the comma-separated list (e.g. the unsupported features
  ! a model reports in one message).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  character(len=*) :: list
  character(len=*) :: item

  if (len_trim(list) > 0) then
    list = trim(list) // ', ' // item
  else
    list = item
  endif

end subroutine TSModelListAdd

! ************************************************************************** !

subroutine TSModelReduceStatus(status,comm)
  !
  ! The worst evaluation status over comm: a status is rank local, a
  ! decision taken on it must be collective.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  PetscInt :: status
  PetscMPIInt :: comm

  PetscErrorCode :: ierr

  call MPI_Allreduce(MPI_IN_PLACE,status,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,comm,ierr);CHKERRQ(ierr)

end subroutine TSModelReduceStatus

! ************************************************************************** !

subroutine TSModelDestroy(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none
  class(ts_model_type) :: this
end subroutine TSModelDestroy

end module TS_Model_class
