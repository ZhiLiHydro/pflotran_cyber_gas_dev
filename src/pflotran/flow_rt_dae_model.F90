module Flow_RT_DAE_Model_class
  !
  ! Richards flow and reactive transport as one DAE model for PETSc TS
  ! (MODE RICHARDS_DAE with MODE GIRT_DAE).
  !
  ! Y = [Y_f; Y_t] on each process: the flow state (Richards_DAE_Model_class,
  ! [P; W] per cell), then the transport state (RT_DAE_Model_class,
  ! [c; N; K] per cell). Transport reads the flow state of Y_f
  ! (TS_Flow_State_module): saturation, porosity, water density, Darcy
  ! fluxes, boundary state, source volume rates. Flow does not depend on
  ! transport, so
  !   J = [ J_ff    0   ]
  !       [ J_tf   J_tt ],
  ! a MATNEST of the flow matrix, the transport pressure derivatives
  ! (RTDAEAddFlowCoupling) and the transport operator. The preconditioner
  ! is block forward substitution (flow solve, then the transport
  ! preconditioner), exact with exact sub-solves; its transpose serves
  ! adjoint solves.
  !
  ! Transport's inventory and fluxes use the flow state of the evaluation
  ! time, as in the classic sequential scheme (RICHARDS, then GIRT); only
  ! the boundary totals are frozen per waypoint interval.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use TS_Model_class
  use TS_Flow_State_module
  use Richards_DAE_Model_class
  use RT_DAE_Model_class
  use RT_DAE_Flow_Coupling_module
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! PCSHELL context
  type :: flow_rt_pc_ctx_type
    class(flow_rt_dae_model_type), pointer :: model => null()
  end type flow_rt_pc_ctx_type

  type, public, extends(ts_model_type) :: flow_rt_dae_model_type
    type(option_type), pointer :: option => null()
    ! the sub-models, owned by their process models
    class(richards_dae_model_type), pointer :: flow => null()
    class(rt_dae_model_type), pointer :: rt => null()
    type(ts_flow_state_type) :: fs
    ! local sizes of the parts, and their global index sets
    PetscInt :: nflow = 0
    PetscInt :: ntran = 0
    IS :: is_flow
    IS :: is_tran
    ! absolute tolerance of the flow model, -ts_model_flow_atol
    PetscReal :: flow_atol = 1.d-8
    ! -ts_model_flow_first (see FlowRTDAEPreCheck); whether the last Newton
    ! update held transport, and the number of held iterations
    PetscBool :: flow_first = PETSC_TRUE
    PetscBool :: transport_held = PETSC_FALSE
    PetscInt :: flow_first_iterations = 0
    ! the preconditioner's sub-solvers, and scratch of the two parts
    type(flow_rt_pc_ctx_type) :: pc_ctx
    KSP :: ksp_flow
    KSP :: ksp_tran
    Vec :: work_flow
    Vec :: work_tran
    PetscBool :: solvers_created = PETSC_FALSE
  contains
    procedure, public :: InitialState => FlowRTDAEInitialState
    procedure, public :: BeginSegment => FlowRTDAEBeginSegment
    procedure, public :: Evaluate => FlowRTDAEEvaluate
    procedure, public :: Linearize => FlowRTDAELinearize
    procedure, public :: CreateVector => FlowRTDAECreateVector
    procedure, public :: CreateJacobian => FlowRTDAECreateJacobian
    procedure, public :: SetupLinearSolver => FlowRTDAESetupLinearSolver
    procedure, public :: SetTolerances => FlowRTDAESetTolerances
    procedure, public :: DefaultTolerances => FlowRTDAEDefaultTolerances
    procedure, public :: ConvergenceTest => FlowRTDAEConvergenceTest
    procedure, public :: Publish => FlowRTDAEPublish
    procedure, public :: PreCheck => FlowRTDAEPreCheck
    procedure, public :: HoldsUpdate => FlowRTDAEHoldsUpdate
    procedure, public :: ProjectInitialGuess => FlowRTDAEProjectInitialGuess
    procedure, public :: NumEvents => FlowRTDAENumEvents
    procedure, public :: EventSetup => FlowRTDAEEventSetup
    procedure, public :: EventIndicator => FlowRTDAEEventIndicator
    procedure, public :: PostEvent => FlowRTDAEPostEvent
    procedure, public :: PrintSummary => FlowRTDAEPrintSummary
    procedure, public :: ForcingTimes => FlowRTDAEForcingTimes
    procedure, public :: RecordForcing => FlowRTDAERecordForcing
    procedure, public :: BeginReplay => FlowRTDAEBeginReplay
    procedure, public :: EndReplay => FlowRTDAEEndReplay
    procedure, public :: InvalidateCache => FlowRTDAEInvalidateCache
    procedure, public :: Destroy => FlowRTDAEDestroy
  end type flow_rt_dae_model_type

  public :: FlowRTDAEModelCreate, &
            FlowRTDAEGetParts, &
            FlowRTDAERestoreParts, &
            FlowRTDAESetFlowState

contains

! ************************************************************************** !

function FlowRTDAEModelCreate(flow,rt)
  !
  ! The coupled model of a flow model on a private DM and a transport model.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type), pointer :: flow
  class(rt_dae_model_type), pointer :: rt
  class(flow_rt_dae_model_type), pointer :: FlowRTDAEModelCreate

  class(flow_rt_dae_model_type), pointer :: this
  PetscInt :: rstart
  PetscBool :: flag, flag2
  Vec :: v
  PetscErrorCode :: ierr

  allocate(this)
  this%name = 'Flow RT DAE'
  this%prefix = rt%prefix
  this%integrates_flow = PETSC_TRUE
  this%integrates_transport = PETSC_TRUE
  ! not laid out cell by cell
  this%ndof = 0
  this%flow => flow
  this%rt => rt
  this%option => rt%option
  if (.not.flow%private_dm) then
    call PrintErrMsg(this%option,'Flow RT DAE: the flow model must have a &
      &private DM.')
  endif

  this%nflow = flow%nlmax*flow%ndof
  this%ntran = rt%nlmax*rt%ndof
  call VecCreateMPI(this%option%mycomm,this%nflow+this%ntran,PETSC_DETERMINE, &
                    v,ierr);CHKERRQ(ierr)
  call VecGetOwnershipRange(v,rstart,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecDestroy(v,ierr);CHKERRQ(ierr)
  call ISCreateStride(this%option%mycomm,this%nflow,rstart,ONE_INTEGER, &
                      this%is_flow,ierr);CHKERRQ(ierr)
  call ISCreateStride(this%option%mycomm,this%ntran,rstart+this%nflow, &
                      ONE_INTEGER,this%is_tran,ierr);CHKERRQ(ierr)

  ! transport reads the flow state's fluxes by its own connection index;
  ! both tables are TSConnectionsBuild's of the same patch
  call TSFlowStateInit(this%fs,flow%ngmax,flow%conn%nconn_int, &
                       flow%conn%nconn_bc,flow%conn%nconn_ss)
  ! no flow state yet
  this%fs%id = -1
  call RTDAESetFlowCoupling(rt,this%fs)
  this%pc_ctx%model => this

  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_flow_atol',this%flow_atol,flag, &
                           ierr);CHKERRQ(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-ts_model_flow_first',this%flow_first,flag, &
                           ierr);CHKERRQ(ierr)
  ! the flow part has no TS of its own: its tolerances would be ignored
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-flow_ts_atol',flag,ierr);CHKERRQ(ierr)
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-flow_ts_rtol',flag2,ierr);CHKERRQ(ierr)
  if (flag .or. flag2) then
    call PrintErrMsg(this%option,'Flow and transport are one TS: the flow &
      &tolerances are -ts_model_flow_atol and -tran_ts_rtol, not &
      &-flow_ts_atol/-flow_ts_rtol.')
  endif

  FlowRTDAEModelCreate => this

end function FlowRTDAEModelCreate

! ************************************************************************** !

subroutine FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  !
  ! The flow and transport parts of Y (VecGetSubVector); return them with
  ! FlowRTDAERestoreParts.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  Vec :: Y
  Vec :: Y_f
  Vec :: Y_t

  PetscErrorCode :: ierr

  call VecGetSubVector(Y,this%is_flow,Y_f,ierr);CHKERRQ(ierr)
  call VecGetSubVector(Y,this%is_tran,Y_t,ierr);CHKERRQ(ierr)

end subroutine FlowRTDAEGetParts

! ************************************************************************** !

subroutine FlowRTDAERestoreParts(this,Y,Y_f,Y_t)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this
  Vec :: Y
  Vec :: Y_f
  Vec :: Y_t

  PetscErrorCode :: ierr

  call VecRestoreSubVector(Y,this%is_flow,Y_f,ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(Y,this%is_tran,Y_t,ierr);CHKERRQ(ierr)

end subroutine FlowRTDAERestoreParts

! ************************************************************************** !

subroutine FlowRTDAESetFlowState(this,time,Y_f,status)
  !
  ! Fills the flow state transport reads from Y_f.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y_f
  PetscInt :: status

  call RichardsDAEFlowState(this%flow,time,Y_f,this%fs,status)

end subroutine FlowRTDAESetFlowState

! ************************************************************************** !

subroutine FlowRTDAEInitialState(this,time,Y)
  !
  ! The flow model's initial state, then the transport model's at the
  ! initial flow state.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y

  Vec :: Y_f, Y_t
  PetscInt :: status

  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call this%flow%InitialState(time,Y_f)
  call FlowRTDAESetFlowState(this,time,Y_f,status)
  ! PrintErrMsg is collective
  call TSModelReduceStatus(status,this%option%mycomm)
  if (status /= TS_MODEL_OK) then
    call PrintErrMsg(this%option,'Flow RT DAE: the initial flow state could &
      &not be evaluated.')
  endif
  call this%rt%InitialState(time,Y_t)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)

end subroutine FlowRTDAEInitialState

! ************************************************************************** !

subroutine FlowRTDAEBeginSegment(this,t_start,t_end,discontinuous)
  !
  ! Both models' forcing of [t_start, t_end].
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: t_start
  PetscReal :: t_end
  PetscBool :: discontinuous

  PetscBool :: flow_jumps, tran_jumps

  call this%flow%BeginSegment(t_start,t_end,flow_jumps)
  call this%rt%BeginSegment(t_start,t_end,tran_jumps)
  discontinuous = flow_jumps .or. tran_jumps
  this%segment_id = this%segment_id + 1

end subroutine FlowRTDAEBeginSegment

! ************************************************************************** !

subroutine FlowRTDAEEvaluate(this,time,Y,Ydot,F,status)
  !
  ! F = [F_f(Y_f); F_t(Y_t; flow state of Y_f)].
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  Vec :: F
  PetscInt :: status

  Vec :: Y_f, Y_t, Ydot_f, Ydot_t, F_f, F_t
  PetscInt :: tran_status
  PetscErrorCode :: ierr

  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call FlowRTDAEGetParts(this,Ydot,Ydot_f,Ydot_t)
  call FlowRTDAEGetParts(this,F,F_f,F_t)
  call this%flow%Evaluate(time,Y_f,Ydot_f,F_f,status)
  ! every rank takes the same branch: the transport evaluation is collective
  call TSModelReduceStatus(status,this%option%mycomm)
  if (status == TS_MODEL_OK) then
    call FlowRTDAESetFlowState(this,time,Y_f,status)
    call this%rt%Evaluate(time,Y_t,Ydot_t,F_t,tran_status)
    status = max(status,tran_status)
  else
    ! no flow state for transport
    call VecSet(F_t,ieee_value(1.d0,ieee_quiet_nan),ierr);CHKERRQ(ierr)
  endif
  call FlowRTDAERestoreParts(this,F,F_f,F_t)
  call FlowRTDAERestoreParts(this,Ydot,Ydot_f,Ydot_t)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)

end subroutine FlowRTDAEEvaluate

! ************************************************************************** !

subroutine FlowRTDAELinearize(this,time,Y,Ydot,shift,A,B,status)
  !
  ! The blocks of the MATNEST A: J_ff and J_tt by the models, J_tf the
  ! transport pressure derivatives. B is not used.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  PetscReal :: shift
  Mat :: A
  Mat :: B

  PetscInt :: status

  Vec :: Y_f, Y_t, Ydot_f, Ydot_t
  Mat :: J_ff, J_tf, J_tt
  PetscInt :: tran_status
  PetscErrorCode :: ierr

  call MatNestGetSubMat(A,ZERO_INTEGER,ZERO_INTEGER,J_ff,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(A,ONE_INTEGER,ZERO_INTEGER,J_tf,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(A,ONE_INTEGER,ONE_INTEGER,J_tt,ierr);CHKERRQ(ierr)
  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call FlowRTDAEGetParts(this,Ydot,Ydot_f,Ydot_t)
  call this%flow%Linearize(time,Y_f,Ydot_f,shift,J_ff,J_ff,status)
  call FlowRTDAESetFlowState(this,time,Y_f,tran_status)
  status = max(status,tran_status)
  call this%rt%Linearize(time,Y_t,Ydot_t,shift,J_tt,J_tt,tran_status)
  status = max(status,tran_status)
  call MatZeroEntries(J_tf,ierr);CHKERRQ(ierr)
  call RTDAEAddFlowCoupling(this%rt,J_tf,this%flow%ndof)
  call MatAssemblyBegin(J_tf,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J_tf,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call FlowRTDAERestoreParts(this,Ydot,Ydot_f,Ydot_t)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine FlowRTDAELinearize

! ************************************************************************** !

subroutine FlowRTDAECreateVector(this,Y,created)
  !
  ! A vector of the model's layout.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  Vec :: Y
  PetscBool :: created

  PetscErrorCode :: ierr

  call VecCreateMPI(this%option%mycomm,this%nflow+this%ntran, &
                    PETSC_DETERMINE,Y,ierr);CHKERRQ(ierr)
  created = PETSC_TRUE

end subroutine FlowRTDAECreateVector

! ************************************************************************** !

subroutine FlowRTDAECreateJacobian(this,J,created)
  !
  ! The MATNEST [[J_ff, 0], [J_tf, J_tt]]; it owns its blocks.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Discretization_module

  implicit none

  class(flow_rt_dae_model_type) :: this
  Mat :: J
  PetscBool :: created

  Mat :: blocks(4)
  IS :: is(2)
  PetscBool :: tran_created
  ! a MatType dummy reads 80 characters; the constant MATBAIJ has 4
  MatType :: mat_type = MATBAIJ
  PetscErrorCode :: ierr

  call RichardsDAECreateMatrix(this%flow,blocks(1))
  blocks(2) = PETSC_NULL_MAT
  call FlowRTDAECreateCoupling(this,blocks(3))
  call this%rt%CreateJacobian(blocks(4),tran_created)
  if (.not.tran_created) then
    call DiscretizationCreateMatrix(this%rt%realization%discretization, &
                                    NDAEDOF,mat_type,blocks(4),PETSC_FALSE, &
                                    this%option)
  endif
  is(1) = this%is_flow
  is(2) = this%is_tran
  call MatCreateNest(this%option%mycomm,TWO_INTEGER,is,TWO_INTEGER,is, &
                     blocks,J,ierr);CHKERRQ(ierr)
  ! the nest holds references
  call MatDestroy(blocks(1),ierr);CHKERRQ(ierr)
  call MatDestroy(blocks(3),ierr);CHKERRQ(ierr)
  call MatDestroy(blocks(4),ierr);CHKERRQ(ierr)
  created = PETSC_TRUE

end subroutine FlowRTDAECreateJacobian

! ************************************************************************** !

subroutine FlowRTDAECreateCoupling(this,J_tf)
  !
  ! J_tf: transport rows, flow columns, in the models' local numbering
  ! (MatSetValuesLocal). Each row of a cell has the pressures of the cell
  ! and its neighbors.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Grid_module
  use Discretization_module

  implicit none

  class(flow_rt_dae_model_type) :: this
  Mat :: J_tf

  type(grid_type), pointer :: grid
  PetscInt, allocatable :: d_nnz(:), o_nnz(:), d_cell(:), o_cell(:)
  PetscInt :: local_id, iconn, lup, ldn, i, ndof
  ISLocalToGlobalMapping :: ltog_row, ltog_col
  PetscErrorCode :: ierr

  grid => this%rt%realization%patch%grid
  ndof = this%rt%ndof
  ! pressure columns per cell: its own, and one per neighbor (owned or not)
  allocate(d_cell(grid%nlmax),o_cell(grid%nlmax))
  d_cell = 1
  o_cell = 0
  do iconn = 1, this%flow%conn%nconn_int
    lup = grid%nG2L(this%flow%conn%int_up(iconn))
    ldn = grid%nG2L(this%flow%conn%int_dn(iconn))
    if (lup > 0) then
      if (ldn > 0) then
        d_cell(lup) = d_cell(lup) + 1
      else
        o_cell(lup) = o_cell(lup) + 1
      endif
    endif
    if (ldn > 0) then
      if (lup > 0) then
        d_cell(ldn) = d_cell(ldn) + 1
      else
        o_cell(ldn) = o_cell(ldn) + 1
      endif
    endif
  enddo
  allocate(d_nnz(this%ntran),o_nnz(this%ntran))
  do local_id = 1, grid%nlmax
    do i = 1, ndof
      d_nnz((local_id-1)*ndof+i) = d_cell(local_id)
      o_nnz((local_id-1)*ndof+i) = o_cell(local_id)
    enddo
  enddo
  call MatCreate(this%option%mycomm,J_tf,ierr);CHKERRQ(ierr)
  call MatSetSizes(J_tf,this%ntran,this%nflow,PETSC_DETERMINE, &
                   PETSC_DETERMINE,ierr);CHKERRQ(ierr)
  call MatSetType(J_tf,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(J_tf,ZERO_INTEGER,d_nnz,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(J_tf,ZERO_INTEGER,d_nnz,ZERO_INTEGER, &
                                 o_nnz,ierr);CHKERRQ(ierr)
  call FlowRTDAELocalToGlobal(this,this%rt%realization%discretization% &
                                dm_ndaedof,ltog_row)
  call FlowRTDAELocalToGlobal(this,this%flow%dm,ltog_col)
  call MatSetLocalToGlobalMapping(J_tf,ltog_row,ltog_col,ierr);CHKERRQ(ierr)
  call MatSetOption(J_tf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE, &
                    ierr);CHKERRQ(ierr)
  call MatSetOption(J_tf,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
  deallocate(d_nnz,o_nnz,d_cell,o_cell)

end subroutine FlowRTDAECreateCoupling

! ************************************************************************** !

subroutine FlowRTDAELocalToGlobal(this,dm_ptr,ltog)
  !
  ! The local-to-global map of the entries of a PFLOTRAN DM's local vectors.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Discretization_module
  use UGDM_Pointer_module
  use Grid_module

  implicit none

  class(flow_rt_dae_model_type) :: this
  type(dm_ptr_type) :: dm_ptr
  ISLocalToGlobalMapping :: ltog

  PetscErrorCode :: ierr

  select case(this%rt%realization%discretization%itype)
    case(STRUCTURED_GRID)
      call DMGetLocalToGlobalMapping(dm_ptr%dm,ltog,ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      ltog = dm_ptr%ugdm%mapping_ltog
  end select

end subroutine FlowRTDAELocalToGlobal

! ************************************************************************** !

subroutine FlowRTDAESetupLinearSolver(this,ksp,J,pc_type,shift_inblocks)
  !
  ! Block forward substitution on ksp: a flow solve (-flow_ksp_*,
  ! -flow_pc_*; preonly LU if pc_type is LU, else ILU on one process and
  ! block Jacobi on several), then the transport model's preconditioner
  ! with pc_type.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use RT_DAE_Linear_module, only : RTDAEOperatorIsCondensed

  implicit none

  class(flow_rt_dae_model_type) :: this
  KSP :: ksp
  Mat :: J
  PCType :: pc_type
  PetscBool :: shift_inblocks

  PC :: pc, pc_sub
  Mat :: J_tt
  PetscMPIInt :: comm_size
  PetscErrorCode :: ierr

  if (.not.this%solvers_created) then
    call KSPCreate(this%option%mycomm,this%ksp_flow,ierr);CHKERRQ(ierr)
    call KSPSetOptionsPrefix(this%ksp_flow,'flow_',ierr);CHKERRQ(ierr)
    call KSPSetType(this%ksp_flow,KSPPREONLY,ierr);CHKERRQ(ierr)
    call KSPCreate(this%option%mycomm,this%ksp_tran,ierr);CHKERRQ(ierr)
    call KSPSetOptionsPrefix(this%ksp_tran,'tran_rt_',ierr);CHKERRQ(ierr)
    call KSPSetType(this%ksp_tran,KSPPREONLY,ierr);CHKERRQ(ierr)
    call VecCreateMPI(this%option%mycomm,this%nflow,PETSC_DETERMINE, &
                      this%work_flow,ierr);CHKERRQ(ierr)
    call VecCreateMPI(this%option%mycomm,this%ntran,PETSC_DETERMINE, &
                      this%work_tran,ierr);CHKERRQ(ierr)
    this%solvers_created = PETSC_TRUE
  endif
  ! the flow PC follows pc_type on every call (the model checks ask for LU)
  call KSPGetPC(this%ksp_flow,pc_sub,ierr);CHKERRQ(ierr)
  call MPI_Comm_size(this%option%mycomm,comm_size,ierr);CHKERRQ(ierr)
  if (pc_type == PCLU) then
    ! a direct deck solver: exact in every block
    call PCSetType(pc_sub,PCLU,ierr);CHKERRQ(ierr)
  else if (comm_size == 1) then
    call PCSetType(pc_sub,PCILU,ierr);CHKERRQ(ierr)
  else
    call PCSetType(pc_sub,PCBJACOBI,ierr);CHKERRQ(ierr)
  endif
  call KSPSetFromOptions(this%ksp_flow,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(J,ONE_INTEGER,ONE_INTEGER,J_tt,ierr);CHKERRQ(ierr)
  call KSPSetOperators(this%ksp_tran,J_tt,J_tt,ierr);CHKERRQ(ierr)
  if (RTDAEOperatorIsCondensed(J_tt)) then
    call this%rt%SetupLinearSolver(this%ksp_tran,J_tt,pc_type,shift_inblocks)
  else
    call KSPGetPC(this%ksp_tran,pc_sub,ierr);CHKERRQ(ierr)
    if (len_trim(pc_type) > 0) then
      call PCSetType(pc_sub,pc_type,ierr);CHKERRQ(ierr)
    endif
    call KSPSetFromOptions(this%ksp_tran,ierr);CHKERRQ(ierr)
  endif

  call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
  call PCSetType(pc,PCSHELL,ierr);CHKERRQ(ierr)
  call PCShellSetContext(pc,this%pc_ctx,ierr);CHKERRQ(ierr)
  call PCShellSetSetUp(pc,FlowRTDAEPCSetUp,ierr);CHKERRQ(ierr)
  call PCShellSetApply(pc,FlowRTDAEPCApply,ierr);CHKERRQ(ierr)
  call PCShellSetApplyTranspose(pc,FlowRTDAEPCApplyTranspose, &
                                ierr);CHKERRQ(ierr)
  call PCShellSetName(pc,'Flow RT DAE block forward substitution', &
                      ierr);CHKERRQ(ierr)

end subroutine FlowRTDAESetupLinearSolver

! ************************************************************************** !

subroutine FlowRTDAEPCContext(pc,model)
  !
  ! The model of a block preconditioner.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

  implicit none

  PC :: pc
  class(flow_rt_dae_model_type), pointer :: model

  type(flow_rt_pc_ctx_type), pointer :: ctx
  type(c_ptr) :: ctx_c
  PetscErrorCode :: ierr

  call PCShellGetContext(pc,ctx_c,ierr);CHKERRQ(ierr)
  call c_f_pointer(ctx_c,ctx)
  model => ctx%model

end subroutine FlowRTDAEPCContext

! ************************************************************************** !

subroutine FlowRTDAEPCSetUp(pc,ierr)
  !
  ! The sub-solvers take the diagonal blocks of the current operator.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  PetscErrorCode :: ierr

  class(flow_rt_dae_model_type), pointer :: model
  Mat :: A, B, J_ff, J_tt

  ! reset the failure reason (see FlowRTDAEPCApply)
  call PCSetFailedReason(pc,PC_NOERROR,ierr);CHKERRQ(ierr)
  call FlowRTDAEPCContext(pc,model)
  call PCGetOperators(pc,A,B,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(A,ZERO_INTEGER,ZERO_INTEGER,J_ff,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(A,ONE_INTEGER,ONE_INTEGER,J_tt,ierr);CHKERRQ(ierr)
  call KSPSetOperators(model%ksp_flow,J_ff,J_ff,ierr);CHKERRQ(ierr)
  call KSPSetOperators(model%ksp_tran,J_tt,J_tt,ierr);CHKERRQ(ierr)
  call KSPSetUp(model%ksp_flow,ierr);CHKERRQ(ierr)
  call KSPSetUp(model%ksp_tran,ierr);CHKERRQ(ierr)
  ierr = 0

end subroutine FlowRTDAEPCSetUp

! ************************************************************************** !

subroutine FlowRTDAEPCApply(pc,r,z,ierr)
  !
  ! z = P^-1 r: z_f = J_ff^-1 r_f; z_t = J_tt^-1 (r_t - J_tf z_f).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  Vec :: r, z
  PetscErrorCode :: ierr

  class(flow_rt_dae_model_type), pointer :: model
  Mat :: A, B, J_tf
  Vec :: r_f, r_t, z_f, z_t

  ! PETSc never clears a PCSHELL's failure reason; reset it per application
  call PCSetFailedReason(pc,PC_NOERROR,ierr);CHKERRQ(ierr)
  call FlowRTDAEPCContext(pc,model)
  call PCGetOperators(pc,A,B,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(A,ONE_INTEGER,ZERO_INTEGER,J_tf,ierr);CHKERRQ(ierr)
  call FlowRTDAEGetParts(model,r,r_f,r_t)
  call FlowRTDAEGetParts(model,z,z_f,z_t)
  call KSPSolve(model%ksp_flow,r_f,z_f,ierr);CHKERRQ(ierr)
  call FlowRTDAECheckSubSolve(pc,model%ksp_flow)
  call MatMult(J_tf,z_f,model%work_tran,ierr);CHKERRQ(ierr)
  call VecAYPX(model%work_tran,-1.d0,r_t,ierr);CHKERRQ(ierr)
  call KSPSolve(model%ksp_tran,model%work_tran,z_t,ierr);CHKERRQ(ierr)
  call FlowRTDAECheckSubSolve(pc,model%ksp_tran)
  call FlowRTDAERestoreParts(model,z,z_f,z_t)
  call FlowRTDAERestoreParts(model,r,r_f,r_t)
  ierr = 0

end subroutine FlowRTDAEPCApply

! ************************************************************************** !

subroutine FlowRTDAEPCApplyTranspose(pc,r,z,ierr)
  !
  ! z = P^-T r: z_t = J_tt^-T r_t; z_f = J_ff^-T (r_f - J_tf^T z_t).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  Vec :: r, z
  PetscErrorCode :: ierr

  class(flow_rt_dae_model_type), pointer :: model
  Mat :: A, B, J_tf
  Vec :: r_f, r_t, z_f, z_t

  ! reset the failure reason (see FlowRTDAEPCApply)
  call PCSetFailedReason(pc,PC_NOERROR,ierr);CHKERRQ(ierr)
  call FlowRTDAEPCContext(pc,model)
  call PCGetOperators(pc,A,B,ierr);CHKERRQ(ierr)
  call MatNestGetSubMat(A,ONE_INTEGER,ZERO_INTEGER,J_tf,ierr);CHKERRQ(ierr)
  call FlowRTDAEGetParts(model,r,r_f,r_t)
  call FlowRTDAEGetParts(model,z,z_f,z_t)
  call KSPSolveTranspose(model%ksp_tran,r_t,z_t,ierr);CHKERRQ(ierr)
  call FlowRTDAECheckSubSolve(pc,model%ksp_tran)
  call MatMultTranspose(J_tf,z_t,model%work_flow,ierr);CHKERRQ(ierr)
  call VecAYPX(model%work_flow,-1.d0,r_f,ierr);CHKERRQ(ierr)
  call KSPSolveTranspose(model%ksp_flow,model%work_flow,z_f,ierr);CHKERRQ(ierr)
  call FlowRTDAECheckSubSolve(pc,model%ksp_flow)
  call FlowRTDAERestoreParts(model,z,z_f,z_t)
  call FlowRTDAERestoreParts(model,r,r_f,r_t)
  ierr = 0

end subroutine FlowRTDAEPCApplyTranspose

! ************************************************************************** !

subroutine FlowRTDAECheckSubSolve(pc,ksp)
  !
  ! A failed sub-solve fails the preconditioner (the step is retried).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  PC :: pc
  KSP :: ksp

  KSPConvergedReason :: reason
  PCFailedReason :: pc_reason
  PC :: pc_sub
  PetscErrorCode :: ierr

  call KSPGetConvergedReason(ksp,reason,ierr);CHKERRQ(ierr)
  call KSPGetPC(ksp,pc_sub,ierr);CHKERRQ(ierr)
  call PCGetFailedReason(pc_sub,pc_reason,ierr);CHKERRQ(ierr)
  ! TODO: here and in RTDAESolveS the output stays finite, so an outer
  ! Krylov KSP never reads this failure while PREONLY does, and the reason
  ! outlives the solve when PCSetUp is skipped; KSPCheckSolve would give
  ! one policy. Only iterative sub-solvers are affected.
  if (reason%v < 0 .or. pc_reason%v /= PC_NOERROR%v) then
    call PCSetFailedReason(pc,PC_SUBPC_ERROR,ierr);CHKERRQ(ierr)
  endif

end subroutine FlowRTDAECheckSubSolve

! ************************************************************************** !

subroutine FlowRTDAEDefaultTolerances(this,atol,rtol)
  !
  ! The transport model's; the flow model's absolute tolerance is
  ! -ts_model_flow_atol (see RichardsDAESetTolerances).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol

  call this%rt%DefaultTolerances(atol,rtol)

end subroutine FlowRTDAEDefaultTolerances

! ************************************************************************** !

subroutine FlowRTDAESetTolerances(this,atol,rtol,vatol,vrtol)
  !
  ! atol and rtol apply to transport, rtol and -ts_model_flow_atol to flow.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol
  Vec :: vatol
  Vec :: vrtol

  Vec :: a_f, a_t, r_f, r_t

  call FlowRTDAEGetParts(this,vatol,a_f,a_t)
  call FlowRTDAEGetParts(this,vrtol,r_f,r_t)
  call this%flow%SetTolerances(this%flow_atol,rtol,a_f,r_f)
  call this%rt%SetTolerances(atol,rtol,a_t,r_t)
  call FlowRTDAERestoreParts(this,vrtol,r_f,r_t)
  call FlowRTDAERestoreParts(this,vatol,a_f,a_t)

end subroutine FlowRTDAESetTolerances

! ************************************************************************** !

subroutine FlowRTDAEConvergenceTest(this,time,X,dX,F,dt,kappa,converged)
  !
  ! Both models' tests; not converged while the last update held transport
  ! (its update test would pass on the zeroed update).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: X
  Vec :: dX
  Vec :: F
  PetscReal :: dt
  PetscReal :: kappa
  PetscBool :: converged

  Vec :: X_f, X_t, dX_f, dX_t, F_f, F_t
  PetscBool :: flow_converged, tran_converged
  PetscInt :: status

  call FlowRTDAEGetParts(this,X,X_f,X_t)
  call FlowRTDAEGetParts(this,dX,dX_f,dX_t)
  call FlowRTDAEGetParts(this,F,F_f,F_t)
  call this%flow%ConvergenceTest(time,X_f,dX_f,F_f,dt,kappa,flow_converged)
  ! transport's test reads the flow state of X (cached after F(X))
  call FlowRTDAESetFlowState(this,time,X_f,status)
  call this%rt%ConvergenceTest(time,X_t,dX_t,F_t,dt,kappa,tran_converged)
  converged = flow_converged .and. tran_converged .and. &
              status == TS_MODEL_OK .and. .not.this%transport_held
  call FlowRTDAERestoreParts(this,F,F_f,F_t)
  call FlowRTDAERestoreParts(this,dX,dX_f,dX_t)
  call FlowRTDAERestoreParts(this,X,X_f,X_t)

end subroutine FlowRTDAEConvergenceTest

! ************************************************************************** !

subroutine FlowRTDAEPublish(this,time,Y)
  !
  ! Flow first: transport output reads the flow fields too.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y

  Vec :: Y_f, Y_t
  PetscInt :: status

  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call this%flow%Publish(time,Y_f)
  call FlowRTDAESetFlowState(this,time,Y_f,status)
  call this%rt%Publish(time,Y_t)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)

end subroutine FlowRTDAEPublish

! ************************************************************************** !

subroutine FlowRTDAEPreCheck(this,Y,dY,changed)
  !
  ! Both models' Newton update strategies, and flow first
  ! (-ts_model_flow_first): while the flow update exceeds the flow's error
  ! tolerance (RichardsDAEUpdateSettled), the transport update is zeroed,
  ! so transport linearizes about a settled flow state, as in the
  ! sequential scheme, instead of extrapolating J_tf far. The converged
  ! state is unchanged.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  Vec :: Y
  Vec :: dY
  PetscBool :: changed

  Vec :: Y_f, Y_t, dY_f, dY_t
  PetscBool :: flow_changed
  PetscErrorCode :: ierr

  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call FlowRTDAEGetParts(this,dY,dY_f,dY_t)
  call this%flow%PreCheck(Y_f,dY_f,flow_changed)
  this%transport_held = PETSC_FALSE
  if (this%flow_first) then
    if (.not.RichardsDAEUpdateSettled(this%flow,dY_f)) then
      call VecZeroEntries(dY_t,ierr);CHKERRQ(ierr)
      flow_changed = PETSC_TRUE
      this%transport_held = PETSC_TRUE
      this%flow_first_iterations = this%flow_first_iterations + 1
    endif
  endif
  call this%rt%PreCheck(Y_t,dY_t,changed)
  changed = changed .or. flow_changed
  call FlowRTDAERestoreParts(this,dY,dY_f,dY_t)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)

end subroutine FlowRTDAEPreCheck

! ************************************************************************** !

function FlowRTDAEHoldsUpdate(this)
  !
  ! Whether the last Newton update held transport (FlowRTDAEPreCheck).
  !
  ! Author: Hansol Suh
  ! Date: 10/13/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscBool :: FlowRTDAEHoldsUpdate

  FlowRTDAEHoldsUpdate = this%transport_held

end function FlowRTDAEHoldsUpdate

! ************************************************************************** !

subroutine FlowRTDAEProjectInitialGuess(this,X,Xstart)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this
  Vec :: X
  Vec :: Xstart

  Vec :: X_f, X_t, S_f, S_t

  call FlowRTDAEGetParts(this,X,X_f,X_t)
  call FlowRTDAEGetParts(this,Xstart,S_f,S_t)
  call this%flow%ProjectInitialGuess(X_f,S_f)
  call this%rt%ProjectInitialGuess(X_t,S_t)
  call FlowRTDAERestoreParts(this,Xstart,S_f,S_t)
  call FlowRTDAERestoreParts(this,X,X_f,X_t)

end subroutine FlowRTDAEProjectInitialGuess

! ************************************************************************** !

function FlowRTDAENumEvents(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscInt :: FlowRTDAENumEvents

  FlowRTDAENumEvents = this%rt%NumEvents()

end function FlowRTDAENumEvents

! ************************************************************************** !

subroutine FlowRTDAEEventSetup(this,direction,vtol)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscInt :: direction(:)
  PetscReal :: vtol(:)

  call this%rt%EventSetup(direction,vtol)

end subroutine FlowRTDAEEventSetup

! ************************************************************************** !

subroutine FlowRTDAEEventIndicator(this,time,Y,fvalue)
  !
  ! The transport model's events, at the flow state of Y.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscReal :: fvalue(:)

  Vec :: Y_f, Y_t
  PetscInt :: status

  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call FlowRTDAESetFlowState(this,time,Y_f,status)
  call this%rt%EventIndicator(time,Y_t,fvalue)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)

end subroutine FlowRTDAEEventIndicator

! ************************************************************************** !

subroutine FlowRTDAEPostEvent(this,time,Y,events,changed)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscInt :: events(:)
  PetscBool :: changed

  Vec :: Y_f, Y_t
  PetscInt :: status

  call FlowRTDAEGetParts(this,Y,Y_f,Y_t)
  call FlowRTDAESetFlowState(this,time,Y_f,status)
  call this%rt%PostEvent(time,Y_t,events,changed)
  call FlowRTDAERestoreParts(this,Y,Y_f,Y_t)

end subroutine FlowRTDAEPostEvent

! ************************************************************************** !

subroutine FlowRTDAEPrintSummary(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this

  call this%flow%PrintSummary()
  call this%rt%PrintSummary()
  write(this%option%io_buffer,'(a,i8,a)') ' Flow RT DAE: transport held &
    &while flow settled in ',this%flow_first_iterations,' Newton iterations'
  call PrintMsg(this%option)

end subroutine FlowRTDAEPrintSummary

! ************************************************************************** !

subroutine FlowRTDAEForcingTimes(this,final_time,times)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this
  PetscReal :: final_time
  PetscReal, allocatable :: times(:)

  PetscReal, allocatable :: flow_times(:), tran_times(:)

  call this%flow%ForcingTimes(final_time,flow_times)
  call this%rt%ForcingTimes(final_time,tran_times)
  times = [flow_times,tran_times]

end subroutine FlowRTDAEForcingTimes

! ************************************************************************** !

subroutine FlowRTDAERecordForcing(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this

  call this%flow%RecordForcing()
  call this%rt%RecordForcing()

end subroutine FlowRTDAERecordForcing

! ************************************************************************** !

subroutine FlowRTDAEBeginReplay(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this

  call this%flow%BeginReplay()
  call this%rt%BeginReplay()

end subroutine FlowRTDAEBeginReplay

! ************************************************************************** !

subroutine FlowRTDAEEndReplay(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this

  call this%flow%EndReplay()
  call this%rt%EndReplay()

end subroutine FlowRTDAEEndReplay

! ************************************************************************** !

subroutine FlowRTDAEInvalidateCache(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(flow_rt_dae_model_type) :: this

  call this%flow%InvalidateCache()
  call this%rt%InvalidateCache()
  ! the flow state is recomputed at the next SetFlowState
  this%fs%id = -1

end subroutine FlowRTDAEInvalidateCache

! ************************************************************************** !

subroutine FlowRTDAEDestroy(this)
  !
  ! The sub-models belong to their process models.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: this

  PetscErrorCode :: ierr

  if (this%solvers_created) then
    call KSPDestroy(this%ksp_flow,ierr);CHKERRQ(ierr)
    call KSPDestroy(this%ksp_tran,ierr);CHKERRQ(ierr)
    call VecDestroy(this%work_flow,ierr);CHKERRQ(ierr)
    call VecDestroy(this%work_tran,ierr);CHKERRQ(ierr)
    this%solvers_created = PETSC_FALSE
  endif
  call ISDestroy(this%is_flow,ierr);CHKERRQ(ierr)
  call ISDestroy(this%is_tran,ierr);CHKERRQ(ierr)
  nullify(this%flow)
  nullify(this%rt)
  nullify(this%pc_ctx%model)

end subroutine FlowRTDAEDestroy

end module Flow_RT_DAE_Model_class
