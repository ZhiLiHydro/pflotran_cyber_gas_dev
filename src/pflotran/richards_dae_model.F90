module Richards_DAE_Model_class
  !
  ! Richards flow as a differential-algebraic equation model for PETSc TS
  ! (MODE RICHARDS_DAE).
  !
  ! State per grid cell, Y = [P; W]:
  !   P   liquid pressure [Pa]
  !   W   water inventory [kmol], the Richards accumulation phi*s*rho*V
  !       (kg, and kg/s below, with USE_MASS_DENSITY)
  ! Equations of an active cell:
  !   F_P = W - M(P)                                  (algebraic)
  !   F_W = dW/dt + sum_faces flux(P) - Q(P)          [kmol/s]
  ! with the fluxes and sources of the classic RICHARDS residual
  ! (Richards_Common_module), which backward Euler reproduces for constant
  ! or STEP forcing (LINEAR data are interpolated at the stage time); W is
  ! conserved for any integrator.
  ! Inactive cells carry F_P = P - P_frozen, F_W = dW/dt.
  !
  ! Evaluation uses private auxiliary variables. The realization supplies
  ! the configuration and, in BeginSegment, the forcing; Publish copies an
  ! accepted state back for output. A coupled transport model reads the
  ! flow state of the last Prepare (RichardsDAEFlowState).
  !
  ! Its connection loops mirror classic RICHARDS's around the shared
  ! Richards_Common_module kernels: the classic loops work on the
  ! realization's auxiliary variables.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use TS_Model_class
  use PFLOTRAN_Constants_module
  use Realization_Subsurface_class
  use Option_module
  use Richards_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use UGDM_Pointer_module
  use Richards_Aniso_module, only : aniso_richards_data_type
  use Richards_Common_module
  use TS_Connections_module

  implicit none

  private

  PetscInt, parameter, public :: RICHARDS_DAE_NDOF = 2

  ! algebraic-row absolute tolerance relative to atol; see
  ! RichardsDAEConvergenceTest
  PetscReal, parameter :: ALG_ABS_FACTOR = 1.d-12

  ! Boundary data and source rates of one waypoint interval
  ! [t_start, t_end], immutable once BeginSegment returns. Data of LINEAR
  ! conditions are interpolated between t_start and t_end; all other data
  ! hold their t_start values (STEP knots are waypoints).
  type :: richards_dae_segment_type
    PetscInt :: id = 0
    PetscReal :: t_start = 0.d0
    PetscReal :: t_end = 0.d0
    PetscBool :: transient = PETSC_FALSE
    ! the boundary couplers' flow_aux_real_var, (naux_bc,nconn_bc)
    PetscReal, allocatable :: bc_aux(:,:)
    PetscReal, allocatable :: bc_aux_end(:,:)
    ! source rates: kg/s for mass rates, m^3/s for volume rates, (nconn_ss)
    PetscReal, allocatable :: ss_rate(:)
    PetscReal, allocatable :: ss_rate_end(:)
  end type richards_dae_segment_type

  type, public, extends(ts_model_type) :: richards_dae_model_type
    class(realization_subsurface_type), pointer :: realization => null()
    type(option_type), pointer :: option => null()
    type(richards_parameter_type), pointer :: richards_parameter => null()
    PetscInt :: nlmax = 0
    PetscInt :: ngmax = 0
    ! state DM: the TS solution's (NDAEDOF) alone, a private one if coupled
    type(dm_ptr_type), pointer :: dm => null()
    PetscBool :: private_dm = PETSC_FALSE
    ! connection table, in the order of the patch's connection lists
    type(ts_connections_type) :: conn
    PetscInt, allocatable :: bc_itype(:)     ! pressure boundary type
    PetscInt :: naux_bc = 0
    PetscInt, allocatable :: ss_itype(:)     ! source type
    ! whether a connection's flow condition is LINEAR, and whether any is
    ! (on any process)
    PetscBool, allocatable :: bc_linear(:)
    PetscBool, allocatable :: ss_linear(:)
    PetscBool :: linear_forcing = PETSC_FALSE
    ! current segment (the record in use during a replay), the recorded
    ! segments (see TS_Model_class), and the forcing at forcing_time of
    ! segment forcing_segment
    type(richards_dae_segment_type) :: seg
    type(richards_dae_segment_type), allocatable :: records(:)
    PetscReal :: forcing_time = -MAX_DOUBLE
    PetscInt :: forcing_segment = -1
    PetscReal, allocatable :: bc_aux(:,:)
    PetscReal, allocatable :: ss_rate(:)
    ! private auxiliary variables of the cells and boundaries
    type(richards_auxvar_type), pointer :: rich(:) => null()
    type(global_auxvar_type), pointer :: gaux(:) => null()
    type(material_auxvar_type), pointer :: maux(:) => null()
    type(richards_auxvar_type), pointer :: rich_bc(:) => null()
    type(global_auxvar_type), pointer :: gaux_bc(:) => null()
    PetscBool, allocatable :: cell_ok(:)
    ! key of the flow state below (see RichardsDAEPrepare)
    PetscBool :: prepared = PETSC_FALSE
    PetscReal :: prepared_time = 0.d0
    PetscInt :: prepared_segment = -1
    PetscBool :: prepared_derivatives = PETSC_FALSE
    PetscInt :: prepared_status = TS_MODEL_OK
    ! number of Prepares that computed a new flow state
    PetscInt :: prepare_count = 0
    Vec :: Y_prepared
    Vec :: Y_loc
    ! flow state at the last Prepare: inventories M [kmol] of the local
    ! cells; Darcy fluxes v [m/s] and water fluxes [kmol/s] of the
    ! connections (internal: from int_up to int_dn; boundary: into the
    ! cell); injections q [kmol/s] and their volume rates [m^3/s]; and the
    ! derivatives with respect to the pressures (internal: of int_up, int_dn)
    PetscReal, allocatable :: acc(:), dacc_dp(:)
    PetscReal, allocatable :: v_int(:), flux_int(:)
    PetscReal, allocatable :: dv_int(:,:), dflux_int(:,:)   ! (2,nconn_int)
    PetscReal, allocatable :: v_bc(:), flux_bc(:), dv_bc(:), dflux_bc(:)
    PetscReal, allocatable :: q_ss(:), dq_ss(:), qvol_ss(:), dqvol_ss(:)
    ! inactive cells
    PetscReal, allocatable :: p_frozen(:)
    ! local error tolerances of the inventories, set by SetTolerances
    PetscReal :: tol_rtol = 0.d0
    PetscReal, allocatable :: tol_atol(:)
  contains
    procedure, public :: InitialState => RichardsDAEInitialState
    procedure, public :: BeginSegment => RichardsDAEBeginSegment
    procedure, public :: Evaluate => RichardsDAEEvaluate
    procedure, public :: Linearize => RichardsDAELinearize
    procedure, public :: SetTolerances => RichardsDAESetTolerances
    procedure, public :: DefaultTolerances => RichardsDAEDefaultTolerances
    procedure, public :: ConvergenceTest => RichardsDAEConvergenceTest
    procedure, public :: PreCheck => RichardsDAEPreCheck
    procedure, public :: Publish => RichardsDAEPublish
    procedure, public :: ForcingTimes => RichardsDAEForcingTimes
    procedure, public :: StoreRecord => RichardsDAEStoreRecord
    procedure, public :: RestoreRecord => RichardsDAERestoreRecord
    procedure, public :: InvalidateCache => RichardsDAEInvalidateCache
    procedure, public :: ResetToInitial => RichardsDAEResetToInitial
    procedure, public :: Destroy => RichardsDAEDestroy
  end type richards_dae_model_type

  public :: RichardsDAEModelCreate, &
            RichardsDAEPrepare, &
            RichardsDAECellInventory, &
            RichardsDAEAddLinearization, &
            RichardsDAECreateMatrix, &
            RichardsDAEFlowState, &
            RichardsDAEPermeabilityDerivative, &
            RichardsDAEUpdateSettled

contains

! ************************************************************************** !

function RichardsDAEModelCreate(realization,private_dm)
  !
  ! Creates the model and checks that the configuration is supported.
  ! private_dm: part of a coupled model, whose TS solution has its own layout.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Discretization_module
  use Grid_module
  use Patch_module

  implicit none

  class(realization_subsurface_type), pointer :: realization
  PetscBool :: private_dm
  class(richards_dae_model_type), pointer :: RichardsDAEModelCreate

  class(richards_dae_model_type), pointer :: this
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: disc
  PetscInt :: ghosted_id, i
  PetscErrorCode :: ierr

  allocate(this)
  this%name = 'Richards DAE'
  this%prefix = 'flow_'
  this%integrates_flow = PETSC_TRUE
  this%realization => realization
  option => realization%option
  patch => realization%patch
  disc => realization%discretization
  this%option => option
  this%richards_parameter => patch%aux%Richards%richards_parameter
  this%ndof = RICHARDS_DAE_NDOF
  this%nlmax = patch%grid%nlmax
  this%ngmax = patch%grid%ngmax

  call RichardsDAECheckCapabilities(this)

  this%private_dm = private_dm
  if (private_dm) then
    allocate(this%dm)
    PetscObjectNullify(this%dm%dm)
    nullify(this%dm%ugdm)
  else
    this%dm => disc%dm_ndaedof
    disc%dm_index_to_ndof(NDAEDOF) = this%ndof
  endif
  call DiscretizationCreateDM(disc,this%dm,this%ndof,disc%stencil_width, &
                              disc%stencil_type,option,'')
  call DMCreateLocalVector(this%dm%dm,this%Y_loc,ierr);CHKERRQ(ierr)
  call DMCreateGlobalVector(this%dm%dm,this%Y_prepared,ierr);CHKERRQ(ierr)

  call RichardsDAEConnectionTable(this)

  allocate(this%rich(this%ngmax))
  allocate(this%gaux(this%ngmax))
  allocate(this%maux(this%ngmax))
  do ghosted_id = 1, this%ngmax
    call RichardsAuxVarInit(this%rich(ghosted_id),option)
    call RichardsAuxVarCopy(patch%aux%Richards%auxvars(ghosted_id), &
                            this%rich(ghosted_id),option)
    call GlobalAuxVarInit(this%gaux(ghosted_id),option)
    call GlobalAuxVarCopy(patch%aux%Global%auxvars(ghosted_id), &
                          this%gaux(ghosted_id),option)
    call MaterialAuxVarInit(this%maux(ghosted_id),option)
  enddo
  call RichardsDAECopyMaterial(this)
  allocate(this%rich_bc(max(this%conn%nconn_bc,1)))
  allocate(this%gaux_bc(max(this%conn%nconn_bc,1)))
  do i = 1, this%conn%nconn_bc
    call RichardsAuxVarInit(this%rich_bc(i),option)
    call GlobalAuxVarInit(this%gaux_bc(i),option)
    call GlobalAuxVarCopy(patch%aux%Global%auxvars_bc(i),this%gaux_bc(i), &
                          option)
  enddo
  allocate(this%cell_ok(this%ngmax))
  this%cell_ok = PETSC_TRUE

  allocate(this%bc_aux(this%naux_bc,max(this%conn%nconn_bc,1)))
  allocate(this%ss_rate(max(this%conn%nconn_ss,1)))
  this%bc_aux = 0.d0
  this%ss_rate = 0.d0
  allocate(this%acc(this%nlmax),this%dacc_dp(this%nlmax))
  allocate(this%v_int(max(this%conn%nconn_int,1)))
  allocate(this%flux_int(max(this%conn%nconn_int,1)))
  allocate(this%dv_int(2,max(this%conn%nconn_int,1)))
  allocate(this%dflux_int(2,max(this%conn%nconn_int,1)))
  allocate(this%v_bc(max(this%conn%nconn_bc,1)))
  allocate(this%flux_bc(max(this%conn%nconn_bc,1)))
  allocate(this%dv_bc(max(this%conn%nconn_bc,1)))
  allocate(this%dflux_bc(max(this%conn%nconn_bc,1)))
  allocate(this%q_ss(max(this%conn%nconn_ss,1)))
  allocate(this%dq_ss(max(this%conn%nconn_ss,1)))
  allocate(this%qvol_ss(max(this%conn%nconn_ss,1)))
  allocate(this%dqvol_ss(max(this%conn%nconn_ss,1)))
  allocate(this%p_frozen(this%nlmax))
  this%acc = 0.d0
  this%dacc_dp = 0.d0
  this%v_int = 0.d0
  this%flux_int = 0.d0
  this%dv_int = 0.d0
  this%dflux_int = 0.d0
  this%v_bc = 0.d0
  this%flux_bc = 0.d0
  this%dv_bc = 0.d0
  this%dflux_bc = 0.d0
  this%q_ss = 0.d0
  this%dq_ss = 0.d0
  this%qvol_ss = 0.d0
  this%dqvol_ss = 0.d0
  this%p_frozen = 0.d0

  RichardsDAEModelCreate => this

end function RichardsDAEModelCreate

! ************************************************************************** !

subroutine RichardsDAECheckCapabilities(this)
  !
  ! Rejects every unsupported configuration, naming all of them.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module
  use Connection_module
  use Coupler_module
  use SrcSink_Sandbox_module, only : ss_sandbox_list

  implicit none

  class(richards_dae_model_type) :: this

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: coupler
  character(len=MAXSTRINGLENGTH) :: unsupported
  PetscBool :: aniso

  option => this%option
  patch => this%realization%patch
  unsupported = ''

  if (option%iflowmode /= RICHARDS_MODE) &
    call TSModelListAdd(unsupported,'a flow mode other than Richards')
  if (option%flow%inline_surface_flow) &
    call TSModelListAdd(unsupported,'inline surface flow')
  if (option%coupled_well .or. richards_well_coupling > 0) &
    call TSModelListAdd(unsupported,'well models')
  if (option%flow%density_depends_on_salinity) &
    call TSModelListAdd(unsupported,'salinity-dependent density')
  if (option%flow%numerical_derivatives) &
    call TSModelListAdd(unsupported,'NUMERICAL_JACOBIAN')
  if (associated(this%richards_parameter%geomech_parameter)) &
    call TSModelListAdd(unsupported,'geomechanics coupling')
  if (option%flow%update_flow_perm_from_pres .or. &
      option%flow%update_flow_perm_from_por) &
    call TSModelListAdd(unsupported,'permeability updates')
  if (associated(ss_sandbox_list)) &
    call TSModelListAdd(unsupported,'source/sink sandboxes')
  if (.not.PetscObjectIsNull(this%realization%field%flow_mass_transfer)) &
    call TSModelListAdd(unsupported,'flow mass transfer')
  if (associated(patch%prescribed_condition_list%first)) &
    call TSModelListAdd(unsupported,'PRESCRIBED_CONDITION')
  if (option%compute_mass_balance_new .or. option%flow%store_fluxes) &
    call TSModelListAdd(unsupported,'mass balance and integral flux output')
  if (option%restart_flag .or. associated(option%checkpoint)) &
    call TSModelListAdd(unsupported,'checkpoint/restart')
  if (option%flow%steady_state) &
    call TSModelListAdd(unsupported,'steady state flow')
  if (option%flow%only_vertical_flow) &
    call TSModelListAdd(unsupported,'ONLY_VERTICAL_FLOW')
  if (option%flow%quasi_3d) call TSModelListAdd(unsupported,'QUASI_3D')

  aniso = PETSC_FALSE
  cur_connection_set => patch%grid%internal_connection_set_list%first
  do while (associated(cur_connection_set))
    if (associated(cur_connection_set%aniso_geom)) aniso = PETSC_TRUE
    cur_connection_set => cur_connection_set%next
  enddo
  if (aniso) call TSModelListAdd(unsupported,'anisotropic flow stencils')
  coupler => patch%boundary_condition_list%first
  do while (associated(coupler))
    if (associated(coupler%flow_condition)) then
      select case(coupler%flow_condition%itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
             HYDROSTATIC_CONDUCTANCE_BC,HET_DIRICHLET_BC, &
             HET_HYDROSTATIC_SEEPAGE_BC,HET_HYDROSTATIC_CONDUCTANCE_BC, &
             NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC)
        case default
          call TSModelListAdd(unsupported,'the flow type of boundary &
                              &condition ' // trim(coupler%name))
      end select
    endif
    coupler => coupler%next
  enddo
  coupler => patch%source_sink_list%first
  do while (associated(coupler))
    if (associated(coupler%flow_condition)) then
      if (associated(coupler%flow_condition%well)) then
        call TSModelListAdd(unsupported,'the well of source/sink ' // &
                            trim(coupler%name))
      else
        select case(coupler%flow_condition%itype(1))
          case(MASS_RATE_SS,SCALED_MASS_RATE_SS,VOLUMETRIC_RATE_SS, &
               SCALED_VOLUMETRIC_RATE_SS,HET_VOL_RATE_SS,HET_MASS_RATE_SS)
          case default
            call TSModelListAdd(unsupported,'the flow type of source/sink ' &
                                // trim(coupler%name))
        end select
      endif
    endif
    coupler => coupler%next
  enddo

  if (len_trim(unsupported) > 0) then
    option%io_buffer = 'MODE RICHARDS_DAE does not yet support ' // &
      trim(unsupported) // '.'
    call PrintErrMsg(option)
  endif

end subroutine RichardsDAECheckCapabilities

! ************************************************************************** !

subroutine RichardsDAEConnectionTable(this)
  !
  ! The patch's connections (TSConnectionsBuild) and the flow condition
  ! type of each boundary and source/sink connection, and whether it is
  ! LINEAR.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Coupler_module

  implicit none

  class(richards_dae_model_type) :: this

  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: coupler
  PetscInt :: iconn, n
  PetscBool :: linear

  patch => this%realization%patch
  call TSConnectionsBuild(this%conn,patch)

  ! at least one entry, as the segment arrays (see BeginSegment)
  allocate(this%bc_itype(max(this%conn%nconn_bc,1)))
  allocate(this%bc_linear(max(this%conn%nconn_bc,1)))
  this%bc_itype = NULL_CONDITION
  this%bc_linear = PETSC_FALSE
  this%naux_bc = 1
  n = 0
  coupler => patch%boundary_condition_list%first
  do while (associated(coupler))
    if (associated(coupler%flow_aux_real_var)) &
      this%naux_bc = max(this%naux_bc,size(coupler%flow_aux_real_var,1))
    linear = RichardsDAELinearCondition(this,coupler%flow_condition)
    this%linear_forcing = this%linear_forcing .or. linear
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      this%bc_linear(n) = linear
      if (associated(coupler%flow_condition)) &
        this%bc_itype(n) = coupler%flow_condition%itype(RICHARDS_PRESSURE_DOF)
    enddo
    coupler => coupler%next
  enddo

  allocate(this%ss_itype(max(this%conn%nconn_ss,1)))
  allocate(this%ss_linear(max(this%conn%nconn_ss,1)))
  this%ss_itype = NULL_CONDITION
  this%ss_linear = PETSC_FALSE
  n = 0
  coupler => patch%source_sink_list%first
  do while (associated(coupler))
    linear = RichardsDAELinearCondition(this,coupler%flow_condition)
    this%linear_forcing = this%linear_forcing .or. linear
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      this%ss_linear(n) = linear
      if (associated(coupler%flow_condition)) &
        this%ss_itype(n) = coupler%flow_condition%itype(1)
    enddo
    coupler => coupler%next
  enddo

end subroutine RichardsDAEConnectionTable

! ************************************************************************** !

subroutine RichardsDAEBeginSegment(this,t_start,t_end,discontinuous)
  !
  ! Captures a new segment's boundary data and source rates at t_start and,
  ! for LINEAR conditions, at t_end, in that order: TimeStorageUpdate cannot
  ! go back in time. Restores option%time but leaves the realization's flow
  ! conditions and couplers at the last time read.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: t_start
  PetscReal :: t_end
  PetscBool :: discontinuous

  type(richards_dae_segment_type) :: new
  PetscReal :: saved_time
  PetscInt :: n
  PetscErrorCode :: ierr

  new%id = this%seg%id + 1
  new%t_start = t_start
  new%t_end = t_end
  allocate(new%bc_aux(this%naux_bc,max(this%conn%nconn_bc,1)))
  allocate(new%ss_rate(max(this%conn%nconn_ss,1)))

  saved_time = this%option%time
  new%transient = this%linear_forcing .and. t_end > t_start
  call RichardsDAECaptureForcing(this,t_start,new%bc_aux,new%ss_rate)
  new%bc_aux_end = new%bc_aux
  new%ss_rate_end = new%ss_rate
  if (new%transient) then
    call RichardsDAECaptureForcing(this,t_end,new%bc_aux_end, &
                                   new%ss_rate_end)
    ! at t_end, STEP data may already be the next piece
    do n = 1, this%conn%nconn_bc
      if (.not.this%bc_linear(n)) new%bc_aux_end(:,n) = new%bc_aux(:,n)
    enddo
    where(.not.this%ss_linear) new%ss_rate_end = new%ss_rate
  endif
  this%option%time = saved_time
  discontinuous = PETSC_TRUE
  if (this%seg%id > 0) then
    discontinuous = any(new%bc_aux /= this%seg%bc_aux_end) .or. &
                    any(new%ss_rate /= this%seg%ss_rate_end)
    call MPI_Allreduce(MPI_IN_PLACE,discontinuous,ONE_INTEGER_MPI, &
                       MPI_C_BOOL,MPI_LOR,this%option%mycomm, &
                       ierr);CHKERRQ(ierr)
  endif

  this%seg = new
  this%segment_id = new%id
  call this%RecordSegment(t_end)
  this%forcing_segment = -1
  this%forcing_time = -MAX_DOUBLE

end subroutine RichardsDAEBeginSegment

! ************************************************************************** !

function RichardsDAELinearCondition(this,condition)
  !
  ! Whether a flow condition's time-varying datasets (values, datums,
  ! gradients) are LINEAR. Mixing LINEAR and STEP is rejected: a segment
  ! interpolates or holds all of a connection's data.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Condition_module
  use Dataset_Base_class

  implicit none

  class(richards_dae_model_type) :: this
  type(flow_condition_type), pointer :: condition
  PetscBool :: RichardsDAELinearCondition

  PetscBool :: linear, step
  PetscInt :: i

  RichardsDAELinearCondition = PETSC_FALSE
  if (.not.associated(condition)) return
  linear = PETSC_FALSE
  step = PETSC_FALSE
  call Classify(condition%datum)
  do i = 1, condition%num_sub_conditions
    if (.not.associated(condition%sub_condition_ptr(i)%ptr)) cycle
    call Classify(condition%sub_condition_ptr(i)%ptr%dataset)
    call Classify(condition%sub_condition_ptr(i)%ptr%gradient)
  enddo
  if (linear .and. step) then
    call PrintErrMsg(this%option,'MODE RICHARDS_DAE does not support &
      &FLOW_CONDITION "' // trim(condition%name) // '": its datasets mix &
      &LINEAR and STEP time interpolation.')
  endif
  RichardsDAELinearCondition = linear

contains

  subroutine Classify(dataset)
    class(dataset_base_type), pointer :: dataset
    if (.not.associated(dataset)) return
    if (.not.associated(dataset%time_storage)) return
    if (dataset%time_storage%max_time_index <= 1) return
    if (dataset%time_storage%time_interpolation_method == &
        INTERPOLATION_LINEAR) then
      linear = PETSC_TRUE
    else
      step = PETSC_TRUE
    endif
  end subroutine Classify

end function RichardsDAELinearCondition

! ************************************************************************** !

subroutine RichardsDAEForcingTimes(this,final_time,times)
  !
  ! Knots up to final_time of the transient datasets of the boundary and
  ! source flow conditions.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Coupler_module
  use Condition_module
  use Dataset_Base_class
  use Time_Storage_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: final_time
  PetscReal, allocatable :: times(:)

  type(coupler_type), pointer :: coupler
  type(flow_condition_type), pointer :: condition
  PetscInt :: i, ilist

  allocate(times(0))
  do ilist = 1, 2
    if (ilist == 1) then
      coupler => this%realization%patch%boundary_condition_list%first
    else
      coupler => this%realization%patch%source_sink_list%first
    endif
    do while (associated(coupler))
      condition => coupler%flow_condition
      if (associated(condition)) then
        call Append(condition%datum)
        do i = 1, condition%num_sub_conditions
          if (.not.associated(condition%sub_condition_ptr(i)%ptr)) cycle
          call Append(condition%sub_condition_ptr(i)%ptr%dataset)
          call Append(condition%sub_condition_ptr(i)%ptr%gradient)
        enddo
      endif
      coupler => coupler%next
    enddo
  enddo

contains

  subroutine Append(dataset)
    class(dataset_base_type), pointer :: dataset
    PetscReal, pointer :: knots(:)
    if (.not.associated(dataset)) return
    if (.not.associated(dataset%time_storage)) return
    if (dataset%time_storage%max_time_index <= 1) return
    call TimeStorageGetTimes(dataset%time_storage,this%option,final_time, &
                             knots)
    if (.not.associated(knots)) return
    times = [times,knots]
    deallocate(knots)
  end subroutine Append

end subroutine RichardsDAEForcingTimes

! ************************************************************************** !

subroutine RichardsDAECaptureForcing(this,time,bc_aux,ss_rate)
  !
  ! The couplers' boundary data and the source rates at time.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Condition_module, only : FlowConditionUpdate
  use Coupler_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  PetscReal :: bc_aux(:,:)
  PetscReal :: ss_rate(:)

  type(coupler_type), pointer :: coupler
  PetscInt :: iconn, n, naux
  PetscReal :: rate, scale

  this%option%time = time
  call FlowConditionUpdate(this%realization%flow_conditions,this%option)
  call RealizUpdateAllCouplerAuxVars(this%realization,PETSC_TRUE)

  bc_aux = 0.d0
  n = 0
  coupler => this%realization%patch%boundary_condition_list%first
  do while (associated(coupler))
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      if (.not.associated(coupler%flow_aux_real_var)) cycle
      naux = size(coupler%flow_aux_real_var,1)
      bc_aux(1:naux,n) = coupler%flow_aux_real_var(:,iconn)
    enddo
    coupler => coupler%next
  enddo

  ss_rate = 0.d0
  n = 0
  coupler => this%realization%patch%source_sink_list%first
  do while (associated(coupler))
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      if (.not.associated(coupler%flow_condition)) cycle
      call RichardsSrcSinkRate(coupler,iconn,rate,scale)
      ss_rate(n) = rate*scale
    enddo
    coupler => coupler%next
  enddo

end subroutine RichardsDAECaptureForcing

! ************************************************************************** !

subroutine RichardsDAEForcingAt(this,time)
  !
  ! The boundary data and source rates at time; recomputed only when the
  ! segment or, for linearly varying data, the time changes.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time

  PetscReal :: w

  if (this%forcing_segment == this%seg%id) then
    if (.not.this%seg%transient) return
    if (time == this%forcing_time) return
  endif
  this%forcing_segment = this%seg%id
  this%forcing_time = time
  w = 0.d0
  if (this%seg%transient) then
    w = (time-this%seg%t_start)/(this%seg%t_end-this%seg%t_start)
    w = max(0.d0,min(1.d0,w))
  endif
  this%bc_aux = (1.d0-w)*this%seg%bc_aux + w*this%seg%bc_aux_end
  this%ss_rate = (1.d0-w)*this%seg%ss_rate + w*this%seg%ss_rate_end

end subroutine RichardsDAEForcingAt

! ************************************************************************** !

subroutine RichardsDAECellAuxVars(this,ghosted_id,pressure,status)
  !
  ! The private auxiliary variables of a cell at a pressure, as the
  ! classic accumulation computes them (porosity updated for soil
  ! compressibility).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: pressure
  PetscInt :: status

  type(patch_type), pointer :: patch
  PetscReal :: x(1)
  PetscInt :: iflag

  patch => this%realization%patch
  status = TS_MODEL_OK
  if (.not.(abs(pressure) < huge(pressure))) then
    status = TS_MODEL_DOMAIN_ERROR
    return
  endif
  x(1) = pressure
  iflag = this%option%iflag
  this%option%iflag = RICHARDS_UPDATE_FOR_ACCUM
  call RichardsAuxVarCompute(x,this%rich(ghosted_id),this%gaux(ghosted_id), &
                             this%maux(ghosted_id),this%richards_parameter, &
                             patch%characteristic_curves_array( &
                               patch%cc_id(ghosted_id))%ptr, &
                             patch%grid%nG2A(ghosted_id),PETSC_TRUE, &
                             this%option)
  this%option%iflag = iflag
  if (.not.(abs(this%gaux(ghosted_id)%den(1)) < huge(1.d0)) .or. &
      this%gaux(ghosted_id)%den(1) <= 0.d0) status = TS_MODEL_DOMAIN_ERROR

end subroutine RichardsDAECellAuxVars

! ************************************************************************** !

subroutine RichardsDAEBoundaryAuxVars(this,iconn,p_cell)
  !
  ! The private auxiliary variables of boundary connection iconn: at the
  ! boundary pressure for prescribed-pressure types, at the cell pressure
  ! p_cell otherwise (as RichardsUpdateAuxVarsPatch).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscInt :: iconn
  PetscReal :: p_cell

  type(patch_type), pointer :: patch
  PetscReal :: x(1)
  PetscInt :: iflag, ghosted_id

  patch => this%realization%patch
  ghosted_id = this%conn%bc_ghost(iconn)
  if (RichardsBCPressureFollowsCell(this%bc_itype(iconn))) then
    x(1) = p_cell
  else
    x(1) = this%bc_aux(RICHARDS_PRESSURE_DOF,iconn)
  endif
  iflag = this%option%iflag
  this%option%iflag = RICHARDS_UPDATE_FOR_BOUNDARY
  call RichardsAuxVarCompute(x,this%rich_bc(iconn),this%gaux_bc(iconn), &
                             this%maux(ghosted_id),this%richards_parameter, &
                             patch%characteristic_curves_array( &
                               patch%cc_id(ghosted_id))%ptr, &
                             -patch%grid%nG2A(ghosted_id),PETSC_FALSE, &
                             this%option)
  this%option%iflag = iflag

end subroutine RichardsDAEBoundaryAuxVars

! ************************************************************************** !

function RichardsDAECellInventory(this,ghosted_id)
  !
  ! The water inventory M [kmol] of a cell at its private auxiliary
  ! variables: the classic accumulation times the step size.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: RichardsDAECellInventory

  RichardsDAECellInventory = &
    RichardsWaterMass(this%rich(ghosted_id),this%gaux(ghosted_id), &
                      this%maux(ghosted_id)%volume)

end function RichardsDAECellInventory

! ************************************************************************** !

subroutine RichardsDAEPrepare(this,time,Y,derivatives,status)
  !
  ! Forcing, ghosted state, auxiliary variables of active cells and
  ! boundaries, and inventories, fluxes and sources (with pressure
  ! derivatives if requested) at (time, Y); reused while neither changed
  ! within a forcing segment. Collective.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscBool :: derivatives
  PetscInt :: status

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: yl_p(:)
  PetscReal :: Res(1), Jup(1,1), Jdn(1,1), v, rho
  PetscInt :: ghosted_id, local_id, iconn, up, dn, cell_status
  PetscBool :: same
  type(aniso_richards_data_type), pointer :: no_aniso
  PetscErrorCode :: ierr

  call this%ReplayAt(time)
  if (this%prepared .and. time == this%prepared_time .and. &
      this%seg%id == this%prepared_segment .and. &
      (this%prepared_derivatives .or. .not.derivatives)) then
    call VecEqual(Y,this%Y_prepared,same,ierr);CHKERRQ(ierr)
    if (same) then
      status = this%prepared_status
      return
    endif
  endif

  patch => this%realization%patch
  grid => patch%grid
  nullify(no_aniso)
  status = TS_MODEL_OK

  call RichardsDAEForcingAt(this,time)
  call DMGlobalToLocalBegin(this%dm%dm,Y,INSERT_VALUES,this%Y_loc, &
                            ierr);CHKERRQ(ierr)
  call DMGlobalToLocalEnd(this%dm%dm,Y,INSERT_VALUES,this%Y_loc, &
                          ierr);CHKERRQ(ierr)
  call VecGetArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)

  ! cells
  do ghosted_id = 1, this%ngmax
    this%cell_ok(ghosted_id) = PETSC_TRUE
    if (grid%nG2L(ghosted_id) < 0) cycle
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsDAECellAuxVars(this,ghosted_id, &
                                yl_p((ghosted_id-1)*this%ndof+1),cell_status)
    if (cell_status /= TS_MODEL_OK) then
      this%cell_ok(ghosted_id) = PETSC_FALSE
      if (grid%nG2L(ghosted_id) > 0) status = TS_MODEL_DOMAIN_ERROR
      cycle
    endif
    local_id = grid%nG2L(ghosted_id)
    if (local_id <= 0) cycle
    this%acc(local_id) = RichardsDAECellInventory(this,ghosted_id)
    this%dacc_dp(local_id) = &
      RichardsWaterMassDerivative(this%rich(ghosted_id), &
                                  this%gaux(ghosted_id), &
                                  this%maux(ghosted_id)%volume)
  enddo

  ! internal connections: flux from int_up to int_dn
  do iconn = 1, this%conn%nconn_int
    up = this%conn%int_up(iconn)
    dn = this%conn%int_dn(iconn)
    this%v_int(iconn) = 0.d0
    this%flux_int(iconn) = 0.d0
    this%dv_int(:,iconn) = 0.d0
    this%dflux_int(:,iconn) = 0.d0
    if (patch%imat(up) <= 0 .or. patch%imat(dn) <= 0) cycle
    if (.not.(this%cell_ok(up) .and. this%cell_ok(dn))) cycle
    call RichardsFlux(this%rich(up),this%gaux(up),this%maux(up), &
                      this%rich(dn),this%gaux(dn),this%maux(dn), &
                      this%conn%int_area(iconn),this%conn%int_dist(:,iconn), &
                      no_aniso,this%option,v,Res)
    this%v_int(iconn) = v
    this%flux_int(iconn) = Res(1)
    if (.not.derivatives) cycle
    call RichardsFluxDerivative(this%rich(up),this%gaux(up),this%maux(up), &
                                this%rich(dn),this%gaux(dn),this%maux(dn), &
                                this%richards_parameter, &
                                this%conn%int_area(iconn), &
                                this%conn%int_dist(:,iconn), &
                                no_aniso,this%option, &
                                patch%characteristic_curves_array( &
                                  patch%cc_id(up))%ptr, &
                                patch%characteristic_curves_array( &
                                  patch%cc_id(dn))%ptr,Jup,Jdn, &
                                this%dv_int(1,iconn),this%dv_int(2,iconn))
    this%dflux_int(1,iconn) = Jup(1,1)
    this%dflux_int(2,iconn) = Jdn(1,1)
  enddo

  ! boundary connections: flux into the cell
  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    this%v_bc(iconn) = 0.d0
    this%flux_bc(iconn) = 0.d0
    this%dv_bc(iconn) = 0.d0
    this%dflux_bc(iconn) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    if (this%bc_itype(iconn) == NULL_CONDITION) cycle
    call RichardsDAEBoundaryAuxVars(this,iconn, &
                                    yl_p((ghosted_id-1)*this%ndof+1))
    call RichardsBCFlux([this%bc_itype(iconn)],this%bc_aux(:,iconn), &
                        this%rich_bc(iconn),this%gaux_bc(iconn), &
                        this%rich(ghosted_id),this%gaux(ghosted_id), &
                        this%maux(ghosted_id),this%conn%bc_area(iconn), &
                        this%conn%bc_dist(:,iconn),this%option,v,Res)
    this%v_bc(iconn) = v
    this%flux_bc(iconn) = Res(1)
    if (.not.derivatives) cycle
    call RichardsBCFluxDerivative([this%bc_itype(iconn)],this%bc_aux(:,iconn), &
                                  this%rich_bc(iconn),this%gaux_bc(iconn), &
                                  this%rich(ghosted_id),this%gaux(ghosted_id), &
                                  this%maux(ghosted_id), &
                                  this%richards_parameter, &
                                  this%conn%bc_area(iconn), &
                                  this%conn%bc_dist(:,iconn), &
                                  this%option, &
                                  patch%characteristic_curves_array( &
                                    patch%cc_id(ghosted_id))%ptr,Jdn, &
                                  this%dv_bc(iconn))
    this%dflux_bc(iconn) = Jdn(1,1)
  enddo

  ! sources: injection [kmol/s] and its volume rate [m^3/s] at the cell
  ! density (the transport source rate, as ss_flow_vol_fluxes)
  do iconn = 1, this%conn%nconn_ss
    ghosted_id = this%conn%ss_ghost(iconn)
    this%q_ss(iconn) = 0.d0
    this%dq_ss(iconn) = 0.d0
    this%qvol_ss(iconn) = 0.d0
    this%dqvol_ss(iconn) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    rho = this%gaux(ghosted_id)%den(1)
    ! (ss_rate holds a SCALED type's factor)
    call RichardsSrcSinkMolarRate(this%ss_itype(iconn),this%ss_rate(iconn), &
                                  1.d0,rho,this%q_ss(iconn), &
                                  this%rich(ghosted_id)%dden_dp, &
                                  this%dq_ss(iconn))
    ! as ss_flow_vol_fluxes
    this%qvol_ss(iconn) = this%q_ss(iconn)/rho
    this%dqvol_ss(iconn) = (this%dq_ss(iconn) - &
                            this%qvol_ss(iconn)* &
                            this%rich(ghosted_id)%dden_dp)/rho
  enddo

  call VecRestoreArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)

  this%prepare_count = this%prepare_count + 1
  call VecCopy(Y,this%Y_prepared,ierr);CHKERRQ(ierr)
  this%prepared = PETSC_TRUE
  this%prepared_time = time
  this%prepared_segment = this%seg%id
  this%prepared_derivatives = derivatives
  this%prepared_status = status

end subroutine RichardsDAEPrepare

! ************************************************************************** !

subroutine RichardsDAEEvaluate(this,time,Y,Ydot,F,status)
  !
  ! F(t,Y,Ydot) of the whole grid.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  Vec :: F
  PetscInt :: status

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: y_p(:), ydot_p(:), f_p(:)
  PetscInt :: local_id, ghosted_id, iconn, lup, ldn, off
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  call RichardsDAEPrepare(this,time,Y,PETSC_TRUE,status)

  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Ydot,ydot_p,ierr);CHKERRQ(ierr)
  call VecGetArray(F,f_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    off = (local_id-1)*this%ndof
    if (patch%imat(ghosted_id) <= 0) then
      f_p(off+1) = y_p(off+1) - this%p_frozen(local_id)
      f_p(off+2) = ydot_p(off+2)
      cycle
    endif
    if (.not.this%cell_ok(ghosted_id)) then
      f_p(off+1:off+2) = ieee_value(1.d0,ieee_quiet_nan)
      cycle
    endif
    f_p(off+1) = y_p(off+2) - this%acc(local_id)
    f_p(off+2) = ydot_p(off+2)
  enddo
  do iconn = 1, this%conn%nconn_int
    lup = grid%nG2L(this%conn%int_up(iconn))
    ldn = grid%nG2L(this%conn%int_dn(iconn))
    if (lup > 0) f_p((lup-1)*this%ndof+2) = f_p((lup-1)*this%ndof+2) + &
                                           this%flux_int(iconn)
    if (ldn > 0) f_p((ldn-1)*this%ndof+2) = f_p((ldn-1)*this%ndof+2) - &
                                           this%flux_int(iconn)
  enddo
  do iconn = 1, this%conn%nconn_bc
    off = (this%conn%bc_local(iconn)-1)*this%ndof
    f_p(off+2) = f_p(off+2) - this%flux_bc(iconn)
  enddo
  do iconn = 1, this%conn%nconn_ss
    off = (this%conn%ss_local(iconn)-1)*this%ndof
    f_p(off+2) = f_p(off+2) - this%q_ss(iconn)
  enddo
  call VecRestoreArray(F,f_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Ydot,ydot_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)

end subroutine RichardsDAEEvaluate

! ************************************************************************** !

subroutine RichardsDAELinearize(this,time,Y,Ydot,shift,A,B,status)
  !
  ! J = dF/dY + shift*dF/dYdot at (time,Y) into the block matrix A (block
  ! size ndof, local numbering of the model's DM). Valid at any shift.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  PetscReal :: shift
  Mat :: A
  Mat :: B
  PetscInt :: status

  PetscErrorCode :: ierr

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  call RichardsDAEPrepare(this,time,Y,PETSC_TRUE,status)
  call RichardsDAEAddLinearization(this,shift,A,0)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine RichardsDAELinearize

! ************************************************************************** !

subroutine RichardsDAEAddLinearization(this,shift,A,offset)
  !
  ! Adds the Jacobian at the last Prepare (with derivatives) to A, whose
  ! local rows and columns number the model's entries from offset in blocks
  ! of ndof per ghosted cell (MatSetValuesBlockedLocal).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module
  use Petsc_Utility_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: shift
  Mat :: A
  PetscInt :: offset

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal :: blk(RICHARDS_DAE_NDOF,RICHARDS_DAE_NDOF)
  PetscInt :: local_id, ghosted_id, iconn, up, dn
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    blk = 0.d0
    if (patch%imat(ghosted_id) <= 0) then
      blk(1,1) = 1.d0
      blk(2,2) = shift
    else if (.not.this%cell_ok(ghosted_id)) then
      blk(1,1) = 1.d0
      blk(2,2) = 1.d0
    else
      ! F_P = W - M(P), F_W = dW/dt + ...
      blk(1,1) = -this%dacc_dp(local_id)
      blk(1,2) = 1.d0
      blk(2,2) = shift
    endif
    call Add(ghosted_id,ghosted_id,blk)
  enddo
  do iconn = 1, this%conn%nconn_int
    up = this%conn%int_up(iconn)
    dn = this%conn%int_dn(iconn)
    if (patch%imat(up) <= 0 .or. patch%imat(dn) <= 0) cycle
    if (grid%nG2L(up) > 0) then
      call AddFlux(up,up,this%dflux_int(1,iconn))
      call AddFlux(up,dn,this%dflux_int(2,iconn))
    endif
    if (grid%nG2L(dn) > 0) then
      call AddFlux(dn,up,-this%dflux_int(1,iconn))
      call AddFlux(dn,dn,-this%dflux_int(2,iconn))
    endif
  enddo
  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    call AddFlux(ghosted_id,ghosted_id,-this%dflux_bc(iconn))
  enddo
  do iconn = 1, this%conn%nconn_ss
    ghosted_id = this%conn%ss_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    call AddFlux(ghosted_id,ghosted_id,-this%dq_ss(iconn))
  enddo

contains

  subroutine AddFlux(row,col,value)
    ! row W of cell row, column P of cell col
    PetscInt :: row, col
    PetscReal :: value
    blk = 0.d0
    blk(2,1) = value
    call Add(row,col,blk)
  end subroutine AddFlux

  subroutine Add(row,col,block)
    PetscInt :: row, col
    PetscReal :: block(:,:)
    call PUMSetValuesBlockedLocal(A,1,offset+row-1,1,offset+col-1,block, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  end subroutine Add

end subroutine RichardsDAEAddLinearization

! ************************************************************************** !

subroutine RichardsDAECreateMatrix(this,A)
  !
  ! A block matrix of ndof per cell on the model's DM, with the stencil of
  ! the grid (as DiscretizationCreateMatrix, for a private DM too).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Grid_module
  use Grid_Unstructured_Aux_module, only : UGridDMCreateMatrix
  use Discretization_module

  implicit none

  class(richards_dae_model_type) :: this
  Mat :: A

  ! a MatType dummy reads 80 characters; the constant MATBAIJ has 4
  MatType :: mat_type = MATBAIJ
  PetscErrorCode :: ierr

  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID)
      call DMSetMatType(this%dm%dm,MATBAIJ,ierr);CHKERRQ(ierr)
      call DMCreateMatrix(this%dm%dm,A,ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      call UGridDMCreateMatrix(this%realization%discretization%grid% &
                                 unstructured_grid, &
                               this%dm%ugdm,mat_type,A,this%option)
  end select
  call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE,ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_ROW_ORIENTED,PETSC_FALSE,ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE,ierr);CHKERRQ(ierr)

end subroutine RichardsDAECreateMatrix

! ************************************************************************** !

subroutine RichardsDAEFlowState(this,time,Y,fs,status)
  !
  ! The flow state at (time,Y) that a coupled transport model reads, with
  ! its pressure derivatives (see TS_Flow_State_module). fs%id changes when
  ! the state does. Collective.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use TS_Flow_State_module
  use Patch_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  type(ts_flow_state_type) :: fs
  PetscInt :: status

  type(patch_type), pointer :: patch
  PetscInt :: ghosted_id, iconn

  patch => this%realization%patch
  call RichardsDAEPrepare(this,time,Y,PETSC_TRUE,status)
  if (fs%id == this%prepare_count) return
  fs%id = this%prepare_count
  do ghosted_id = 1, this%ngmax
    fs%sat(ghosted_id) = this%gaux(ghosted_id)%sat(1)
    fs%den(ghosted_id) = this%gaux(ghosted_id)%den(1)
    fs%den_kg(ghosted_id) = this%gaux(ghosted_id)%den_kg(1)
    fs%pres(ghosted_id) = this%gaux(ghosted_id)%pres(1)
    fs%por(ghosted_id) = this%maux(ghosted_id)%porosity
    fs%dsat_dp(ghosted_id) = 0.d0
    fs%dden_kg_dp(ghosted_id) = 0.d0
    fs%dpor_dp(ghosted_id) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    if (patch%grid%nG2L(ghosted_id) < 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    fs%por(ghosted_id) = this%rich(ghosted_id)%effective_porosity
    fs%dsat_dp(ghosted_id) = this%rich(ghosted_id)%dsat_dp
    fs%dpor_dp(ghosted_id) = this%rich(ghosted_id)%dpor_dp
    ! den and den_kg are proportional
    fs%dden_kg_dp(ghosted_id) = this%gaux(ghosted_id)%den_kg(1)* &
                                this%rich(ghosted_id)%dden_dp/ &
                                this%gaux(ghosted_id)%den(1)
  enddo
  if (this%conn%nconn_int > 0) then
    fs%v_int(1:this%conn%nconn_int) = this%v_int(1:this%conn%nconn_int)
    fs%dv_int(:,1:this%conn%nconn_int) = this%dv_int(:,1:this%conn%nconn_int)
  endif
  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    fs%v_bc(iconn) = this%v_bc(iconn)
    fs%dv_bc(iconn) = this%dv_bc(iconn)
    fs%sat_bc(iconn) = this%gaux_bc(iconn)%sat(1)
    fs%den_bc(iconn) = this%gaux_bc(iconn)%den(1)
    fs%den_kg_bc(iconn) = this%gaux_bc(iconn)%den_kg(1)
    fs%pres_bc(iconn) = this%gaux_bc(iconn)%pres(1)
    fs%dsat_bc_dp(iconn) = 0.d0
    fs%dden_kg_bc_dp(iconn) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    if (this%bc_itype(iconn) == NULL_CONDITION) cycle
    if (.not.RichardsBCPressureFollowsCell(this%bc_itype(iconn))) cycle
    fs%dsat_bc_dp(iconn) = this%rich_bc(iconn)%dsat_dp
    fs%dden_kg_bc_dp(iconn) = this%gaux_bc(iconn)%den_kg(1)* &
                              this%rich_bc(iconn)%dden_dp/ &
                              this%gaux_bc(iconn)%den(1)
  enddo
  if (this%conn%nconn_ss > 0) then
    fs%qvol_ss(1:this%conn%nconn_ss) = this%qvol_ss(1:this%conn%nconn_ss)
    fs%dqvol_ss(1:this%conn%nconn_ss) = this%dqvol_ss(1:this%conn%nconn_ss)
  endif

end subroutine RichardsDAEFlowState

! ************************************************************************** !

subroutine RichardsDAEPermeabilityDerivative(this,material_id,k,F,dv_int, &
                                             dv_bc)
  !
  ! d/dk at the last Prepare, k the isotropic permeability of material
  ! material_id (its cells' diagonal entries scale with it): of the water
  ! rows, added to F, and of the Darcy fluxes, into dv_int and dv_bc. With
  ! S = d_up k_dn + d_dn k_up, a flux is linear in k_up k_dn/S (RichardsFlux's
  ! harmonic average), so
  !   dflux/dk = flux/k (w_up + w_dn),  w_up = d_up k_dn/S, w_dn = d_dn k_up/S
  ! counting w_up (w_dn) if the up (down) cell is of the material. A
  ! boundary flux is linear in the cell's k unless Neumann or conductance.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module
  use Connection_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscInt :: material_id
  PetscReal :: k
  Vec :: F
  PetscReal :: dv_int(:)
  PetscReal :: dv_bc(:)

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: f_p(:)
  PetscReal :: dd_up, dd_dn, dist_gravity, upweight, perm_up, perm_dn, w
  PetscInt :: iconn, up, dn, lup, ldn, ghosted_id
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  dv_int = 0.d0
  dv_bc = 0.d0
  call VecGetArray(F,f_p,ierr);CHKERRQ(ierr)
  do iconn = 1, this%conn%nconn_int
    up = this%conn%int_up(iconn)
    dn = this%conn%int_dn(iconn)
    if (patch%imat(up) <= 0 .or. patch%imat(dn) <= 0) cycle
    if (patch%imat(up) /= material_id .and. &
        patch%imat(dn) /= material_id) cycle
    if (.not.(this%cell_ok(up) .and. this%cell_ok(dn))) cycle
    call ConnectionCalculateDistances(this%conn%int_dist(:,iconn), &
                                      this%option%gravity,dd_up,dd_dn, &
                                      dist_gravity,upweight)
    call PermeabilityTensorToScalar(this%maux(up),this%conn%int_dist(:,iconn), &
                                    perm_up)
    call PermeabilityTensorToScalar(this%maux(dn),this%conn%int_dist(:,iconn), &
                                    perm_dn)
    w = 0.d0
    if (patch%imat(up) == material_id) &
      w = w + dd_up*perm_dn/(dd_up*perm_dn + dd_dn*perm_up)
    if (patch%imat(dn) == material_id) &
      w = w + dd_dn*perm_up/(dd_up*perm_dn + dd_dn*perm_up)
    dv_int(iconn) = this%v_int(iconn)*w/k
    lup = grid%nG2L(up)
    ldn = grid%nG2L(dn)
    if (lup > 0) f_p((lup-1)*this%ndof+2) = f_p((lup-1)*this%ndof+2) + &
                                            this%flux_int(iconn)*w/k
    if (ldn > 0) f_p((ldn-1)*this%ndof+2) = f_p((ldn-1)*this%ndof+2) - &
                                            this%flux_int(iconn)*w/k
  enddo
  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    if (patch%imat(ghosted_id) /= material_id) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    select case(this%bc_itype(iconn))
      case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,HYDROSTATIC_BC, &
           HYDROSTATIC_SEEPAGE_BC,HET_DIRICHLET_BC, &
           HET_HYDROSTATIC_SEEPAGE_BC,UNIT_GRADIENT_BC)
        dv_bc(iconn) = this%v_bc(iconn)/k
        f_p((this%conn%bc_local(iconn)-1)*this%ndof+2) = &
          f_p((this%conn%bc_local(iconn)-1)*this%ndof+2) - this%flux_bc(iconn)/k
    end select
  enddo
  call VecRestoreArray(F,f_p,ierr);CHKERRQ(ierr)

end subroutine RichardsDAEPermeabilityDerivative

! ************************************************************************** !

subroutine RichardsDAEInitialState(this,time,Y)
  !
  ! P from the flow initial condition (field%flow_xx) and W = M(P), with
  ! the material of the initialized realization.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: y_p(:), xx_p(:)
  PetscInt :: local_id, ghosted_id, off, status
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  call RichardsDAECopyMaterial(this)
  this%prepared = PETSC_FALSE
  call VecGetArray(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(this%realization%field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    off = (local_id-1)*this%ndof
    y_p(off+1) = xx_p(local_id)
    this%p_frozen(local_id) = xx_p(local_id)
    y_p(off+2) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    call RichardsDAECellAuxVars(this,ghosted_id,y_p(off+1),status)
    if (status /= TS_MODEL_OK) then
      write(this%option%io_buffer,'(a,i9,a)') 'RICHARDS_DAE: the flow state &
        &of cell ',grid%nG2A(ghosted_id),' could not be evaluated at the &
        &initial condition.'
      call PrintErrMsgByRank(this%option)
    endif
    y_p(off+2) = RichardsDAECellInventory(this,ghosted_id)
  enddo
  call VecRestoreArrayRead(this%realization%field%flow_xx,xx_p, &
                           ierr);CHKERRQ(ierr)
  call VecRestoreArray(Y,y_p,ierr);CHKERRQ(ierr)

end subroutine RichardsDAEInitialState

! ************************************************************************** !

subroutine RichardsDAEDefaultTolerances(this,atol,rtol)
  !
  ! atol is a fraction of a cell's pore water; see RichardsDAESetTolerances.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol

  atol = 1.d-8
  rtol = 1.d-5

end subroutine RichardsDAEDefaultTolerances

! ************************************************************************** !

subroutine RichardsDAESetTolerances(this,atol,rtol,vatol,vrtol)
  !
  ! Pressure entries are algebraic: excluded from the local error estimate.
  ! An inventory's absolute tolerance is atol times the cell's pore water
  ! at full saturation.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol
  Vec :: vatol
  Vec :: vrtol

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: a_p(:), r_p(:)
  PetscInt :: local_id, ghosted_id, off
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  this%tol_rtol = rtol
  if (.not.allocated(this%tol_atol)) allocate(this%tol_atol(this%nlmax))
  call VecGetArray(vatol,a_p,ierr);CHKERRQ(ierr)
  call VecGetArray(vrtol,r_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    off = (local_id-1)*this%ndof
    ! pore water at full saturation [kmol, or kg], at the density of pure
    ! water
    this%tol_atol(local_id) = atol*max(patch%aux%Material%auxvars(ghosted_id)% &
                                         porosity,1.d-3)* &
                              patch%aux%Material%auxvars(ghosted_id)%volume* &
                              1000.d0/richards_density_kmol_to_kg
    a_p(off+1) = 0.d0
    r_p(off+1) = 0.d0
    a_p(off+2) = this%tol_atol(local_id)
    r_p(off+2) = rtol
  enddo
  call VecRestoreArray(vatol,a_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(vrtol,r_p,ierr);CHKERRQ(ierr)

end subroutine RichardsDAESetTolerances

! ************************************************************************** !

subroutine RichardsDAEConvergenceTest(this,time,X,dX,F,dt,kappa,converged)
  !
  ! In inventory units [kmol], with tol = rtol*|W| + atol:
  !   rate rows       dt*|F_W| <= kappa*tol
  !   updates         |dW|     <= kappa*tol
  !   algebraic rows  |W - M(P)| <= kappa*rtol*max(|W|,M) + ALG_ABS_FACTOR*atol
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: X
  Vec :: dX
  Vec :: F
  PetscReal :: dt
  PetscReal :: kappa
  PetscBool :: converged

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: x_p(:), dx_p(:), f_p(:)
  PetscReal :: tol, w
  PetscInt :: local_id, ghosted_id, off
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  converged = PETSC_TRUE
  call VecGetArrayRead(X,x_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(dX,dx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(F,f_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    off = (local_id-1)*this%ndof
    w = x_p(off+2)
    tol = this%tol_rtol*abs(w) + this%tol_atol(local_id)
    if (.not.(dt*abs(f_p(off+2)) <= kappa*tol)) converged = PETSC_FALSE
    if (.not.(abs(dx_p(off+2)) <= kappa*tol)) converged = PETSC_FALSE
    if (.not.(abs(f_p(off+1)) <= kappa*this%tol_rtol* &
              max(abs(w),abs(w-f_p(off+1))) + &
              ALG_ABS_FACTOR*this%tol_atol(local_id))) converged = PETSC_FALSE
    if (.not.converged) exit
  enddo
  call VecRestoreArrayRead(F,f_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(dX,dx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X,x_p,ierr);CHKERRQ(ierr)

end subroutine RichardsDAEConvergenceTest

! ************************************************************************** !

subroutine RichardsDAEPreCheck(this,Y,dY,changed)
  !
  ! Sets the inventory update so the new iterate is consistent,
  ! W = M(P - dP) (SNES applies Y - dY); Newton on P then matches classic
  ! RICHARDS. A pressure outside the domain of M leaves the update as is.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  Vec :: Y
  Vec :: dY
  PetscBool :: changed

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: y_p(:), dy_p(:)
  PetscInt :: local_id, ghosted_id, off, status
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  changed = PETSC_FALSE
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArray(dY,dy_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    off = (local_id-1)*this%ndof
    call RichardsDAECellAuxVars(this,ghosted_id,y_p(off+1)-dy_p(off+1),status)
    if (status /= TS_MODEL_OK) cycle
    dy_p(off+2) = y_p(off+2) - RichardsDAECellInventory(this,ghosted_id)
    changed = PETSC_TRUE
  enddo
  call VecRestoreArray(dY,dy_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  ! the private auxiliary variables now hold the trial pressures
  this%prepared = PETSC_FALSE

end subroutine RichardsDAEPreCheck

! ************************************************************************** !

function RichardsDAEUpdateSettled(this,dY)
  !
  ! Whether a Newton update, after RichardsDAEPreCheck, changes no
  ! inventory by more than its local error tolerance. Collective.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(richards_dae_model_type) :: this
  Vec :: dY
  PetscBool :: RichardsDAEUpdateSettled

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: dy_p(:)
  PetscInt :: local_id, off
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  RichardsDAEUpdateSettled = PETSC_TRUE
  call VecGetArrayRead(dY,dy_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    off = (local_id-1)*this%ndof
    if (.not.(abs(dy_p(off+2)) <= this%tol_rtol*abs(this%acc(local_id)) + &
                                  this%tol_atol(local_id))) then
      RichardsDAEUpdateSettled = PETSC_FALSE
      exit
    endif
  enddo
  call VecRestoreArrayRead(dY,dy_p,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,RichardsDAEUpdateSettled,ONE_INTEGER_MPI, &
                     MPI_C_BOOL,MPI_LAND,this%option%mycomm,ierr);CHKERRQ(ierr)

end function RichardsDAEUpdateSettled

! ************************************************************************** !

subroutine RichardsDAEPublish(this,time,Y)
  !
  ! Copies an accepted state to the realization for output only: flow_xx,
  ! the Richards and global auxiliary variables of cells and boundaries,
  ! porosity, and the Darcy and source fluxes.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module

  implicit none

  class(richards_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  PetscReal, pointer :: y_p(:), xx_p(:)
  PetscInt :: local_id, ghosted_id, iconn, status
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  field => this%realization%field
  call RichardsDAEPrepare(this,time,Y,PETSC_FALSE,status)
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    xx_p(local_id) = y_p((local_id-1)*this%ndof+1)
  enddo
  call VecRestoreArray(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   field%flow_xx,field%flow_xx_loc,NFLOWDOF)
  do ghosted_id = 1, this%ngmax
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    if (grid%nG2L(ghosted_id) < 0) cycle
    call RichardsAuxVarCopy(this%rich(ghosted_id), &
                            patch%aux%Richards%auxvars(ghosted_id),this%option)
    call GlobalAuxVarCopy(this%gaux(ghosted_id), &
                          patch%aux%Global%auxvars(ghosted_id),this%option)
    patch%aux%Material%auxvars(ghosted_id)%porosity = &
      this%rich(ghosted_id)%effective_porosity
  enddo
  do iconn = 1, this%conn%nconn_bc
    if (patch%imat(this%conn%bc_ghost(iconn)) <= 0) cycle
    call RichardsAuxVarCopy(this%rich_bc(iconn), &
                            patch%aux%Richards%auxvars_bc(iconn),this%option)
    call GlobalAuxVarCopy(this%gaux_bc(iconn), &
                          patch%aux%Global%auxvars_bc(iconn),this%option)
  enddo
  if (this%conn%nconn_int > 0) &
    patch%internal_velocities(1,1:this%conn%nconn_int) = &
      this%v_int(1:this%conn%nconn_int)
  if (this%conn%nconn_bc > 0) &
    patch%boundary_velocities(1,1:this%conn%nconn_bc) = &
      this%v_bc(1:this%conn%nconn_bc)
  if (associated(patch%ss_flow_vol_fluxes) .and. this%conn%nconn_ss > 0) &
    patch%ss_flow_vol_fluxes(1,1:this%conn%nconn_ss) = &
      this%qvol_ss(1:this%conn%nconn_ss)

end subroutine RichardsDAEPublish

! ************************************************************************** !

subroutine RichardsDAEStoreRecord(this,k)
  !
  ! Keeps the current segment as record k (see TS_Model_class).
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscInt :: k

  type(richards_dae_segment_type), allocatable :: grown(:)

  if (.not.allocated(this%records)) allocate(this%records(16))
  if (k > size(this%records)) then
    allocate(grown(2*size(this%records)))
    grown(1:size(this%records)) = this%records
    call move_alloc(grown,this%records)
  endif
  this%records(k) = this%seg

end subroutine RichardsDAEStoreRecord

! ************************************************************************** !

subroutine RichardsDAERestoreRecord(this,k)
  !
  ! Makes record k the current segment.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this
  PetscInt :: k

  this%seg = this%records(k)
  this%segment_id = this%seg%id
  this%forcing_segment = -1
  this%forcing_time = -MAX_DOUBLE
  this%prepared = PETSC_FALSE

end subroutine RichardsDAERestoreRecord

! ************************************************************************** !

subroutine RichardsDAECopyMaterial(this)
  !
  ! Copies the realization's material auxiliary variables. Repeated after
  ! creation: PM initialization completes them (soil reference pressures,
  ! porosity from mineral volume fractions) and an adjoint's parameter
  ! check changes them in place.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this

  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id

  material_auxvars => this%realization%patch%aux%Material%auxvars
  do ghosted_id = 1, this%ngmax
    call MaterialAuxVarCopy(material_auxvars(ghosted_id), &
                            this%maux(ghosted_id),this%option)
  enddo

end subroutine RichardsDAECopyMaterial

! ************************************************************************** !

subroutine RichardsDAEInvalidateCache(this)
  !
  ! Material parameters changed in place: the private material auxiliary
  ! variables take the realization's, and the prepared flow state expires.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(richards_dae_model_type) :: this

  call RichardsDAECopyMaterial(this)
  this%prepared = PETSC_FALSE

end subroutine RichardsDAEInvalidateCache

! ************************************************************************** !

subroutine RichardsDAEResetToInitial(this)
  !
  ! The forcing and flow-state caches expire, for a rerun.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(richards_dae_model_type) :: this

  this%forcing_time = -MAX_DOUBLE
  this%forcing_segment = -1
  call this%InvalidateCache()

end subroutine RichardsDAEResetToInitial

! ************************************************************************** !

subroutine RichardsDAEDestroy(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  use Petsc_Utility_module, only : PUVecDestroy, PUDMDestroy

  implicit none

  class(richards_dae_model_type) :: this

  PetscInt :: i

  ! (the Richards auxiliary variables own nothing without a well model)
  if (associated(this%rich)) then
    do i = 1, size(this%rich)
      call GlobalAuxVarStrip(this%gaux(i))
      call MaterialAuxVarStrip(this%maux(i))
    enddo
    deallocate(this%rich,this%gaux,this%maux)
  endif
  if (associated(this%rich_bc)) then
    do i = 1, this%conn%nconn_bc
      call GlobalAuxVarStrip(this%gaux_bc(i))
    enddo
    deallocate(this%rich_bc,this%gaux_bc)
  endif
  call PUVecDestroy(this%Y_loc)
  call PUVecDestroy(this%Y_prepared)
  if (this%private_dm .and. associated(this%dm)) then
    call PUDMDestroy(this%dm%dm)
    if (associated(this%dm%ugdm)) call UGridDMDestroy(this%dm%ugdm)
    deallocate(this%dm)
  endif
  nullify(this%dm)

end subroutine RichardsDAEDestroy

end module Richards_DAE_Model_class
