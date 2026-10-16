module RT_DAE_Model_class
  !
  ! Reactive transport as a DAE model for PETSc TS (MODE GIRT_DAE).
  !
  ! State per grid cell, Y = [c; N; K]:
  !   c(1:ncomp)   free-ion molality of the aqueous primaries [mol/kg] and
  !                immobile concentrations [mol/m^3 bulk]; their logarithms
  !                under LOG_FORMULATION
  !   N(1:ncomp)   component inventory [mol]
  !   K(1:nkin)    kinetic inventories [mol]; minerals first
  !
  ! Equations of an active cell:
  !   F_c = N - A(t,c,K)                              (algebraic)
  !   F_N = dN/dt + sum_faces flux(t,c) + R(t,c,K) - Q(t,c)
  !   F_K = dK/dt - j(t,c,K)
  ! A: equilibrium inventory, V*(phi*s*1000*T(c) + S(c,K)) for aqueous
  ! components (S: equilibrium sorbed) and V*c for immobile ones. R: reaction
  ! sink of the transported components (precipitation consumes). j: kinetic
  ! transfer rates, each distributed to both balances through its
  ! stoichiometry B, so N + B K is conserved for any integrator.
  ! Inactive cells: F_c = c - c_frozen, F_N = dN/dt, F_K = dK/dt.
  !
  ! The Jacobian is a condensed operator (RT_DAE_Linear_module) with linear
  ! solves of ncomp per cell; -rt_dae_assembled_jacobian assembles a block
  ! matrix of ndof per cell instead.
  !
  ! Evaluation reads the realization only for configuration (geometry,
  ! connectivity, material and reaction parameters); the forcing
  ! (RT_DAE_Forcing_module) copies each segment's flow state, fluxes and
  ! boundary and source chemistry from it.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use TS_Model_class
  use PFLOTRAN_Constants_module
  use Realization_Subsurface_class
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use RT_DAE_Chem_module
  use RT_DAE_Linear_module
  use RT_DAE_Forcing_module
  use TS_Flow_State_module
  use TS_Connections_module

  implicit none

  private

  ! regime of a kinetic mineral in a cell: present (dissolves and
  ! precipitates) or absent (zero inventory, precipitation only)
  PetscInt, parameter, public :: RT_DAE_MNRL_PRESENT = 1
  PetscInt, parameter, public :: RT_DAE_MNRL_ABSENT = 2

  ! smallest fraction of a concentration a Newton update may leave
  PetscReal, parameter :: PRECHECK_FLOOR = 1.d-10
  ! LOG_FORMULATION: largest change of ln(c) per Newton iteration, and the
  ! smallest ln(c) of an iterate or initial guess (exp would underflow
  ! near -708)
  PetscReal, parameter :: LOG_MAX_UPDATE = 25.d0
  PetscReal, parameter :: LOG_MIN = -690.d0
  ! absolute part of the algebraic convergence test, relative to the
  ! inventory's atol; smaller residuals are not resolvable against the
  ! other rows of the linear system
  PetscReal, parameter :: ALG_ABS_FACTOR = 1.d-12
  ! value tolerance of the mineral event functions (K/mnrl_scale): a
  ! mineral within it of zero is absent, or its depletion would go unseen
  PetscReal, parameter :: RT_DAE_EVENT_VTOL = 1.d-9

  type, public, extends(ts_model_type) :: rt_dae_model_type
    class(realization_subsurface_type), pointer :: realization => null()
    type(option_type), pointer :: option => null()
    class(reaction_rt_type), pointer :: reaction => null()
    PetscInt :: ncomp = 0
    PetscInt :: naqcomp = 0
    PetscInt :: nimcomp = 0
    PetscInt :: nkin = 0
    PetscInt :: nmnrl = 0
    ! multirate surface complexation: kinetic sorbed inventories [mol] in K
    ! after the minerals (see RTDAEMRIndex)
    PetscInt :: nmr = 0
    PetscInt, allocatable :: mr_ncomp(:)
    PetscInt, allocatable :: mr_comp(:,:)      ! (naqcomp,nmr)
    PetscInt, allocatable :: mr_offset(:)
    ! equilibrium surface complexation reactions in the accumulation
    PetscBool :: eq_srfcplx = PETSC_FALSE
    PetscReal, allocatable :: sorb_work(:,:)    ! (naqcomp,naqcomp)
    ! dS/dVF of the equilibrium sorbed totals on each kinetic mineral
    PetscReal, allocatable :: sorb_dvf(:,:)     ! (naqcomp,nmnrl)
    ! coupled with soil compressibility: the sorbed totals on rock sites,
    ! whose density follows the solid fraction 1 - phi(P)
    PetscReal, allocatable :: sorb_rock(:,:)    ! (naqcomp,ngmax)
    PetscInt :: nlmax = 0
    PetscInt :: ngmax = 0
    ! transport phases (liquid only)
    PetscInt :: nphase = 0
    ! connection table, in the order of the patch's connection lists (and
    ! of its velocity and transport coefficient arrays)
    type(ts_connections_type) :: conn
    PetscInt, allocatable :: bc_itype(:)     ! (nconn_bc) transport BC type
    PetscInt, allocatable :: ss_itype(:)     ! (nconn_ss) transport SS type
    PetscBool :: log_formulation = PETSC_FALSE
    ! -rt_dae_debug: report failed local closures
    PetscBool :: debug = PETSC_FALSE
    ! -rt_dae_fail_closure n (a test of the recovery): the n-th preparation
    ! fails the closure of every local cell
    PetscInt :: fail_closure = 0
    PetscInt :: num_preparations = 0
    ! -rt_dae_assembled_jacobian: block matrix, not the condensed operator
    PetscBool :: assembled_jacobian = PETSC_FALSE
    ! local error tolerance of the concentrations (relative; absolute in
    ! log coordinates), -rt_dae_conc_tol; 0 excludes them
    PetscReal :: conc_tol = 0.d0
    ! components whose every species has a positive stoichiometric
    ! coefficient: their inventory cannot be negative
    PetscBool, allocatable :: positive_comp(:)  ! (ncomp)
    ! local error tolerances of the inventories, set by SetTolerances
    PetscReal :: tol_rtol = 0.d0
    PetscReal, allocatable :: tol_atol(:,:)   ! (ncomp+nkin,nlmax) [mol]
    ! the forcing segments, the flow state and face fluxes at the forcing
    ! time, and the records (see TS_Model_class)
    type(rt_dae_forcing_type) :: forcing
    ! coupled (Flow_RT_DAE_Model_class): fs, set by the coupled model before
    ! every evaluation, replaces the segment's flow state
    type(ts_flow_state_type), pointer :: fs => null()
    ! private material auxiliary variables; the forcing sets the porosity
    type(material_auxvar_type), pointer :: maux(:) => null()
    ! dispersion coefficients over distance at the forcing time; pointers,
    ! as RTComputeTransportCoefs takes them
    PetscReal, pointer :: tran_coefs_int(:,:,:) => null()
    PetscReal, pointer :: tran_coefs_bc(:,:,:) => null()
    ! transport coefficients [L/s] at the forcing time
    PetscReal, allocatable :: coef_up(:,:)
    PetscReal, allocatable :: coef_dn(:,:)
    PetscReal, allocatable :: coef_bc_up(:,:)
    PetscReal, allocatable :: coef_bc_dn(:,:)
    PetscReal, allocatable :: coef_ss_in(:)
    PetscReal, allocatable :: coef_ss_out(:)
    ! cache key of the per-cell arrays below (see RTDAEPrepare)
    PetscBool :: prepared = PETSC_FALSE
    PetscReal :: prepared_time = 0.d0
    PetscInt :: prepared_segment = -1
    PetscBool :: prepared_derivatives = PETSC_FALSE
    PetscInt :: prepared_status = TS_MODEL_OK
    PetscInt :: prepared_flow_id = -1
    Vec :: Y_prepared
    ! private chemistry scratch
    type(reactive_transport_auxvar_type), pointer :: scratch => null()
    type(rt_dae_gamma_type) :: gam
    ! per ghosted cell: aqueous totals [mol/L] and dT/dY_c [kg/L, or
    ! mol/L under LOG_FORMULATION]
    PetscReal, allocatable :: tot(:,:)
    PetscReal, allocatable :: dtot(:,:,:)
    PetscBool, allocatable :: cell_ok(:)
    ! after RTDAEKeepActivities: the activity coefficients of each cell's
    ! last closure
    PetscReal, allocatable :: pri_act_coef(:,:)  ! (naqcomp,ngmax)
    ! per local cell at the last Prepare: accumulation, reaction sink and
    ! kinetic rates, and their derivatives with respect to the cell's own
    ! c coordinates (log chain applied) and kinetic inventories
    PetscReal, allocatable :: acc(:,:)          ! (ncomp,nlmax)
    ! gross equilibrium inventory: species amounts summed with |nu|, the
    ! scale of the roundoff and of the Newton error in acc
    PetscReal, allocatable :: acc_gross(:,:)    ! (ncomp,nlmax)
    PetscReal, allocatable :: rn(:,:)           ! (ncomp,nlmax)
    PetscReal, allocatable :: jk(:,:)           ! (nkin,nlmax)
    PetscReal, allocatable :: dacc_dc(:,:,:)    ! (ncomp,ncomp,nlmax)
    PetscReal, allocatable :: drn_dc(:,:,:)     ! (ncomp,ncomp,nlmax)
    PetscReal, allocatable :: drn_dk(:,:,:)     ! (ncomp,nkin,nlmax)
    PetscReal, allocatable :: djk_dc(:,:,:)     ! (nkin,ncomp,nlmax)
    PetscReal, allocatable :: djk_dk(:,:,:)     ! (nkin,nkin,nlmax)
    PetscReal, allocatable :: mnrl_rate(:,:)    ! (nmnrl,nlmax)
    PetscReal, allocatable :: mnrl_area(:,:)    ! (nmnrl,nlmax)
    PetscReal, allocatable :: mnrl_lnQK(:,:)    ! (nmnrl,nlmax)
    PetscReal, allocatable :: dacc_dk(:,:,:)    ! (ncomp,nkin,nlmax)
    ! discrete state of the mineral regimes and its transition record
    PetscInt, allocatable :: mnrl_regime(:,:)   ! (nmnrl,nlmax)
    ! the regimes InitialState set, for ResetToInitial
    PetscInt, allocatable :: mnrl_regime0(:,:)  ! (nmnrl,nlmax)
    PetscReal, allocatable :: mnrl_scale(:,:)   ! (nmnrl,nlmax) [mol]
    PetscInt :: num_depletions = 0
    PetscInt :: num_reappearances = 0
    ! inventory moved to the dissolved components at depletions [mol]
    PetscReal :: depletion_residual = 0.d0
    ! immutable per-cell data
    PetscReal, allocatable :: c_frozen(:,:)     ! (ncomp,nlmax) inactive cells
    PetscReal, allocatable :: mnrl_vf0(:,:)     ! (nmnrl,ngmax)
    PetscReal, allocatable :: mnrl_area0(:,:)   ! (nmnrl,ngmax)
    Vec :: Y_loc
    ! work blocks of an assembled Jacobian
    PetscReal, allocatable :: blk(:,:)
    PetscReal, allocatable :: blk2(:,:)
  contains
    procedure, public :: InitialState => RTDAEInitialState
    procedure, public :: BeginSegment => RTDAEBeginSegment
    procedure, public :: Evaluate => RTDAEEvaluate
    procedure, public :: Linearize => RTDAELinearize
    procedure, public :: CreateJacobian => RTDAECreateJacobian
    procedure, public :: SetupLinearSolver => RTDAESetupLinearSolver
    procedure, public :: SetTolerances => RTDAESetTolerances
    procedure, public :: ConvergenceTest => RTDAEConvergenceTest
    procedure, public :: DefaultTolerances => RTDAEDefaultTolerances
    procedure, public :: Publish => RTDAEPublish
    procedure, public :: PreCheck => RTDAEPreCheck
    procedure, public :: ProjectInitialGuess => RTDAEProjectInitialGuess
    procedure, public :: NumEvents => RTDAENumEvents
    procedure, public :: EventSetup => RTDAEEventSetup
    procedure, public :: EventIndicator => RTDAEEventIndicator
    procedure, public :: PostEvent => RTDAEPostEvent
    procedure, public :: PrintSummary => RTDAEPrintSummary
    procedure, public :: ForcingTimes => RTDAEForcingTimes
    procedure, public :: StoreRecord => RTDAEStoreRecord
    procedure, public :: RestoreRecord => RTDAERestoreRecord
    procedure, public :: InvalidateCache => RTDAEInvalidateCache
    procedure, public :: ResetToInitial => RTDAEResetToInitial
    procedure, public :: Destroy => RTDAEDestroy
  end type rt_dae_model_type

  public :: RTDAEModelCreate, &
            RTDAEPrepareReactions, &
            RTDAELoadCell, &
            RTDAELoadCellRates, &
            RTDAEPoreLiters, &
            RTDAEInitialSensitivity, &
            RTDAECellObservables, &
            RTDAEInitialCellState

contains

! ************************************************************************** !

function RTDAEModelCreate(realization,flow_coupled)
  !
  ! Creates the model, rejects unsupported configurations and allocates its
  ! scratch. flow_coupled: the flow state will be part of the solution
  ! (RTDAESetFlowCoupling) rather than forcing.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Discretization_module
  use Grid_module
  use Patch_module
  use Reactive_Transport_Aux_module, only : RTAuxVarInit
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(realization_subsurface_type), pointer :: realization
  PetscBool :: flow_coupled
  class(rt_dae_model_type), pointer :: RTDAEModelCreate

  class(rt_dae_model_type), pointer :: this
  type(option_type), pointer :: option
  class(reaction_rt_type), pointer :: reaction
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(surface_complexation_type), pointer :: srfcplx
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id, icplx, i, j, ieqrxn, irxn
  PetscBool :: flag
  PetscErrorCode :: ierr

  allocate(this)
  this%name = 'RT DAE'
  this%integrates_transport = PETSC_TRUE
  this%realization => realization
  option => realization%option
  reaction => realization%reaction
  patch => realization%patch
  grid => patch%grid
  this%option => option
  this%reaction => reaction

  call RTDAECheckCapabilities(this,flow_coupled)
  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) &
    call PrintMsg(option,'  GIRT_DAE converges the activity coefficients &
      &at every evaluation: the update frequency and algorithm of &
      &ACTIVITY_COEFFICIENTS do not apply.')
  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-rt_dae_debug',this%debug,ierr);CHKERRQ(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                          '-rt_dae_fail_closure',this%fail_closure,flag, &
                          ierr);CHKERRQ(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-rt_dae_assembled_jacobian', &
                           this%assembled_jacobian,flag,ierr);CHKERRQ(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-rt_dae_conc_tol',this%conc_tol,flag, &
                           ierr);CHKERRQ(ierr)

  this%ncomp = reaction%ncomp
  this%naqcomp = reaction%naqcomp
  this%nimcomp = reaction%ncomp - reaction%naqcomp
  this%nmnrl = reaction%mineral%nkinmnrl
  call RTDAEMultirateLayout(this)
  this%nkin = this%nmnrl + &
    sum(this%mr_ncomp(1:this%nmr)* &
        reaction%surface_complexation%kinmr_nrate(1:this%nmr))
  this%ndof = 2*this%ncomp + this%nkin
  this%eq_srfcplx = (reaction%surface_complexation%neqsrfcplxrxn > 0)

  ! the DM of the TS solution: ndof per grid cell
  call DiscretizationCreateDM(realization%discretization, &
                              realization%discretization%dm_ndaedof, &
                              this%ndof, &
                              realization%discretization%stencil_width, &
                              realization%discretization%stencil_type, &
                              option,'')
  realization%discretization%dm_index_to_ndof(NDAEDOF) = this%ndof
  this%log_formulation = reaction%use_log_formulation
  this%nlmax = grid%nlmax
  this%ngmax = grid%ngmax

  call RTDAEConnectionTable(this)
  call RTDAEForcingInit(this%forcing,realization,this%conn)

  ! private copy of the material auxiliary variables; the porosity is
  ! overwritten at each forcing time
  material_auxvars => patch%aux%Material%auxvars
  allocate(this%maux(this%ngmax))
  do ghosted_id = 1, this%ngmax
    call MaterialAuxVarInit(this%maux(ghosted_id),option)
    call MaterialAuxVarCopy(material_auxvars(ghosted_id), &
                            this%maux(ghosted_id),option)
  enddo

  allocate(this%coef_up(this%naqcomp,max(this%conn%nconn_int,1)))
  allocate(this%coef_dn(this%naqcomp,max(this%conn%nconn_int,1)))
  allocate(this%coef_bc_up(this%naqcomp,max(this%conn%nconn_bc,1)))
  allocate(this%coef_bc_dn(this%naqcomp,max(this%conn%nconn_bc,1)))
  allocate(this%coef_ss_in(max(this%conn%nconn_ss,1)))
  allocate(this%coef_ss_out(max(this%conn%nconn_ss,1)))
  this%nphase = patch%aux%RT%rt_parameter%nphase
  allocate(this%tran_coefs_int(patch%aux%RT%rt_parameter%ndiffcoef, &
                               this%nphase,max(this%conn%nconn_int,1)))
  allocate(this%tran_coefs_bc(patch%aux%RT%rt_parameter%ndiffcoef, &
                              this%nphase,max(this%conn%nconn_bc,1)))
  this%tran_coefs_int = 0.d0
  this%tran_coefs_bc = 0.d0
  this%coef_up = 0.d0
  this%coef_dn = 0.d0
  this%coef_bc_up = 0.d0
  this%coef_bc_dn = 0.d0
  this%coef_ss_in = 0.d0
  this%coef_ss_out = 0.d0

  allocate(this%scratch)
  call RTAuxVarInit(this%scratch,reaction,PETSC_FALSE,option)
  call RTDAEGammaInit(this%gam,reaction)

  allocate(this%sorb_work(this%naqcomp,this%naqcomp))
  allocate(this%sorb_dvf(this%naqcomp,max(this%nmnrl,1)))
  this%sorb_dvf = 0.d0
  allocate(this%tot(this%naqcomp,this%ngmax))
  allocate(this%dtot(this%naqcomp,this%naqcomp,this%ngmax))
  allocate(this%cell_ok(this%ngmax))
  this%tot = 0.d0
  this%dtot = 0.d0
  this%cell_ok = PETSC_TRUE

  allocate(this%acc(this%ncomp,this%nlmax))
  allocate(this%acc_gross(this%ncomp,this%nlmax))
  this%acc_gross = 0.d0
  allocate(this%rn(this%ncomp,this%nlmax))
  allocate(this%jk(max(this%nkin,1),this%nlmax))
  allocate(this%dacc_dc(this%ncomp,this%ncomp,this%nlmax))
  allocate(this%drn_dc(this%ncomp,this%ncomp,this%nlmax))
  allocate(this%drn_dk(this%ncomp,max(this%nkin,1),this%nlmax))
  allocate(this%djk_dc(max(this%nkin,1),this%ncomp,this%nlmax))
  allocate(this%djk_dk(max(this%nkin,1),max(this%nkin,1),this%nlmax))
  allocate(this%dacc_dk(this%ncomp,max(this%nkin,1),this%nlmax))
  this%dacc_dk = 0.d0
  allocate(this%mnrl_rate(max(this%nmnrl,1),this%nlmax))
  allocate(this%mnrl_area(max(this%nmnrl,1),this%nlmax))
  allocate(this%mnrl_lnQK(max(this%nmnrl,1),this%nlmax))
  allocate(this%mnrl_regime(max(this%nmnrl,1),this%nlmax))
  allocate(this%mnrl_scale(max(this%nmnrl,1),this%nlmax))
  this%mnrl_lnQK = 0.d0
  this%mnrl_regime = RT_DAE_MNRL_PRESENT
  this%mnrl_scale = 1.d0
  this%acc = 0.d0
  this%rn = 0.d0
  this%jk = 0.d0
  this%dacc_dc = 0.d0
  this%drn_dc = 0.d0
  this%drn_dk = 0.d0
  this%djk_dc = 0.d0
  this%djk_dk = 0.d0
  this%mnrl_rate = 0.d0
  this%mnrl_area = 0.d0
  allocate(this%c_frozen(this%ncomp,this%nlmax))
  this%c_frozen = 0.d0
  allocate(this%mnrl_vf0(max(this%nmnrl,1),this%ngmax))
  allocate(this%mnrl_area0(max(this%nmnrl,1),this%ngmax))
  this%mnrl_vf0 = 0.d0
  this%mnrl_area0 = 0.d0

  ! components whose inventory cannot be negative: no aqueous or
  ! equilibrium surface complex has a negative coefficient in them
  allocate(this%positive_comp(this%ncomp))
  this%positive_comp = PETSC_TRUE
  do icplx = 1, reaction%neqcplx
    do i = 1, reaction%eqcplxspecid(0,icplx)
      if (reaction%eqcplxstoich(i,icplx) < 0.d0) &
        this%positive_comp(reaction%eqcplxspecid(i,icplx)) = PETSC_FALSE
    enddo
  enddo
  srfcplx => reaction%surface_complexation
  do ieqrxn = 1, srfcplx%neqsrfcplxrxn
    irxn = srfcplx%eqsrfcplxrxn_to_srfcplxrxn(ieqrxn)
    do j = 1, srfcplx%srfcplxrxn_to_complex(0,irxn)
      icplx = srfcplx%srfcplxrxn_to_complex(j,irxn)
      do i = 1, srfcplx%srfcplxspecid(0,icplx)
        if (srfcplx%srfcplxstoich(i,icplx) < 0.d0) &
          this%positive_comp(srfcplx%srfcplxspecid(i,icplx)) = PETSC_FALSE
      enddo
    enddo
  enddo

  call DiscretizationCreateVector(realization%discretization,NDAEDOF, &
                                  this%Y_loc,LOCAL,option)
  call DiscretizationCreateVector(realization%discretization,NDAEDOF, &
                                  this%Y_prepared,GLOBAL,option)

  RTDAEModelCreate => this

end function RTDAEModelCreate

! ************************************************************************** !

subroutine RTDAEConnectionTable(this)
  !
  ! The patch's connections (TSConnectionsBuild) and the transport
  ! condition type of each boundary and source/sink connection. Immutable:
  ! geometry and transport condition types do not change during a run.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Coupler_module

  implicit none

  class(rt_dae_model_type) :: this

  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: coupler
  PetscInt :: iconn, n

  patch => this%realization%patch
  call TSConnectionsBuild(this%conn,patch)

  allocate(this%bc_itype(this%conn%nconn_bc))
  n = 0
  coupler => patch%boundary_condition_list%first
  do while (associated(coupler))
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      this%bc_itype(n) = coupler%tran_condition%itype
    enddo
    coupler => coupler%next
  enddo

  allocate(this%ss_itype(this%conn%nconn_ss))
  n = 0
  coupler => patch%source_sink_list%first
  do while (associated(coupler))
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      this%ss_itype(n) = coupler%tran_condition%itype
    enddo
    coupler => coupler%next
  enddo

end subroutine RTDAEConnectionTable

! ************************************************************************** !

subroutine RTDAEMultirateLayout(this)
  !
  ! Sorbing components of each multirate surface complexation reaction and
  ! the offsets of their kinetic inventories in K.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(rt_dae_model_type) :: this

  type(surface_complexation_type), pointer :: srfcplx
  PetscInt :: ikinmr, irxn, j, i, icplx, icomp, offset
  PetscBool :: sorbs(this%realization%reaction%naqcomp)

  srfcplx => this%realization%reaction%surface_complexation
  this%nmr = srfcplx%nkinmrsrfcplxrxn
  allocate(this%mr_ncomp(max(this%nmr,1)))
  allocate(this%mr_comp(this%naqcomp,max(this%nmr,1)))
  allocate(this%mr_offset(max(this%nmr,1)))
  this%mr_ncomp = 0
  this%mr_comp = 0
  this%mr_offset = 0
  offset = 0
  do ikinmr = 1, this%nmr
    irxn = srfcplx%kinmrsrfcplxrxn_to_srfcplxrxn(ikinmr)
    sorbs = PETSC_FALSE
    do j = 1, srfcplx%srfcplxrxn_to_complex(0,irxn)
      icplx = srfcplx%srfcplxrxn_to_complex(j,irxn)
      do i = 1, srfcplx%srfcplxspecid(0,icplx)
        icomp = srfcplx%srfcplxspecid(i,icplx)
        if (srfcplx%srfcplxstoich(i,icplx) /= 0.d0) sorbs(icomp) = PETSC_TRUE
      enddo
    enddo
    do icomp = 1, this%naqcomp
      if (sorbs(icomp)) then
        this%mr_ncomp(ikinmr) = this%mr_ncomp(ikinmr) + 1
        this%mr_comp(this%mr_ncomp(ikinmr),ikinmr) = icomp
      endif
    enddo
    this%mr_offset(ikinmr) = offset
    offset = offset + this%mr_ncomp(ikinmr)*srfcplx%kinmr_nrate(ikinmr)
  enddo

end subroutine RTDAEMultirateLayout

! ************************************************************************** !

function RTDAEMRIndex(this,ikinmr,irate,iq)
  !
  ! Position in K of the kinetic sorbed inventory of multirate reaction
  ! ikinmr, rate irate and sorbing component iq (see RTDAEMultirateLayout).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ikinmr, irate, iq
  PetscInt :: RTDAEMRIndex

  RTDAEMRIndex = this%nmnrl + this%mr_offset(ikinmr) + &
                 (irate-1)*this%mr_ncomp(ikinmr) + iq

end function RTDAEMRIndex

! ************************************************************************** !

function RTDAEMineralVolFrac(this,imnrl,inventory,volume)
  !
  ! Volume fraction of kinetic mineral imnrl from its inventory [mol] in a
  ! cell of the given bulk volume [m^3].
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: imnrl
  PetscReal :: inventory, volume
  PetscReal :: RTDAEMineralVolFrac

  RTDAEMineralVolFrac = &
    this%reaction%mineral%kinmnrl_molar_vol(imnrl)*inventory/volume

end function RTDAEMineralVolFrac

! ************************************************************************** !

function RTDAEMineralInventory(this,imnrl,volfrac,volume)
  !
  ! Inventory [mol] of kinetic mineral imnrl at a volume fraction in a cell
  ! of the given bulk volume [m^3]; the inverse of RTDAEMineralVolFrac.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: imnrl
  PetscReal :: volfrac, volume
  PetscReal :: RTDAEMineralInventory

  RTDAEMineralInventory = &
    volfrac*volume/this%reaction%mineral%kinmnrl_molar_vol(imnrl)

end function RTDAEMineralInventory

! ************************************************************************** !

function RTDAEMineralAbove(this,imnrl,local_id,K)
  !
  ! Whether inventory K of kinetic mineral imnrl in local cell local_id
  ! exceeds the event tolerance; below it the event function could not
  ! see the mineral deplete.
  !
  ! Author: Hansol Suh
  ! Date: 10/13/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: imnrl, local_id
  PetscReal :: K
  PetscBool :: RTDAEMineralAbove

  RTDAEMineralAbove = K > RT_DAE_EVENT_VTOL*this%mnrl_scale(imnrl,local_id)

end function RTDAEMineralAbove

! ************************************************************************** !

subroutine RTDAEMineralToSolution(this,imnrl,K,N)
  !
  ! Moves inventory K of kinetic mineral imnrl to the components it
  ! dissolves into, N_i += nu_i K: the map of a depletion.
  !
  ! Author: Hansol Suh
  ! Date: 10/13/26
  !
  use Reaction_Mineral_Aux_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: imnrl
  PetscReal :: K
  PetscReal :: N(:)

  type(mineral_type), pointer :: mineral
  PetscInt :: i, icomp

  mineral => this%reaction%mineral
  do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
    icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
    N(icomp) = N(icomp) + mineral%kinmnrlstoich_in_residual(i,imnrl)*K
  enddo

end subroutine RTDAEMineralToSolution

! ************************************************************************** !

function RTDAEPoreLiters(sat_por,volume)
  !
  ! Liters of pore water in a cell: liquid saturation times porosity [-]
  ! and bulk volume [m^3].
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  PetscReal :: sat_por, volume
  PetscReal :: RTDAEPoreLiters

  RTDAEPoreLiters = sat_por*1000.d0*volume

end function RTDAEPoreLiters

! ************************************************************************** !

subroutine RTDAESorbedEq1(this,ghosted_id,irxn,free_site_conc,srfcplx_conc, &
                          S,dS_dm)
  !
  ! Adds the equilibrium sorbed totals of reaction irxn [mol/m^3 bulk] and
  ! their derivative (activity coefficients fixed) to S and dS_dm, at the
  ! chemistry in scratch; complex concentrations go to srfcplx_conc if
  ! associated. With unit free-site stoichiometry the free-site guess
  ! cancels: the result depends on the chemistry alone.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_Surface_Complexation_module, only : ReactionSrfCplxTotalSorbEq1

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id, irxn
  PetscReal :: free_site_conc
  PetscReal, pointer :: srfcplx_conc(:)
  PetscReal :: S(:), dS_dm(:,:)

  free_site_conc = 1.d-9
  call ReactionSrfCplxTotalSorbEq1(this%scratch,this%forcing%gaux(ghosted_id), &
                                   this%maux(ghosted_id),this%reaction, &
                                   this%option,irxn,free_site_conc, &
                                   srfcplx_conc,S,dS_dm)

end subroutine RTDAESorbedEq1

! ************************************************************************** !

subroutine RTDAESurfaceComplexation(this,ghosted_id,irxn,S,dS_dm,dS_dvf, &
                                    imnrl)
  !
  ! Equilibrium sorbed totals S [mol/m^3 bulk] of surface complexation
  ! reaction irxn at the chemistry in scratch, with the complete dS/dm
  ! (activity coefficient sensitivity included) and, on a kinetic mineral's
  ! surface, the kinetic mineral index imnrl (0 otherwise) and dS/dVF.
  ! With unit free-site stoichiometry S is linear in the site density:
  ! dS/dVF = S/VF, or where S vanishes (VF <= 0, or too few sites) the
  ! slope from the positive side, S at a unit volume fraction.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id, irxn
  PetscReal :: S(:), dS_dm(:,:), dS_dvf(:)
  PetscInt :: imnrl

  type(surface_complexation_type), pointer :: srfcplx
  PetscReal, pointer :: null_array(:)
  PetscReal :: volfrac

  srfcplx => this%reaction%surface_complexation
  nullify(null_array)
  S = 0.d0
  dS_dm = 0.d0
  dS_dvf = 0.d0
  imnrl = 0
  call RTDAESorbedEq1(this,ghosted_id,irxn, &
                      this%scratch%srfcplxrxn_free_site_conc(irxn), &
                      null_array,S,dS_dm)
  call RTDAEGammaChain(dS_dm,this%naqcomp,this%naqcomp, &
                       this%scratch%pri_molal,this%gam)
  if (srfcplx%srfcplxrxn_surf_type(irxn) /= MINERAL_SURFACE) return
  imnrl = srfcplx%srfcplxrxn_to_surf(irxn)
  volfrac = this%scratch%mnrl_volfrac(imnrl)
  if (volfrac > 0.d0 .and. any(S /= 0.d0)) then
    dS_dvf = S/volfrac
  else
    this%scratch%mnrl_volfrac(imnrl) = 1.d0
    this%sorb_work = 0.d0
    call RTDAESorbedEq1(this,ghosted_id,irxn, &
                        this%scratch%srfcplxrxn_free_site_conc(irxn), &
                        null_array,dS_dvf,this%sorb_work)
    this%scratch%mnrl_volfrac(imnrl) = volfrac
  endif

end subroutine RTDAESurfaceComplexation

! ************************************************************************** !

subroutine RTDAECheckCapabilities(this,flow_coupled)
  !
  ! Rejects every configuration whose physics this model does not define,
  ! naming all of them. flow_coupled: see RTDAEModelCreate.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_Sandbox_module, only : rxn_sandbox_list
  use Carbon_Sandbox_module, only : carbon_sandbox_list
  use CLM_Rxn_module, only : clmrxn_list
  use Reaction_Mineral_Aux_module
  use Reaction_Surface_Complexation_Aux_module
  use Patch_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscBool :: flow_coupled

  type(option_type), pointer :: option
  class(reaction_rt_type), pointer :: reaction
  type(patch_type), pointer :: patch
  character(len=MAXSTRINGLENGTH) :: unsupported

  option => this%realization%option
  reaction => this%realization%reaction
  patch => this%realization%patch
  unsupported = ''

  if (option%transport%nphase > 1 .or. reaction%gas%nactive_gas > 0) &
    call TSModelListAdd(unsupported,'ACTIVE_GAS_SPECIES')
  if (option%transport%couple_co2) &
    call TSModelListAdd(unsupported,'CO2 coupling')
  if (reaction%neqionxrxn > 0) call TSModelListAdd(unsupported,'ion exchange')
  if (reaction%neqdynamickdrxn > 0 .or. reaction%isotherm%neqkdrxn > 0) &
    call TSModelListAdd(unsupported,'sorption isotherms')
  if (reaction%surface_complexation%nkinsrfcplxrxn > 0) &
    call TSModelListAdd(unsupported,'kinetic surface complexation')
  ! unit free-site stoichiometry gives the free sites in closed form (else
  ! a Newton iteration from RTDAESorbedEq1's fixed guess) and dS/dVF = S/VF
  ! on a mineral (see RTDAESurfaceComplexation)
  associate(srfcplx => reaction%surface_complexation)
    ! (srfcplxrxn_stoich_flag marks only stoichiometries above one)
    if (srfcplx%nsrfcplx > 0) then
      if (any(srfcplx%srfcplx_free_site_stoich(1:srfcplx%nsrfcplx) /= &
              1.d0)) &
        call TSModelListAdd(unsupported,'surface complexation with a &
                            &free-site stoichiometry other than one')
    endif
  end associate
  if (associated(rxn_sandbox_list) .or. associated(carbon_sandbox_list) .or. &
      associated(clmrxn_list)) &
    call TSModelListAdd(unsupported,'reaction sandboxes')
  if (reaction%update_porosity .or. reaction%update_tortuosity .or. &
      reaction%update_permeability .or. &
      ReactionMnrlAnyUpdatePorosity(reaction%mineral)) &
    call TSModelListAdd(unsupported,'UPDATE_POROSITY/TORTUOSITY/PERMEABILITY')
  if (reaction%update_armor_mineral_surface) &
    call TSModelListAdd(unsupported,'mineral surface armoring')
  if (reaction%use_activity_h2o) call TSModelListAdd(unsupported,'ACTIVITY_H2O')
  if (reaction%calculate_water_age .or. reaction%calculate_tracer_age) &
    call TSModelListAdd(unsupported,'water/tracer age')
  if (reaction%nauxiliary > 0) &
    call TSModelListAdd(unsupported,'auxiliary reaction data')
  if (associated(reaction%mineral%nucleation_array)) &
    call TSModelListAdd(unsupported,'mineral nucleation')
  if (reaction%mineral%nkinmnrl > 0) then
    if (any(reaction%mineral%kinmnrl_tst_itype(1:reaction%mineral%nkinmnrl) &
            /= MINERAL_KINETICS_TST_SIMPLE)) &
      call TSModelListAdd(unsupported,'mineral rate laws with prefactors, &
        &Temkin constants, activation energies, thresholds or rate limiters')
  endif
  if (.not.option%transport%isothermal_reaction) &
    call TSModelListAdd(unsupported,'temperature-dependent reaction &
                                    &coefficients')
  if (option%use_sc) call TSModelListAdd(unsupported,'SECONDARY_CONTINUUM')
  if (option%transport%use_np) &
    call TSModelListAdd(unsupported,'Nernst-Planck transport')
  ! coupled, porosity is flow state; otherwise it is the initial one
  if (option%flow%transient_porosity .and. .not.flow_coupled) &
    call TSModelListAdd(unsupported,'transient porosity (soil &
      &compressibility, evolving STRATA) without MODE RICHARDS_DAE')
  if (associated(patch%prescribed_condition_list%first)) &
    call TSModelListAdd(unsupported,'PRESCRIBED_CONDITION')
  if (.not.PetscObjectIsNull(this%realization%field%tran_mass_transfer)) &
    call TSModelListAdd(unsupported,'transport mass transfer')
  if (option%compute_mass_balance_new .or. option%transport%store_fluxes) &
    call TSModelListAdd(unsupported,'mass balance and integral flux output')
  if (option%restart_flag .or. associated(option%checkpoint)) &
    call TSModelListAdd(unsupported,'checkpoint/restart')
  if (option%transport%steady_state) &
    call TSModelListAdd(unsupported,'steady state transport')

  if (len_trim(unsupported) > 0) then
    option%io_buffer = 'MODE GIRT_DAE does not yet support ' // &
      trim(unsupported) // '.'
    call PrintErrMsg(option)
  endif

end subroutine RTDAECheckCapabilities

! ************************************************************************** !

subroutine RTDAEBeginSegment(this,t_start,t_end,discontinuous)
  !
  ! The forcing of [t_start, t_end] (RTDAEForcingBegin), recorded after
  ! RecordForcing.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: t_start
  PetscReal :: t_end
  PetscBool :: discontinuous

  PetscBool :: coupled

  coupled = associated(this%fs)
  call RTDAEForcingBegin(this%forcing,coupled,t_start,t_end,discontinuous)
  this%segment_id = this%forcing%seg%id
  call this%RecordSegment(t_end)

end subroutine RTDAEBeginSegment

! ************************************************************************** !

subroutine RTDAEStoreRecord(this,k)
  !
  ! Keeps the current forcing segment as record k (see TS_Model_class).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: k

  call RTDAEForcingStore(this%forcing,k)

end subroutine RTDAEStoreRecord

! ************************************************************************** !

subroutine RTDAERestoreRecord(this,k)
  !
  ! Makes record k the current forcing segment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: k

  call RTDAEForcingRestore(this%forcing,k)
  this%segment_id = this%forcing%seg%id
  this%prepared = PETSC_FALSE

end subroutine RTDAERestoreRecord

! ************************************************************************** !

subroutine RTDAEInvalidateCache(this)
  !
  ! Reaction parameters changed in place: the prepared chemistry expires.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this

  this%prepared = PETSC_FALSE

end subroutine RTDAEInvalidateCache

! ************************************************************************** !

subroutine RTDAEResetToInitial(this)
  !
  ! The mineral regimes and transition counts of InitialState, and no
  ! preparations yet (-rt_dae_fail_closure), for a rerun.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  class(rt_dae_model_type) :: this

  if (allocated(this%mnrl_regime0)) this%mnrl_regime = this%mnrl_regime0
  this%num_depletions = 0
  this%num_reappearances = 0
  this%num_preparations = 0
  this%depletion_residual = 0.d0
  this%forcing%forcing_time = -MAX_DOUBLE
  this%forcing%forcing_segment = -1
  this%forcing%forcing_flow_id = -1
  call this%InvalidateCache()

end subroutine RTDAEResetToInitial

! ************************************************************************** !

subroutine RTDAEForcingAt(this,time)
  !
  ! The forcing at time, from the segment or, coupled, from fs, and the
  ! transport coefficients, recomputed only when the forcing changes.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reactive_Transport_module, only : RTComputeTransportCoefs

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time

  PetscBool :: changed

  call this%ReplayAt(time)
  if (associated(this%fs)) then
    call RTDAEForcingFromFlowState(this%forcing,this%fs,this%maux,changed)
  else
    call RTDAEForcingInterpolate(this%forcing,time,this%maux,changed)
  endif
  if (.not.changed) return

  call RTComputeTransportCoefs(this%realization%patch,this%option, &
                               this%forcing%gaux,this%forcing%gaux_bc, &
                               this%maux,this%forcing%vel_int, &
                               this%forcing%vel_bc,this%forcing%ccv, &
                               this%tran_coefs_int,this%tran_coefs_bc)
  call RTDAEFluxCoefficients(this)

end subroutine RTDAEForcingAt


! ************************************************************************** !

subroutine RTDAEFluxCoefficients(this)
  !
  ! Coefficients [L/s] of the connection fluxes at the forcing in gaux,
  ! vel_* and tran_coefs_*; computed once per forcing time.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Transport_module

  implicit none

  class(rt_dae_model_type) :: this

  type(patch_type), pointer :: patch
  type(reactive_transport_param_type), pointer :: rt_parameter
  PetscInt :: iconn, up, dn, ghosted_id
  PetscReal :: qsrc(2), coef_in(2), coef_out(2)
  PetscReal :: T_up(this%naqcomp,this%nphase)
  PetscReal :: T_dn(this%naqcomp,this%nphase)

  patch => this%realization%patch
  rt_parameter => patch%aux%RT%rt_parameter

  do iconn = 1, this%conn%nconn_int
    up = this%conn%int_up(iconn)
    dn = this%conn%int_dn(iconn)
    this%coef_up(:,iconn) = 0.d0
    this%coef_dn(:,iconn) = 0.d0
    if (patch%imat(up) <= 0 .or. patch%imat(dn) <= 0) cycle
    call TFluxCoef(rt_parameter,this%forcing%gaux(up), &
                   this%forcing%gaux(dn),this%option, &
                   this%conn%int_area(iconn),this%forcing%vel_int(:,iconn), &
                   this%tran_coefs_int(:,:,iconn), &
                   this%conn%int_dist(-1,iconn), &
                   PETSC_TRUE,T_up,T_dn)
    this%coef_up(:,iconn) = T_up(:,1)
    this%coef_dn(:,iconn) = T_dn(:,1)
  enddo

  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    this%coef_bc_up(:,iconn) = 0.d0
    this%coef_bc_dn(:,iconn) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    call TFluxCoefBC(this%bc_itype(iconn),rt_parameter, &
                     this%forcing%gaux_bc(iconn), &
                     this%forcing%gaux(ghosted_id),this%option, &
                     this%conn%bc_area(iconn),this%forcing%vel_bc(:,iconn), &
                     this%tran_coefs_bc(:,:,iconn),0.5d0,T_up,T_dn)
    this%coef_bc_up(:,iconn) = T_up(:,1)
    this%coef_bc_dn(:,iconn) = T_dn(:,1)
  enddo

  do iconn = 1, this%conn%nconn_ss
    ghosted_id = this%conn%ss_ghost(iconn)
    this%coef_ss_in(iconn) = 0.d0
    this%coef_ss_out(iconn) = 0.d0
    if (patch%imat(ghosted_id) <= 0) cycle
    ! uncoupled there is no flow: the sources carry no water
    qsrc = 0.d0
    if (associated(this%fs)) qsrc(1) = this%fs%qvol_ss(iconn)
    call TSrcSinkCoef(rt_parameter,this%forcing%gaux(ghosted_id),qsrc, &
                      this%ss_itype(iconn),coef_in,coef_out)
    this%coef_ss_in(iconn) = coef_in(1)
    this%coef_ss_out(iconn) = coef_out(1)
  enddo

end subroutine RTDAEFluxCoefficients

! ************************************************************************** !

subroutine RTDAEInternalFlux(this,iconn,flux,dflux_up,dflux_dn)
  !
  ! Aqueous flux [mol/s] across internal connection iconn, positive from
  ! int_up to int_dn, at the last Prepare, and its derivatives with respect
  ! to the two cells' c (leading naqcomp x naqcomp blocks).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: iconn
  PetscReal, optional :: flux(:)
  PetscReal, optional :: dflux_up(:,:), dflux_dn(:,:)

  PetscInt :: up, dn, i

  up = this%conn%int_up(iconn)
  dn = this%conn%int_dn(iconn)
  if (present(flux)) &
    flux = this%coef_up(:,iconn)*this%tot(:,up) + &
           this%coef_dn(:,iconn)*this%tot(:,dn)
  if (present(dflux_up)) then
    do i = 1, this%naqcomp
      dflux_up(i,1:this%naqcomp) = this%coef_up(i,iconn)*this%dtot(i,:,up)
    enddo
  endif
  if (present(dflux_dn)) then
    do i = 1, this%naqcomp
      dflux_dn(i,1:this%naqcomp) = this%coef_dn(i,iconn)*this%dtot(i,:,dn)
    enddo
  endif

end subroutine RTDAEInternalFlux

! ************************************************************************** !

subroutine RTDAEBoundaryFlux(this,iconn,flux,dflux)
  !
  ! Aqueous flux [mol/s] into the cell of boundary connection iconn, and
  ! its derivative with respect to the cell's c (leading naqcomp block).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: iconn
  PetscReal, optional :: flux(:)
  PetscReal, optional :: dflux(:,:)

  PetscInt :: cell, i

  cell = this%conn%bc_ghost(iconn)
  if (present(flux)) &
    flux = this%coef_bc_up(:,iconn)*this%forcing%seg%total_bc(:,iconn) + &
           this%coef_bc_dn(:,iconn)*this%tot(:,cell)
  if (present(dflux)) then
    do i = 1, this%naqcomp
      dflux(i,1:this%naqcomp) = this%coef_bc_dn(i,iconn)*this%dtot(i,:,cell)
    enddo
  endif

end subroutine RTDAEBoundaryFlux

! ************************************************************************** !

subroutine RTDAESourceFlux(this,iconn,flux,dflux)
  !
  ! Aqueous flux [mol/s] out of the cell of source/sink connection iconn
  ! (injection negative), and its derivative with respect to the cell's c
  ! (leading naqcomp block).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: iconn
  PetscReal, optional :: flux(:)
  PetscReal, optional :: dflux(:,:)

  PetscInt :: cell, i

  cell = this%conn%ss_ghost(iconn)
  if (present(flux)) &
    flux = this%coef_ss_in(iconn)*this%tot(:,cell) + &
           this%coef_ss_out(iconn)*this%forcing%seg%total_ss(:,iconn)
  if (present(dflux)) then
    do i = 1, this%naqcomp
      dflux(i,1:this%naqcomp) = this%coef_ss_in(iconn)*this%dtot(i,:,cell)
    enddo
  endif

end subroutine RTDAESourceFlux

! ************************************************************************** !

subroutine RTDAECellChemistry(this,ghosted_id,y_c,derivatives,status)
  !
  ! Converged local chemistry of one cell from its c coordinates; fills
  ! scratch, tot(:,ghosted_id) and, if requested, dtot(:,:,ghosted_id)
  ! with respect to the c coordinates (log chain applied).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: y_c(:)
  PetscBool :: derivatives
  PetscInt :: status

  PetscInt :: k
  PetscReal :: c(this%ncomp)

  if (this%log_formulation) then
    c = exp(y_c(1:this%ncomp))
    ! an underflow would leave dc/dy = 0: a singular Jacobian column
    if (any(c(1:this%naqcomp) < tiny(1.d0))) then
      status = RT_DAE_CHEM_FAILED
      return
    endif
  else
    c = y_c(1:this%ncomp)
  endif
  this%scratch%pri_molal = c(1:this%naqcomp)
  if (this%nimcomp > 0) &
    this%scratch%immobile = c(this%naqcomp+1:this%ncomp)
  call RTDAEAqueousClosure(this%scratch,this%forcing%gaux(ghosted_id), &
                           this%reaction,this%option,this%gam,status)
  if (status /= RT_DAE_CHEM_OK) return
  this%tot(:,ghosted_id) = this%scratch%total(:,1)
  if (allocated(this%pri_act_coef)) &
    this%pri_act_coef(:,ghosted_id) = this%scratch%pri_act_coef
  if (derivatives) then
    this%dtot(:,:,ghosted_id) = this%scratch%aqueous%dtotal(:,:,1)
    if (this%log_formulation) then
      do k = 1, this%naqcomp
        this%dtot(:,k,ghosted_id) = this%dtot(:,k,ghosted_id)*c(k)
      enddo
    endif
  endif

end subroutine RTDAECellChemistry

! ************************************************************************** !

subroutine RTDAECellReactions(this,ghosted_id,Kin,regime,Rn,jK,derivatives, &
                              dRn_dc,dRn_dK,djK_dc,djK_dK,mnrl_rate, &
                              mnrl_area,mnrl_lnQK)
  !
  ! Reaction terms of one active cell at the chemistry in scratch.
  !   Rn(ncomp)  sink of the transported components [mol/s]
  !   jK(nkin)   kinetic transfer rates into the kinetic inventories [mol/s]
  ! Derivatives are with respect to the physical c (molality or immobile
  ! concentration); the caller applies the log chain.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_Mineral_Aux_module
  use Reaction_Mineral_module, only : ReactionMnrlSpecSurfArea, &
                                      ReactionMnrlTSTSimpleRate
  use Reaction_module, only : RRadioactiveDecay, RGeneral
  use Reaction_Microbial_module, only : ReactionMicrobRate
  use Reaction_Microbial_Aux_module, only : MICROBIAL_ACTIVITY
  use Reaction_Immobile_module, only : ReactionImDecay

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: Kin(:)
  PetscInt :: regime(:)
  PetscReal :: Rn(:), jK(:)
  PetscBool :: derivatives
  PetscReal :: dRn_dc(:,:), dRn_dK(:,:), djK_dc(:,:), djK_dK(:,:)
  PetscReal :: mnrl_rate(:), mnrl_area(:), mnrl_lnQK(:)

  type(mineral_type), pointer :: mineral
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: imnrl, i, icomp, naq
  PetscReal :: volume, volfrac, area, darea_dvf
  PetscReal :: ln_act(this%naqcomp)
  PetscBool :: may_dissolve
  PetscReal :: lnQK, QK, rate_constant, rate, j_rate
  PetscReal :: dj_dm(1,this%naqcomp), dj_dK
  PetscReal :: Jleg(this%ncomp,this%ncomp)
  PetscReal :: Seq(this%naqcomp), dSeq_dm(this%naqcomp,this%naqcomp)
  PetscReal :: dSeq_dvf(this%naqcomp)
  PetscReal :: kr, fr, dsum_dK
  PetscInt :: ikinmr, irxn, irate, iq, iK, isurf, jcomp

  reaction => this%reaction
  mineral => reaction%mineral
  naq = this%naqcomp
  volume = this%maux(ghosted_id)%volume

  Rn = 0.d0
  jK = 0.d0
  if (derivatives) then
    dRn_dc = 0.d0
    dRn_dK = 0.d0
    djK_dc = 0.d0
    djK_dK = 0.d0
  endif

  ! kinetic minerals: transition-state rate, precipitation positive
  if (this%nmnrl > 0) &
    ln_act = log(this%scratch%pri_molal(1:naq)* &
                 this%scratch%pri_act_coef(1:naq))
  do imnrl = 1, this%nmnrl
    volfrac = RTDAEMineralVolFrac(this,imnrl,Kin(imnrl),volume)
    area = this%mnrl_area0(imnrl,ghosted_id)
    darea_dvf = 0.d0
    ! porosity factor 1: GIRT_DAE does not update porosity
    if (mineral%update_surface_area) &
      call ReactionMnrlSpecSurfArea(mineral,imnrl, &
                                    this%mnrl_area0(imnrl,ghosted_id), &
                                    this%mnrl_vf0(imnrl,ghosted_id), &
                                    volfrac,1.d0,area,darea_dvf)
    ! an absent mineral cannot dissolve
    may_dissolve = regime(imnrl) /= RT_DAE_MNRL_ABSENT
    ! mol/m^3 bulk/s
    call ReactionMnrlTSTSimpleRate(mineral,imnrl,ln_act, &
                                   this%scratch%ln_act_h2o,area, &
                                   may_dissolve,lnQK,QK,rate_constant,rate)
    mnrl_lnQK(imnrl) = lnQK
    mnrl_rate(imnrl) = rate
    mnrl_area(imnrl) = area
    j_rate = volume*rate
    jK(imnrl) = jK(imnrl) + j_rate
    do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
      icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
      Rn(icomp) = Rn(icomp) + &
        mineral%kinmnrlstoich_in_residual(i,imnrl)*j_rate
    enddo
    if (.not.derivatives) cycle
    dj_dm = 0.d0
    do i = 1, mineral%kinmnrlspecid(0,imnrl)
      icomp = mineral%kinmnrlspecid(i,imnrl)
      dj_dm(1,icomp) = dj_dm(1,icomp) + volume*area*rate_constant*QK* &
        mineral%kinmnrlstoich(i,imnrl)/this%scratch%pri_molal(icomp)
    enddo
    call RTDAEGammaChain(dj_dm,1,naq,this%scratch%pri_molal,this%gam)
    ! d(V r)/dK = V k (QK-1) darea/dvf * Vbar/V
    dj_dK = rate_constant*(QK-1.d0)*darea_dvf* &
            mineral%kinmnrl_molar_vol(imnrl)
    djK_dc(imnrl,1:naq) = djK_dc(imnrl,1:naq) + dj_dm(1,:)
    djK_dK(imnrl,imnrl) = djK_dK(imnrl,imnrl) + dj_dK
    do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
      icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
      dRn_dc(icomp,1:naq) = dRn_dc(icomp,1:naq) + &
        mineral%kinmnrlstoich_in_residual(i,imnrl)*dj_dm(1,:)
      dRn_dK(icomp,imnrl) = dRn_dK(icomp,imnrl) + &
        mineral%kinmnrlstoich_in_residual(i,imnrl)*dj_dK
    enddo
  enddo

  ! multirate surface complexation: for rate r and sorbing component q,
  ! K = V*S_rq relaxes to its share of the equilibrium sorbed amount,
  !   dK/dt = k_r (f_r V S_eq,q(c,VF) - K),
  ! and the same transfer leaves the transported component
  do ikinmr = 1, this%nmr
    irxn = reaction%surface_complexation%kinmrsrfcplxrxn_to_srfcplxrxn(ikinmr)
    call RTDAESurfaceComplexation(this,ghosted_id,irxn,Seq,dSeq_dm, &
                                  dSeq_dvf,isurf)
    do irate = 1, reaction%surface_complexation%kinmr_nrate(ikinmr)
      kr = reaction%surface_complexation%kinmr_rate(irate,ikinmr)
      fr = reaction%surface_complexation%kinmr_frac(irate,ikinmr)
      do iq = 1, this%mr_ncomp(ikinmr)
        icomp = this%mr_comp(iq,ikinmr)
        iK = RTDAEMRIndex(this,ikinmr,irate,iq)
        j_rate = kr*(fr*volume*Seq(icomp)-Kin(iK))
        jK(iK) = jK(iK) + j_rate
        Rn(icomp) = Rn(icomp) + j_rate
        if (.not.derivatives) cycle
        djK_dc(iK,1:naq) = djK_dc(iK,1:naq) + kr*fr*volume*dSeq_dm(icomp,:)
        dRn_dc(icomp,1:naq) = dRn_dc(icomp,1:naq) + &
                              kr*fr*volume*dSeq_dm(icomp,:)
        djK_dK(iK,iK) = djK_dK(iK,iK) - kr
        dRn_dK(icomp,iK) = dRn_dK(icomp,iK) - kr
        if (isurf > 0) then
          ! dVF/dK_mineral = Vbar/V
          djK_dK(iK,isurf) = djK_dK(iK,isurf) + kr*fr*dSeq_dvf(icomp)* &
                             mineral%kinmnrl_molar_vol(isurf)
          dRn_dK(icomp,isurf) = dRn_dK(icomp,isurf) + kr*fr* &
                                dSeq_dvf(icomp)* &
                                mineral%kinmnrl_molar_vol(isurf)
        endif
      enddo
    enddo
  enddo

  ! instantaneous kinetics without kinetic state, at the chemistry in
  ! scratch (RTDAEAccumulation leaves the sorbed totals there) and the
  ! private flow/material state
  if (reaction%nradiodecay_rxn > 0) then
    ! the aqueous and sorbed totals, with their complete dT/dm
    Jleg = 0.d0
    call RRadioactiveDecay(Rn,Jleg,derivatives,this%scratch, &
                           this%forcing%gaux(ghosted_id), &
                           this%maux(ghosted_id), &
                           reaction,this%option)
    if (derivatives) then
      dRn_dc = dRn_dc + Jleg
      ! a parent sorbed on a kinetic mineral's surface: dVF/dK = Vbar/V
      do irxn = 1, reaction%nradiodecay_rxn
        jcomp = reaction%radiodecayforwardspecid(irxn)
        do imnrl = 1, this%nmnrl
          dsum_dK = this%sorb_dvf(jcomp,imnrl)* &
                    reaction%mineral%kinmnrl_molar_vol(imnrl)
          if (dsum_dK == 0.d0) cycle
          do i = 1, reaction%radiodecayspecid(0,irxn)
            icomp = reaction%radiodecayspecid(i,irxn)
            dRn_dK(icomp,imnrl) = dRn_dK(icomp,imnrl) - &
              reaction%radiodecaystoich(i,irxn)* &
              reaction%radiodecay_kf(irxn)*dsum_dK
          enddo
        enddo
      enddo
    endif
  endif
  if (reaction%ngeneral_rxn > 0) then
    ! activity based: add the activity-coefficient sensitivity
    Jleg = 0.d0
    call RGeneral(Rn,Jleg,derivatives,this%scratch, &
                  this%forcing%gaux(ghosted_id),this%maux(ghosted_id), &
                  reaction,this%option)
    if (derivatives) then
      call RTDAEGammaChain(Jleg,this%ncomp,naq,this%scratch%pri_molal, &
                           this%gam)
      dRn_dc = dRn_dc + Jleg
    endif
  endif
  if (reaction%microbial%nrxn > 0) then
    Jleg = 0.d0
    call ReactionMicrobRate(Rn,Jleg,derivatives,this%scratch, &
                            this%forcing%gaux(ghosted_id), &
                            this%maux(ghosted_id), &
                            reaction,this%option)
    if (derivatives) then
      ! in activity units the rate reads the primary activities
      if (reaction%microbial%concentration_units == MICROBIAL_ACTIVITY) then
        call RTDAEGammaChain(Jleg,this%ncomp,naq,this%scratch%pri_molal, &
                             this%gam)
      endif
      dRn_dc = dRn_dc + Jleg
    endif
  endif
  if (reaction%immobile%ndecay_rxn > 0) then
    Jleg = 0.d0
    call ReactionImDecay(Rn,Jleg,derivatives,this%scratch, &
                         this%forcing%gaux(ghosted_id),this%maux(ghosted_id), &
                         reaction,this%option)
    if (derivatives) dRn_dc = dRn_dc + Jleg
  endif

end subroutine RTDAECellReactions

! ************************************************************************** !

subroutine RTDAEAccumulation(this,ghosted_id,A,gross,dA_dc,dA_dk)
  !
  ! Equilibrium inventory A [mol] of one cell, RTAccumulation's at the
  ! chemistry in scratch, whose sorbed totals and dS/dm this fills
  ! (RRadioactiveDecay reads them too; dS/dVF goes to sorb_dvf); its gross
  ! magnitude (species amounts summed with |nu|), and, when present, dA_dc
  ! (log chain applied) and dA_dk.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_module, only : RTAccumulation
  use Reaction_Surface_Complexation_Aux_module, only : ROCK_SURFACE

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: A(:)
  PetscReal :: gross(:)
  PetscReal, optional :: dA_dc(:,:)
  PetscReal, optional :: dA_dk(:,:)

  PetscReal :: psv, volume, den_kg_per_L
  PetscReal :: S(this%naqcomp), dS_dm(this%naqcomp,this%naqcomp)
  PetscReal :: dS_dvf(this%naqcomp)
  PetscInt :: i, icplx, icomp, ieqrxn, irxn, imnrl

  volume = this%maux(ghosted_id)%volume
  ! L water per cell
  psv = RTDAEPoreLiters(this%maux(ghosted_id)%porosity* &
                        this%forcing%gaux(ghosted_id)%sat(1),volume)
  gross(1:this%naqcomp) = this%scratch%pri_molal
  do icplx = 1, this%reaction%neqcplx
    do i = 1, this%reaction%eqcplxspecid(0,icplx)
      icomp = this%reaction%eqcplxspecid(i,icplx)
      gross(icomp) = gross(icomp) + &
        abs(this%reaction%eqcplxstoich(i,icplx))*this%scratch%sec_molal(icplx)
    enddo
  enddo
  den_kg_per_L = this%forcing%gaux(ghosted_id)%den_kg(1)*1.d-3
  gross(1:this%naqcomp) = psv*den_kg_per_L*gross(1:this%naqcomp)
  if (this%nimcomp > 0) &
    gross(this%naqcomp+1:this%ncomp) = abs(volume*this%scratch%immobile(:))
  if (present(dA_dc)) then
    dA_dc = 0.d0
    dA_dc(1:this%naqcomp,1:this%naqcomp) = psv*this%dtot(:,:,ghosted_id)
    do i = 1, this%nimcomp
      if (this%log_formulation) then
        dA_dc(this%naqcomp+i,this%naqcomp+i) = &
          volume*this%scratch%immobile(i)
      else
        dA_dc(this%naqcomp+i,this%naqcomp+i) = volume
      endif
    enddo
  endif
  if (present(dA_dk)) dA_dk = 0.d0

  ! equilibrium surface complexation: the sorbed totals S [mol/m^3 bulk]
  if (this%nmnrl > 0) this%sorb_dvf = 0.d0
  if (allocated(this%sorb_rock)) this%sorb_rock(:,ghosted_id) = 0.d0
  if (this%eq_srfcplx) then
    this%scratch%total_sorb_eq = 0.d0
    this%scratch%dtotal_sorb_eq = 0.d0
    do ieqrxn = 1, this%reaction%surface_complexation%neqsrfcplxrxn
      irxn = this%reaction%surface_complexation% &
               eqsrfcplxrxn_to_srfcplxrxn(ieqrxn)
      call RTDAESurfaceComplexation(this,ghosted_id,irxn,S,dS_dm,dS_dvf, &
                                    imnrl)
      this%scratch%total_sorb_eq = this%scratch%total_sorb_eq + S
      this%scratch%dtotal_sorb_eq = this%scratch%dtotal_sorb_eq + dS_dm
      if (allocated(this%sorb_rock)) then
        if (this%reaction%surface_complexation% &
              srfcplxrxn_surf_type(irxn) == ROCK_SURFACE) &
          this%sorb_rock(:,ghosted_id) = this%sorb_rock(:,ghosted_id) + S
      endif
      gross(1:this%naqcomp) = gross(1:this%naqcomp) + volume*abs(S)
      if (present(dA_dc)) then
        if (this%log_formulation) then
          do i = 1, this%naqcomp
            dS_dm(:,i) = dS_dm(:,i)*this%scratch%pri_molal(i)
          enddo
        endif
        dA_dc(1:this%naqcomp,1:this%naqcomp) = &
          dA_dc(1:this%naqcomp,1:this%naqcomp) + volume*dS_dm
      endif
      if (imnrl > 0) then
        this%sorb_dvf(:,imnrl) = this%sorb_dvf(:,imnrl) + dS_dvf
        ! dVF/dK = Vbar/V
        if (present(dA_dk)) &
          dA_dk(1:this%naqcomp,imnrl) = dA_dk(1:this%naqcomp,imnrl) + &
            dS_dvf*this%reaction%mineral%kinmnrl_molar_vol(imnrl)
      endif
    enddo
  endif

  call RTAccumulation(this%scratch,this%forcing%gaux(ghosted_id), &
                      this%maux(ghosted_id),this%reaction,this%option,A)

end subroutine RTDAEAccumulation

! ************************************************************************** !

subroutine RTDAEPrepare(this,time,Y,derivatives,status)
  !
  ! Forcing at time, ghosted state, the closure of every ghosted active
  ! cell, and for local active cells the accumulation, reaction and kinetic
  ! terms (with derivatives if requested). status reports a failed closure
  ! in a local cell; a ghost cell's failure is reported by its owner.
  !
  ! Reused while (time, Y) is unchanged within a segment and mineral
  ! regime: Newton linearizes and tests convergence at the state it just
  ! evaluated. The key compares the values of Y: SNES copies its trial
  ! vector into the iterate. Collective (ghost update).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Discretization_module
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscBool :: derivatives
  PetscInt :: status

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: yl_p(:)
  PetscInt :: ghosted_id, local_id, offset, cell_status, nc, nk, k, imnrl
  PetscReal :: c(this%ncomp)
  PetscBool :: same
  PetscErrorCode :: ierr

  call this%ReplayAt(time)
  ! every rank takes the same branch: VecEqual is collective
  same = PETSC_TRUE
  if (associated(this%fs)) same = (this%prepared_flow_id == this%fs%id)
  if (this%prepared .and. same .and. time == this%prepared_time .and. &
      this%forcing%seg%id == this%prepared_segment .and. &
      (this%prepared_derivatives .or. .not.derivatives)) then
    call VecEqual(Y,this%Y_prepared,same,ierr);CHKERRQ(ierr)
    if (same) then
      status = this%prepared_status
      return
    endif
  endif

  patch => this%realization%patch
  grid => patch%grid
  nc = this%ncomp
  nk = this%nkin
  status = TS_MODEL_OK
  this%num_preparations = this%num_preparations + 1

  call RTDAEForcingAt(this,time)
  call DiscretizationGlobalToLocal(this%realization%discretization,Y, &
                                   this%Y_loc,NDAEDOF)
  call VecGetArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)
  if (this%debug .and. derivatives) call RTDAEDebugClosure(this,yl_p)
  do ghosted_id = 1, this%ngmax
    this%cell_ok(ghosted_id) = PETSC_TRUE
    if (grid%nG2L(ghosted_id) < 0) cycle
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (ghosted_id-1)*this%ndof
    call RTDAECellChemistry(this,ghosted_id,yl_p(offset+1:offset+nc), &
                            derivatives,cell_status)
    if (this%num_preparations == this%fail_closure) &
      cell_status = RT_DAE_CHEM_FAILED
    if (cell_status /= RT_DAE_CHEM_OK) then
      if (this%debug) then
        write(*,'(a,i8,a,es12.4,a,*(es12.4))') ' GIRT_DAE closure failed: &
          &cell ',grid%nG2A(ghosted_id),' t=',time,' c=', &
          yl_p(offset+1:offset+nc)
      endif
      ! finite scratch keeps the remaining assembly well defined
      this%cell_ok(ghosted_id) = PETSC_FALSE
      this%tot(:,ghosted_id) = 0.d0
      if (derivatives) this%dtot(:,:,ghosted_id) = 0.d0
      if (grid%nG2L(ghosted_id) > 0) status = TS_MODEL_DOMAIN_ERROR
      cycle
    endif
    local_id = grid%nG2L(ghosted_id)
    if (local_id <= 0) cycle
    ! mineral volume fractions (site densities) from the inventories
    do imnrl = 1, this%nmnrl
      this%scratch%mnrl_volfrac(imnrl) = &
        RTDAEMineralVolFrac(this,imnrl,yl_p(offset+2*nc+imnrl), &
                            this%maux(ghosted_id)%volume)
    enddo
    if (derivatives) then
      call RTDAEAccumulation(this,ghosted_id,this%acc(:,local_id), &
                             this%acc_gross(:,local_id), &
                             this%dacc_dc(:,:,local_id), &
                             this%dacc_dk(:,:,local_id))
    else
      call RTDAEAccumulation(this,ghosted_id,this%acc(:,local_id), &
                             this%acc_gross(:,local_id))
    endif
    call RTDAECellReactions(this,ghosted_id, &
                            yl_p(offset+2*nc+1:offset+2*nc+nk), &
                            this%mnrl_regime(:,local_id), &
                            this%rn(:,local_id),this%jk(:,local_id), &
                            derivatives, &
                            this%drn_dc(:,:,local_id), &
                            this%drn_dk(:,:,local_id), &
                            this%djk_dc(:,:,local_id), &
                            this%djk_dk(:,:,local_id), &
                            this%mnrl_rate(:,local_id), &
                            this%mnrl_area(:,local_id), &
                            this%mnrl_lnQK(:,local_id))
    if (derivatives .and. this%log_formulation) then
      c(1:this%naqcomp) = this%scratch%pri_molal
      if (this%nimcomp > 0) c(this%naqcomp+1:nc) = this%scratch%immobile
      do k = 1, nc
        this%drn_dc(:,k,local_id) = this%drn_dc(:,k,local_id)*c(k)
        if (nk > 0) this%djk_dc(:,k,local_id) = &
                      this%djk_dc(:,k,local_id)*c(k)
      enddo
    endif
  enddo
  call VecRestoreArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)

  call VecCopy(Y,this%Y_prepared,ierr);CHKERRQ(ierr)
  this%prepared = PETSC_TRUE
  this%prepared_time = time
  this%prepared_segment = this%forcing%seg%id
  this%prepared_derivatives = derivatives
  this%prepared_status = status
  if (associated(this%fs)) this%prepared_flow_id = this%fs%id

end subroutine RTDAEPrepare

! ************************************************************************** !

subroutine RTDAEDebugClosure(this,yl_p)
  !
  ! -rt_dae_debug: central differences of the closure totals of the first
  ! active cell against their analytical derivative.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: yl_p(:)

  PetscReal :: y(this%ncomp), yp(this%ncomp), tp(this%naqcomp)
  PetscReal :: tm(this%naqcomp), dt_an(this%naqcomp,this%naqcomp)
  PetscReal :: h, err
  PetscInt :: k, status, ghosted_id

  ghosted_id = this%realization%patch%grid%nL2G(1)
  y = yl_p((ghosted_id-1)*this%ndof+1:(ghosted_id-1)*this%ndof+this%ncomp)
  call RTDAECellChemistry(this,ghosted_id,y,PETSC_TRUE,status)
  dt_an = this%dtot(:,:,ghosted_id)
  do k = 1, this%naqcomp
    h = 1.d-6*max(abs(y(k)),1.d0)
    if (.not.this%log_formulation) h = 1.d-6*abs(y(k))
    yp = y
    yp(k) = y(k) + h
    call RTDAECellChemistry(this,ghosted_id,yp,PETSC_FALSE,status)
    tp = this%tot(:,ghosted_id)
    yp(k) = y(k) - h
    call RTDAECellChemistry(this,ghosted_id,yp,PETSC_FALSE,status)
    tm = this%tot(:,ghosted_id)
    err = maxval(abs((tp-tm)/(2.d0*h)-dt_an(:,k)))/ &
          max(maxval(abs(dt_an(:,k))),1.d-300)
    write(*,'(a,i4,a,es10.3,a,i4)') ' GIRT_DAE closure FD check: column ', &
      k,' rel. error ',err,' gamma its ',this%gam%iterations
  enddo
  call RTDAECellChemistry(this,ghosted_id,y,PETSC_TRUE,status)

end subroutine RTDAEDebugClosure

! ************************************************************************** !

subroutine RTDAEEvaluate(this,time,Y,Ydot,F,status)
  !
  ! F(t,Y,Ydot) of the whole grid. Prepares derivatives too, for Newton's
  ! linearization at this state to reuse (see RTDAEPrepare).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  Vec :: F
  PetscInt :: status

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: y_p(:), ydot_p(:), f_p(:)
  PetscInt :: local_id, ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: iconn, off, offN, offK, nc, nk, naq
  PetscReal :: flux(this%naqcomp)
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  nc = this%ncomp
  nk = this%nkin
  naq = this%naqcomp

  call RTDAEPrepare(this,time,Y,PETSC_TRUE,status)

  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Ydot,ydot_p,ierr);CHKERRQ(ierr)
  call VecGetArray(F,f_p,ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    off = (local_id-1)*this%ndof
    offN = off + nc
    offK = off + 2*nc
    if (patch%imat(ghosted_id) <= 0) then
      f_p(off+1:off+nc) = y_p(off+1:off+nc) - this%c_frozen(:,local_id)
      f_p(offN+1:off+this%ndof) = ydot_p(offN+1:off+this%ndof)
      cycle
    endif
    if (.not.this%cell_ok(ghosted_id)) then
      ! non-finite for callers other than SNES (which sets Inf itself on a
      ! domain error); a quiet NaN made without an invalid operation, which
      ! would trap under -fp_trap
      f_p(off+1:off+this%ndof) = ieee_value(1.d0,ieee_quiet_nan)
      cycle
    endif
    f_p(off+1:off+nc) = y_p(offN+1:offN+nc) - this%acc(:,local_id)
    f_p(offN+1:offN+nc) = ydot_p(offN+1:offN+nc) + this%rn(:,local_id)
    if (nk > 0) f_p(offK+1:offK+nk) = ydot_p(offK+1:offK+nk) - &
                                      this%jk(1:nk,local_id)
  enddo

  ! internal faces
  do iconn = 1, this%conn%nconn_int
    ghosted_id_up = this%conn%int_up(iconn)
    ghosted_id_dn = this%conn%int_dn(iconn)
    if (patch%imat(ghosted_id_up) <= 0 .or. &
        patch%imat(ghosted_id_dn) <= 0) cycle
    call RTDAEInternalFlux(this,iconn,flux)
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    if (local_id_up > 0) then
      offN = (local_id_up-1)*this%ndof + nc
      f_p(offN+1:offN+naq) = f_p(offN+1:offN+naq) + flux
    endif
    if (local_id_dn > 0) then
      offN = (local_id_dn-1)*this%ndof + nc
      f_p(offN+1:offN+naq) = f_p(offN+1:offN+naq) - flux
    endif
  enddo

  ! boundary faces
  do iconn = 1, this%conn%nconn_bc
    if (patch%imat(this%conn%bc_ghost(iconn)) <= 0) cycle
    call RTDAEBoundaryFlux(this,iconn,flux)
    offN = (this%conn%bc_local(iconn)-1)*this%ndof + nc
    f_p(offN+1:offN+naq) = f_p(offN+1:offN+naq) - flux
  enddo

  ! sources and sinks
  do iconn = 1, this%conn%nconn_ss
    if (patch%imat(this%conn%ss_ghost(iconn)) <= 0) cycle
    call RTDAESourceFlux(this,iconn,flux)
    offN = (this%conn%ss_local(iconn)-1)*this%ndof + nc
    f_p(offN+1:offN+naq) = f_p(offN+1:offN+naq) + flux
  enddo

  call VecRestoreArray(F,f_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Ydot,ydot_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)

end subroutine RTDAEEvaluate

! ************************************************************************** !

subroutine RTDAELinearize(this,time,Y,Ydot,shift,A,B,status)
  !
  ! J = dF/dY + shift*dF/dYdot at (time,Y), for any shift including zero;
  ! reuses an Evaluate at the same state but needs none. A is the condensed
  ! operator or a block matrix of ndof per cell; B is unused (the
  ! preconditioner is built from A).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module
  use Petsc_Utility_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  PetscReal :: shift
  Mat :: A
  Mat :: B
  PetscInt :: status

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(rt_dae_operator_type), pointer :: op
  PetscReal :: fb(this%ncomp,this%ncomp), fb2(this%ncomp,this%ncomp)
  PetscReal :: diag(this%ndof)
  PetscBool :: condensed
  PetscInt :: local_id, ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: iconn, nc, nk, i
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  nc = this%ncomp
  nk = this%nkin

  condensed = RTDAEOperatorIsCondensed(A)
  if (condensed) then
    call RTDAEOperatorGet(A,op)
    call RTDAEOperatorBeginLinearization(op,shift)
  else
    call MatZeroEntries(A,ierr);CHKERRQ(ierr)
    if (.not.allocated(this%blk)) then
      allocate(this%blk(this%ndof,this%ndof))
      allocate(this%blk2(this%ndof,this%ndof))
    endif
  endif

  call RTDAEPrepare(this,time,Y,PETSC_TRUE,status)

  ! cell blocks
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      diag(1:nc) = 1.d0
      diag(nc+1:this%ndof) = shift
      if (condensed) then
        call RTDAEOperatorSetDiagonalCell(op,local_id,diag)
        fb = 0.d0
        call AddFluxBlock(ghosted_id,ghosted_id,fb)
      else
        this%blk = 0.d0
        do i = 1, this%ndof
          this%blk(i,i) = diag(i)
        enddo
        call AddCellBlock(ghosted_id)
      endif
      cycle
    endif
    if (.not.this%cell_ok(ghosted_id)) then
      ! the residual is not finite here; any nonsingular block will do
      if (condensed) then
        diag = 1.d0
        call RTDAEOperatorSetDiagonalCell(op,local_id,diag)
        fb = 0.d0
        call AddFluxBlock(ghosted_id,ghosted_id,fb)
      else
        this%blk = 0.d0
        do i = 1, this%ndof
          this%blk(i,i) = 1.d0
        enddo
        call AddCellBlock(ghosted_id)
      endif
      cycle
    endif
    if (condensed) then
      call RTDAEOperatorSetCell(op,local_id,this%dacc_dc(:,:,local_id), &
                                this%dacc_dk(:,1:nk,local_id), &
                                this%drn_dk(:,1:nk,local_id), &
                                this%djk_dc(1:nk,:,local_id), &
                                this%djk_dk(1:nk,1:nk,local_id))
      call AddFluxBlock(ghosted_id,ghosted_id,this%drn_dc(:,:,local_id))
      cycle
    endif
    this%blk = 0.d0
    ! F_c = N - A
    this%blk(1:nc,1:nc) = -this%dacc_dc(:,:,local_id)
    do i = 1, nc
      this%blk(i,nc+i) = 1.d0
    enddo
    ! F_N = dN/dt + R (+ flux, added below)
    this%blk(nc+1:2*nc,1:nc) = this%drn_dc(:,:,local_id)
    do i = 1, nc
      this%blk(nc+i,nc+i) = shift
    enddo
    if (nk > 0) then
      this%blk(1:nc,2*nc+1:2*nc+nk) = -this%dacc_dk(:,1:nk,local_id)
      this%blk(nc+1:2*nc,2*nc+1:2*nc+nk) = this%drn_dk(:,1:nk,local_id)
      ! F_K = dK/dt - j
      this%blk(2*nc+1:2*nc+nk,1:nc) = -this%djk_dc(1:nk,:,local_id)
      this%blk(2*nc+1:2*nc+nk,2*nc+1:2*nc+nk) = -this%djk_dk(1:nk,1:nk,local_id)
      do i = 1, nk
        this%blk(2*nc+i,2*nc+i) = this%blk(2*nc+i,2*nc+i) + shift
      enddo
    endif
    call AddCellBlock(ghosted_id)
  enddo

  ! internal faces: flux rows N, columns c of both cells
  do iconn = 1, this%conn%nconn_int
    ghosted_id_up = this%conn%int_up(iconn)
    ghosted_id_dn = this%conn%int_dn(iconn)
    if (patch%imat(ghosted_id_up) <= 0 .or. &
        patch%imat(ghosted_id_dn) <= 0) cycle
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    fb = 0.d0
    fb2 = 0.d0
    call RTDAEInternalFlux(this,iconn,dflux_up=fb,dflux_dn=fb2)
    if (local_id_up > 0) then
      call AddFluxBlock(ghosted_id_up,ghosted_id_up,fb)
      call AddFluxBlock(ghosted_id_up,ghosted_id_dn,fb2)
    endif
    if (local_id_dn > 0) then
      call AddFluxBlock(ghosted_id_dn,ghosted_id_up,-fb)
      call AddFluxBlock(ghosted_id_dn,ghosted_id_dn,-fb2)
    endif
  enddo

  ! boundary faces (flux into the cell)
  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    fb = 0.d0
    call RTDAEBoundaryFlux(this,iconn,dflux=fb)
    call AddFluxBlock(ghosted_id,ghosted_id,-fb)
  enddo

  ! sources and sinks (flux out of the cell)
  do iconn = 1, this%conn%nconn_ss
    ghosted_id = this%conn%ss_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    fb = 0.d0
    call RTDAESourceFlux(this,iconn,dflux=fb)
    call AddFluxBlock(ghosted_id,ghosted_id,fb)
  enddo

  if (condensed) call RTDAEOperatorEndLinearization(op)
  ! the shell operator too: its state changes with the new values
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

contains

  subroutine AddCellBlock(ghosted)
    ! the full ndof block of a cell, in this%blk
    PetscInt :: ghosted
    call PUMSetValuesBlockedLocal(A,1,ghosted-1,1,ghosted-1,this%blk, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  end subroutine AddCellBlock

  subroutine AddFluxBlock(ghosted_row,ghosted_col,block)
    ! a coupling of rows N (of ghosted_row) to columns c (of ghosted_col)
    PetscInt :: ghosted_row, ghosted_col
    PetscReal :: block(:,:)
    if (condensed) then
      call PUMSetValuesBlockedLocal(op%G,1,ghosted_row-1,1,ghosted_col-1, &
                                    block,ADD_VALUES,ierr);CHKERRQ(ierr)
    else
      this%blk2 = 0.d0
      this%blk2(nc+1:2*nc,1:nc) = block
      call PUMSetValuesBlockedLocal(A,1,ghosted_row-1,1,ghosted_col-1, &
                                    this%blk2,ADD_VALUES,ierr);CHKERRQ(ierr)
    endif
  end subroutine AddFluxBlock

end subroutine RTDAELinearize

! ************************************************************************** !

subroutine RTDAECreateJacobian(this,J,created)
  !
  ! The condensed operator, unless -rt_dae_assembled_jacobian.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Discretization_module

  implicit none

  class(rt_dae_model_type) :: this
  Mat :: J
  PetscBool :: created

  Mat :: G
  MatType :: mat_type

  created = PETSC_FALSE
  if (this%assembled_jacobian) return
  ! the couplings between cells: block size ncomp, the flux stencil
  mat_type = MATBAIJ
  call DiscretizationCreateMatrix(this%realization%discretization,NTRANDOF, &
                                  mat_type,G,PETSC_FALSE,this%option)
  call RTDAEOperatorCreate(G,this%nkin,this%nlmax,this%option%mycomm,J)
  created = PETSC_TRUE

end subroutine RTDAECreateJacobian

! ************************************************************************** !

subroutine RTDAESetupLinearSolver(this,ksp,J,pc_type,shift_inblocks)
  !
  ! Elimination of N and K; pc_type preconditions the condensed system.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  KSP :: ksp
  Mat :: J
  PCType :: pc_type
  PetscBool :: shift_inblocks

  call RTDAEOperatorSetupKSP(ksp,J,pc_type,shift_inblocks)

end subroutine RTDAESetupLinearSolver

! ************************************************************************** !

subroutine RTDAEInitialState(this,time,Y)
  !
  ! c0 from the initial constraint, kinetic inventories from the initial
  ! mineral volume fractions, and N0 = A(c0,K0) with this model's closure.
  ! Must follow BeginSegment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module
  use Field_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal, pointer :: y_p(:), xx_p(:)
  PetscInt :: local_id, ghosted_id, off, nc, imnrl, status
  PetscReal :: Kres(max(this%nmnrl,1))
  PetscInt :: ikinmr, irate, iq
  PetscReal :: A(this%ncomp), gross(this%ncomp)
  PetscReal :: yc(this%ncomp)
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  field => this%realization%field
  rt_auxvars => patch%aux%RT%auxvars
  nc = this%ncomp

  ! immutable mineral reference data
  do ghosted_id = 1, this%ngmax
    do imnrl = 1, this%nmnrl
      this%mnrl_vf0(imnrl,ghosted_id) = rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl)
      this%mnrl_area0(imnrl,ghosted_id) = rt_auxvars(ghosted_id)%mnrl_area0(imnrl)
    enddo
  enddo

  call RTDAEForcingAt(this,time)
  ! the closures below overwrite the per-cell chemistry
  this%prepared = PETSC_FALSE
  call VecGetArray(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    off = (local_id-1)*this%ndof
    if (this%log_formulation) then
      yc = log(xx_p((local_id-1)*nc+1:local_id*nc))
    else
      yc = xx_p((local_id-1)*nc+1:local_id*nc)
    endif
    y_p(off+1:off+nc) = yc
    y_p(off+nc+1:off+this%ndof) = 0.d0
    do imnrl = 1, this%nmnrl
      y_p(off+2*nc+imnrl) = &
        RTDAEMineralInventory(this,imnrl, &
                              rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl), &
                              this%maux(ghosted_id)%volume)
      ! event functions are inventories relative to that of a volume
      ! fraction of max(initial, 1e-6)
      this%mnrl_scale(imnrl,local_id) = &
        RTDAEMineralInventory(this,imnrl, &
          max(rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl),1.d-6), &
          this%maux(ghosted_id)%volume)
      ! within the event tolerance of zero the mineral is absent: its
      ! inventory goes to solution below, as at a depletion
      Kres(imnrl) = 0.d0
      if (RTDAEMineralAbove(this,imnrl,local_id,y_p(off+2*nc+imnrl))) then
        this%mnrl_regime(imnrl,local_id) = RT_DAE_MNRL_PRESENT
      else
        Kres(imnrl) = max(y_p(off+2*nc+imnrl),0.d0)
        y_p(off+2*nc+imnrl) = 0.d0
        this%mnrl_regime(imnrl,local_id) = RT_DAE_MNRL_ABSENT
      endif
    enddo
    do ikinmr = 1, this%nmr
      do irate = 1, this%reaction%surface_complexation%kinmr_nrate(ikinmr)
        do iq = 1, this%mr_ncomp(ikinmr)
          y_p(off+2*nc+RTDAEMRIndex(this,ikinmr,irate,iq)) = &
            this%maux(ghosted_id)%volume* &
            rt_auxvars(ghosted_id)%kinmr_total_sorb( &
              this%mr_comp(iq,ikinmr),irate,ikinmr)
        enddo
      enddo
    enddo
    if (patch%imat(ghosted_id) <= 0) then
      this%c_frozen(:,local_id) = yc
      cycle
    endif
    call RTDAECellChemistry(this,ghosted_id,yc,PETSC_FALSE,status)
    do imnrl = 1, this%nmnrl
      this%scratch%mnrl_volfrac(imnrl) = &
        RTDAEMineralVolFrac(this,imnrl,y_p(off+2*nc+imnrl), &
                            this%maux(ghosted_id)%volume)
    enddo
    if (status /= RT_DAE_CHEM_OK) then
      write(this%option%io_buffer,'(a,i9,a)') 'GIRT_DAE: the chemistry of &
        &cell ', grid%nG2A(ghosted_id), ' could not be evaluated at the &
        &initial condition.'
      call PrintErrMsgByRank(this%option)
    endif
    call RTDAEAccumulation(this,ghosted_id,A,gross)
    y_p(off+nc+1:off+2*nc) = A
    do imnrl = 1, this%nmnrl
      if (Kres(imnrl) > 0.d0) &
        call RTDAEMineralToSolution(this,imnrl,Kres(imnrl), &
                                    y_p(off+nc+1:off+2*nc))
    enddo
  enddo
  call VecRestoreArrayRead(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(Y,y_p,ierr);CHKERRQ(ierr)
  this%mnrl_regime0 = this%mnrl_regime

end subroutine RTDAEInitialState

! ************************************************************************** !

subroutine RTDAEDefaultTolerances(this,atol,rtol)
  !
  ! atol is a concentration [mol/L pore water]; see RTDAESetTolerances.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol

  atol = 1.d-12
  rtol = 1.d-5

end subroutine RTDAEDefaultTolerances

! ************************************************************************** !

subroutine RTDAESetTolerances(this,atol,rtol,vatol,vrtol)
  !
  ! c entries are algebraic, excluded from the local error estimate unless
  ! -rt_dae_conc_tol. atol is a concentration [mol/L]; an inventory's
  ! absolute tolerance is atol times the cell's initial pore water [L]
  ! (bulk volume [m^3] for immobile components).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: atol
  PetscReal :: rtol
  Vec :: vatol
  Vec :: vrtol

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: a_p(:), r_p(:)
  PetscReal :: liters, volume, sat_por
  PetscInt :: local_id, ghosted_id, off, nc, naq
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  nc = this%ncomp
  naq = this%naqcomp
  this%tol_rtol = rtol
  if (.not.allocated(this%tol_atol)) &
    allocate(this%tol_atol(nc+this%nkin,this%nlmax))
  call VecGetArray(vatol,a_p,ierr);CHKERRQ(ierr)
  call VecGetArray(vrtol,r_p,ierr);CHKERRQ(ierr)
  do local_id = 1, this%nlmax
    ghosted_id = grid%nL2G(local_id)
    off = (local_id-1)*this%ndof
    volume = patch%aux%Material%auxvars(ghosted_id)%volume
    sat_por = max(patch%aux%Material%auxvars(ghosted_id)%porosity* &
                  patch%aux%Global%auxvars(ghosted_id)%sat(1),1.d-3)
    liters = RTDAEPoreLiters(sat_por,volume)
    a_p(off+1:off+nc) = 0.d0
    r_p(off+1:off+nc) = 0.d0
    if (this%conc_tol > 0.d0) then
      if (this%log_formulation) then
        a_p(off+1:off+nc) = this%conc_tol
      else
        r_p(off+1:off+nc) = this%conc_tol
      endif
    endif
    a_p(off+nc+1:off+nc+naq) = atol*liters
    a_p(off+nc+naq+1:off+2*nc) = atol*volume
    r_p(off+nc+1:off+2*nc) = rtol
    if (this%nkin > 0) then
      a_p(off+2*nc+1:off+this%ndof) = atol*liters
      r_p(off+2*nc+1:off+this%ndof) = rtol
    endif
    this%tol_atol(:,local_id) = a_p(off+nc+1:off+this%ndof)
  enddo
  call VecRestoreArray(vatol,a_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(vrtol,r_p,ierr);CHKERRQ(ierr)

end subroutine RTDAESetTolerances

! ************************************************************************** !

subroutine RTDAEConvergenceTest(this,time,X,dX,F,dt,kappa,converged)
  !
  ! Converged when, in inventory units [mol], with the local error
  ! tolerance tol_i = rtol*|N_i| + atol_i:
  !   rate rows   dt*|F_N,i|, dt*|F_K,m|  <= kappa*tol_i
  !   updates     |dN_i|, |dK_m|          <= kappa*tol_i
  !   algebraic   |N_i - A_i(c)| <= kappa*rtol*max(|N_i|,G_i)
  !                                 + ALG_ABS_FACTOR*atol_i
  ! with G_i the gross equilibrium inventory at X (prepared here if SNES
  ! did not just evaluate X). LOG_FORMULATION has no absolute part, which
  ! would accept the negative N a large drop of ln(c) can leave. A
  ! component whose inventory cannot be negative needs
  ! N_i >= -ALG_ABS_FACTOR*atol_i: a stage that drives it negative fails.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
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
  PetscReal :: tol, tol_alg, f_alg
  PetscInt :: local_id, i, off, nc, nk, status
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  nc = this%ncomp
  nk = this%nkin
  call RTDAEPrepare(this,time,X,PETSC_FALSE,status)
  ! an iterate outside the domain has converged to nothing
  converged = PETSC_FALSE
  if (status /= TS_MODEL_OK) return
  converged = PETSC_TRUE
  call VecGetArrayRead(X,x_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(dX,dx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(F,f_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    off = (local_id-1)*this%ndof
    do i = 1, nc
      tol = kappa*(this%tol_rtol*abs(x_p(off+nc+i)) + this%tol_atol(i,local_id))
      tol_alg = kappa*this%tol_rtol* &
                max(abs(x_p(off+nc+i)),this%acc_gross(i,local_id))
      f_alg = abs(f_p(off+i))
      if (this%log_formulation) then
        if (x_p(off+i) <= LOG_MIN) f_alg = f_p(off+i)
      else
        tol_alg = tol_alg + ALG_ABS_FACTOR*this%tol_atol(i,local_id)
      endif
      if (.not.(f_alg <= tol_alg .and. &
                dt*abs(f_p(off+nc+i)) <= tol .and. &
                abs(dx_p(off+nc+i)) <= tol) .or. &
          (this%positive_comp(i) .and. &
           x_p(off+nc+i) < -ALG_ABS_FACTOR*this%tol_atol(i,local_id))) then
        if (this%debug) then
          write(*,'(a,2i4,*(es11.3))') ' GIRT_DAE not converged: cell,comp,&
            &F_c,tol_alg,dt*F_N,dN,tol,c,dc: ',local_id,i,f_p(off+i), &
            tol_alg,dt*f_p(off+nc+i),dx_p(off+nc+i),tol,x_p(off+i), &
            dx_p(off+i)
        endif
        converged = PETSC_FALSE
        exit
      endif
    enddo
    if (.not.converged) exit
    do i = 1, nk
      tol = kappa*(this%tol_rtol*abs(x_p(off+2*nc+i)) + &
                   this%tol_atol(nc+i,local_id))
      if (.not.(dt*abs(f_p(off+2*nc+i)) <= tol .and. &
                abs(dx_p(off+2*nc+i)) <= tol)) then
        converged = PETSC_FALSE
        exit
      endif
    enddo
    if (.not.converged) exit
  enddo
  call VecRestoreArrayRead(F,f_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(dX,dx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X,x_p,ierr);CHKERRQ(ierr)

end subroutine RTDAEConvergenceTest

! ************************************************************************** !

subroutine RTDAEPreCheck(this,Y,dY,changed)
  !
  ! Newton update limiter (SNES applies Y - dY). Without LOG_FORMULATION an
  ! update that would take a concentration below PRECHECK_FLOOR of its
  ! current value is truncated to land there, other entries take the full
  ! step; with it, a change of ln(c) is capped at LOG_MAX_UPDATE and ln(c)
  ! floored at LOG_MIN. A solver strategy only: the convergence test
  ! decides.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  Vec :: Y
  Vec :: dY
  PetscBool :: changed

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: y_p(:), dy_p(:)
  PetscInt :: local_id, i, idx
  PetscBool :: truncate
  PetscErrorCode :: ierr

  changed = PETSC_FALSE
  patch => this%realization%patch
  grid => patch%grid

  if (this%log_formulation) then
    ! no positivity to protect, but exp(ln(c)) must stay a normal number;
    ! a vanishing species' Newton step in ln(c) can be arbitrarily large
    call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
    call VecGetArray(dY,dy_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
      do i = 1, this%ncomp
        idx = (local_id-1)*this%ndof + i
        if (abs(dy_p(idx)) > LOG_MAX_UPDATE) then
          dy_p(idx) = sign(LOG_MAX_UPDATE,dy_p(idx))
          changed = PETSC_TRUE
        endif
        if (dy_p(idx) > 0.d0 .and. y_p(idx)-dy_p(idx) < LOG_MIN) then
          dy_p(idx) = max(y_p(idx)-LOG_MIN,0.d0)
          changed = PETSC_TRUE
        endif
      enddo
    enddo
    call VecRestoreArray(dY,dy_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
    return
  endif

  truncate = PETSC_FALSE
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(dY,dy_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    do i = 1, this%ncomp
      idx = (local_id-1)*this%ndof + i
      if (dy_p(idx) > (1.d0-PRECHECK_FLOOR)*y_p(idx)) truncate = PETSC_TRUE
    enddo
  enddo
  call VecRestoreArrayRead(dY,dy_p,ierr);CHKERRQ(ierr)
  if (truncate) then
    call VecGetArray(dY,dy_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
      do i = 1, this%ncomp
        idx = (local_id-1)*this%ndof + i
        if (dy_p(idx) > (1.d0-PRECHECK_FLOOR)*y_p(idx)) &
          dy_p(idx) = (1.d0-PRECHECK_FLOOR)*y_p(idx)
      enddo
    enddo
    call VecRestoreArray(dY,dy_p,ierr);CHKERRQ(ierr)
  endif
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  changed = truncate

end subroutine RTDAEPreCheck

! ************************************************************************** !

subroutine RTDAEProjectInitialGuess(this,X,Xstart)
  !
  ! A multistep extrapolation of a vanishing species can leave the
  ! closure's domain (a negative concentration, or a logarithm whose
  ! exponential underflows); such an entry starts from its step-start value.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  Vec :: X
  Vec :: Xstart

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: x_p(:), x0_p(:)
  PetscInt :: local_id, i, idx
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  call VecGetArray(X,x_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Xstart,x0_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    do i = 1, this%ncomp
      idx = (local_id-1)*this%ndof + i
      if (this%log_formulation) then
        if (.not.(x_p(idx) >= LOG_MIN .and. x_p(idx) < 700.d0)) &
          x_p(idx) = x0_p(idx)
      else
        if (.not.(x_p(idx) >= 0.d0)) x_p(idx) = x0_p(idx)
      endif
    enddo
  enddo
  call VecRestoreArrayRead(Xstart,x0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(X,x_p,ierr);CHKERRQ(ierr)

end subroutine RTDAEProjectInitialGuess

! ************************************************************************** !

subroutine RTDAEPublish(this,time,Y)
  !
  ! Derives the legacy output fields from an accepted state: tran_xx, the
  ! realization's RT auxiliary variables, the mineral volume fractions,
  ! areas and rates and, uncoupled, the patch's face fluxes. One way only:
  ! evaluation never reads these.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal, pointer :: yl_p(:), xx_p(:)
  PetscReal :: c(this%ncomp)
  PetscInt :: local_id, ghosted_id, off, nc, status, imnrl
  PetscInt :: ikinmr, irate, iq
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  field => this%realization%field
  rt_auxvars => patch%aux%RT%auxvars
  nc = this%ncomp

  call RTDAEPrepare(this,time,Y,PETSC_FALSE,status)
  ! uncoupled, the face fluxes at time, for output: BeginSegment leaves the
  ! segment end's in the patch for a LINEAR velocity dataset (coupled, the
  ! flow model publishes them)
  if (.not.associated(this%fs)) then
    call RTDAEForcingAt(this,time)
    if (this%conn%nconn_int > 0) &
      patch%internal_velocities(1:this%nphase,1:this%conn%nconn_int) = &
        this%forcing%vel_int(:,1:this%conn%nconn_int)
    if (this%conn%nconn_bc > 0) &
      patch%boundary_velocities(1:this%nphase,1:this%conn%nconn_bc) = &
        this%forcing%vel_bc(:,1:this%conn%nconn_bc)
  endif
  call VecGetArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, this%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle
    off = (ghosted_id-1)*this%ndof
    if (this%log_formulation) then
      c = exp(yl_p(off+1:off+nc))
    else
      c = yl_p(off+1:off+nc)
    endif
    local_id = grid%nG2L(ghosted_id)
    if (local_id > 0) xx_p((local_id-1)*nc+1:local_id*nc) = c
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    call RTDAECellChemistry(this,ghosted_id,yl_p(off+1:off+nc),PETSC_FALSE, &
                            status)
    rt_auxvars(ghosted_id)%pri_molal = this%scratch%pri_molal
    rt_auxvars(ghosted_id)%pri_act_coef = this%scratch%pri_act_coef
    rt_auxvars(ghosted_id)%total = this%scratch%total
    rt_auxvars(ghosted_id)%aqueous%dtotal = this%scratch%aqueous%dtotal
    if (this%reaction%neqcplx > 0) then
      rt_auxvars(ghosted_id)%sec_molal = this%scratch%sec_molal
      rt_auxvars(ghosted_id)%sec_act_coef = this%scratch%sec_act_coef
    endif
    if (this%nimcomp > 0) &
      rt_auxvars(ghosted_id)%immobile = this%scratch%immobile
    ! free site and surface complex concentrations (output only)
    if (this%reaction%surface_complexation%nsrfcplxrxn > 0) then
      do imnrl = 1, this%nmnrl
        this%scratch%mnrl_volfrac(imnrl) = &
          RTDAEMineralVolFrac(this,imnrl,yl_p(off+2*nc+imnrl), &
                              this%maux(ghosted_id)%volume)
      enddo
      call RTDAEPublishSurfaceComplexes(this,ghosted_id,rt_auxvars(ghosted_id))
    endif
    do ikinmr = 1, this%nmr
      do irate = 1, this%reaction%surface_complexation%kinmr_nrate(ikinmr)
        do iq = 1, this%mr_ncomp(ikinmr)
          rt_auxvars(ghosted_id)%kinmr_total_sorb( &
            this%mr_comp(iq,ikinmr),irate,ikinmr) = &
            yl_p(off+2*nc+RTDAEMRIndex(this,ikinmr,irate,iq))/ &
            this%maux(ghosted_id)%volume
        enddo
      enddo
    enddo
    if (this%eq_srfcplx .and. local_id > 0) then
      ! sorbed totals V*S of the accumulation, per bulk volume
      rt_auxvars(ghosted_id)%total_sorb_eq = &
        (this%acc(1:this%naqcomp,local_id) - &
         RTDAEPoreLiters(this%maux(ghosted_id)%porosity* &
                         this%forcing%gaux(ghosted_id)%sat(1), &
                         this%maux(ghosted_id)%volume)* &
         this%scratch%total(:,1))/this%maux(ghosted_id)%volume
    endif
    do imnrl = 1, this%nmnrl
      rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl) = &
        RTDAEMineralVolFrac(this,imnrl,yl_p(off+2*nc+imnrl), &
                            this%maux(ghosted_id)%volume)
      if (local_id > 0) then
        rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
          this%mnrl_area(imnrl,local_id)
        rt_auxvars(ghosted_id)%mnrl_rate(imnrl) = &
          this%mnrl_rate(imnrl,local_id)
      endif
    enddo
  enddo
  call VecRestoreArray(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   field%tran_xx,field%tran_xx_loc,NTRANDOF)
  if (this%log_formulation) then
    call VecCopy(field%tran_xx,field%tran_log_xx,ierr);CHKERRQ(ierr)
    call VecLog(field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif

end subroutine RTDAEPublish

! ************************************************************************** !

function RTDAENumEvents(this)
  !
  ! One event function per kinetic mineral and local cell.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: RTDAENumEvents

  RTDAENumEvents = this%nmnrl*this%nlmax

end function RTDAENumEvents

! ************************************************************************** !

subroutine RTDAEEventSetup(this,direction,vtol)
  !
  ! Present mineral: K/mnrl_scale falls through zero (depletion). Absent
  ! mineral: ln(Q/K) falls through zero as the solution turns
  ! undersaturated (PostEvent decides whether dissolution resumes). Both
  ! decrease.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: direction(:)
  PetscReal :: vtol(:)

  direction = -1
  vtol = RT_DAE_EVENT_VTOL

end subroutine RTDAEEventSetup

! ************************************************************************** !

subroutine RTDAEEventIndicator(this,time,Y,fvalue)
  !
  ! K/mnrl_scale of present minerals, ln(Q/K) of absent ones.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscReal :: fvalue(:)

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: y_p(:)
  PetscInt :: local_id, imnrl, ievent, status
  PetscBool :: absent
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  fvalue = 1.d0
  if (this%nmnrl == 0) return
  ! ln(Q/K) of absent minerals needs the chemistry at Y (after an accepted
  ! step, reused from the last Newton iterate). Prepare is collective:
  ! every rank takes the same branch.
  absent = any(this%mnrl_regime == RT_DAE_MNRL_ABSENT)
  call MPI_Allreduce(MPI_IN_PLACE,absent,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,this%option%mycomm,ierr);CHKERRQ(ierr)
  if (absent) then
    call RTDAEPrepare(this,time,Y,PETSC_FALSE,status)
  endif
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    do imnrl = 1, this%nmnrl
      ievent = (local_id-1)*this%nmnrl + imnrl
      if (this%mnrl_regime(imnrl,local_id) == RT_DAE_MNRL_PRESENT) then
        fvalue(ievent) = y_p((local_id-1)*this%ndof+2*this%ncomp+imnrl)/ &
                         this%mnrl_scale(imnrl,local_id)
      else if (this%cell_ok(grid%nL2G(local_id))) then
        fvalue(ievent) = this%mnrl_lnQK(imnrl,local_id)
      endif
    enddo
  enddo
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)

end subroutine RTDAEEventIndicator

! ************************************************************************** !

subroutine RTDAEPostEvent(this,time,Y,events,changed)
  !
  ! Regime transitions at located events.
  !   depletion (present -> absent): the remaining K, zero to the event
  !     tolerance, dissolves: N_i += nu_i K, K := 0 (N_i + nu_i K kept);
  !   reappearance (absent -> present): precipitated mineral exists and
  !     the solution turns undersaturated; dissolution resumes, state
  !     unchanged. Within the event tolerance of zero the mineral instead
  !     dissolves at once, as at a depletion, and stays absent.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscInt :: events(:)
  PetscBool :: changed

  PetscReal, pointer :: y_p(:)
  PetscInt :: ie, ievent, local_id, imnrl, offK, offN
  PetscBool :: write_y
  PetscErrorCode :: ierr

  changed = PETSC_FALSE
  ! the regimes enter the reaction terms: the prepared chemistry expires
  ! on every rank (the driver calls this collectively)
  this%prepared = PETSC_FALSE
  if (size(events) == 0) return

  ! write access only to change the state: PETSc restarts the integrator
  ! whenever the state of Y changes in the post-event handler
  write_y = PETSC_FALSE
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  do ie = 1, size(events)
    ievent = events(ie) + 1
    local_id = (ievent-1)/this%nmnrl + 1
    imnrl = ievent - (local_id-1)*this%nmnrl
    offK = (local_id-1)*this%ndof + 2*this%ncomp
    if (this%mnrl_regime(imnrl,local_id) == RT_DAE_MNRL_PRESENT) then
      write_y = PETSC_TRUE
    else if (y_p(offK+imnrl) > 0.d0) then
      if (.not.RTDAEMineralAbove(this,imnrl,local_id,y_p(offK+imnrl))) &
        write_y = PETSC_TRUE
    endif
  enddo
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  if (write_y) then
    call VecGetArray(Y,y_p,ierr);CHKERRQ(ierr)
  else
    call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  endif
  do ie = 1, size(events)
    ! PETSc numbers events from zero
    ievent = events(ie) + 1
    local_id = (ievent-1)/this%nmnrl + 1
    imnrl = ievent - (local_id-1)*this%nmnrl
    offN = (local_id-1)*this%ndof + this%ncomp
    offK = (local_id-1)*this%ndof + 2*this%ncomp
    if (this%mnrl_regime(imnrl,local_id) == RT_DAE_MNRL_PRESENT .or. &
        (y_p(offK+imnrl) > 0.d0 .and. &
         .not.RTDAEMineralAbove(this,imnrl,local_id,y_p(offK+imnrl)))) then
      this%depletion_residual = this%depletion_residual + &
                                abs(y_p(offK+imnrl))
      call RTDAEMineralToSolution(this,imnrl,y_p(offK+imnrl), &
                                  y_p(offN+1:offN+this%ncomp))
      y_p(offK+imnrl) = 0.d0
      this%mnrl_regime(imnrl,local_id) = RT_DAE_MNRL_ABSENT
      this%num_depletions = this%num_depletions + 1
      changed = PETSC_TRUE
    else if (y_p(offK+imnrl) > 0.d0) then
      this%mnrl_regime(imnrl,local_id) = RT_DAE_MNRL_PRESENT
      this%num_reappearances = this%num_reappearances + 1
      changed = PETSC_TRUE
    endif
  enddo
  if (write_y) then
    call VecRestoreArray(Y,y_p,ierr);CHKERRQ(ierr)
  else
    call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  endif

end subroutine RTDAEPostEvent

! ************************************************************************** !

subroutine RTDAEPrintSummary(this)
  !
  ! Mineral regime transitions of the run, summed over ranks.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this

  PetscInt :: counts(2)
  PetscReal :: residual
  PetscErrorCode :: ierr

  if (this%nmnrl == 0) return
  counts(1) = this%num_depletions
  counts(2) = this%num_reappearances
  call MPI_Allreduce(MPI_IN_PLACE,counts,TWO_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,this%option%mycomm,ierr);CHKERRQ(ierr)
  residual = this%depletion_residual
  call MPI_Allreduce(MPI_IN_PLACE,residual,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,this%option%mycomm, &
                     ierr);CHKERRQ(ierr)
  write(this%option%io_buffer,'(a,i8,a,i8,a,es10.3,a)') &
    ' GIRT_DAE mineral regime transitions: ',counts(1),' depletions, ', &
    counts(2),' reappearances; ',residual, &
    ' mol moved to solution at depletion events'
  call PrintMsg(this%option)

end subroutine RTDAEPrintSummary

! ************************************************************************** !

subroutine RTDAEForcingTimes(this,final_time,times)
  !
  ! Knots up to final_time of the uniform velocity dataset, every cycle of
  ! a CYCLIC one (PFLOTRAN makes waypoints of the first cycle only).
  !
  ! Author: Hansol Suh
  ! Date: 10/06/26
  !
  use Time_Storage_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: final_time
  PetscReal, allocatable :: times(:)

  PetscReal, pointer :: knots(:)

  allocate(times(0))
  if (.not.associated(this%realization%uniform_velocity_dataset)) return
  if (.not.associated(this%realization%uniform_velocity_dataset% &
                      time_storage)) return
  call TimeStorageGetTimes(this%realization%uniform_velocity_dataset% &
                             time_storage,this%option,final_time,knots)
  if (.not.associated(knots)) return
  times = knots
  deallocate(knots)

end subroutine RTDAEForcingTimes

! ************************************************************************** !

subroutine RTDAEPrepareReactions(this,time,Y,status)
  !
  ! The forcing at time and ghosted state Y for RTDAELoadCell, and every
  ! local cell's mnrl_area and mnrl_lnQK. Collective.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  PetscInt :: status

  call RTDAEPrepare(this,time,Y,PETSC_FALSE,status)

end subroutine RTDAEPrepareReactions

! ************************************************************************** !

subroutine RTDAELoadCell(this,local_id)
  !
  ! Loads into scratch the chemistry of active local cell local_id at the
  ! last RTDAEPrepareReactions (see RTDAELoadCellRates) and its mineral
  ! volume fractions, for the rate kernels with the cell's gaux and maux.
  ! Requires RTDAEKeepActivities.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: local_id

  PetscReal, pointer :: yl_p(:)
  PetscInt :: ghosted_id, offset, imnrl
  PetscErrorCode :: ierr

  if (.not.allocated(this%pri_act_coef)) then
    call PrintErrMsg(this%option,'RTDAELoadCell requires &
                     &RTDAEKeepActivities.')
  endif
  ghosted_id = this%realization%patch%grid%nL2G(local_id)
  offset = (ghosted_id-1)*this%ndof
  call VecGetArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)
  call RTDAELoadCellRates(this,ghosted_id,yl_p)
  do imnrl = 1, this%nmnrl
    this%scratch%mnrl_volfrac(imnrl) = &
      RTDAEMineralVolFrac(this,imnrl,yl_p(offset+2*this%ncomp+imnrl), &
                          this%maux(ghosted_id)%volume)
  enddo
  call VecRestoreArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)

end subroutine RTDAELoadCell

! ************************************************************************** !

subroutine RTDAELoadCellRates(this,ghosted_id,yl_p)
  !
  ! Loads into scratch what the microbial and general rate kernels read of
  ! active cell ghosted_id at the last Prepare: c from yl_p and the kept
  ! activity coefficients and ln water activity (RTDAEKeepActivities).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: yl_p(:)

  PetscReal :: c(this%ncomp)
  PetscInt :: offset

  offset = (ghosted_id-1)*this%ndof
  c = yl_p(offset+1:offset+this%ncomp)
  if (this%log_formulation) c = exp(c)
  this%scratch%pri_molal = c(1:this%naqcomp)
  if (this%nimcomp > 0) &
    this%scratch%immobile = c(this%naqcomp+1:this%ncomp)
  this%scratch%pri_act_coef = this%pri_act_coef(:,ghosted_id)
  ! water activity 1 (ACTIVITY_H2O is not supported)
  this%scratch%ln_act_h2o = 0.d0

end subroutine RTDAELoadCellRates

! ************************************************************************** !

subroutine RTDAEInitialSensitivity(this,time,Y,lambda,dJ_dc,dJ_dK)
  !
  ! Pulls lambda = dJ/dY back through Y0 = [c0; A(c0,K0); K0]
  ! (RTDAEInitialState) to each local cell's c0 (physical units, also under
  ! LOG_FORMULATION) and K0 [mol]:
  !   dJ/dc0 = lambda_c + (dA/dc)^T lambda_N
  !   dJ/dK0 = lambda_K + (dA/dK)^T lambda_N
  ! c0 and K0 are independent: dJ/dc0 holds K0 fixed, also where the input
  ! initializes K0 from c0 (multirate sorption at equilibrium). Zero in
  ! inactive cells. Collective.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  Vec :: Y
  Vec :: lambda
  PetscReal :: dJ_dc(:,:)
  PetscReal :: dJ_dK(:,:)

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: lam_p(:), y_p(:)
  PetscInt :: local_id, ghosted_id, offset, nc, nk, status
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  nc = this%ncomp
  nk = this%nkin
  dJ_dc = 0.d0
  dJ_dK = 0.d0
  call RTDAEPrepare(this,time,Y,PETSC_TRUE,status)
  if (status /= TS_MODEL_OK) then
    call PrintErrMsgByRank(this%option,'GIRT_DAE: the chemistry of the &
      &initial state could not be evaluated for the adjoint.')
  endif
  call VecGetArrayRead(lambda,lam_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*this%ndof
    dJ_dc(:,local_id) = lam_p(offset+1:offset+nc) + &
      matmul(lam_p(offset+nc+1:offset+2*nc),this%dacc_dc(:,:,local_id))
    ! under LOG_FORMULATION the c coordinates are ln(c)
    if (this%log_formulation) &
      dJ_dc(:,local_id) = dJ_dc(:,local_id)/exp(y_p(offset+1:offset+nc))
    if (nk > 0) dJ_dK(1:nk,local_id) = lam_p(offset+2*nc+1:offset+2*nc+nk) + &
      matmul(lam_p(offset+nc+1:offset+2*nc),this%dacc_dk(:,1:nk,local_id))
  enddo
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(lambda,lam_p,ierr);CHKERRQ(ierr)

end subroutine RTDAEInitialSensitivity

! ************************************************************************** !

subroutine RTDAECellObservables(this,time,local_id,y_c,total,dtotal,ph, &
                                dph,status)
  !
  ! At active local cell local_id with c coordinates y_c, under the forcing
  ! at time (coupled: of the flow state set): each aqueous component's
  ! total in solution [mol/kg] and, with H+ primary, the pH, -log10 of its
  ! activity (else 0), and their derivatives with respect to y_c (log
  ! chain applied), activity-coefficient sensitivity included. status:
  ! TS_MODEL_OK, or TS_MODEL_DOMAIN_ERROR where the closure fails.
  !
  ! Author: Hansol Suh
  ! Date: 10/16/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  PetscInt :: local_id
  PetscReal :: y_c(:)
  PetscReal :: total(:)
  PetscReal :: dtotal(:,:)
  PetscReal :: ph
  PetscReal :: dph(:)
  PetscInt :: status

  PetscReal :: den_kg_per_L, m
  PetscReal :: J(1,this%naqcomp)
  PetscInt :: ghosted_id, ih, naq

  naq = this%naqcomp
  ghosted_id = this%realization%patch%grid%nL2G(local_id)
  call RTDAEForcingAt(this,time)
  call RTDAECellChemistry(this,ghosted_id,y_c,PETSC_TRUE,status)
  if (status /= RT_DAE_CHEM_OK) then
    status = TS_MODEL_DOMAIN_ERROR
    return
  endif
  status = TS_MODEL_OK
  ! the closure's totals are per liter at this density
  den_kg_per_L = this%forcing%gaux(ghosted_id)%den_kg(1)*1.d-3
  total = this%tot(:,ghosted_id)/den_kg_per_L
  dtotal = 0.d0
  dtotal(:,1:naq) = this%dtot(:,:,ghosted_id)/den_kg_per_L
  ph = 0.d0
  dph = 0.d0
  ih = this%reaction%species_idx%h_ion_id
  if (ih <= 0) return
  m = this%scratch%pri_molal(ih)
  ph = -log10(this%scratch%pri_act_coef(ih)*m)
  ! at fixed activity coefficients, then their sensitivity
  J = 0.d0
  J(1,ih) = -1.d0/(log(10.d0)*m)
  call RTDAEGammaChain(J,1,naq,this%scratch%pri_molal,this%gam)
  dph(1:naq) = J(1,:)
  if (this%log_formulation) &
    dph(1:naq) = dph(1:naq)*this%scratch%pri_molal(1:naq)

end subroutine RTDAECellObservables

! ************************************************************************** !

subroutine RTDAEInitialCellState(this,time,local_id,y_ref,c,y,status)
  !
  ! Block y of active local cell local_id of an initial state, from block
  ! y_ref of InitialState's with the physical c coordinates c: K0 stays,
  ! N0 = N0_ref + A(c,K0) - A(c_ref,K0) with A as in InitialState, so
  ! N0 = N0_ref where c = c_ref (also with the inventory of a mineral
  ! absent at the start), under the forcing at time (coupled: of the flow
  ! state set). status: TS_MODEL_OK, or TS_MODEL_DOMAIN_ERROR where the
  ! closure fails (y is then y_ref with c).
  !
  ! Author: Hansol Suh
  ! Date: 10/16/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: time
  PetscInt :: local_id
  PetscReal :: y_ref(:)
  PetscReal :: c(:)
  PetscReal :: y(:)
  PetscInt :: status

  PetscReal :: A_ref(this%ncomp), A(this%ncomp), gross(this%ncomp)
  PetscReal :: yc(this%ncomp)
  PetscInt :: ghosted_id, nc

  nc = this%ncomp
  ghosted_id = this%realization%patch%grid%nL2G(local_id)
  if (this%log_formulation) then
    yc = log(c)
  else
    yc = c
  endif
  y = y_ref
  y(1:nc) = yc
  call RTDAEForcingAt(this,time)
  call CellInventory(y_ref(1:nc),A_ref)
  if (status == RT_DAE_CHEM_OK) call CellInventory(yc,A)
  if (status /= RT_DAE_CHEM_OK) then
    status = TS_MODEL_DOMAIN_ERROR
    return
  endif
  y(nc+1:2*nc) = y_ref(nc+1:2*nc) + (A - A_ref)
  status = TS_MODEL_OK

contains

subroutine CellInventory(y_c,inventory)

  PetscReal :: y_c(:)
  PetscReal :: inventory(:)

  PetscInt :: imnrl

  call RTDAECellChemistry(this,ghosted_id,y_c,PETSC_FALSE,status)
  if (status /= RT_DAE_CHEM_OK) return
  do imnrl = 1, this%nmnrl
    this%scratch%mnrl_volfrac(imnrl) = &
      RTDAEMineralVolFrac(this,imnrl,y_ref(2*nc+imnrl), &
                          this%maux(ghosted_id)%volume)
  enddo
  call RTDAEAccumulation(this,ghosted_id,inventory,gross)

end subroutine CellInventory

end subroutine RTDAEInitialCellState

! ************************************************************************** !

subroutine RTDAEPublishSurfaceComplexes(this,ghosted_id,rt_auxvar)
  !
  ! Free site and equilibrium surface complex concentrations at the
  ! chemistry in scratch, for output.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id
  type(reactive_transport_auxvar_type) :: rt_auxvar

  type(surface_complexation_type), pointer :: srfcplx
  PetscReal, pointer :: srfcplx_conc(:)
  PetscReal :: S(this%naqcomp)
  PetscInt :: irxn

  srfcplx => this%reaction%surface_complexation
  nullify(srfcplx_conc)
  if (associated(rt_auxvar%eqsrfcplx_conc)) then
    rt_auxvar%eqsrfcplx_conc = 0.d0
    srfcplx_conc => rt_auxvar%eqsrfcplx_conc
  endif
  do irxn = 1, srfcplx%nsrfcplxrxn
    S = 0.d0
    this%sorb_work = 0.d0
    call RTDAESorbedEq1(this,ghosted_id,irxn, &
                        rt_auxvar%srfcplxrxn_free_site_conc(irxn), &
                        srfcplx_conc,S,this%sorb_work)
  enddo

end subroutine RTDAEPublishSurfaceComplexes

! ************************************************************************** !

subroutine RTDAEDestroy(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  use Reactive_Transport_Aux_module, only : RTAuxVarStrip

  implicit none

  class(rt_dae_model_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  call RTDAEForcingDestroy(this%forcing)
  this%recorded%n = 0
  if (associated(this%tran_coefs_int)) deallocate(this%tran_coefs_int)
  if (associated(this%tran_coefs_bc)) deallocate(this%tran_coefs_bc)
  if (associated(this%maux)) then
    do i = 1, size(this%maux)
      call MaterialAuxVarStrip(this%maux(i))
    enddo
    deallocate(this%maux)
  endif
  if (associated(this%scratch)) then
    call RTAuxVarStrip(this%scratch)
    deallocate(this%scratch)
  endif
  call RTDAEGammaDestroy(this%gam)
  if (.not.PetscObjectIsNull(this%Y_loc)) then
    call VecDestroy(this%Y_loc,ierr);CHKERRQ(ierr)
  endif
  if (.not.PetscObjectIsNull(this%Y_prepared)) then
    call VecDestroy(this%Y_prepared,ierr);CHKERRQ(ierr)
  endif

end subroutine RTDAEDestroy

end module RT_DAE_Model_class
