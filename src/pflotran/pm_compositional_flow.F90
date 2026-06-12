module PM_Compositional_Flow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module
  use Compositional_Flow_Aux_module

  implicit none

  private

  PetscInt, parameter :: MAX_CHANGE_LIQ_PRES_NI = 1
  PetscInt, parameter :: MAX_CHANGE_CONC_NI = 2
  PetscInt, parameter :: MAX_RES_LIQ_EQ = 3
  PetscInt, parameter :: MAX_RES_SOL_EQ = 4

  type, public, extends(pm_subsurface_flow_type) :: pm_compflow_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: max_allow_liq_pres_change_ni
    PetscReal :: liq_pres_change_ts_governor
    PetscReal :: liq_sat_change_ts_governor
    PetscInt :: convergence_flags(MAX_RES_SOL_EQ)
    PetscReal :: convergence_reals(MAX_RES_SOL_EQ)
    PetscReal :: sat_update_trunc_ni
    PetscReal :: unsat_to_sat_pres_damping_ni
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMCompFlowReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMCompFlowReadTSSelectCase
    procedure, public :: ReadNewtonBlock => PMCompFlowReadNewtonSelectCase
    procedure, public :: Setup => PMCompFlowSetup
    procedure, public :: InitializeRun => PMCompFlowInitializeRun
    procedure, public :: InitializeTimestep => PMCompFlowInitializeTimestep
    procedure, public :: Residual => PMCompFlowResidual
    procedure, public :: Jacobian => PMCompFlowJacobian
    procedure, public :: UpdateTimestep => PMCompFlowUpdateTimestep
    procedure, public :: FinalizeTimestep => PMCompFlowFinalizeTimestep
    procedure, public :: PreSolve => PMCompFlowPreSolve
    procedure, public :: PostSolve => PMCompFlowPostSolve
    procedure, public :: CheckUpdatePre => PMCompFlowCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMCompFlowCheckUpdatePost
    procedure, public :: CheckConvergence => PMCompFlowCheckConvergence
    procedure, public :: TimeCut => PMCompFlowTimeCut
    procedure, public :: UpdateSolution => PMCompFlowUpdateSolution
    procedure, public :: UpdateAuxVars => PMCompFlowUpdateAuxVars
    procedure, public :: MaxChange => PMCompFlowMaxChange
    procedure, public :: InputRecord => PMCompFlowInputRecord
    procedure, public :: CheckpointBinary => PMCompFlowCheckpointBinary
    procedure, public :: RestartBinary => PMCompFlowRestartBinary
    procedure, public :: Destroy => PMCompFlowDestroy
  end type pm_compflow_type

  public :: PMCompFlowCreate, &
            PMCompFlowInitObject, &
            PMCompFlowInitializeRun, &
            PMCompFlowFinalizeTimestep, &
            PMCompFlowCheckUpdatePre, &
            PMCompFlowDestroy

contains

! ************************************************************************** !

function PMCompFlowCreate()
  !
  ! Creates CompFlow process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  implicit none

  class(pm_compflow_type), pointer :: PMCompFlowCreate

  class(pm_compflow_type), pointer :: compflow_pm

  allocate(compflow_pm)
  call PMCompFlowInitObject(compflow_pm)

  PMCompFlowCreate => compflow_pm

end function PMCompFlowCreate

! ************************************************************************** !

subroutine PMCompFlowInitObject(this)
  !
  ! Creates CompFlow process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Option_module
  use String_module
  use Variables_module

  implicit none

  class(pm_compflow_type) :: this

  call PMSubsurfaceFlowInit(this)
  this%name = 'Z Flow'
  this%header = 'Z FLOW'

  nullify(this%max_change_ivar)

  ! set to UNINITIALIZED_DOUBLE and report error below is set from input
  this%pressure_change_governor = UNINITIALIZED_DOUBLE
  this%temperature_change_governor = UNINITIALIZED_DOUBLE
  this%saturation_change_governor = UNINITIALIZED_DOUBLE
  this%xmol_change_governor = UNINITIALIZED_DOUBLE

  this%check_post_convergence = PETSC_TRUE

  this%max_allow_liq_pres_change_ni = UNINITIALIZED_DOUBLE
  this%liq_pres_change_ts_governor = 5.d5    ! [Pa]
  this%liq_sat_change_ts_governor = 1.d0
  this%sat_update_trunc_ni = UNINITIALIZED_DOUBLE
  this%unsat_to_sat_pres_damping_ni = UNINITIALIZED_DOUBLE

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

  compflow_debug_cell_id = UNINITIALIZED_INTEGER

end subroutine PMCompFlowInitObject

! ************************************************************************** !

subroutine PMCompFlowReadSimOptionsBlock(this,input)
  !
  ! Read ZFLOW options input block
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Compositional_Flow_module
  use Compositional_Flow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module
  use EOS_Water_module

  implicit none

  class(pm_compflow_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH) :: local_error_string
  PetscBool :: found
  PetscReal :: array(1)

  option => this%option

  error_string = 'CompFlow Options'

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('PROCESSES')
        local_error_string = trim(error_string) // ',' // keyword
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword',local_error_string)
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('LIQUID_FLOW')
              compflow_liq_flow_eq = 1     ! these serve as flags until updated
            case('HEAT_TRANSFER')
              compflow_heat_tran_eq = 1
            case('SOLUTE_TRANSPORT')
              compflow_sol_tran_eq = 1
            case default
              call InputKeywordUnrecognized(input,keyword,'CompFlow Mode',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('NO_ACCUMULATION')
        compflow_calc_accum = PETSC_FALSE
      case('NO_FLUX')
        compflow_calc_flux = PETSC_FALSE
      case('NO_BCFLUX')
        compflow_calc_bcflux = PETSC_FALSE
      case('LIQUID_DENSITY')
        call InputReadDouble(input,option,compflow_density_kg)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,compflow_density_kg,'kg/m^3', &
                                      trim(error_string)//','//keyword,option)
      case('LIQUID_VISCOSITY')
        call InputReadDouble(input,option,compflow_viscosity)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,compflow_viscosity,'Pa-s', &
                                      trim(error_string)//','//keyword,option)
      case('DEBUG_CELL_ID')
        call InputReadInt(input,option,compflow_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,'CompFlow Mode',option)
    end select
  enddo
  call InputPopBlock(input,option)

  compflow_density_kmol = compflow_density_kg / FMWH2O
  array(1) = compflow_density_kg ! dist is the aux array
  call EOSWaterSetDensity('CONSTANT',array)

end subroutine PMCompFlowReadSimOptionsBlock

! ************************************************************************** !

subroutine PMCompFlowReadTSSelectCase(this,input,keyword,found, &
                                   error_string,option)
  !
  ! Read timestepper settings specific to the ZFLOW process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_compflow_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  found = PETSC_TRUE
  call PMSubsurfaceFlowReadTSSelectCase(this,input,keyword,found, &
                                        error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('LIQ_SAT_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_sat_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('LIQ_PRES_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_pres_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%liq_pres_change_ts_governor, &
                                    'Pa',keyword,option)
    case default
      found = PETSC_FALSE
  end select

end subroutine PMCompFlowReadTSSelectCase

! ************************************************************************** !

subroutine PMCompFlowReadNewtonSelectCase(this,input,keyword,found, &
                                       error_string,option)
  !
  ! Reads input file parameters associated with the ZFLOW process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Compositional_Flow_Aux_module

  implicit none

  class(pm_compflow_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  error_string = 'ZFLOW Newton Solver'

  select case(trim(keyword))
    case('ITOL_UPDATE')
      option%io_buffer = 'ITOL_UPDATE not supported with ZFLOW. Please &
        &use MAX_ALLOW_LIQ_PRES_CHANGE_NI.'
      call PrintErrMsg(option)
  end select

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('REL_PERTURBATION')
      call InputReadDouble(input,option,compflow_rel_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      ! no units conversion since it is relative
    case('MIN_LIQ_PRESSURE_PERTURBATION')
      call InputReadDouble(input,option,compflow_pres_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      call InputReadAndConvertUnits(input,compflow_pres_min_pert, &
                                    'Pa',keyword,option)
    case('MAX_ALLOW_LIQ_PRES_CHANGE_NI')
      call InputReadDouble(input,option,this%max_allow_liq_pres_change_ni)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%max_allow_liq_pres_change_ni, &
                                    'Pa',keyword,option)
    case('SATURATION_UPDATE_TRUNCATION_NI')
      call InputReadDouble(input,option,this%sat_update_trunc_ni)
      call InputErrorMsg(input,option,keyword,error_string)
    case('UNSAT_TO_SAT_PRESSURE_DAMPING_NI')
      call InputReadDouble(input,option,this%unsat_to_sat_pres_damping_ni)
      call InputErrorMsg(input,option,keyword,error_string)
    case default
      found = PETSC_FALSE

  end select

end subroutine PMCompFlowReadNewtonSelectCase

! ************************************************************************** !

subroutine PMCompFlowSetup(this)
  !
  ! Sets up auxvars and parameters
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Compositional_Flow_module

  implicit none

  class(pm_compflow_type) :: this

  call this%SetRealization()
  call CompFlowSetup(this%realization)
  call PMSubsurfaceFlowSetup(this)

end subroutine PMCompFlowSetup

! ************************************************************************** !

recursive subroutine PMCompFlowInitializeRun(this)
  !
  ! Initializes the ZFLOW mode run.
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Realization_Base_class
  use Patch_module
  use Field_module
  use Material_Aux_module
  use Option_module
  use Variables_module

  implicit none

  class(pm_compflow_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscInt :: ivar

  patch => this%realization%patch
  field => this%realization%field
  option => this%option

  if (this%steady_state) compflow_calc_accum = PETSC_FALSE

  ivar = 0
  if (compflow_liq_flow_eq > 0) then
    ivar = ivar + 2
  endif
  if (compflow_sol_tran_eq > 0) then
    ivar = ivar + 1
  endif
  allocate(this%max_change_ivar(ivar))
  ivar = 1
  if (compflow_liq_flow_eq > 0) then
    this%max_change_ivar(ivar) = LIQUID_PRESSURE
    ivar = ivar + 1
    this%max_change_ivar(ivar) = LIQUID_SATURATION
    ivar = ivar + 1
  endif
  if (compflow_sol_tran_eq > 0) then
    this%max_change_ivar(ivar) = SOLUTE_CONCENTRATION
    ivar = ivar + 1
  endif

  ! need to allocate vectors for max change
  i = size(this%max_change_ivar)
  call VecDuplicateVecs(field%work,i,field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, size(field%max_change_vecs)
    call RealizationGetVariable(this%realization,field%max_change_vecs(i), &
                                this%max_change_ivar(i),ZERO_INTEGER)
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

  if (Initialized(this%temperature_change_governor)) then
    option%io_buffer = 'TEMPERATURE_CHANGE_GOVERNOR &
      &may not be used with ZFLOW.'
    call PrintErrMsg(option)
  endif

  if (soil_compressibility_index == 0 .and. associated(option%inversion)) then
    option%io_buffer = 'Soil compressibility must be employed for CompFlow &
      &when used for inversion.'
    call PrintErrMsg(option)
  endif

end subroutine PMCompFlowInitializeRun

! ************************************************************************** !

subroutine PMCompFlowInitializeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Compositional_Flow_module, only : CompFlowInitializeTimestep

  implicit none

  class(pm_compflow_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call CompFlowInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMCompFlowInitializeTimestep

! ************************************************************************** !

subroutine PMCompFlowFinalizeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  implicit none

  class(pm_compflow_type) :: this

  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMCompFlowFinalizeTimestep

! ************************************************************************** !

subroutine PMCompFlowPreSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  implicit none

  class(pm_compflow_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMCompFlowPreSolve

! ************************************************************************** !

subroutine PMCompFlowPostSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_compflow_type) :: this

end subroutine PMCompFlowPostSolve

! ************************************************************************** !

subroutine PMCompFlowUpdateTimestep(this,update_dt, &
                                 dt,dt_min,dt_max,iacceleration, &
                                 num_newton_iterations,tfac, &
                                 time_step_max_growth_factor)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Option_module
  use Utility_module, only : Equal

  implicit none

  class(pm_compflow_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: sat_ratio, pres_ratio
  PetscReal :: dt_prev

  if (update_dt .and. iacceleration /= 0) then
    dt_prev = dt
    ! calculate the time step ramping factor
    sat_ratio = (2.d0*this%liq_sat_change_ts_governor)/ &
                (this%liq_sat_change_ts_governor+this%max_saturation_change)
    pres_ratio = (2.d0*this%liq_pres_change_ts_governor)/ &
                (this%liq_pres_change_ts_governor+this%max_pressure_change)
    ! pick minimum time step from calc'd ramping factor or maximum ramping factor
    dt = min(min(sat_ratio,pres_ratio)*dt,time_step_max_growth_factor*dt)
    ! make sure time step is within bounds given in the input deck
    dt = min(dt,dt_max)
    if (this%logging_verbosity > 0) then
      if (Equal(dt,dt_max)) then
        string = 'maximum time step size'
      else if (min(sat_ratio,pres_ratio) > time_step_max_growth_factor) then
        string = 'maximum time step growth factor'
      else if (sat_ratio < pres_ratio) then
        string = 'liquid saturation governor'
      else
        string = 'liquid pressure governor'
      endif
      string = 'TS update: ' // trim(string)
      call PrintMsg(this%option,string)
    endif
  endif

  if (Initialized(this%cfl_governor)) then
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMCompFlowUpdateTimestep

! ************************************************************************** !

subroutine PMCompFlowResidual(this,snes,xx,r,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Compositional_Flow_module, only : CompFlowResidual
  use Debug_module
  use Grid_module

  implicit none

  class(pm_compflow_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  Mat :: M

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  ! calculate residual
  if (compflow_simult_function_evals) then
    call SNESGetJacobian(snes,M,PETSC_NULL_MAT,PETSC_NULL_FUNCTION, &
                         PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
    call CompFlowResidual(snes,xx,r,M,this%realization,this%debug,ierr)
  else
    call CompFlowResidual(snes,xx,r,PETSC_NULL_MAT,this%realization, &
                       this%debug,ierr)
  endif

  if (this%debug%vecview_residual) then
    call DebugVecView(this%debug,r,'ZFresidual.out',this%option)
  endif
  if (this%debug%vecview_solution) then
    call DebugVecView(this%debug,xx,'ZFxx.out',this%option)
  endif

end subroutine PMCompFlowResidual

! ************************************************************************** !

subroutine PMCompFlowJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Debug_module
  use Option_module

  implicit none

  class(pm_compflow_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  PetscReal :: norm

  ! the Jacobian was already calculated in PMCompFlowResidual

  if (this%debug%matview_Matrix) then
    call DebugMatView(this%debug,A,'ZFjacobian','', &
                      compflow_ts_count,compflow_ts_cut_count, &
                      compflow_ni_count,this%option)
  endif

  if (this%debug%norm_Matrix) then
    call MatNorm(A,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("1 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("2 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("inf norm: ",es11.4)') norm
    call PrintMsg(this%option)
  endif

  compflow_ni_count = compflow_ni_count + 1

end subroutine PMCompFlowJacobian

! ************************************************************************** !

subroutine PMCompFlowCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Grid_module
  use Material_Aux_module
  use Option_module
  use Patch_module

  implicit none

  class(pm_compflow_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: p_index
  PetscReal :: p_ref, pc, p_target
  PetscReal :: p0, p1, dp
  PetscReal :: sl
  PetscReal :: tempreal
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscBool :: unsat_to_sat_damping_flag
  PetscBool :: sat_update_trunc_flag
  type(compflow_auxvar_type), pointer :: compflow_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option

  compflow_auxvars => patch%aux%CompFlow%auxvars
  material_auxvars => patch%aux%Material%auxvars

  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  changed = PETSC_FALSE

  p_ref = option%flow%reference_pressure
  unsat_to_sat_damping_flag = Initialized(this%unsat_to_sat_pres_damping_ni)
  sat_update_trunc_flag = Initialized(this%sat_update_trunc_ni)

  call VecGetArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(X,X_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    offset = (local_id-1)*option%nflowdof
    if (compflow_liq_flow_eq > 0) then
      p_index = offset+compflow_liq_flow_eq
      dp = -dX_p(p_index)
      p0 = X_p(p_index)
      p1 = p0+dp
      if (unsat_to_sat_damping_flag) then
        sl = compflow_auxvars(ZERO_INTEGER,ghosted_id)%sat + &
            sign(this%sat_update_trunc_ni,dp)
        call patch%characteristic_curves_array( &
              patch%cc_id(ghosted_id))%ptr%saturation_function% &
                CapillaryPressure(sl,pc,tempreal,option)
        if (pc > 0.d0) then
          p_target = p_ref-pc
          if ((dp >= 0.d0 .and. p1 > p_target) .or. &
              (dp < 0.d0 .and. p1 < p_target)) then
            dX_p(p_index) = p0-p_target ! p1 = p0 - dX_p()
            changed = PETSC_TRUE
          endif
        endif
        ! update these incase used below
        dp = -dX_p(p_index)
        p1 = p0+dp
      endif
      if (unsat_to_sat_damping_flag) then
        ! the following initiate damping when transitioning from
        ! unsaturated to saturated state
        if (p0 < p_ref .and. p1 > p_ref) then
          dX_p(p_index) = this%unsat_to_sat_pres_damping_ni*dX_p(p_index)
          changed = PETSC_TRUE
        endif
      endif
    endif
  enddo
  call VecRestoreArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X,X_p,ierr);CHKERRQ(ierr)

end subroutine PMCompFlowCheckUpdatePre

! ************************************************************************** !

subroutine PMCompFlowCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                  X1_changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_module
  use Compositional_Flow_Aux_module

  implicit none

  class(pm_compflow_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscReal :: tempreal

  PetscBool :: converged_liquid_pressure
  PetscReal :: max_abs_pressure_change_NI
  PetscInt :: max_abs_pressure_change_NI_cell
  PetscBool :: converged_concentration
  PetscReal :: max_abs_conc_change_NI
  PetscInt :: max_abs_conc_change_NI_cell

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  ! If these are changed from true, we must add a global reduction on both
  ! variables to ensure that their values match across all processes. Otherwise
  ! PETSc will throw an error in debug mode or ignore the error in optimized.
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecGetArray(X1,X1_p,ierr);CHKERRQ(ierr)
  converged_liquid_pressure = PETSC_TRUE
  converged_concentration = PETSC_FALSE
  max_abs_pressure_change_NI = 0.d0
  max_abs_pressure_change_NI_cell = 0
  max_abs_conc_change_NI = 0.d0
  max_abs_conc_change_NI_cell = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    offset = (local_id-1)*option%nflowdof
    if (compflow_liq_flow_eq > 0) then
      ! maximum absolute change in liquid pressure over Newton iteration
      tempreal = dabs(dX_p(offset+compflow_liq_flow_eq))
      if (tempreal > dabs(max_abs_pressure_change_NI)) then
        max_abs_pressure_change_NI_cell = grid%nG2A(ghosted_id)
        max_abs_pressure_change_NI = tempreal
      endif
    endif
    if (compflow_sol_tran_eq > 0) then
      ! maximum absolute change in liquid pressure over Newton iteration
      tempreal = dabs(dX_p(offset+compflow_sol_tran_eq))
      if (tempreal > dabs(max_abs_conc_change_NI)) then
        max_abs_conc_change_NI_cell = grid%nG2A(ghosted_id)
        max_abs_conc_change_NI = tempreal
      endif
    endif
  enddo

  if (Initialized(this%max_allow_liq_pres_change_ni) .and. &
      max_abs_pressure_change_NI > this%max_allow_liq_pres_change_ni) then
    converged_liquid_pressure = PETSC_FALSE
  endif

  ! the following flags are used in detemining convergence
  if (.not.converged_liquid_pressure) then
    this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = &
      max_abs_pressure_change_NI_cell
  endif

  ! the following flags are for REPORTING purposes only
  this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI) = max_abs_pressure_change_NI
  this%convergence_reals(MAX_CHANGE_CONC_NI) = max_abs_conc_change_NI

  call VecRestoreArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(X1,X1_p,ierr);CHKERRQ(ierr)

end subroutine PMCompFlowCheckUpdatePost

! ************************************************************************** !

subroutine PMCompFlowCheckConvergence(this,snes,it,xnorm,unorm, &
                                   fnorm,reason,ierr)
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_module
  use Compositional_Flow_Aux_module
  use Convergence_module

  implicit none

  class(pm_compflow_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  Vec :: residual_vec
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_tpdt_p(:)
  PetscReal, pointer :: X1_p(:)
  character(len=10) :: reason_string

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(compflow_auxvar_type), pointer :: compflow_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: converged_flag

  PetscReal :: max_abs_res_liq_
  PetscInt :: max_abs_res_liq_cell
  PetscReal :: max_abs_res_sol_
  PetscInt :: max_abs_res_sol_cell
  PetscInt :: offset
  PetscMPIInt :: int_mpi

  PetscReal :: accumulation
  PetscReal :: residual
  PetscReal :: tempreal

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  compflow_auxvars => patch%aux%CompFlow%auxvars
  material_auxvars => patch%aux%Material%auxvars

  residual_vec = field%flow_r
  ! check residual terms
  call VecGetArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_accum_tpdt,accum_tpdt_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (compflow_liq_flow_eq > 0) then
      residual = r_p(offset+compflow_liq_flow_eq)
      accumulation = accum_tpdt_p(offset+compflow_liq_flow_eq)
      ! residual
      tempreal = dabs(residual)
      if (tempreal > max_abs_res_liq_) then
        max_abs_res_liq_ = tempreal
        max_abs_res_liq_cell = grid%nG2A(ghosted_id)
      endif
    endif
    if (compflow_sol_tran_eq > 0) then
      residual = r_p(offset+compflow_sol_tran_eq)
      accumulation = accum_tpdt_p(offset+compflow_sol_tran_eq)
      ! residual
      tempreal = dabs(residual)
      if (tempreal > max_abs_res_sol_) then
        max_abs_res_sol_ = tempreal
        max_abs_res_sol_cell = grid%nG2A(ghosted_id)
      endif
    endif
  enddo

  ! the following flags are used in detemining convergence
  ! currently none

  ! the following flags are for REPORTING purposes only
  this%convergence_flags(MAX_RES_LIQ_EQ) = max_abs_res_liq_cell
  this%convergence_reals(MAX_RES_LIQ_EQ) = max_abs_res_liq_
  this%convergence_flags(MAX_RES_SOL_EQ) = max_abs_res_sol_cell
  this%convergence_reals(MAX_RES_SOL_EQ) = max_abs_res_sol_

  if (this%logging_verbosity >= 10) then
    print *, option%myrank, &
      this%convergence_flags(MAX_RES_LIQ_EQ), &
      this%convergence_reals(MAX_RES_LIQ_EQ), &
      this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI), &
      this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI)
  endif

  int_mpi = size(this%convergence_flags)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_flags,int_mpi,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  int_mpi = size(this%convergence_reals)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_reals,int_mpi, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)

  ! these conditionals cannot change order
  reason_string = '---| '
  converged_flag = CONVERGENCE_CONVERGED
  if (this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) > 0) then
    reason_string(1:1) = 'P'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_CONC_NI) > 0) then
    reason_string(2:2) = 'C'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif

  if (this%logging_verbosity > 0 .and. &
      OptionPrintToScreen(option)) then
    if (option%comm%size > 1) then
      write(*,'(4x,"Rsn: ",a10,2es10.2)') reason_string, &
        this%convergence_reals(MAX_RES_LIQ_EQ), &
        this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI)
    else if (grid%nmax > 9999) then
      write(*,'(4x,"Rsn: ",a10,2(i8,es10.2))') reason_string, &
        this%convergence_flags(MAX_RES_LIQ_EQ), &
        this%convergence_reals(MAX_RES_LIQ_EQ), &
        this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI), &
        this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI)
    else
      write(*,'(4x,"Rsn: ",a10,2(i5,es10.2))') reason_string, &
        this%convergence_flags(MAX_RES_LIQ_EQ), &
        this%convergence_reals(MAX_RES_LIQ_EQ), &
        this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI), &
        this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI)
    endif
  endif

  if (Initialized(this%max_allow_liq_pres_change_ni)) then
    option%convergence = converged_flag
  else
    ! forced standard 2 norms
    option%convergence = CONVERGENCE_OFF
  endif

  call VecRestoreArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%flow_accum_tpdt,accum_tpdt_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMCompFlowCheckConvergence

! ************************************************************************** !

subroutine PMCompFlowTimeCut(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Compositional_Flow_module, only : CompFlowTimeCut

  implicit none

  class(pm_compflow_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call CompFlowTimeCut(this%realization)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMCompFlowTimeCut

! ************************************************************************** !

subroutine PMCompFlowUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Compositional_Flow_module, only : CompFlowUpdateSolution, &
                           CompFlowMapBCAuxVarsToGlobal

  implicit none

  class(pm_compflow_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call CompFlowUpdateSolution(this%realization)
  call CompFlowMapBCAuxVarsToGlobal(this%realization)

end subroutine PMCompFlowUpdateSolution

! ************************************************************************** !

subroutine PMCompFlowUpdateAuxVars(this)
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  use Compositional_Flow_module, only : CompFlowUpdateAuxVars

  implicit none

  class(pm_compflow_type) :: this

  call CompFlowUpdateAuxVars(this%realization)

end subroutine PMCompFlowUpdateAuxVars

! ************************************************************************** !

subroutine PMCompFlowMaxChange(this)
  !
  ! Not needed given CompFlowMaxChange is called in PostSolve
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use String_module
  use Compositional_Flow_Aux_module

  implicit none

  class(pm_compflow_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal, allocatable :: max_change_global(:)
  PetscReal :: max_change, change
  PetscInt :: i, j
  PetscInt :: ivar
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  allocate(max_change_global(size(this%max_change_ivar)))
  max_change_global = 0.d0

  do i = 1, size(this%max_change_ivar)
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i),ZERO_INTEGER)
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArray(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecGetArray(field%max_change_vecs(i),vec_old_ptr, &
                        ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      change = dabs(vec_new_ptr(j)-vec_old_ptr(j))
      max_change = max(max_change,change)
    enddo
    max_change_global(i) = max_change
    call VecRestoreArray(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArray(field%max_change_vecs(i),vec_old_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  i = size(max_change_global)
  call MPI_Allreduce(MPI_IN_PLACE,max_change_global,i,MPI_DOUBLE_PRECISION, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)

  ivar = 1
  if (compflow_liq_flow_eq > 0) then
    write(option%io_buffer,'("  --> max change: dpl= ",1pe12.4, " dsl= ",&
                           &1pe12.4)') &
      max_change_global(ivar:ivar+1)
    this%max_pressure_change = max_change_global(ivar)
    this%max_saturation_change = max_change_global(ivar+1)
    ivar = ivar+2
    call PrintMsg(option)
  endif
  if (compflow_sol_tran_eq > 0) then
    write(option%io_buffer,'(19x,"dc= ",1pe12.4)') max_change_global(ivar)
    ! hijacking xmol_change
    this%max_xmol_change = max_change_global(ivar)
    ivar = ivar+1
    call PrintMsg(option)
  endif
  deallocate(max_change_global)

end subroutine PMCompFlowMaxChange

! ************************************************************************** !

subroutine PMCompFlowInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 06/04/26
  !

  implicit none

  class(pm_compflow_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'brine'

end subroutine PMCompFlowInputRecord

! ************************************************************************** !

subroutine PMCompFlowCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with CompFlow PM
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Checkpoint_module
  use Global_module

  implicit none

  class(pm_compflow_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMCompFlowCheckpointBinary

! ************************************************************************** !

subroutine PMCompFlowRestartBinary(this,viewer)
  !
  ! Restarts data associated with CompFlow PM
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Checkpoint_module
  use Global_module

  implicit none

  class(pm_compflow_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMCompFlowRestartBinary

! ************************************************************************** !

subroutine PMCompFlowDestroy(this)
  !
  ! Destroys CompFlow process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Compositional_Flow_module, only : CompFlowDestroy
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_compflow_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call DeallocateArray(this%max_change_ivar)
  call CompFlowDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMCompFlowDestroy

end module PM_Compositional_Flow_class
