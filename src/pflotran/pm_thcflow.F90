module PM_THCFlow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module
  use THCFlow_Aux_module

  implicit none

  private

  PetscInt, parameter :: MAX_CHANGE_LIQ_PRES_NI = 1
  PetscInt, parameter :: MAX_CHANGE_TEMP_NI = 2
  PetscInt, parameter :: MAX_CHANGE_CONC_NI = 3
  PetscInt, parameter :: MAX_RES_LIQ_EQ = 4
  PetscInt, parameter :: MAX_RES_ENERGY_EQ = 5
  PetscInt, parameter :: MAX_RES_SOL_EQ = 6

  type, public, extends(pm_subsurface_flow_type) :: pm_thcflow_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: max_allow_liq_pres_change_ni
    PetscReal :: liq_pres_change_ts_governor
    PetscReal :: liq_sat_change_ts_governor
    PetscReal :: temp_change_ts_governor
    PetscReal :: conc_change_ts_governor
    PetscInt :: convergence_flags(MAX_RES_SOL_EQ)
    PetscReal :: convergence_reals(MAX_RES_SOL_EQ)
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMTHCFlowReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMTHCFlowReadTSSelectCase
    procedure, public :: ReadNewtonBlock => PMTHCFlowReadNewtonSelectCase
    procedure, public :: Setup => PMTHCFlowSetup
    procedure, public :: InitializeRun => PMTHCFlowInitializeRun
    procedure, public :: InitializeTimestep => PMTHCFlowInitializeTimestep
    procedure, public :: Residual => PMTHCFlowResidual
    procedure, public :: Jacobian => PMTHCFlowJacobian
    procedure, public :: UpdateTimestep => PMTHCFlowUpdateTimestep
    procedure, public :: FinalizeTimestep => PMTHCFlowFinalizeTimestep
    procedure, public :: PreSolve => PMTHCFlowPreSolve
    procedure, public :: PostSolve => PMTHCFlowPostSolve
    procedure, public :: CheckUpdatePost => PMTHCFlowCheckUpdatePost
    procedure, public :: CheckConvergence => PMTHCFlowCheckConvergence
    procedure, public :: TimeCut => PMTHCFlowTimeCut
    procedure, public :: UpdateSolution => PMTHCFlowUpdateSolution
    procedure, public :: UpdateAuxVars => PMTHCFlowUpdateAuxVars
    procedure, public :: MaxChange => PMTHCFlowMaxChange
    procedure, public :: ComputeMassBalance => PMTHCFlowComputeMassBalance
    procedure, public :: InputRecord => PMTHCFlowInputRecord
    procedure, public :: CheckpointBinary => PMTHCFlowCheckpointBinary
    procedure, public :: RestartBinary => PMTHCFlowRestartBinary
    procedure, public :: Destroy => PMTHCFlowDestroy
  end type pm_thcflow_type

  public :: PMTHCFlowCreate, &
            PMTHCFlowInitObject, &
            PMTHCFlowInitializeRun, &
            PMTHCFlowFinalizeTimestep, &
            PMTHCFlowDestroy

contains

! ************************************************************************** !

function PMTHCFlowCreate()
  !
  ! Creates THCFlow process model shell
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  implicit none

  class(pm_thcflow_type), pointer :: PMTHCFlowCreate

  class(pm_thcflow_type), pointer :: thcflow_pm

  allocate(thcflow_pm)
  call PMTHCFlowInitObject(thcflow_pm)

  PMTHCFlowCreate => thcflow_pm

end function PMTHCFlowCreate

! ************************************************************************** !

subroutine PMTHCFlowInitObject(this)
  !
  ! Initializes THCFlow process model shell
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Option_module
  use String_module
  use Variables_module

  implicit none

  class(pm_thcflow_type) :: this

  call PMSubsurfaceFlowInit(this)
  this%name = 'THC Flow'
  this%header = 'THC FLOW'

  nullify(this%max_change_ivar)

  ! set to UNINITIALIZED_DOUBLE and report error below if set from input
  this%pressure_change_governor = UNINITIALIZED_DOUBLE
  this%temperature_change_governor = UNINITIALIZED_DOUBLE
  this%saturation_change_governor = UNINITIALIZED_DOUBLE
  this%xmol_change_governor = UNINITIALIZED_DOUBLE

  this%check_post_convergence = PETSC_TRUE

  this%max_allow_liq_pres_change_ni = UNINITIALIZED_DOUBLE
  this%liq_pres_change_ts_governor = 5.d5    ! [Pa]
  this%liq_sat_change_ts_governor = 1.d0
  this%temp_change_ts_governor = 5.d0        ! [C]
  this%conc_change_ts_governor = 1.d0        ! [mol/L]

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

  thcflow_debug_cell_id = UNINITIALIZED_INTEGER

end subroutine PMTHCFlowInitObject

! ************************************************************************** !

subroutine PMTHCFlowReadSimOptionsBlock(this,input)
  !
  ! Read THCFLOW options input block
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use THCFlow_module
  use THCFlow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module

  implicit none

  class(pm_thcflow_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'THCFlow Options'

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
      case('NO_ACCUMULATION')
        thcflow_calc_accum = PETSC_FALSE
      case('NO_FLUX')
        thcflow_calc_flux = PETSC_FALSE
      case('NO_BCFLUX')
        thcflow_calc_bcflux = PETSC_FALSE
      case('TENSORIAL_RELATIVE_PERMEABILITY')
        thcflow_tensorial_rel_perm = PETSC_TRUE
      case('SOLID_DENSITY','SOLID_GRAIN_DENSITY','ROCK_DENSITY')
        call InputReadDouble(input,option,thcflow_density_solid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thcflow_density_solid,'kg/m^3', &
                                      trim(error_string)//','//keyword,option)
      case('SPECIFIC_HEAT_LIQUID')
        call InputReadDouble(input,option,thcflow_specific_heat_liquid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thcflow_specific_heat_liquid, &
                                      'J/kg-C',trim(error_string)//','// &
                                      keyword,option)
      case('SPECIFIC_HEAT_SOLID')
        call InputReadDouble(input,option,thcflow_specific_heat_solid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thcflow_specific_heat_solid, &
                                      'J/kg-C',trim(error_string)//','// &
                                      keyword,option)
      case('THERMAL_CONDUCTIVITY_SOLID','SOLID_THERMAL_CONDUCTIVITY')
        call InputReadDouble(input,option,thcflow_kappa_solid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thcflow_kappa_solid,'W/m-C', &
                                      trim(error_string)//','//keyword,option)
      case('THERMAL_CONDUCTIVITY_DRY')
        call InputReadDouble(input,option,thcflow_kappa_dry)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thcflow_kappa_dry,'W/m-C', &
                                      trim(error_string)//','//keyword,option)
      case('THERMAL_CONDUCTIVITY_WET')
        call InputReadDouble(input,option,thcflow_kappa_wet)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thcflow_kappa_wet,'W/m-C', &
                                      trim(error_string)//','//keyword,option)
      case('ENERGY_FORMULATION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(trim(word))
          case('RHO_CP_T','RHOCPT')
            thcflow_energy_mode = THCFLOW_ENERGY_RHO_CP_T
          case('FULL_EOS','FULLEOS','EOS')
            thcflow_energy_mode = THCFLOW_ENERGY_FULL_EOS
          case default
            call InputKeywordUnrecognized(input,word, &
                   'THCFlow Mode,ENERGY_FORMULATION',option)
        end select
      case('DEBUG_CELL_ID')
        call InputReadInt(input,option,thcflow_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,'THCFlow Mode',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMTHCFlowReadSimOptionsBlock

! ************************************************************************** !

subroutine PMTHCFlowReadTSSelectCase(this,input,keyword,found, &
                                     error_string,option)
  !
  ! Read timestepper settings specific to the THCFLOW process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_thcflow_type) :: this
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
    case('TEMP_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%temp_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('CONC_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%conc_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case default
      found = PETSC_FALSE
  end select

end subroutine PMTHCFlowReadTSSelectCase

! ************************************************************************** !

subroutine PMTHCFlowReadNewtonSelectCase(this,input,keyword,found, &
                                         error_string,option)
  !
  ! Reads input file parameters associated with the THCFLOW process model
  ! Newton solver convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use THCFlow_Aux_module

  implicit none

  class(pm_thcflow_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  error_string = 'THCFLOW Newton Solver'

  select case(trim(keyword))
    case('ITOL_UPDATE')
      option%io_buffer = 'ITOL_UPDATE not supported with THCFLOW. Please &
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
      call InputReadDouble(input,option,thcflow_rel_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      ! no units conversion since it is relative
    case('MIN_LIQ_PRESSURE_PERTURBATION')
      call InputReadDouble(input,option,thcflow_pres_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      call InputReadAndConvertUnits(input,thcflow_pres_min_pert, &
                                    'Pa',keyword,option)
    case('MIN_TEMPERATURE_PERTURBATION')
      call InputReadDouble(input,option,thcflow_temp_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
    case('MIN_CONCENTRATION_PERTURBATION')
      call InputReadDouble(input,option,thcflow_conc_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
    case('MAX_ALLOW_LIQ_PRES_CHANGE_NI')
      call InputReadDouble(input,option,this%max_allow_liq_pres_change_ni)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%max_allow_liq_pres_change_ni, &
                                    'Pa',keyword,option)
    case default
      found = PETSC_FALSE

  end select

end subroutine PMTHCFlowReadNewtonSelectCase

! ************************************************************************** !

subroutine PMTHCFlowSetup(this)
  !
  ! Sets up auxvars and parameters
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use THCFlow_module

  implicit none

  class(pm_thcflow_type) :: this

  call this%SetRealization()
  call THCFlowSetup(this%realization)
  call PMSubsurfaceFlowSetup(this)

end subroutine PMTHCFlowSetup

! ************************************************************************** !

recursive subroutine PMTHCFlowInitializeRun(this)
  !
  ! Initializes the THCFLOW mode run.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Realization_Base_class
  use Patch_module
  use Field_module
  use Material_Aux_module
  use Option_module
  use Variables_module

  implicit none

  class(pm_thcflow_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option

  patch => this%realization%patch
  field => this%realization%field
  option => this%option

  if (this%steady_state) thcflow_calc_accum = PETSC_FALSE

  ! THCFLOW always solves all three equations; the max change variables
  ! are pressure, saturation, temperature, and concentration.
  allocate(this%max_change_ivar(4))
  this%max_change_ivar(1) = LIQUID_PRESSURE
  this%max_change_ivar(2) = LIQUID_SATURATION
  this%max_change_ivar(3) = TEMPERATURE
  this%max_change_ivar(4) = SOLUTE_CONCENTRATION

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

  if (soil_compressibility_index == 0 .and. associated(option%inversion)) then
    option%io_buffer = 'Soil compressibility must be employed for THCFlow &
      &when used for inversion.'
    call PrintErrMsg(option)
  endif

end subroutine PMTHCFlowInitializeRun

! ************************************************************************** !

subroutine PMTHCFlowInitializeTimestep(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use THCFlow_module, only : THCFlowInitializeTimestep

  implicit none

  class(pm_thcflow_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call THCFlowInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMTHCFlowInitializeTimestep

! ************************************************************************** !

subroutine PMTHCFlowFinalizeTimestep(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  implicit none

  class(pm_thcflow_type) :: this

  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMTHCFlowFinalizeTimestep

! ************************************************************************** !

subroutine PMTHCFlowPreSolve(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  implicit none

  class(pm_thcflow_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMTHCFlowPreSolve

! ************************************************************************** !

subroutine PMTHCFlowPostSolve(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  implicit none

  class(pm_thcflow_type) :: this

end subroutine PMTHCFlowPostSolve

! ************************************************************************** !

subroutine PMTHCFlowUpdateTimestep(this,update_dt, &
                                   dt,dt_min,dt_max,iacceleration, &
                                   num_newton_iterations,tfac, &
                                   time_step_max_growth_factor)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Option_module
  use Utility_module, only : Equal

  implicit none

  class(pm_thcflow_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: sat_ratio, pres_ratio, temp_ratio, conc_ratio
  PetscReal :: min_ratio
  PetscReal :: dt_prev

  if (update_dt .and. iacceleration /= 0) then
    dt_prev = dt
    ! calculate the time step ramping factor
    sat_ratio = (2.d0*this%liq_sat_change_ts_governor)/ &
                (this%liq_sat_change_ts_governor+this%max_saturation_change)
    pres_ratio = (2.d0*this%liq_pres_change_ts_governor)/ &
                (this%liq_pres_change_ts_governor+this%max_pressure_change)
    temp_ratio = (2.d0*this%temp_change_ts_governor)/ &
                (this%temp_change_ts_governor+this%max_temperature_change)
    conc_ratio = (2.d0*this%conc_change_ts_governor)/ &
                (this%conc_change_ts_governor+this%max_xmol_change)
    min_ratio = min(sat_ratio,pres_ratio,temp_ratio,conc_ratio)
    ! pick minimum time step from calc'd ramping factor or maximum ramping factor
    dt = min(min_ratio*dt,time_step_max_growth_factor*dt)
    ! make sure time step is within bounds given in the input deck
    dt = min(dt,dt_max)
    if (this%logging_verbosity > 0) then
      if (Equal(dt,dt_max)) then
        string = 'maximum time step size'
      else if (min_ratio > time_step_max_growth_factor) then
        string = 'maximum time step growth factor'
      else if (Equal(min_ratio,sat_ratio)) then
        string = 'liquid saturation governor'
      else if (Equal(min_ratio,pres_ratio)) then
        string = 'liquid pressure governor'
      else if (Equal(min_ratio,temp_ratio)) then
        string = 'temperature governor'
      else
        string = 'concentration governor'
      endif
      string = 'TS update: ' // trim(string)
      call PrintMsg(this%option,string)
    endif
  endif

  if (Initialized(this%cfl_governor)) then
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMTHCFlowUpdateTimestep

! ************************************************************************** !

subroutine PMTHCFlowResidual(this,snes,xx,r,ierr)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use THCFlow_module, only : THCFlowResidual
  use Debug_module
  use Grid_module

  implicit none

  class(pm_thcflow_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  Mat :: M

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  ! calculate residual
  if (thcflow_simult_function_evals) then
    call SNESGetJacobian(snes,M,PETSC_NULL_MAT,PETSC_NULL_FUNCTION, &
                         PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
    call THCFlowResidual(snes,xx,r,M,this%realization,this%debug,ierr)
  else
    call THCFlowResidual(snes,xx,r,PETSC_NULL_MAT,this%realization, &
                         this%debug,ierr)
  endif

  if (this%debug%vecview_residual) then
    call DebugVecView(this%debug,r,'THCFresidual.out',this%option)
  endif
  if (this%debug%vecview_solution) then
    call DebugVecView(this%debug,xx,'THCFxx.out',this%option)
  endif

end subroutine PMTHCFlowResidual

! ************************************************************************** !

subroutine PMTHCFlowJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Debug_module
  use Option_module

  implicit none

  class(pm_thcflow_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  PetscReal :: norm

  ! the Jacobian was already calculated in PMTHCFlowResidual

  if (this%debug%matview_Matrix) then
    call DebugMatView(this%debug,A,'THCFjacobian','', &
                      thcflow_ts_count,thcflow_ts_cut_count, &
                      thcflow_ni_count,this%option)
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

  thcflow_ni_count = thcflow_ni_count + 1

end subroutine PMTHCFlowJacobian

! ************************************************************************** !

subroutine PMTHCFlowCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  !
  ! Tracks the maximum absolute change in each primary variable (P, T, C)
  ! over a Newton iteration so that PMTHCFlowCheckConvergence can enforce a
  ! solution-change convergence criterion.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Material_Aux_module
  use THCFlow_Aux_module

  implicit none

  class(pm_thcflow_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch

  PetscReal, pointer :: dX_p(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscReal :: tempreal

  PetscReal :: max_abs_pressure_change_NI
  PetscInt :: max_abs_pressure_change_NI_cell
  PetscReal :: max_abs_temp_change_NI
  PetscInt :: max_abs_temp_change_NI_cell
  PetscReal :: max_abs_conc_change_NI
  PetscInt :: max_abs_conc_change_NI_cell

  grid => this%realization%patch%grid
  option => this%realization%option
  patch => this%realization%patch

  ! If these are changed from true, we must add a global reduction on both
  ! variables to ensure that their values match across all processes. Otherwise
  ! PETSc will throw an error in debug mode or ignore the error in optimized.
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  ! reset the solution-change flags each Newton iteration (THCFLOW has no
  ! pre-check, which is where ZFLOW performs this reset)
  this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = 0
  this%convergence_flags(MAX_CHANGE_TEMP_NI) = 0
  this%convergence_flags(MAX_CHANGE_CONC_NI) = 0

  call VecGetArray(dX,dX_p,ierr);CHKERRQ(ierr)
  max_abs_pressure_change_NI = 0.d0
  max_abs_pressure_change_NI_cell = 0
  max_abs_temp_change_NI = 0.d0
  max_abs_temp_change_NI_cell = 0
  max_abs_conc_change_NI = 0.d0
  max_abs_conc_change_NI_cell = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    offset = (local_id-1)*option%nflowdof
    tempreal = dabs(dX_p(offset+thcflow_pressure_dof))
    if (tempreal > dabs(max_abs_pressure_change_NI)) then
      max_abs_pressure_change_NI_cell = grid%nG2A(ghosted_id)
      max_abs_pressure_change_NI = tempreal
    endif
    tempreal = dabs(dX_p(offset+thcflow_temperature_dof))
    if (tempreal > dabs(max_abs_temp_change_NI)) then
      max_abs_temp_change_NI_cell = grid%nG2A(ghosted_id)
      max_abs_temp_change_NI = tempreal
    endif
    tempreal = dabs(dX_p(offset+thcflow_concentration_dof))
    if (tempreal > dabs(max_abs_conc_change_NI)) then
      max_abs_conc_change_NI_cell = grid%nG2A(ghosted_id)
      max_abs_conc_change_NI = tempreal
    endif
  enddo

  ! flag non-convergence in liquid pressure when a tolerance has been set
  if (Initialized(this%max_allow_liq_pres_change_ni) .and. &
      max_abs_pressure_change_NI > this%max_allow_liq_pres_change_ni) then
    this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = &
      max_abs_pressure_change_NI_cell
  endif

  ! the following reals are for REPORTING purposes only
  this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI) = max_abs_pressure_change_NI
  this%convergence_reals(MAX_CHANGE_TEMP_NI) = max_abs_temp_change_NI
  this%convergence_reals(MAX_CHANGE_CONC_NI) = max_abs_conc_change_NI

  call VecRestoreArray(dX,dX_p,ierr);CHKERRQ(ierr)

end subroutine PMTHCFlowCheckUpdatePost

! ************************************************************************** !

subroutine PMTHCFlowCheckConvergence(this,snes,it,xnorm,unorm, &
                                     fnorm,reason,ierr)
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Material_Aux_module
  use THCFlow_Aux_module
  use Convergence_module

  implicit none

  class(pm_thcflow_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  Vec :: residual_vec
  PetscReal, pointer :: r_p(:)
  character(len=10) :: reason_string

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(thcflow_auxvar_type), pointer :: thcflow_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: converged_flag

  PetscReal :: max_abs_res_liq_
  PetscInt :: max_abs_res_liq_cell
  PetscReal :: max_abs_res_energy_
  PetscInt :: max_abs_res_energy_cell
  PetscReal :: max_abs_res_sol_
  PetscInt :: max_abs_res_sol_cell
  PetscInt :: offset
  PetscMPIInt :: int_mpi

  PetscReal :: residual
  PetscReal :: tempreal

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  thcflow_auxvars => patch%aux%THCFlow%auxvars
  material_auxvars => patch%aux%Material%auxvars

  max_abs_res_liq_ = 0.d0
  max_abs_res_liq_cell = 0
  max_abs_res_energy_ = 0.d0
  max_abs_res_energy_cell = 0
  max_abs_res_sol_ = 0.d0
  max_abs_res_sol_cell = 0

  residual_vec = field%flow_r
  ! check residual terms
  call VecGetArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    residual = r_p(offset+thcflow_pressure_dof)
    tempreal = dabs(residual)
    if (tempreal > max_abs_res_liq_) then
      max_abs_res_liq_ = tempreal
      max_abs_res_liq_cell = grid%nG2A(ghosted_id)
    endif
    residual = r_p(offset+thcflow_temperature_dof)
    tempreal = dabs(residual)
    if (tempreal > max_abs_res_energy_) then
      max_abs_res_energy_ = tempreal
      max_abs_res_energy_cell = grid%nG2A(ghosted_id)
    endif
    residual = r_p(offset+thcflow_concentration_dof)
    tempreal = dabs(residual)
    if (tempreal > max_abs_res_sol_) then
      max_abs_res_sol_ = tempreal
      max_abs_res_sol_cell = grid%nG2A(ghosted_id)
    endif
  enddo

  ! the following flags are for REPORTING purposes only
  this%convergence_flags(MAX_RES_LIQ_EQ) = max_abs_res_liq_cell
  this%convergence_reals(MAX_RES_LIQ_EQ) = max_abs_res_liq_
  this%convergence_flags(MAX_RES_ENERGY_EQ) = max_abs_res_energy_cell
  this%convergence_reals(MAX_RES_ENERGY_EQ) = max_abs_res_energy_
  this%convergence_flags(MAX_RES_SOL_EQ) = max_abs_res_sol_cell
  this%convergence_reals(MAX_RES_SOL_EQ) = max_abs_res_sol_

  int_mpi = size(this%convergence_flags)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_flags,int_mpi,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  int_mpi = size(this%convergence_reals)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_reals,int_mpi, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)

  ! these conditionals cannot change order
  reason_string = '----| '
  converged_flag = CONVERGENCE_CONVERGED
  if (this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) > 0) then
    reason_string(1:1) = 'P'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_TEMP_NI) > 0) then
    reason_string(2:2) = 'T'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_CONC_NI) > 0) then
    reason_string(3:3) = 'C'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif

  if (Initialized(this%max_allow_liq_pres_change_ni)) then
    option%convergence = converged_flag
  else
    ! forced standard 2 norms
    option%convergence = CONVERGENCE_OFF
  endif

  call VecRestoreArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMTHCFlowCheckConvergence

! ************************************************************************** !

subroutine PMTHCFlowTimeCut(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THCFlow_module, only : THCFlowTimeCut

  implicit none

  class(pm_thcflow_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call THCFlowTimeCut(this%realization)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMTHCFlowTimeCut

! ************************************************************************** !

subroutine PMTHCFlowUpdateSolution(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THCFlow_module, only : THCFlowUpdateSolution, &
                             THCFlowMapBCAuxVarsToGlobal

  implicit none

  class(pm_thcflow_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call THCFlowUpdateSolution(this%realization)
  call THCFlowMapBCAuxVarsToGlobal(this%realization)

end subroutine PMTHCFlowUpdateSolution

! ************************************************************************** !

subroutine PMTHCFlowUpdateAuxVars(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  use THCFlow_module, only : THCFlowUpdateAuxVars

  implicit none

  class(pm_thcflow_type) :: this

  call THCFlowUpdateAuxVars(this%realization)

end subroutine PMTHCFlowUpdateAuxVars

! ************************************************************************** !

subroutine PMTHCFlowMaxChange(this)
  !
  ! Computes the maximum change in the solution for the timestep governors.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use String_module
  use THCFlow_Aux_module

  implicit none

  class(pm_thcflow_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal, allocatable :: max_change_global(:)
  PetscReal :: max_change, change
  PetscInt :: i, j
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

  write(option%io_buffer,'("  --> max change: dpl= ",1pe12.4," dsl= ",&
                         &1pe12.4," dt= ",1pe12.4," dc= ",1pe12.4)') &
    max_change_global(1:4)
  call PrintMsg(option)
  this%max_pressure_change = max_change_global(1)
  this%max_saturation_change = max_change_global(2)
  this%max_temperature_change = max_change_global(3)
  ! hijacking xmol_change for concentration
  this%max_xmol_change = max_change_global(4)

  deallocate(max_change_global)

end subroutine PMTHCFlowMaxChange

! ************************************************************************** !

subroutine PMTHCFlowComputeMassBalance(this,mass_balance_array)
  !
  ! Computes the current solute/water/energy mass in the domain
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THCFlow_module, only : THCFlowComputeMassBalance

  implicit none

  class(pm_thcflow_type) :: this
  PetscReal :: mass_balance_array(:)

  call THCFlowComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTHCFlowComputeMassBalance

! ************************************************************************** !

subroutine PMTHCFlowInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  class(pm_thcflow_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'thcflow'

end subroutine PMTHCFlowInputRecord

! ************************************************************************** !

subroutine PMTHCFlowCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with THCFlow PM
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Checkpoint_module
  use Global_module

  implicit none

  class(pm_thcflow_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMTHCFlowCheckpointBinary

! ************************************************************************** !

subroutine PMTHCFlowRestartBinary(this,viewer)
  !
  ! Restarts data associated with THCFlow PM
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Checkpoint_module
  use Global_module

  implicit none

  class(pm_thcflow_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMTHCFlowRestartBinary

! ************************************************************************** !

subroutine PMTHCFlowDestroy(this)
  !
  ! Destroys THCFlow process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THCFlow_module, only : THCFlowDestroy
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_thcflow_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call DeallocateArray(this%max_change_ivar)
  call THCFlowDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMTHCFlowDestroy

end module PM_THCFlow_class
