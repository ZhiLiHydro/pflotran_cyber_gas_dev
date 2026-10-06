module PM_RT_DAE_class
  !
  ! Process model of MODE GIRT_DAE: reactive transport integrated by PETSc
  ! TS through an equation model (RT_DAE_Model_class) and the TS driver
  ! coupler (PMC_TS_Driver_class). Reuses the input, setup and
  ! initialization of PM_RT_class. The TS driver calls none of the classic
  ! per-step bindings; those overridden here stop the run if called.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
#include "petsc/finclude/petscsnes.h"
  use petscsnes

  use PM_Base_class
  use PM_RT_class
  use TS_Model_class
  use RT_DAE_Model_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_rt_type) :: pm_rt_dae_type
    class(rt_dae_model_type), pointer :: model => null()
  contains
    procedure, public :: Setup => PMRTDAESetup
    procedure, public :: InitializeTimestep => PMRTDAENotInThisMode
    procedure, public :: PreSolve => PMRTDAENotInThisMode
    procedure, public :: UpdateSolution => PMRTDAENotInThisMode
    procedure, public :: TimeCut => PMRTDAENotInThisMode
    procedure, public :: FinalizeTimestep => PMRTDAENotInThisMode
    procedure, public :: Residual => PMRTDAEResidual
    procedure, public :: Jacobian => PMRTDAEJacobian
    procedure, public :: Destroy => PMRTDAEDestroy
  end type pm_rt_dae_type

  public :: PMRTDAECreate, &
            PMRTDAEGetModel

contains

! ************************************************************************** !

function PMRTDAECreate()
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  class(pm_rt_dae_type), pointer :: PMRTDAECreate

  class(pm_rt_dae_type), pointer :: pm

  allocate(pm)
  call PMRTInit(pm)
  pm%name = 'Reactive Transport DAE'
  pm%header = 'REACTIVE TRANSPORT DAE'
  PMRTDAECreate => pm

end function PMRTDAECreate

! ************************************************************************** !

subroutine PMRTDAESetup(this)
  !
  ! Classic reactive transport setup, then the equation model.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Realization_Subsurface_class

  implicit none

  class(pm_rt_dae_type) :: this

  class(realization_subsurface_type), pointer :: realization

  call this%pm_rt_type%Setup()
  call PMRTDAECheckOptions(this)

  ! the model defines the state layout and creates the TS solution's DM
  realization => this%realization
  this%model => RTDAEModelCreate(realization,PETSC_FALSE)

end subroutine PMRTDAESetup

! ************************************************************************** !

subroutine PMRTDAECheckOptions(this)
  !
  ! Rejects classic GIRT time step and Newton options (NUMERICAL_METHODS
  ! TRANSPORT, transport SIMULATION solver diagnostics), which PETSc TS
  ! replaces and which would otherwise be silently ignored.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Option_module

  implicit none

  class(pm_rt_dae_type) :: this

  character(len=MAXSTRINGLENGTH) :: unsupported

  unsupported = ''
  if (Initialized(this%cfl_governor)) &
    call TSModelListAdd(unsupported,'CFL_GOVERNOR')
  if (this%volfrac_change_governor < 1.d0) &
    call TSModelListAdd(unsupported,'VOLUME_FRACTION_CHANGE_GOVERNOR')
  if (this%option%transport%numerical_derivatives) &
    call TSModelListAdd(unsupported,'NUMERICAL_JACOBIAN')
  if (this%refactored_convergence) then
    call TSModelListAdd(unsupported,'REFACTORED_CONVERGENCE')
  else if (this%check_post_convergence) then
    call TSModelListAdd(unsupported, &
                        'ITOL_ABSOLUTE_UPDATE/ITOL_RELATIVE_UPDATE')
  endif
  if (this%dampen_oscillatory_updates .or. &
      Initialized(this%dampen_start_iteration) .or. &
      Initialized(this%dampening_factor_from_input)) &
    call TSModelListAdd(unsupported,'DAMPEN_OSCILLATION/&
                        &DAMPEN_START_ITERATION/DAMPENING_FACTOR')
  if (this%scale_linear_system) &
    call TSModelListAdd(unsupported,'SCALE_LINEAR_SYSTEM')
  if (this%debug_update) call TSModelListAdd(unsupported,'DEBUG_UPDATE')
  ! (written by the classic time steppers and convergence tests only)
  if (this%print_ekg) call TSModelListAdd(unsupported,'PRINT_EKG')
  if (this%option%transport%debug_derivatives) &
    call TSModelListAdd(unsupported,'DEBUG_DERIVATIVES')
  ! the update limits of PMRTCheckUpdatePre (RTDAEPreCheck has its own)
  if (this%realization%reaction%max_dlnC /= 5.d0) &   ! (default 5)
    call TSModelListAdd(unsupported,'MAX_DLNC')
  if (Initialized(this%realization%reaction%truncated_concentration)) &
    call TSModelListAdd(unsupported,'TRUNCATE_CONCENTRATION')
  if (len_trim(unsupported) > 0) then
    this%option%io_buffer = 'MODE GIRT_DAE does not support the classic &
      &GIRT solver options ' // trim(unsupported) // ' (PETSc TS controls &
      &the time steps; see -tran_ts_* and -ts_model_* options).'
    call PrintErrMsg(this%option)
  endif

end subroutine PMRTDAECheckOptions

! ************************************************************************** !

subroutine PMRTDAEGetModel(pm,model)
  !
  ! The TS driver's equation model: the transport model.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(pm_base_type) :: pm
  class(ts_model_type), pointer :: model

  nullify(model)
  select type(pm)
    class is(pm_rt_dae_type)
      model => pm%model
  end select

end subroutine PMRTDAEGetModel

! ************************************************************************** !

subroutine PMRTDAENotInThisMode(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  class(pm_rt_dae_type) :: this

  call PrintErrMsg(this%option,'A classic GIRT time step routine was called &
    &under MODE GIRT_DAE, where PETSc TS owns the time steps.')

end subroutine PMRTDAENotInThisMode

! ************************************************************************** !

subroutine PMRTDAEResidual(this,snes,xx,r,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  class(pm_rt_dae_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMRTDAENotInThisMode(this)

end subroutine PMRTDAEResidual

! ************************************************************************** !

subroutine PMRTDAEJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  class(pm_rt_dae_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call PMRTDAENotInThisMode(this)

end subroutine PMRTDAEJacobian

! ************************************************************************** !

subroutine PMRTDAEDestroy(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  class(pm_rt_dae_type) :: this

  if (associated(this%model)) then
    call this%model%Destroy()
    deallocate(this%model)
    nullify(this%model)
  endif
  call this%pm_rt_type%Destroy()

end subroutine PMRTDAEDestroy

end module PM_RT_DAE_class
