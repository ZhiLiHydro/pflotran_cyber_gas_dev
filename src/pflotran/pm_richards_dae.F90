module PM_Richards_DAE_class
  !
  ! Process model of MODE RICHARDS_DAE: Richards flow integrated by PETSc TS
  ! through an equation model (Richards_DAE_Model_class) and
  ! PMC_TS_Driver_class. With GIRT_DAE, flow and transport are one model
  ! (Flow_RT_DAE_Model_class).
  !
  ! Reuses the RICHARDS input, setup and initialization and owns the
  ! equation model. The classic per-timestep bindings stop the run.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscsnes.h"
  use petscsnes

  use PM_Base_class
  use PM_Richards_class
  use TS_Model_class
  use Richards_DAE_Model_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_richards_type) :: pm_richards_dae_type
    class(richards_dae_model_type), pointer :: model => null()
    ! part of a coupled flow and transport model, set by the factory
    PetscBool :: coupled = PETSC_FALSE
  contains
    procedure, public :: Setup => PMRichardsDAESetup
    procedure, public :: InitializeTimestep => PMRichardsDAENotInThisMode
    procedure, public :: TimeCut => PMRichardsDAENotInThisMode
    procedure, public :: FinalizeTimestep => PMRichardsDAENotInThisMode
    procedure, public :: Residual => PMRichardsDAEResidual
    procedure, public :: Jacobian => PMRichardsDAEJacobian
    procedure, public :: Destroy => PMRichardsDAEDestroy
  end type pm_richards_dae_type

  public :: PMRichardsDAECreate, &
            PMRichardsDAEGetModel

contains

! ************************************************************************** !

function PMRichardsDAECreate()
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(pm_richards_dae_type), pointer :: PMRichardsDAECreate

  class(pm_richards_dae_type), pointer :: pm

  allocate(pm)
  call PMRichardsInit(pm)
  pm%name = 'Richards Flow DAE'
  pm%header = 'RICHARDS FLOW DAE'
  ! the Newton tolerances and time step governors of PMRichardsCreate do
  ! not apply; unset for PMRichardsDAECheckOptions
  pm%pressure_change_governor = UNINITIALIZED_DOUBLE
  pm%saturation_change_governor = UNINITIALIZED_DOUBLE
  pm%residual_abs_inf_tol = UNINITIALIZED_DOUBLE
  pm%residual_scaled_inf_tol = UNINITIALIZED_DOUBLE
  pm%abs_update_inf_tol = UNINITIALIZED_DOUBLE
  pm%rel_update_inf_tol = UNINITIALIZED_DOUBLE
  nullify(pm%pmwell_ptr)
  PMRichardsDAECreate => pm

end function PMRichardsDAECreate

! ************************************************************************** !

subroutine PMRichardsDAESetup(this)
  !
  ! Classic Richards setup (auxiliary data, initial condition), then the
  ! equation model.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Realization_Subsurface_class

  implicit none

  class(pm_richards_dae_type) :: this

  class(realization_subsurface_type), pointer :: realization

  call this%pm_richards_type%Setup()
  call PMRichardsDAECheckOptions(this)
  realization => this%realization
  this%model => RichardsDAEModelCreate(realization,this%coupled)

end subroutine PMRichardsDAESetup

! ************************************************************************** !

subroutine PMRichardsDAECheckOptions(this)
  !
  ! Rejects classic RICHARDS time step governors and Newton options, which
  ! PETSc TS and the model replace and would otherwise ignore silently.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(pm_richards_dae_type) :: this

  character(len=MAXSTRINGLENGTH) :: unsupported

  if (Initialized(this%pressure_change_governor) .or. &
      Initialized(this%saturation_change_governor) .or. &
      Initialized(this%cfl_governor)) then
    call PrintErrMsg(this%option,'MODE RICHARDS_DAE does not support the &
      &classic RICHARDS time step governors (PETSc TS controls the time &
      &steps; see the -flow_ts_* and -ts_model_* options).')
  endif

  ! unset by PMSubsurfaceFlowInit and PMRichardsDAECreate
  unsupported = ''
  if (Initialized(this%pressure_change_limit)) &
    call TSModelListAdd(unsupported,'PRESSURE_CHANGE_LIMIT')
  if (Initialized(this%saturation_change_limit)) &
    call TSModelListAdd(unsupported,'SATURATION_CHANGE_LIMIT')
  if (Initialized(this%temperature_change_limit)) &
    call TSModelListAdd(unsupported,'TEMPERATURE_CHANGE_LIMIT')
  if (Initialized(this%pressure_dampening_factor)) &
    call TSModelListAdd(unsupported,'PRESSURE_DAMPENING_FACTOR')
  if (this%check_post_convergence) &
    call TSModelListAdd(unsupported,'USE_INFINITY_NORM_CONVERGENCE')
  if (Initialized(this%residual_abs_inf_tol) .or. &
      Initialized(this%residual_scaled_inf_tol)) &
    call TSModelListAdd(unsupported,'RESIDUAL_INF_TOL/RESIDUAL_ABS_INF_TOL/&
                        &ITOL_SCALED_RESIDUAL')
  if (Initialized(this%abs_update_inf_tol) .or. &
      Initialized(this%rel_update_inf_tol)) &
    call TSModelListAdd(unsupported,'UPDATE_INF_TOL/ABS_UPDATE_INF_TOL/&
                        &ITOL_RELATIVE_UPDATE')
  if (this%scale_linear_system) &
    call TSModelListAdd(unsupported,'SCALE_LINEAR_SYSTEM')
  ! (written by the classic time steppers and convergence tests only)
  if (this%print_ekg) call TSModelListAdd(unsupported,'PRINT_EKG')
  if (len_trim(unsupported) > 0) then
    this%option%io_buffer = 'MODE RICHARDS_DAE does not support the classic &
      &RICHARDS Newton solver options ' // trim(unsupported) // ' (PETSc &
      &TS and the model''s Newton test replace them; see the -ts_model_* &
      &options).'
    call PrintErrMsg(this%option)
  endif

end subroutine PMRichardsDAECheckOptions

! ************************************************************************** !

subroutine PMRichardsDAEGetModel(pm,model)
  !
  ! Hands the TS driver the equation model of a RICHARDS_DAE process model.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(pm_base_type) :: pm
  class(ts_model_type), pointer :: model

  nullify(model)
  select type(pm)
    class is(pm_richards_dae_type)
      model => pm%model
  end select

end subroutine PMRichardsDAEGetModel

! ************************************************************************** !

subroutine PMRichardsDAENotInThisMode(this)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(pm_richards_dae_type) :: this

  call PrintErrMsg(this%option,'A classic RICHARDS time step routine was &
    &called under MODE RICHARDS_DAE, where PETSc TS owns the time steps.')

end subroutine PMRichardsDAENotInThisMode

! ************************************************************************** !

subroutine PMRichardsDAEResidual(this,snes,xx,r,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(pm_richards_dae_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMRichardsDAENotInThisMode(this)

end subroutine PMRichardsDAEResidual

! ************************************************************************** !

subroutine PMRichardsDAEJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  class(pm_richards_dae_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call PMRichardsDAENotInThisMode(this)

end subroutine PMRichardsDAEJacobian

! ************************************************************************** !

subroutine PMRichardsDAEDestroy(this)
  !
  ! Coupled, the transport PM is this%next, which the coupler does not
  ! free; free it here, before PMRichardsDestroy would destroy it.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(pm_richards_dae_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
    deallocate(this%next)
    nullify(this%next)
  endif
  if (associated(this%model)) then
    call this%model%Destroy()
    deallocate(this%model)
    nullify(this%model)
  endif
  call this%pm_richards_type%Destroy()

end subroutine PMRichardsDAEDestroy

end module PM_Richards_DAE_class
