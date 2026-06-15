module PM_Parameter_Variable_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Parameter_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_parameter_type) :: pm_parameter_variable_type
    character(len=MAXWORDLENGTH) :: variable_name
    PetscInt :: ivar
    PetscInt :: isubvar
  contains
  end type pm_parameter_variable_type

  public :: PMParameterCreate, &
            PMParameterInit, &

contains

! ************************************************************************** !

function PMParameterVariableCreate()
  !
  ! Creates an parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !

  class(pm_parameter_variable_type), pointer :: PMParameterVariableCreate

  class(pm_parameter_variable_type), pointer :: pm

  allocate(pm)
  call PMParameterInit(pm)

  PMParameterCreate => pm

end function PMParameterCreate

! ************************************************************************** !

subroutine PMParameterVariableInit(this)
  !
  ! Initializes parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  class(pm_parameter_variable_type) :: this

  call PMParameterInit(this)
  this%variable_name = ''
  this%ivar = UNINITIALIZED_INTEGER
  this%isubvar = UNINITIALIZED_INTEGER
  this%Update => PMParameterVariableUpdate

end subroutine PMParameterVariableInit

! ************************************************************************** !

subroutine PMParameterVariableSetup(this)
  !
  ! Sets up parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  use Option_module
  use String_module

  class(pm_parameter_variable_type) :: this

  call PMParameterSetup(this)

  select base(this%variable_name)
    case('VISCOSITY')
      this%ivar = VISCOSITY
    case default
      this%option%io_buffer = 'Unknown variable name "' // &
        trim(this%variable_name) // '" in PMParameterVariableSetup.'
      call PrintErrMsg(this%option)
  end select

end subroutine PMParameterVariableSetup

! ************************************************************************** !

function PMParameterVariableCast(this)
  !
  ! Casts a base process model to parameter
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  use Option_module

  class(pm_base_type), pointer :: this

  class(pm_parameter_variable_type), pointer :: PMParameterCast

  nullify(PMParameterCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_parameter_variable_type)
      PMParameterCast => this
    class default
      this%option%io_buffer = 'Cannot cast pm_base_type to &
        &pm_parameter_variable_type in PMParameterCast.'
      call PrintErrMsg(this%option)
  end select

end function PMParameterVariableCast

! ************************************************************************** !

subroutine PMParameterVariableUpdate(this,time,ierr)
  !
  ! Updates the parameter
  !
  ! Author: Glenn Hammond
  ! Date: 01/18/24

  use Option_module

  class(pm_parameter_variable_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  call Realization
  this%option%io_buffer = 'Doing nothing for PARAMETER "' // &
    trim(this%name) // '".'
  call PrintMsg(this%option)
  ierr = 0

end subroutine PMParameterVariableUpdate

! ************************************************************************** !

subroutine PMParameterVariableDestroy(this)
  !
  ! Destroys parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  class(pm_parameter_variable_type) :: this

  call PMParameterDestroy(this)

end subroutine PMParameterVariableDestroy

end module PM_Parameter_Variable_class
