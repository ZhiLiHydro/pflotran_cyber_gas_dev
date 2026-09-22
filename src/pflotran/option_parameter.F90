module Option_Parameter_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: parameter_option_type
    PetscInt :: num_parameters
    character(len=MAXWORDLENGTH), pointer :: parameter_names(:)
  end type parameter_option_type

  public :: OptionParameterCreate, &
            OptionParameterInit, &
            OptionParameterDestroy, &
            OptionParameterGetIDFromName

contains

! ************************************************************************** !

function OptionParameterCreate()
  !
  ! Allocates and initializes a new OptionParameter object
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/24

  type(parameter_option_type), pointer :: OptionParameterCreate

  type(parameter_option_type), pointer :: option

  allocate(option)

  call OptionParameterInit(option)
  OptionParameterCreate => option

end function OptionParameterCreate

! ************************************************************************** !

subroutine OptionParameterInit(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/24

  type(parameter_option_type) :: option

  option%num_parameters = 0
  nullify(option%parameter_names)

end subroutine OptionParameterInit

! ************************************************************************** !

function OptionParameterGetIDFromName(parameter_name,parameter_option,driver)
  !
  ! Returns the parameter id/index given a name. Lives here (not Parameter_module)
  ! so chemistry/pflotran_rxn can resolve names without linking dataset-aware
  ! parameter.o.
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/23 / 07/11/26

  use Driver_class
  use String_module

  character(len=*) :: parameter_name
  type(parameter_option_type) :: parameter_option
  class(driver_type) :: driver

  PetscInt :: OptionParameterGetIDFromName

  PetscInt :: i
  character(len=MAXSTRINGLENGTH) :: string

  OptionParameterGetIDFromName = 0
  do i = 1, parameter_option%num_parameters
    if (StringCompare(parameter_name,parameter_option%parameter_names(i))) then
      OptionParameterGetIDFromName = i
      return
    endif
  enddo
  if (parameter_option%num_parameters > 0) then
    string = 'Parameter "' // trim(parameter_name) // '" was not &
      &found among available parameters: ' // &
      StringsMerge(parameter_option%parameter_names,',')
  else
    string = 'Parameter "' // trim(parameter_name) // '" was &
      &requested in OptionParameterGetIDFromName(), but no PARAMETERs &
      &were included in the input deck.'
  endif
  call driver%PrintErrMsg(string)

end function OptionParameterGetIDFromName

! ************************************************************************** !

subroutine OptionParameterDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/24

  type(parameter_option_type), pointer :: option

  if (.not.associated(option)) return

  ! cannot use Utility_module:DeallocateArray due to circular dependency
  if (associated(option%parameter_names)) deallocate(option%parameter_names)
  nullify(option%parameter_names)

  deallocate(option)
  nullify(option)

end subroutine OptionParameterDestroy

end module Option_Parameter_module
