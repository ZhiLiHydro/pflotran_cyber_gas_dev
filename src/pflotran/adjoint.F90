module Adjoint_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Adjoint_Cost_module
  use Adjoint_Param_module
  use PFLOTRAN_Constants_module

  ! Reader for the ADJOINT deck block (MODE GIRT_DAE only, enforced by
  ! factory_subsurface_read).
  !
  !   ADJOINT
  !     COST_FUNCTION <class>
  !       ...
  !     /
  !     PARAMETER <class>
  !       ...
  !     /
  !   /
  !
  ! COST_FUNCTION and PARAMETER may be repeated.
  !
  ! Author: Hansol Suh
  ! Date: 2026-05-27

  implicit none

  private

  ! OUTPUT FULL_FIELD: also write the initial sensitivities as binary Vecs
  PetscBool, public :: adjoint_output_full_field = PETSC_FALSE

  public :: AdjointReadInput, &
            AdjointIsActive, &
            AdjointDestroy

contains

! ************************************************************************** !

subroutine AdjointInit()
  !
  ! Initializes the cost and parameter lists before the deck is read.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  call AdjCostInit()
  call AdjParamInit()
  adjoint_output_full_field = PETSC_FALSE

end subroutine AdjointInit

! ************************************************************************** !

subroutine AdjointReadInput(input,option)
  !
  ! Reads the ADJOINT block.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SUBSURFACE,ADJOINT'

  call AdjointInit()

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('COST_FUNCTION')
        call AdjCostRead(input,option)
      case('PARAMETER')
        call AdjParamRead(input,option)
      case('OUTPUT')
        call AdjointReadOutput(input,option)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (AdjCostCount() == 0) then
    option%io_buffer = 'ADJOINT requires at least one COST_FUNCTION; ' // &
      'there is no objective to differentiate.'
    call PrintErrMsg(option)
  endif
  if (AdjParamCount() == 0) then
    option%io_buffer = 'ADJOINT requires at least one PARAMETER because ' // &
      'zero-parameter adjoint output is not implemented.'
    call PrintErrMsg(option)
  endif

end subroutine AdjointReadInput

! ************************************************************************** !

function AdjointIsActive()
  !
  ! True when the deck registered at least one cost and one parameter.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  PetscBool :: AdjointIsActive

  AdjointIsActive = associated(adj_cost_list) .and. &
                    associated(adj_param_list)

end function AdjointIsActive

! ************************************************************************** !

subroutine AdjointReadOutput(input,option)
  !
  ! Reads the OUTPUT block; the summary is always printed.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call StringToUpper(word)
    select case(trim(word))
      case('FULL_FIELD')
        adjoint_output_full_field = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,word,'ADJOINT,OUTPUT',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine AdjointReadOutput

! ************************************************************************** !

subroutine AdjointDestroy()
  !
  ! Destroys the cost and parameter lists.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  call AdjCostDestroyList()
  call AdjParamDestroyList()
  adjoint_output_full_field = PETSC_FALSE

end subroutine AdjointDestroy

! ************************************************************************** !

end module Adjoint_module
