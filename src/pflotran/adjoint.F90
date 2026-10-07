module Adjoint_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Adjoint_Cost_module
  use Adjoint_Param_module
  use TS_Optimizer_class
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
  !     OPTIMIZATION                         [optional]
  !       ...
  !     /
  !   /
  !
  ! COST_FUNCTION and PARAMETER may be repeated. OPTIMIZATION minimizes the
  ! sum of the costs over the parameters (see AdjointReadOptimization).
  !
  ! Author: Hansol Suh
  ! Date: 2026-05-27

  implicit none

  private

  ! OUTPUT FULL_FIELD: also write the initial sensitivities as binary Vecs
  PetscBool, public :: adjoint_output_full_field = PETSC_FALSE
  ! OPTIMIZATION: the optimizer, until the factory hands it to the TS driver
  type(ts_optimizer_type), pointer, public :: adjoint_optimizer => null()

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
  call AdjointDestroyOptimizer()

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
      case('OPTIMIZATION')
        call AdjointReadOptimization(input,option)
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
  call AdjointDestroyOptimizer()

end subroutine AdjointDestroy

! ************************************************************************** !

subroutine AdjointReadOptimization(input,option)
  !
  ! Reads the OPTIMIZATION block, which minimizes the sum of the costs over
  ! the parameters with PETSc TAO from their deck values:
  !   OPTIMIZATION
  !     TYPE <TAO type>                     [default BQNLS]
  !     MAXIMUM_ITERATIONS <int > 0>        [TAO's default]
  !     GRADIENT_TOLERANCE <abs> <rel>      [TAO's defaults; >= 0]
  !     PARAMETER <number>                  (in ADJOINT PARAMETER order)
  !       TRANSFORM LOG10 | NONE            [LOG10; NONE for a smoothstep
  !                                          interval]
  !       BOUNDS <lower> <upper>            [none; in the parameter's units]
  !     /
  !   /
  ! TAO options take the prefix adj_ (e.g. -adj_tao_monitor).
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: iparam, transform
  PetscReal :: lower, upper

  error_string = 'SUBSURFACE,ADJOINT,OPTIMIZATION'
  if (associated(adjoint_optimizer)) then
    option%io_buffer = 'ADJOINT takes one OPTIMIZATION block.'
    call PrintErrMsg(option)
  endif
  adjoint_optimizer => TSOptimizerCreate()

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToLower(word)
        adjoint_optimizer%tao_type = word
      case('MAXIMUM_ITERATIONS')
        call InputReadInt(input,option,adjoint_optimizer%max_iterations)
        call InputErrorMsg(input,option,keyword,error_string)
        ! TAO reads 0 and negative values as errors or as its defaults
        if (adjoint_optimizer%max_iterations < 1) then
          call PrintErrMsg(option,trim(error_string) // ',' // &
            trim(keyword) // ' must be positive.')
        endif
      case('GRADIENT_TOLERANCE')
        call InputReadDouble(input,option,adjoint_optimizer%gatol)
        call InputErrorMsg(input,option,'absolute gradient tolerance', &
                           error_string)
        call InputReadDouble(input,option,adjoint_optimizer%grtol)
        call InputErrorMsg(input,option,'relative gradient tolerance', &
                           error_string)
        if (adjoint_optimizer%gatol < 0.d0 .or. &
            adjoint_optimizer%grtol < 0.d0) then
          call PrintErrMsg(option,trim(error_string) // ',' // &
            trim(keyword) // ' must not be negative.')
        endif
      case('PARAMETER')
        call InputReadInt(input,option,iparam)
        call InputErrorMsg(input,option,'parameter number',error_string)
        transform = TS_OPT_TRANSFORM_DEFAULT
        lower = UNINITIALIZED_DOUBLE
        upper = UNINITIALIZED_DOUBLE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword', &
                             trim(error_string) // ',PARAMETER')
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('TRANSFORM')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,keyword, &
                                 trim(error_string) // ',PARAMETER')
              call StringToUpper(word)
              select case(trim(word))
                case('LOG10')
                  transform = TS_OPT_TRANSFORM_LOG10
                case('NONE')
                  transform = TS_OPT_TRANSFORM_NONE
                case default
                  call InputKeywordUnrecognized(input,word, &
                    trim(error_string) // ',PARAMETER,TRANSFORM',option)
              end select
            case('BOUNDS')
              call InputReadDouble(input,option,lower)
              call InputErrorMsg(input,option,'lower bound', &
                                 trim(error_string) // ',PARAMETER')
              call InputReadDouble(input,option,upper)
              call InputErrorMsg(input,option,'upper bound', &
                                 trim(error_string) // ',PARAMETER')
            case default
              call InputKeywordUnrecognized(input,keyword, &
                trim(error_string) // ',PARAMETER',option)
          end select
        enddo
        call InputPopBlock(input,option)
        call adjoint_optimizer%AddParameterSettings(iparam,transform, &
                                                    lower,upper)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine AdjointReadOptimization

! ************************************************************************** !

subroutine AdjointDestroyOptimizer()
  !
  ! An optimizer the factory did not hand over.
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  implicit none

  if (.not.associated(adjoint_optimizer)) return
  call adjoint_optimizer%Destroy()
  deallocate(adjoint_optimizer)
  nullify(adjoint_optimizer)

end subroutine AdjointDestroyOptimizer

! ************************************************************************** !

end module Adjoint_module
