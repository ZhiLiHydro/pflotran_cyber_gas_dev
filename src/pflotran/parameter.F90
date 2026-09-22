module Parameter_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Dataset_Base_class

  implicit none

  ! source of parameter values
  PetscInt, parameter, public :: PARAMETER_SCALAR = 1
  PetscInt, parameter, public :: PARAMETER_VARIABLE = 2
  PetscInt, parameter, public :: PARAMETER_DATASET = 3

  ! numeric type of the parameter payload
  PetscInt, parameter, public :: PARAMETER_DATA_REAL = 1
  PetscInt, parameter, public :: PARAMETER_DATA_INTEGER = 2

  ! when the parameter is re-applied on the PMC tree (INITIAL = never)
  PetscInt, parameter, public :: UPDATE_NEVER = 0
  PetscInt, parameter, public :: UPDATE_AFTER_LAST_PM = 1
  PetscInt, parameter, public :: UPDATE_BEFORE_STEP = 2
  PetscInt, parameter, public :: UPDATE_AFTER_STEP = 3
  PetscInt, parameter, public :: UPDATE_BEFORE_FIRST_PM = 4

  PetscInt, parameter, public :: TARGET_PMC_NULL = 0
  PetscInt, parameter, public :: TARGET_PMC_FLOW = 1
  PetscInt, parameter, public :: TARGET_PMC_TRANSPORT = 2
  PetscInt, parameter, public :: TARGET_PMC_GEOMECH = 3
  PetscInt, parameter, public :: TARGET_PMC_GEOPHYSICS = 4

  private

  type, public :: parameter_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    ! SCALAR value (real payload; integers stored as nint(value) when needed)
    PetscReal :: value
    PetscInt :: ivalue
    PetscInt :: itype          ! SCALAR, VARIABLE, or DATASET
    PetscInt :: data_type      ! REAL or INTEGER
    ! named DATASET linkage (not SOURCE/TARGET variables)
    character(len=MAXWORDLENGTH) :: dataset_name
    ! optional get/set linkage to realization variables
    character(len=MAXWORDLENGTH) :: source_variable_name
    character(len=MAXWORDLENGTH) :: target_variable_name
    PetscInt :: source_ivar
    PetscInt :: source_isubvar
    PetscInt :: target_ivar
    PetscInt :: target_isubvar
    ! schedule: INITIAL (UPDATE_NEVER) or BEFORE/AFTER a named PMC
    PetscInt :: when_to_update
    PetscInt :: target_pmc
    PetscInt :: update_order
    ! dataset pointer: named DATASET (shared) or LIST (owned ascii)
    class(dataset_base_type), pointer :: dataset
    PetscBool :: dataset_is_owned
    type(parameter_type), pointer :: next
  end type parameter_type

  public :: ParameterCreate, ParameterDestroy, &
            ParameterRead, ParameterAddToList, &
            ParameterSetup, &
            ParameterQualityCheck, &
            ParameterLinkDatasets, &
            ParameterMapVariables, &
            ParameterMapVariableList, &
            ParameterOrderUpdates

contains

! ************************************************************************** !

function ParameterCreate()
  !
  ! Creates a parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  type(parameter_type), pointer :: ParameterCreate

  type(parameter_type), pointer :: parameter

  allocate(parameter)
  parameter%name = ''
  parameter%id = UNINITIALIZED_INTEGER
  parameter%value = UNINITIALIZED_DOUBLE
  parameter%ivalue = UNINITIALIZED_INTEGER
  parameter%itype = UNINITIALIZED_INTEGER
  parameter%data_type = PARAMETER_DATA_REAL
  parameter%dataset_name = ''
  parameter%source_variable_name = ''
  parameter%target_variable_name = ''
  parameter%source_ivar = UNINITIALIZED_INTEGER
  parameter%source_isubvar = UNINITIALIZED_INTEGER
  parameter%target_ivar = UNINITIALIZED_INTEGER
  parameter%target_isubvar = UNINITIALIZED_INTEGER
  parameter%when_to_update = UPDATE_NEVER
  parameter%target_pmc = TARGET_PMC_NULL
  parameter%update_order = 0
  parameter%dataset_is_owned = PETSC_FALSE
  nullify(parameter%dataset)
  nullify(parameter%next)
  ParameterCreate => parameter

end function ParameterCreate

! ************************************************************************** !

subroutine ParameterRead(parameter,input,option)
  !
  ! Reads a PARAMETER block: SCALAR, INTEGER, LIST, DATASET,
  ! SOURCE_VARIABLE, TARGET_VARIABLE, WHEN_TO_UPDATE.
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  use Option_module
  use Input_Aux_module
  use String_module
  use Dataset_Ascii_class
  use Dataset_module, only : DatasetVerify
  use Time_Storage_module

  type(parameter_type) :: parameter
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: error_str
  character(len=MAXSTRINGLENGTH) :: string
  class(dataset_ascii_type), pointer :: dataset_ascii
  class(dataset_base_type), pointer :: dataset_base
  type(time_storage_type), pointer :: default_time_storage

  error_str = 'PARAMETER'
  nullify(default_time_storage)
  nullify(dataset_base)

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_str)
    call StringToUpper(keyword)

    select case(trim(keyword))

      case('NAME')
        call InputReadWord(input,option,parameter%name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)

      case('DATA_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        call StringToUpper(word)
        select case(trim(word))
          case('REAL')
            parameter%data_type = PARAMETER_DATA_REAL
          case('INTEGER','INT')
            parameter%data_type = PARAMETER_DATA_INTEGER
          case default
            call InputKeywordUnrecognized(input,word, &
                   trim(error_str)//',DATA_TYPE',option)
        end select

      case('SCALAR')
        parameter%itype = PARAMETER_SCALAR
        if (parameter%data_type == PARAMETER_DATA_INTEGER) then
          call InputReadInt(input,option,parameter%ivalue)
          call InputErrorMsg(input,option,keyword,error_str)
          parameter%value = dble(parameter%ivalue)
        else
          call InputReadDouble(input,option,parameter%value)
          call InputErrorMsg(input,option,keyword,error_str)
        endif

      case('INTEGER')
        ! SCALAR integer shorthand: INTEGER <ival>
        parameter%itype = PARAMETER_SCALAR
        parameter%data_type = PARAMETER_DATA_INTEGER
        call InputReadInt(input,option,parameter%ivalue)
        call InputErrorMsg(input,option,keyword,error_str)
        parameter%value = dble(parameter%ivalue)

      case('LIST','FILE')
        ! Time-constant or time-varying ASCII scalar list (optionally
        ! multi-column with time in column 1). Stored as owned ascii dataset.
        parameter%itype = PARAMETER_DATASET
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 1
        if (parameter%data_type == PARAMETER_DATA_INTEGER) then
          dataset_ascii%data_type = DATASET_INTEGER
        else
          dataset_ascii%data_type = DATASET_REAL
        endif
        dataset_ascii%name = parameter%name
        string = trim(error_str) // ' "' // trim(parameter%name) // '",LIST'
        if (StringCompare(keyword,'FILE')) then
          call InputReadFilename(input,option,filename)
          call InputErrorMsg(input,option,'filename',string)
          call DatasetAsciiReadFile(dataset_ascii,filename,'','unitless', &
                                    string,option)
        else
          call DatasetAsciiReadList(dataset_ascii,input,'','unitless', &
                                    string,option)
        endif
        dataset_base => dataset_ascii
        call DatasetVerify(dataset_base,default_time_storage,string,option)
        parameter%dataset => dataset_ascii
        parameter%dataset_is_owned = PETSC_TRUE
        nullify(dataset_ascii)
        nullify(dataset_base)

      case('DATASET')
        parameter%itype = PARAMETER_DATASET
        call InputReadWord(input,option,parameter%dataset_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)

      case('SOURCE_VARIABLE','VARIABLE')
        ! VARIABLE is retained as a synonym for SOURCE_VARIABLE
        parameter%itype = PARAMETER_VARIABLE
        call InputReadWord(input,option,parameter%source_variable_name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)

      case('TARGET_VARIABLE')
        call InputReadWord(input,option,parameter%target_variable_name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)

      case('WHEN_TO_UPDATE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        call StringToUpper(word)
        select case(trim(word))
          case('INITIAL')
            parameter%when_to_update = UPDATE_NEVER
          case('BEFORE_STEP')
            parameter%when_to_update = UPDATE_BEFORE_STEP
            call ParameterReadTargetPMC(parameter,input,option,error_str)
          case('AFTER_STEP')
            parameter%when_to_update = UPDATE_AFTER_STEP
            call ParameterReadTargetPMC(parameter,input,option,error_str)
          case('BEFORE_FIRST_PM')
            parameter%when_to_update = UPDATE_BEFORE_FIRST_PM
          case('AFTER_LAST_PM')
            parameter%when_to_update = UPDATE_AFTER_LAST_PM
          case('AFTER_FLOW')
            parameter%when_to_update = UPDATE_AFTER_STEP
            parameter%target_pmc = TARGET_PMC_FLOW
          case default
            call InputKeywordUnrecognized(input,word, &
                   trim(error_str)//',WHEN_TO_UPDATE',option)
        end select

      case default
        call InputKeywordUnrecognized(input,keyword,error_str,option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (len_trim(parameter%name) == 0) then
    option%io_buffer = 'No NAME specified for a PARAMETER in the input file.'
    call PrintErrMsg(option)
  endif

  ! TARGET_VARIABLE may accompany SCALAR/DATASET; alone it is not a source
  if (Uninitialized(parameter%itype)) then
    if (len_trim(parameter%source_variable_name) > 0) then
      parameter%itype = PARAMETER_VARIABLE
    else if (len_trim(parameter%target_variable_name) > 0 .and. &
             (Initialized(parameter%value) .or. &
              Initialized(parameter%ivalue))) then
      parameter%itype = PARAMETER_SCALAR
    else
      option%io_buffer = 'No SCALAR, INTEGER, LIST, DATASET, or &
        &SOURCE_VARIABLE specified for PARAMETER "' // &
        trim(parameter%name) // '".'
      call PrintErrMsg(option)
    endif
  endif

  if (parameter%itype == PARAMETER_SCALAR) then
    if (parameter%data_type == PARAMETER_DATA_INTEGER) then
      if (Uninitialized(parameter%ivalue)) then
        option%io_buffer = 'No integer value for PARAMETER "' // &
          trim(parameter%name) // '".'
        call PrintErrMsg(option)
      endif
    else if (Uninitialized(parameter%value)) then
      option%io_buffer = 'No scalar value for PARAMETER "' // &
        trim(parameter%name) // '".'
      call PrintErrMsg(option)
    endif
  else if (parameter%itype == PARAMETER_DATASET) then
    if (.not.associated(parameter%dataset) .and. &
        len_trim(parameter%dataset_name) == 0) then
      option%io_buffer = 'No DATASET name or LIST for PARAMETER "' // &
        trim(parameter%name) // '".'
      call PrintErrMsg(option)
    endif
  else if (parameter%itype == PARAMETER_VARIABLE) then
    if (len_trim(parameter%source_variable_name) == 0) then
      option%io_buffer = 'No SOURCE_VARIABLE name for PARAMETER "' // &
        trim(parameter%name) // '".'
      call PrintErrMsg(option)
    endif
  endif

end subroutine ParameterRead

! ************************************************************************** !

subroutine ParameterReadTargetPMC(parameter,input,option,error_str)
  !
  ! Reads FLOW / TRANSPORT / GEOMECHANICS / GEOPHYSICS after BEFORE_STEP
  ! or AFTER_STEP.
  !
  ! Author: Glenn Hammond
  ! Date: 09/21/26
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(parameter_type) :: parameter
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=*) :: error_str

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option, &
         'FLOW, TRANSPORT, GEOMECHANICS, or GEOPHYSICS', &
         trim(error_str)//',WHEN_TO_UPDATE')
  call StringToUpper(word)
  select case(trim(word))
    case('FLOW')
      parameter%target_pmc = TARGET_PMC_FLOW
    case('TRANSPORT')
      parameter%target_pmc = TARGET_PMC_TRANSPORT
    case('GEOMECHANICS')
      parameter%target_pmc = TARGET_PMC_GEOMECH
    case('GEOPHYSICS')
      parameter%target_pmc = TARGET_PMC_GEOPHYSICS
    case default
      call InputKeywordUnrecognized(input,word, &
             trim(error_str)//',WHEN_TO_UPDATE',option)
  end select

end subroutine ParameterReadTargetPMC

! ************************************************************************** !

subroutine ParameterOrderUpdates(parameter_list,option)
  !
  ! Assigns update_order for scheduled PARAMETERs so that a SOURCE_VARIABLE
  ! is applied after another PARAMETER that TARGET_VARIABLEs the same
  ! field in the same WHEN_TO_UPDATE slot. Errors on cycles and duplicate
  ! targets. Independent parameters keep input order.
  !
  ! Author: Glenn Hammond
  ! Date: 09/21/26
  !
  use Option_module
  use String_module

  implicit none

  type(parameter_type), pointer :: parameter_list
  type(option_type) :: option

  type :: plist_type
    type(parameter_type), pointer :: p
  end type plist_type

  type(plist_type), allocatable :: sched(:)
  type(parameter_type), pointer :: cur
  PetscBool, allocatable :: adj(:,:)
  PetscBool, allocatable :: done(:)
  PetscInt, allocatable :: indeg(:)
  PetscInt, allocatable :: nprod(:)
  PetscInt :: n
  PetscInt :: i, j, k
  PetscInt :: pick
  PetscInt :: seq
  PetscInt :: nleft
  PetscBool :: reordered
  PetscBool :: overlap
  character(len=MAXWORDLENGTH), allocatable :: prodname(:,:)
  character(len=MAXSTRINGLENGTH) :: string

  n = 0
  cur => parameter_list
  do
    if (.not.associated(cur)) exit
    if (cur%when_to_update /= UPDATE_NEVER) n = n + 1
    cur => cur%next
  enddo
  if (n == 0) return

  allocate(sched(n))
  allocate(adj(n,n))
  allocate(done(n))
  allocate(indeg(n))
  allocate(nprod(n))
  allocate(prodname(4,n))
  adj = PETSC_FALSE
  done = PETSC_FALSE
  indeg = 0
  nprod = 0
  prodname = ''

  i = 0
  cur => parameter_list
  do
    if (.not.associated(cur)) exit
    if (cur%when_to_update /= UPDATE_NEVER) then
      i = i + 1
      sched(i)%p => cur
      nprod(i) = 0
      if (len_trim(cur%target_variable_name) > 0) then
        nprod(i) = 1
        prodname(1,i) = cur%target_variable_name
        if (StringCompare(cur%target_variable_name,'POROSITY')) then
          nprod(i) = 3
          prodname(2,i) = 'BASE_POROSITY'
          prodname(3,i) = 'INITIAL_POROSITY'
        endif
      endif
    endif
    cur => cur%next
  enddo

  do i = 1, n
    do j = 1, n
      if (i == j) cycle
      if (sched(i)%p%when_to_update /= &
          sched(j)%p%when_to_update) cycle
      select case(sched(i)%p%when_to_update)
        case(UPDATE_BEFORE_STEP,UPDATE_AFTER_STEP)
          if (sched(i)%p%target_pmc /= sched(j)%p%target_pmc) cycle
      end select
      if (i < j) then
        overlap = PETSC_FALSE
        do k = 1, nprod(i)
          if (StringCompare(prodname(k,i),prodname(1,j)) .or. &
              (nprod(j) >= 2 .and. &
               StringCompare(prodname(k,i),prodname(2,j))) .or. &
              (nprod(j) >= 3 .and. &
               StringCompare(prodname(k,i),prodname(3,j)))) then
            overlap = PETSC_TRUE
            exit
          endif
        enddo
        if (overlap) then
          option%io_buffer = 'PARAMETERs "' // &
            trim(sched(i)%p%name) // '" and "' // &
            trim(sched(j)%p%name) // '" TARGET the same variable in &
            &the same WHEN_TO_UPDATE slot.'
          call PrintErrMsg(option)
        endif
      endif
      if (len_trim(sched(j)%p%source_variable_name) > 0) then
        do k = 1, nprod(i)
          if (StringCompare(sched(j)%p%source_variable_name, &
                            prodname(k,i))) then
            if (.not.adj(i,j)) then
              adj(i,j) = PETSC_TRUE
              indeg(j) = indeg(j) + 1
            endif
            exit
          endif
        enddo
      endif
    enddo
  enddo

  seq = 0
  nleft = n
  do
    if (nleft == 0) exit
    pick = 0
    do i = 1, n
      if (.not.done(i) .and. indeg(i) == 0) then
        pick = i
        exit
      endif
    enddo
    if (pick == 0) then
      string = ''
      do i = 1, n
        if (.not.done(i)) then
          string = trim(string) // ' "' // trim(sched(i)%p%name) // '"'
        endif
      enddo
      option%io_buffer = 'Circular PARAMETER SOURCE_VARIABLE / &
        &TARGET_VARIABLE dependence in the same WHEN_TO_UPDATE slot:' &
        // trim(string) // '.'
      call PrintErrMsg(option)
    endif
    seq = seq + 1
    sched(pick)%p%update_order = seq
    done(pick) = PETSC_TRUE
    nleft = nleft - 1
    do j = 1, n
      if (adj(pick,j)) indeg(j) = indeg(j) - 1
    enddo
  enddo

  reordered = PETSC_FALSE
  do i = 1, n
    if (sched(i)%p%update_order /= i) reordered = PETSC_TRUE
  enddo
  if (reordered) then
    option%io_buffer = 'PARAMETER update order adjusted so that a &
      &SOURCE_VARIABLE is applied after the PARAMETER that &
      &TARGET_VARIABLEs the same field in the same WHEN_TO_UPDATE slot.'
    call PrintMsg(option)
  endif

  deallocate(sched)
  deallocate(adj)
  deallocate(done)
  deallocate(indeg)
  deallocate(nprod)
  deallocate(prodname)

end subroutine ParameterOrderUpdates

! ************************************************************************** !

subroutine ParameterAddToList(parameter,list)
  !
  ! Adds a parameter to linked list
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  type(parameter_type), pointer :: parameter
  type(parameter_type), pointer :: list

  type(parameter_type), pointer :: cur_parameter

  if (associated(list)) then
    cur_parameter => list
    do
      if (.not.associated(cur_parameter%next)) exit
      cur_parameter => cur_parameter%next
    enddo
    cur_parameter%next => parameter
  else
    list => parameter
  endif

end subroutine ParameterAddToList

! ************************************************************************** !

subroutine ParameterSetup(parameter_list,option)
  !
  ! Assigns ids and registers parameter names on option
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/23

  use Option_module
  use Option_Parameter_module

  type(parameter_type), pointer :: parameter_list
  type(option_type) :: option

  type(parameter_option_type), pointer :: parameter_option
  type(parameter_type), pointer :: cur_parameter
  PetscInt :: iparameter

  parameter_option => option%parameter

  iparameter = 0
  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    iparameter = iparameter + 1
    cur_parameter%id = iparameter
    cur_parameter => cur_parameter%next
  enddo
  parameter_option%num_parameters = iparameter
  allocate(parameter_option%parameter_names(iparameter))
  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    parameter_option%parameter_names(cur_parameter%id) = cur_parameter%name
    cur_parameter => cur_parameter%next
  enddo

end subroutine ParameterSetup

! ************************************************************************** !

subroutine ParameterLinkDatasets(parameter_list,datasets,option)
  !
  ! Binds named DATASET parameters to realization dataset list entries.
  !
  ! Author: Glenn Hammond
  ! Date: 07/10/26

  use Option_module

  type(parameter_type), pointer :: parameter_list
  class(dataset_base_type), pointer :: datasets
  type(option_type) :: option

  type(parameter_type), pointer :: cur_parameter
  character(len=MAXSTRINGLENGTH) :: err_string

  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    if (cur_parameter%itype == PARAMETER_DATASET .and. &
        .not.associated(cur_parameter%dataset) .and. &
        len_trim(cur_parameter%dataset_name) > 0) then
      err_string = 'PARAMETER "' // trim(cur_parameter%name) // '"'
      cur_parameter%dataset => &
        DatasetBaseGetPointer(datasets,cur_parameter%dataset_name, &
                              err_string,option)
      cur_parameter%dataset_is_owned = PETSC_FALSE
    endif
    cur_parameter => cur_parameter%next
  enddo

end subroutine ParameterLinkDatasets

! ************************************************************************** !

subroutine ParameterMapVariables(parameter,option)
  !
  ! Maps SOURCE_VARIABLE / TARGET_VARIABLE names to ivar/isubvar and rejects
  ! targets that PatchSetVariable cannot write.
  !
  ! Author: Glenn Hammond
  ! Date: 07/10/26

  use Option_module
  use Output_Aux_module
  use Variables_module

  type(parameter_type) :: parameter
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: units
  character(len=MAXWORDLENGTH) :: output_name
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscInt :: category, ivar, isubvar, isubsubvar

  err_string = 'PARAMETER "' // trim(parameter%name) // '"'

  if (len_trim(parameter%source_variable_name) > 0) then
    call OutputVariableToID(parameter%source_variable_name, &
                            output_name,units,category, &
                            ivar,isubvar,isubsubvar,option)
    parameter%source_ivar = ivar
    parameter%source_isubvar = isubvar
  endif

  if (len_trim(parameter%target_variable_name) > 0) then
    call OutputVariableToID(parameter%target_variable_name, &
                            output_name,units,category, &
                            ivar,isubvar,isubsubvar,option)
    select case(ivar)
      case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y,PERMEABILITY_Z, &
           PERMEABILITY_XY,PERMEABILITY_XZ,PERMEABILITY_YZ, &
           GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
           GAS_PERMEABILITY_Z)
        option%io_buffer = trim(err_string) // &
          ': TARGET_VARIABLE "' // trim(parameter%target_variable_name) // &
          '" (permeability) cannot be set through PatchSetVariable. Use a &
          &MATERIAL_PROPERTY permeability dataset or a dedicated update path.'
        call PrintErrMsg(option)
      case(MATERIAL_ID)
        option%io_buffer = trim(err_string) // &
          ': TARGET_VARIABLE Material_ID cannot be set through &
          &PatchSetVariable.'
        call PrintErrMsg(option)
      case(VOLUME,TORTUOSITY,TORTUOSITY_Y,TORTUOSITY_Z,SOIL_COMPRESSIBILITY, &
           SOIL_REFERENCE_PRESSURE)
        option%io_buffer = trim(err_string) // &
          ': TARGET_VARIABLE "' // trim(parameter%target_variable_name) // &
          '" cannot be set through PatchSetVariable.'
        call PrintErrMsg(option)
      case default
        ! remaining targets: PatchSetVariable will error if unsupported
    end select
    parameter%target_ivar = ivar
    parameter%target_isubvar = isubvar
  endif

end subroutine ParameterMapVariables

! ************************************************************************** !

subroutine ParameterMapVariableList(parameter_list,option)
  !
  ! Maps SOURCE/TARGET variable names for all parameters in the list.
  !
  ! Author: Glenn Hammond
  ! Date: 07/10/26

  use Option_module

  type(parameter_type), pointer :: parameter_list
  type(option_type), pointer :: option

  type(parameter_type), pointer :: cur_parameter

  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    call ParameterMapVariables(cur_parameter,option)
    cur_parameter => cur_parameter%next
  enddo

end subroutine ParameterMapVariableList

! ************************************************************************** !

subroutine ParameterQualityCheck(parameter_list,option)
  !
  ! Verifies parameter list size matches option registration
  !
  ! Author: Glenn Hammond
  ! Date: 01/05/23

  use Option_module
  use String_module

  type(parameter_type), pointer :: parameter_list
  type(option_type) :: option

  type(parameter_type), pointer :: cur_parameter
  PetscInt :: iparameter

  iparameter = 0
  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    iparameter = iparameter + 1
    cur_parameter => cur_parameter%next
  enddo
  if (iparameter /= option%parameter%num_parameters) then
    option%io_buffer = &
      StringWrite(iparameter-option%parameter%num_parameters) // &
      ' parameters were added to the parameter list after it was "setup". &
      &Please send your input deck to pflotran-dev@googlegroups.com'
    call PrintErrMsg(option)
  endif

end subroutine ParameterQualityCheck

! ************************************************************************** !

recursive subroutine ParameterDestroy(parameter)
  !
  ! Destroys a parameter
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  use Dataset_module, only : DatasetDestroy

  type(parameter_type), pointer :: parameter

  class(dataset_base_type), pointer :: dataset

  if (.not.associated(parameter)) return

  call ParameterDestroy(parameter%next)

  if (parameter%dataset_is_owned .and. associated(parameter%dataset)) then
    dataset => parameter%dataset
    call DatasetDestroy(dataset)
  endif
  nullify(parameter%dataset)

  deallocate(parameter)
  nullify(parameter)

end subroutine ParameterDestroy

end module Parameter_module
