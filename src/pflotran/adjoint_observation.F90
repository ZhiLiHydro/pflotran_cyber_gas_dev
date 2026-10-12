module Adjoint_Observation_module
  !
  ! The ADJOINT OBSERVATIONS block: point observations at any time of the
  ! run, and their misfit
  !   phi = 1/2 sum_i ((H_i - d_i)/sigma_i)^2,
  ! H_i the simulated value at observation i's cell and time. The forward
  ! run records H_i and dH_i/dY at the end of the waypoint interval ending
  ! at t_i (observation times are waypoints); the backward sweep adds the
  ! kick (H_i - d_i)/sigma_i^2 dH_i/dY to lambda at that step.
  !
  !   OBSERVATIONS
  !     TIME_SERIES
  !       CELL_ID <natural id> | COORDINATE <x> <y> <z>
  !       OBSERVED_VARIABLE FREE_ION <primary species> |
  !                         IMMOBILE <immobile species>
  !       ERROR ABSOLUTE <a> RELATIVE <r> | STANDARD_DEVIATION <sigma>
  !       FILE <path>   ! TIME_UNITS and DATA_UNITS lines, then: time value
  !     /
  !     MEASUREMENT
  !       TIME <t> <units>
  !       CELL_ID ... | COORDINATE ...
  !       OBSERVED_VARIABLE ...
  !       VALUE <d> <units>
  !       ERROR ... | STANDARD_DEVIATION ...
  !     /
  !   /
  !
  ! Nothing has a default. sigma = a + r|d|, a and STANDARD_DEVIATION in
  ! the data units (mol/kg aqueous, mol/m^3 bulk immobile). Times are after
  ! the start and at most the final time.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
#include "petsc/finclude/petscvec.h"
  use petscvec

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use PFLOTRAN_Constants_module
  use Geometry_module
  use Option_module

  implicit none

  private

  PetscInt, parameter :: ADJ_OBS_FREE_ION = 1
  PetscInt, parameter :: ADJ_OBS_IMMOBILE = 2

  type, public :: adjoint_observations_type
    PetscInt :: n = 0
    ! per observation: the deck's
    PetscReal, allocatable :: time(:)
    PetscReal, allocatable :: value(:)
    PetscReal, allocatable :: sigma(:)
    PetscInt, allocatable :: variable(:)
    character(len=MAXWORDLENGTH), allocatable :: species_name(:)
    ! CELL_ID, or 0 for a COORDINATE
    PetscInt, allocatable :: natural_id(:)
    type(point3d_type), allocatable :: coordinate(:)
    ! the first observation of its TIME_SERIES or MEASUREMENT
    PetscInt, allocatable :: first(:)
    ! resolved at Setup: the entry in a cell's block of the transport state,
    ! the local cell (0 on other processes) and the cell's natural id
    PetscInt, allocatable :: species_id(:)
    PetscInt, allocatable :: local_id(:)
    PetscInt, allocatable :: cell(:)
    ! of the last forward run: the step at time(i) (-1: not reached), H_i,
    ! and dH_i/dY on the owning process
    PetscInt, allocatable :: step(:)
    PetscReal, allocatable :: simulated(:)
    PetscReal, allocatable :: dH_dY(:)
  end type adjoint_observations_type

  type(adjoint_observations_type), pointer, public :: adj_observations => &
    null()

  public :: AdjObsRead, &
            AdjObsCount, &
            AdjObsSetup, &
            AdjObsTimes, &
            AdjObsReset, &
            AdjObsRecord, &
            AdjObsKick, &
            AdjObsObjective, &
            AdjObsWrite, &
            AdjObsDestroy

contains

! ************************************************************************** !

function AdjObsCount()
  !
  ! The number of observations read.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  PetscInt :: AdjObsCount

  AdjObsCount = 0
  if (associated(adj_observations)) AdjObsCount = adj_observations%n

end function AdjObsCount

! ************************************************************************** !

subroutine AdjObsRead(input,option)
  !
  ! Reads the OBSERVATIONS block (see the module header).
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  use Input_Aux_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SUBSURFACE,ADJOINT,OBSERVATIONS'
  if (associated(adj_observations)) then
    call PrintErrMsg(option,'ADJOINT takes one OBSERVATIONS block.')
  endif
  allocate(adj_observations)
  allocate(adj_observations%time(0),adj_observations%value(0), &
           adj_observations%sigma(0),adj_observations%variable(0), &
           adj_observations%species_name(0),adj_observations%natural_id(0), &
           adj_observations%coordinate(0),adj_observations%first(0))

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(trim(keyword))
      case('TIME_SERIES','MEASUREMENT')
        call AdjObsReadEntry(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (adj_observations%n == 0) then
    call PrintErrMsg(option,trim(error_string) // ' has no TIME_SERIES or &
      &MEASUREMENT.')
  endif

end subroutine AdjObsRead

! ************************************************************************** !

subroutine AdjObsReadEntry(input,option,kind,block_string)
  !
  ! One TIME_SERIES or MEASUREMENT: every card is required.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  use Input_Aux_module
  use String_module
  use Units_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=*) :: kind
  character(len=*) :: block_string

  character(len=MAXWORDLENGTH) :: keyword, word, species, time_units
  character(len=MAXWORDLENGTH) :: data_units, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string, filename
  type(point3d_type) :: coordinate
  PetscReal :: time, value, sigma_abs, sigma_rel
  PetscReal :: time_conversion, data_conversion
  PetscReal, allocatable :: times(:), values(:), sigmas(:)
  PetscInt :: natural_id, variable, i
  PetscBool :: has_location, has_error, has_time, has_value, has_file

  error_string = trim(block_string) // ',' // trim(kind)
  natural_id = 0
  coordinate%x = UNINITIALIZED_DOUBLE
  coordinate%y = UNINITIALIZED_DOUBLE
  coordinate%z = UNINITIALIZED_DOUBLE
  variable = 0
  species = ''
  sigma_abs = UNINITIALIZED_DOUBLE
  sigma_rel = UNINITIALIZED_DOUBLE
  time_units = ''
  data_units = ''
  filename = ''
  has_location = PETSC_FALSE
  has_error = PETSC_FALSE
  has_time = PETSC_FALSE
  has_value = PETSC_FALSE
  has_file = PETSC_FALSE

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(trim(keyword))
      case('CELL_ID')
        if (has_location) call AdjObsTwice(option,error_string,'location')
        call InputReadInt(input,option,natural_id)
        call InputErrorMsg(input,option,keyword,error_string)
        if (natural_id < 1) then
          call PrintErrMsg(option,trim(error_string) // ',CELL_ID must be &
            &positive.')
        endif
        has_location = PETSC_TRUE
      case('COORDINATE')
        if (has_location) call AdjObsTwice(option,error_string,'location')
        call GeometryReadCoordinate(input,option,coordinate,error_string)
        has_location = PETSC_TRUE
      case('REGION')
        call PrintErrMsg(option,trim(error_string) // ': an observation is &
          &at one cell (CELL_ID or COORDINATE); a REGION would need an &
          &operator such as a volume average, which is not supported.')
      case('OBSERVED_VARIABLE')
        if (variable /= 0) call AdjObsTwice(option,error_string,keyword)
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(trim(word))
          case('FREE_ION')
            variable = ADJ_OBS_FREE_ION
          case('IMMOBILE')
            variable = ADJ_OBS_IMMOBILE
          case('TOTAL')
            call PrintErrMsg(option,trim(error_string) // ',OBSERVED_VARIABLE &
              &TOTAL is not supported yet; use FREE_ION (equal to the total &
              &for a species without aqueous complexes).')
          case default
            call InputKeywordUnrecognized(input,word, &
              trim(error_string) // ',OBSERVED_VARIABLE',option)
        end select
        call InputReadWord(input,option,species,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name', &
                           trim(error_string) // ',OBSERVED_VARIABLE')
      case('ERROR')
        if (has_error) call AdjObsTwice(option,error_string,'error')
        do i = 1, 2
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'ABSOLUTE and RELATIVE', &
                             trim(error_string) // ',ERROR')
          call StringToUpper(word)
          select case(trim(word))
            case('ABSOLUTE')
              call InputReadDouble(input,option,sigma_abs)
            case('RELATIVE')
              call InputReadDouble(input,option,sigma_rel)
            case default
              call InputKeywordUnrecognized(input,word, &
                trim(error_string) // ',ERROR',option)
          end select
          call InputErrorMsg(input,option,word,trim(error_string) // ',ERROR')
        enddo
        if (Uninitialized(sigma_abs) .or. Uninitialized(sigma_rel)) then
          call PrintErrMsg(option,trim(error_string) // ',ERROR requires &
            &ABSOLUTE and RELATIVE.')
        endif
        if (.not.(ieee_is_finite(sigma_abs) .and. ieee_is_finite(sigma_rel) &
                  .and. sigma_abs >= 0.d0 .and. sigma_rel >= 0.d0)) then
          call PrintErrMsg(option,trim(error_string) // ',ERROR must be &
            &finite and not negative.')
        endif
        has_error = PETSC_TRUE
      case('STANDARD_DEVIATION')
        if (has_error) call AdjObsTwice(option,error_string,'error')
        call InputReadDouble(input,option,sigma_abs)
        call InputErrorMsg(input,option,keyword,error_string)
        if (.not.(ieee_is_finite(sigma_abs) .and. sigma_abs > 0.d0)) then
          call PrintErrMsg(option,trim(error_string) // ',STANDARD_DEVIATION &
            &must be positive and finite.')
        endif
        sigma_rel = 0.d0
        has_error = PETSC_TRUE
      case('TIME','VALUE')
        if (trim(kind) /= 'MEASUREMENT') &
          call InputKeywordUnrecognized(input,keyword,error_string,option)
        if (trim(keyword) == 'TIME') then
          if (has_time) call AdjObsTwice(option,error_string,keyword)
          call InputReadDouble(input,option,time)
          call InputErrorMsg(input,option,keyword,error_string)
          call InputReadWord(input,option,time_units,PETSC_TRUE)
          call InputErrorMsg(input,option,'TIME units',error_string)
          has_time = PETSC_TRUE
        else
          if (has_value) call AdjObsTwice(option,error_string,keyword)
          call InputReadDouble(input,option,value)
          call InputErrorMsg(input,option,keyword,error_string)
          call InputReadWord(input,option,data_units,PETSC_TRUE)
          call InputErrorMsg(input,option,'VALUE units',error_string)
          has_value = PETSC_TRUE
        endif
      case('FILE')
        if (trim(kind) /= 'TIME_SERIES') &
          call InputKeywordUnrecognized(input,keyword,error_string,option)
        if (has_file) call AdjObsTwice(option,error_string,keyword)
        call InputReadFilename(input,option,filename)
        call InputErrorMsg(input,option,keyword,error_string)
        has_file = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.has_location) &
    call AdjObsMissing(option,error_string,'CELL_ID or COORDINATE')
  if (variable == 0) call AdjObsMissing(option,error_string,'OBSERVED_VARIABLE')
  if (.not.has_error) &
    call AdjObsMissing(option,error_string,'ERROR or STANDARD_DEVIATION')
  if (trim(kind) == 'MEASUREMENT') then
    if (.not.has_time) call AdjObsMissing(option,error_string,'TIME')
    if (.not.has_value) call AdjObsMissing(option,error_string,'VALUE')
    allocate(times(1),values(1))
    times(1) = time
    values(1) = value
  else
    if (.not.has_file) call AdjObsMissing(option,error_string,'FILE')
    call AdjObsReadFile(filename,option,error_string,time_units,data_units, &
                        times,values)
  endif

  internal_units = 'sec'
  time_conversion = UnitsConvertToInternal(time_units,internal_units, &
                                           trim(error_string) // &
                                           ',time units',option)
  internal_units = 'mol/kg'
  if (variable == ADJ_OBS_IMMOBILE) internal_units = 'mol/m^3'
  data_conversion = UnitsConvertToInternal(data_units,internal_units, &
                                           trim(error_string) // &
                                           ',data units',option)
  allocate(sigmas(size(times)))
  do i = 1, size(times)
    if (.not.ieee_is_finite(values(i))) then
      write(option%io_buffer,'(a,es12.5,a)') trim(error_string) // &
        ': the value at t = ',times(i),' ' // trim(time_units) // &
        ' is not finite.'
      call PrintErrMsg(option)
    endif
    sigmas(i) = sigma_abs*data_conversion + &
                sigma_rel*abs(values(i)*data_conversion)
    if (.not.(sigmas(i) > 0.d0)) then
      write(option%io_buffer,'(a,es12.5,a)') trim(error_string) // &
        ': the error of the observation at t = ',times(i),' ' // &
        trim(time_units) // ' is zero (a zero value needs ERROR ABSOLUTE &
        &> 0).'
      call PrintErrMsg(option)
    endif
  enddo
  call AdjObsAppend(adj_observations,times*time_conversion, &
                    values*data_conversion,sigmas,variable,species, &
                    natural_id,coordinate)

end subroutine AdjObsReadEntry

! ************************************************************************** !

subroutine AdjObsReadFile(filename,option,error_string,time_units, &
                          data_units,times,values)
  !
  ! A TIME_SERIES file: TIME_UNITS and DATA_UNITS lines (both required),
  ! then rows of time and value with increasing times.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  use Input_Aux_module
  use String_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  character(len=*) :: error_string
  character(len=MAXWORDLENGTH) :: time_units
  character(len=MAXWORDLENGTH) :: data_units
  PetscReal, allocatable :: times(:)
  PetscReal, allocatable :: values(:)

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: file_string, line
  character(len=MAXWORDLENGTH) :: word
  PetscReal, allocatable :: grown(:)
  PetscReal :: time, value
  PetscInt :: n
  PetscErrorCode :: ierr

  file_string = trim(error_string) // ',FILE ' // trim(filename)
  allocate(times(16),values(16))
  n = 0
  input => InputCreate(IUNIT_TEMP,filename,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    ! the first word, without consuming it
    line = input%buf
    ierr = INPUT_ERROR_NONE
    call InputReadWord(line,word,PETSC_TRUE,ierr)
    call StringToUpper(word)
    select case(trim(word))
      case('TIME_UNITS','DATA_UNITS')
        if (n > 0) call PrintErrMsg(option,trim(file_string) // ': ' // &
          trim(word) // ' must come before the data.')
        call InputReadCard(input,option,word)
        call StringToUpper(word)
        if (trim(word) == 'TIME_UNITS') then
          call InputReadWord(input,option,time_units,PETSC_TRUE)
        else
          call InputReadWord(input,option,data_units,PETSC_TRUE)
        endif
        call InputErrorMsg(input,option,word,file_string)
      case default
        call InputReadDouble(input,option,time)
        call InputErrorMsg(input,option,'time',file_string)
        call InputReadDouble(input,option,value)
        call InputErrorMsg(input,option,'value',file_string)
        if (len_trim(input%buf) > 0) then
          call PrintErrMsg(option,trim(file_string) // ': a row is a time &
            &and one value; "' // trim(adjustl(input%buf)) // '" follows.')
        endif
        n = n + 1
        if (n > size(times)) then
          allocate(grown(2*size(times)))
          grown(1:n-1) = times(1:n-1)
          call move_alloc(grown,times)
          allocate(grown(2*size(values)))
          grown(1:n-1) = values(1:n-1)
          call move_alloc(grown,values)
        endif
        times(n) = time
        values(n) = value
        if (n > 1) then
          if (.not.(times(n) > times(n-1))) then
            call PrintErrMsg(option,trim(file_string) // ': times must &
              &increase.')
          endif
        endif
    end select
  enddo
  call InputDestroy(input)

  if (len_trim(time_units) == 0) call AdjObsMissing(option,file_string, &
                                                    'TIME_UNITS')
  if (len_trim(data_units) == 0) call AdjObsMissing(option,file_string, &
                                                    'DATA_UNITS')
  if (n == 0) call PrintErrMsg(option,trim(file_string) // ' has no data.')
  times = times(1:n)
  values = values(1:n)

end subroutine AdjObsReadFile

! ************************************************************************** !

subroutine AdjObsMissing(option,error_string,card)
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  type(option_type) :: option
  character(len=*) :: error_string
  character(len=*) :: card

  call PrintErrMsg(option,trim(error_string) // ' requires ' // trim(card) // &
    ' (observations have no defaults).')

end subroutine AdjObsMissing

! ************************************************************************** !

subroutine AdjObsTwice(option,error_string,what)
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  type(option_type) :: option
  character(len=*) :: error_string
  character(len=*) :: what

  call PrintErrMsg(option,trim(error_string) // ' gives the ' // trim(what) // &
    ' twice.')

end subroutine AdjObsTwice

! ************************************************************************** !

subroutine AdjObsAppend(obs,times,values,sigmas,variable,species,natural_id, &
                        coordinate)
  !
  ! The observations of one TIME_SERIES or MEASUREMENT.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  type(adjoint_observations_type) :: obs
  PetscReal :: times(:)
  PetscReal :: values(:)
  PetscReal :: sigmas(:)
  PetscInt :: variable
  character(len=MAXWORDLENGTH) :: species
  PetscInt :: natural_id
  type(point3d_type) :: coordinate

  PetscInt :: k, n

  n = size(times)
  obs%time = [obs%time, times]
  obs%value = [obs%value, values]
  obs%sigma = [obs%sigma, sigmas]
  obs%variable = [obs%variable, (variable, k = 1, n)]
  obs%species_name = [character(len=MAXWORDLENGTH) :: obs%species_name, &
                      (species, k = 1, n)]
  obs%natural_id = [obs%natural_id, (natural_id, k = 1, n)]
  obs%coordinate = [obs%coordinate, (coordinate, k = 1, n)]
  obs%first = [obs%first, (obs%n+1, k = 1, n)]
  obs%n = obs%n + n

end subroutine AdjObsAppend

! ************************************************************************** !

subroutine AdjObsSetup(realization,final_time,option)
  !
  ! Resolves each observation's species and cell, and checks its time.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  use Realization_Subsurface_class
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module, only : ReactionImGetSpeciesIDFromName
  use Grid_module
  use Adjoint_Cost_module, only : AdjCostLocateCell

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: final_time
  type(option_type) :: option

  type(adjoint_observations_type), pointer :: obs
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: i, j, id, local_id, cell
  PetscErrorCode :: ierr

  obs => adj_observations
  grid => realization%patch%grid
  allocate(obs%species_id(obs%n),obs%local_id(obs%n),obs%cell(obs%n))
  allocate(obs%step(obs%n),obs%simulated(obs%n),obs%dH_dY(obs%n))
  obs%local_id = 0
  obs%dH_dY = 0.d0
  call AdjObsReset()

  do i = 1, obs%n
    write(error_string,'(a,i0,a,es12.5,a)') 'ADJOINT,OBSERVATIONS &
      &observation ',i,' (at ',obs%time(i),' s)'
    if (.not.(obs%time(i) > 0.d0)) then
      call PrintErrMsg(option,trim(error_string) // ': its time is not &
        &after the start.')
    else if (.not.(obs%time(i) <= &
                   final_time*(1.d0+64.d0*epsilon(1.d0)))) then
      call PrintErrMsg(option,trim(error_string) // ': its time is after &
        &the final time.')
    endif
    ! the rows of a TIME_SERIES share the first's species and cell
    j = obs%first(i)
    if (j < i) then
      obs%species_id(i) = obs%species_id(j)
      obs%local_id(i) = obs%local_id(j)
      obs%cell(i) = obs%cell(j)
      cycle
    endif
    word = obs%species_name(i)
    select case(obs%variable(i))
      case(ADJ_OBS_FREE_ION)
        id = ReactionAuxGetPriSpecIDFromName(word,realization%reaction, &
                                             PETSC_FALSE,option)
        if (id <= 0) then
          call PrintErrMsg(option,trim(error_string) // ': FREE_ION "' // &
            trim(word) // '" is not a primary species.')
        endif
      case(ADJ_OBS_IMMOBILE)
        id = ReactionImGetSpeciesIDFromName(word,realization%reaction% &
                                            immobile,PETSC_FALSE,option)
        if (id <= 0) then
          call PrintErrMsg(option,trim(error_string) // ': IMMOBILE "' // &
            trim(word) // '" is not an immobile species.')
        endif
        id = realization%reaction%offset_immobile + id
    end select
    obs%species_id(i) = id

    ! the cell: on one process, its natural id on all
    call AdjCostLocateCell(grid,obs%natural_id(i),obs%coordinate(i), &
                           error_string,option,local_id)
    cell = 0
    if (local_id > 0) then
      if (realization%patch%imat(grid%nL2G(local_id)) > 0) &
        cell = grid%nG2A(grid%nL2G(local_id))
    endif
    call MPI_Allreduce(MPI_IN_PLACE,cell,ONE_INTEGER_MPI,MPIU_INTEGER, &
                       MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
    if (cell == 0) then
      call PrintErrMsg(option,trim(error_string) // ': its cell is &
        &inactive.')
    endif
    obs%local_id(i) = local_id
    obs%cell(i) = cell
  enddo

end subroutine AdjObsSetup

! ************************************************************************** !

subroutine AdjObsTimes(times)
  !
  ! The observation times (repeated ones too), for the waypoints of the
  ! run.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  PetscReal, allocatable :: times(:)

  if (.not.associated(adj_observations)) then
    allocate(times(0))
    return
  endif
  times = adj_observations%time(1:adj_observations%n)

end subroutine AdjObsTimes

! ************************************************************************** !

subroutine AdjObsReset()
  !
  ! Before a forward run: nothing is observed yet.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  if (.not.associated(adj_observations)) return
  if (.not.allocated(adj_observations%step)) return
  adj_observations%step = -1
  adj_observations%simulated = 0.d0

end subroutine AdjObsReset

! ************************************************************************** !

subroutine AdjObsRecord(Y,ndof,use_log,time,step,option,recorded)
  !
  ! At the end of a waypoint interval, state Y (transport, ndof entries per
  ! local cell): the observations at this time take H_i and dH_i/dY.
  ! recorded: any did.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  Vec :: Y
  PetscInt :: ndof
  PetscBool :: use_log
  PetscReal :: time
  PetscInt :: step
  type(option_type) :: option
  PetscBool :: recorded

  type(adjoint_observations_type), pointer :: obs
  PetscReal, pointer :: y_p(:)
  PetscReal, allocatable :: values(:)
  PetscInt, allocatable :: list(:)
  PetscInt :: i, k, n, offset
  PetscErrorCode :: ierr

  recorded = PETSC_FALSE
  obs => adj_observations
  if (.not.associated(obs)) return
  allocate(list(obs%n),values(obs%n))
  n = 0
  do i = 1, obs%n
    if (abs(obs%time(i)-time) <= 64.d0*epsilon(1.d0)*max(1.d0,abs(time))) then
      n = n + 1
      list(n) = i
    endif
  enddo
  if (n == 0) return

  values = 0.d0
  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  do k = 1, n
    i = list(k)
    if (obs%local_id(i) == 0) cycle
    offset = (obs%local_id(i)-1)*ndof + obs%species_id(i)
    if (use_log) then
      values(k) = exp(y_p(offset))
      obs%dH_dY(i) = values(k)
    else
      values(k) = y_p(offset)
      obs%dH_dY(i) = 1.d0
    endif
  enddo
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,values,n,MPI_DOUBLE_PRECISION,MPI_SUM, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  do k = 1, n
    obs%simulated(list(k)) = values(k)
    obs%step(list(k)) = step
  enddo
  recorded = PETSC_TRUE

end subroutine AdjObsRecord

! ************************************************************************** !

subroutine AdjObsKick(lambda,ndof,step)
  !
  ! Adds the kicks of the observations at step to lambda (transport, ndof
  ! entries per local cell): dphi/dY of those observations.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  Vec :: lambda
  PetscInt :: ndof
  PetscInt :: step

  type(adjoint_observations_type), pointer :: obs
  PetscReal, pointer :: lam_p(:)
  PetscInt :: i, offset
  PetscErrorCode :: ierr

  obs => adj_observations
  if (.not.associated(obs)) return
  call VecGetArray(lambda,lam_p,ierr);CHKERRQ(ierr)
  do i = 1, obs%n
    if (obs%step(i) /= step .or. obs%local_id(i) == 0) cycle
    offset = (obs%local_id(i)-1)*ndof + obs%species_id(i)
    lam_p(offset) = lam_p(offset) + &
      (obs%simulated(i)-obs%value(i))/obs%sigma(i)**2*obs%dH_dY(i)
  enddo
  call VecRestoreArray(lambda,lam_p,ierr);CHKERRQ(ierr)

end subroutine AdjObsKick

! ************************************************************************** !

subroutine AdjObsObjective(phi,complete)
  !
  ! phi of the last forward run; complete: every observation was recorded.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  PetscReal :: phi
  PetscBool :: complete

  type(adjoint_observations_type), pointer :: obs

  phi = 0.d0
  complete = PETSC_TRUE
  obs => adj_observations
  if (.not.associated(obs)) return
  complete = all(obs%step(1:obs%n) >= 0)
  phi = 0.5d0*sum(((obs%simulated(1:obs%n)-obs%value(1:obs%n))/ &
                   obs%sigma(1:obs%n))**2)

end subroutine AdjObsObjective

! ************************************************************************** !

subroutine AdjObsWrite(option)
  !
  ! <prefix>-observations.dat: one row per observation of the last forward
  ! run.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  type(option_type) :: option

  type(adjoint_observations_type), pointer :: obs
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: variable
  PetscInt :: i, unit

  obs => adj_observations
  if (.not.associated(obs)) return
  if (.not.OptionIsIORank(option)) return
  filename = trim(option%global_prefix) // trim(option%group_prefix) // &
             '-observations.dat'
  open(newunit=unit,file=filename,action='write',status='replace')
  write(unit,'(a)') '# ADJOINT OBSERVATIONS of the last forward run: id, &
    &time [s], cell (natural id), variable, species, observed, simulated, &
    &sigma, residual (simulated - observed)/sigma; concentrations in mol/kg &
    &(FREE_ION) or mol/m^3 bulk (IMMOBILE)'
  do i = 1, obs%n
    select case(obs%variable(i))
      case(ADJ_OBS_FREE_ION)
        variable = 'FREE_ION'
      case(ADJ_OBS_IMMOBILE)
        variable = 'IMMOBILE'
    end select
    write(unit,'(i8,1x,es24.16,1x,i10,1x,a,1x,a,4(1x,es24.16))') i, &
      obs%time(i),obs%cell(i),trim(variable),trim(obs%species_name(i)), &
      obs%value(i),obs%simulated(i),obs%sigma(i), &
      (obs%simulated(i)-obs%value(i))/obs%sigma(i)
  enddo
  close(unit)

end subroutine AdjObsWrite

! ************************************************************************** !

subroutine AdjObsDestroy()
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  implicit none

  if (.not.associated(adj_observations)) return
  deallocate(adj_observations)
  nullify(adj_observations)

end subroutine AdjObsDestroy

end module Adjoint_Observation_module
