module Adjoint_Cost_module

#include "petsc/finclude/petscts.h"
  use petscts

  use Geometry_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module

  ! The ADJOINT block's cost functions of the final state, in deck order:
  !   TERMINAL_SPECIES_CONCENTRATION  J = WEIGHT * sum over the selected
  !                                       cells of C,
  ! C the final free-ion molality (aqueous) or concentration [mol/m^3
  ! bulk] (immobile) of a primary species. Cells are selected by
  ! COORDINATE, REGION or CELL_IDS, as in SrcSink_Sandbox_Base_class.
  ! Misfits of data are OBSERVATIONS (Adjoint_Observation_module).
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26

  implicit none

  private

  type, public :: adjoint_cost_type
    PetscInt :: id = UNINITIALIZED_INTEGER        ! position in list
    character(len=MAXWORDLENGTH) :: name = ''     ! label for output
    PetscReal :: weight = 1.d0
    character(len=MAXWORDLENGTH) :: species_name = ''
    PetscInt :: species_id = UNINITIALIZED_INTEGER  ! resolved at Setup
    ! cell selector
    PetscInt, pointer :: local_cell_ids(:) => null()
    PetscInt, pointer :: natural_cell_ids(:) => null()
    type(point3d_type) :: coordinate
    character(len=MAXWORDLENGTH) :: region_name = ''
    PetscReal :: J_value = 0.d0                   ! J(T) (FillTerminalLambda)
    class(realization_subsurface_type), pointer :: realization => null()
    type(adjoint_cost_type), pointer :: next => null()
  contains
    procedure, public :: Setup => AdjCostSetup
    procedure, public :: FillTerminalLambda => AdjCostFillTerminalLambda
  end type adjoint_cost_type

  type(adjoint_cost_type), pointer, public :: adj_cost_list => null()

  public :: AdjCostInit, &
            AdjCostRead, &
            AdjCostCount, &
            AdjCostLocateCell, &
            AdjCostDestroyList

contains

! ************************************************************************** !

subroutine AdjCostInit()
  !
  ! Initializes the cost list, destroying any existing one.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  if (associated(adj_cost_list)) then
    call AdjCostDestroyList()
  endif
  nullify(adj_cost_list)

end subroutine AdjCostInit

! ************************************************************************** !

subroutine AdjCostRead(input,option)
  !
  ! Reads one COST_FUNCTION <class> block and appends the cost to the list:
  !   COST_FUNCTION TERMINAL_SPECIES_CONCENTRATION
  !     SPECIES_NAME <string>              [required]
  !     COORDINATE x y z | CELL_IDS i [j ...] | REGION <name>   [one]
  !     WEIGHT <double>                    [optional, default 1]
  !   /
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: class_keyword, keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(adjoint_cost_type), pointer :: new_cost, cur_cost
  PetscBool :: has_location

  error_string = 'SUBSURFACE,ADJOINT,COST_FUNCTION'
  has_location = PETSC_FALSE

  call InputReadCard(input,option,class_keyword)
  call InputErrorMsg(input,option,'class keyword',error_string)
  call StringToUpper(class_keyword)

  allocate(new_cost)
  new_cost%coordinate%x = UNINITIALIZED_DOUBLE
  new_cost%coordinate%y = UNINITIALIZED_DOUBLE
  new_cost%coordinate%z = UNINITIALIZED_DOUBLE
  select case(trim(class_keyword))
    case('TERMINAL_SPECIES_CONCENTRATION')
    case default
      deallocate(new_cost)
      call InputKeywordUnrecognized(input,class_keyword,error_string,option)
      return
  end select
  error_string = trim(error_string) // ',' // trim(class_keyword)

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    ! one cell selector
    select case(trim(keyword))
      case('COORDINATE','CELL_IDS','REGION')
        if (has_location) then
          option%io_buffer = trim(error_string) // ' takes one of &
            &COORDINATE, CELL_IDS or REGION.'
          call PrintErrMsg(option)
        endif
        has_location = PETSC_TRUE
    end select
    select case(trim(keyword))
      case('SPECIES_NAME')
        call InputReadWord(input,option,new_cost%species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
      case('COORDINATE')
        call GeometryReadCoordinate(input,option,new_cost%coordinate, &
                                    'ADJOINT,COST_FUNCTION')
      case('CELL_IDS')
        call UtilityReadArray(new_cost%natural_cell_ids, &
                              UNINITIALIZED_INTEGER, &
                              'ADJOINT,COST_FUNCTION,' // trim(keyword), &
                              input,option)
      case('REGION')
        call InputReadWord(input,option,new_cost%region_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,'ADJOINT,COST_FUNCTION')
      case('WEIGHT')
        call InputReadDouble(input,option,new_cost%weight)
        call InputErrorMsg(input,option,keyword,'ADJOINT,COST_FUNCTION')
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (len_trim(new_cost%species_name) == 0) then
    option%io_buffer = 'SPECIES_NAME must be specified for ' // &
      trim(class_keyword) // ' cost function.'
    call PrintErrMsg(option)
  endif
  ! TODO: name holds MAXWORDLENGTH characters, so a species name longer
  ! than 19 characters stops this write ("End of record"); give name
  ! MAXSTRINGLENGTH if that matters
  write(new_cost%name,'(a,"_",a)') 'TerminalConc', &
    trim(new_cost%species_name)

  if (.not.associated(adj_cost_list)) then
    adj_cost_list => new_cost
    new_cost%id = 1
  else
    cur_cost => adj_cost_list
    do
      if (.not.associated(cur_cost%next)) exit
      cur_cost => cur_cost%next
    enddo
    cur_cost%next => new_cost
    new_cost%id = cur_cost%id + 1
  endif

end subroutine AdjCostRead

! ************************************************************************** !

subroutine AdjCostSetup(this,realization,option)
  !
  ! Resolves the species and the cell selector.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module, only : ReactionImGetSpeciesIDFromName

  implicit none

  class(adjoint_cost_type) :: this
  class(realization_subsurface_type), target :: realization
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  ! a primary species, or an immobile one (after them in the state)
  word = this%species_name
  this%species_id = &
    ReactionAuxGetPriSpecIDFromName(word,realization%reaction,PETSC_FALSE, &
                                    option)
  if (this%species_id <= 0) then
    this%species_id = &
      ReactionImGetSpeciesIDFromName(word,realization%reaction%immobile, &
                                     PETSC_FALSE,option)
    if (this%species_id <= 0) then
      option%io_buffer = 'ADJOINT COST_FUNCTION SPECIES_NAME "' // &
        trim(word) // '" is neither a primary nor an immobile species.'
      call PrintErrMsg(option)
    endif
    this%species_id = realization%reaction%offset_immobile + this%species_id
  endif

  call AdjCostResolveCells(this,realization,option)

  this%realization => realization

end subroutine AdjCostSetup

! ************************************************************************** !

subroutine AdjCostResolveCells(this,realization,option)
  !
  ! Resolves the cell selector into this%local_cell_ids.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module
  use Grid_module
  use Region_module
  use Realization_Subsurface_class

  implicit none

  class(adjoint_cost_type) :: this
  class(realization_subsurface_type) :: realization
  type(option_type) :: option

  type(grid_type), pointer :: grid
  type(region_type), pointer :: region
  PetscInt :: local_id, natural_id
  PetscInt, allocatable :: local_cell_ids(:)
  PetscInt :: icell, num_local, num_global, num_active
  PetscInt :: max_natural, min_natural
  PetscErrorCode :: ierr

  grid => realization%patch%grid

  num_local = 0
  num_global = 0
  if (Initialized(this%coordinate%x)) then
    call AdjCostLocateCell(grid,ZERO_INTEGER,this%coordinate, &
                           'ADJOINT COST_FUNCTION COORDINATE',option,local_id)
    if (local_id > 0) then
      num_local = 1
      allocate(local_cell_ids(num_local))
      local_cell_ids(num_local) = local_id
    endif
    num_global = 1
  else if (len_trim(this%region_name) > 0) then
    region => RegionGetPtrFromList(this%region_name, &
                                   realization%patch%region_list)
    if (.not.associated(region)) then
      option%io_buffer = 'Region "' // trim(this%region_name) // &
                '" in ADJOINT cost function not found in region list.'
      call PrintErrMsg(option)
    endif
    num_local = region%num_cells
    if (region%num_cells > 0) then
      allocate(this%local_cell_ids(region%num_cells))
      this%local_cell_ids(:) = region%cell_ids(:)
    endif
    ! a region has no expected count: only the no-cells check applies
  else if (associated(this%natural_cell_ids)) then
    allocate(local_cell_ids(size(this%natural_cell_ids)))
    local_cell_ids = UNINITIALIZED_INTEGER
    num_global = size(this%natural_cell_ids)
    max_natural = maxval(this%natural_cell_ids)
    min_natural = minval(this%natural_cell_ids)
    do local_id = 1, grid%nlmax
      natural_id = grid%nG2A(grid%nL2G(local_id))
      if (natural_id >= min_natural .and. natural_id <= max_natural) then
        do icell = 1, num_global
          if (natural_id == this%natural_cell_ids(icell)) then
            num_local = num_local + 1
            local_cell_ids(num_local) = local_id
            exit
          endif
        enddo
      endif
    enddo
  else
    option%io_buffer = 'ADJOINT cost function not associated with the &
      &domain through CELL_IDS, COORDINATE, or REGION.'
    call PrintErrMsg(option)
  endif

  if (allocated(local_cell_ids)) then
    if (num_local > 0) then
      allocate(this%local_cell_ids(num_local))
      this%local_cell_ids = local_cell_ids(1:num_local)
    endif
    deallocate(local_cell_ids)
  endif

  ! verify each cell maps once across ranks
  call MPI_Allreduce(num_local,icell,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)

  if (icell == 0) then
    call PrintErrMsg(option,'ADJOINT cost function: no grid cells &
      &selected.')
  else if (icell > num_global .and. num_global > 0) then
    call PrintErrMsg(option,'ADJOINT cost function: more grid cells &
      &selected than listed.')
  else if (icell < num_global) then
    call PrintErrMsg(option,'ADJOINT cost function: CELL_IDS lists a cell &
      &twice or one not in the grid.')
  endif

  ! the cost skips inactive cells
  num_active = 0
  if (associated(this%local_cell_ids)) num_active = &
    count(realization%patch%imat(grid%nL2G(this%local_cell_ids)) > 0)
  call MPI_Allreduce(MPI_IN_PLACE,num_active,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  if (num_active == 0) then
    call PrintErrMsg(option,'ADJOINT cost function: its grid cells are all &
      &inactive.')
  endif

end subroutine AdjCostResolveCells

! ************************************************************************** !

subroutine AdjCostLocateCell(grid,natural_id,coordinate,error_string, &
                             option,local_id)
  !
  ! The cell of a natural id (natural_id > 0) or else of a coordinate: its
  ! local id on the process that has it, 0 on the others. Stops unless one
  ! process has it: a coordinate on a face between two processes' cells of
  ! an unstructured grid is in both, as for a REGION's.
  !
  ! Author: Hansol Suh
  ! Date: 10/12/26
  !
  use Option_module
  use Grid_module

  implicit none

  type(grid_type) :: grid
  PetscInt :: natural_id
  type(point3d_type) :: coordinate
  character(len=*) :: error_string
  type(option_type) :: option
  PetscInt :: local_id

  PetscInt :: id, count
  PetscErrorCode :: ierr

  local_id = 0
  if (natural_id > 0) then
    do id = 1, grid%nlmax
      if (grid%nG2A(grid%nL2G(id)) == natural_id) then
        local_id = id
        exit
      endif
    enddo
  else
    ! negative on the processes that do not have it
    call GridGetLocalIDFromCoordinate(grid,coordinate,option,id)
    if (id > 0) local_id = id
  endif
  count = 0
  if (local_id > 0) count = 1
  call MPI_Allreduce(MPI_IN_PLACE,count,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  if (count == 0) then
    if (natural_id > 0) then
      call PrintErrMsg(option,trim(error_string) // ': CELL_ID is not a &
        &cell of the grid.')
    else
      call PrintErrMsg(option,trim(error_string) // ': its COORDINATE is &
        &not in the grid.')
    endif
  else if (count > 1) then
    write(option%io_buffer,'(a,i0,a)') trim(error_string) // ': its &
      &COORDINATE is in cells of ',count,' processes (on a face between &
      &them); move it into one cell.'
    call PrintErrMsg(option)
  endif

end subroutine AdjCostLocateCell

! ************************************************************************** !

subroutine AdjCostFillTerminalLambda(this,Y,lambda,ndof,option)
  !
  ! Adds lambda(T) = dJ/dY(T) at the final transport state Y (ndof entries
  ! per local cell, primary species first), nonzero only at the cost
  ! species of the selected cells, and sets J(T). Under LOG_FORMULATION Y
  ! holds ln(C), so dJ/dY = C dJ/dC.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Option_module

  implicit none

  class(adjoint_cost_type) :: this
  Vec :: Y
  Vec :: lambda
  PetscInt :: ndof
  type(option_type) :: option

  PetscErrorCode :: ierr
  PetscScalar, pointer :: lam_p(:), y_p(:)
  PetscInt :: i, local_id, ghosted_id, offset
  PetscReal :: local_J, global_J, c_val, cell_J, dJ_dC
  PetscBool :: use_log

  use_log = this%realization%reaction%use_log_formulation

  call VecGetArrayRead(Y,y_p,ierr);CHKERRQ(ierr)
  call VecGetArray(lambda,lam_p,ierr);CHKERRQ(ierr)
  local_J = 0.d0
  if (associated(this%local_cell_ids)) then
    do i = 1, size(this%local_cell_ids)
      local_id = this%local_cell_ids(i)
      ! inactive cells are not part of the forward model
      ghosted_id = this%realization%patch%grid%nL2G(local_id)
      if (this%realization%patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*ndof + this%species_id
      if (use_log) then
        c_val = exp(y_p(offset))
      else
        c_val = y_p(offset)
      endif
      cell_J = this%weight * c_val
      dJ_dC = this%weight
      if (use_log) dJ_dC = dJ_dC * c_val
      lam_p(offset) = lam_p(offset) + dJ_dC
      local_J = local_J + cell_J
    enddo
  endif
  call VecRestoreArray(lambda,lam_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(Y,y_p,ierr);CHKERRQ(ierr)

  call MPI_Allreduce(local_J,global_J,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  this%J_value = global_J

end subroutine AdjCostFillTerminalLambda

! ************************************************************************** !

function AdjCostCount()
  !
  ! Returns the number of cost functions currently in the list.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  implicit none

  PetscInt :: AdjCostCount
  type(adjoint_cost_type), pointer :: cur

  AdjCostCount = 0
  cur => adj_cost_list
  do
    if (.not.associated(cur)) exit
    AdjCostCount = AdjCostCount + 1
    cur => cur%next
  enddo

end function AdjCostCount

! ************************************************************************** !

subroutine AdjCostDestroyList()
  !
  ! Destroys every cost in the list.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Utility_module

  implicit none

  type(adjoint_cost_type), pointer :: cur, prev

  cur => adj_cost_list
  do
    if (.not.associated(cur)) exit
    prev => cur%next
    call DeallocateArray(cur%local_cell_ids)
    call DeallocateArray(cur%natural_cell_ids)
    deallocate(cur)
    cur => prev
  enddo
  nullify(adj_cost_list)

end subroutine AdjCostDestroyList

! ************************************************************************** !

end module Adjoint_Cost_module
