module Stokes_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none

  private

  !--------------------------------------------------------------------
  ! DOF indices for Stokes (4 DOFs per cell: u, v, w, p)
  !--------------------------------------------------------------------
  integer, parameter, public :: STOKES_U_DOF = 1
  integer, parameter, public :: STOKES_V_DOF = 2
  integer, parameter, public :: STOKES_W_DOF = 3
  integer, parameter, public :: STOKES_P_DOF = 4
  integer, parameter, public :: STOKES_NUM_DOF = 4

  ! Equation indices if for tracking BCs per equation
  integer, parameter, public :: STOKES_MOMENTUM_EQUATION_INDEX = 1
  integer, parameter, public :: STOKES_PRESSURE_EQUATION_INDEX = 2
  integer, parameter, public :: STOKES_MAX_EQUATION_INDEX      = 2

  ! Debugging counters
  PetscInt, public :: stokes_ts_cut_count
  PetscInt, public :: stokes_ts_count

  !--------------------------------------------------------------------
  ! Per-cell auxiliary variables for Stokes
  !--------------------------------------------------------------------
  type, public :: stokes_auxvar_type
    PetscReal :: pressure
    PetscReal :: u
    PetscReal :: v
    PetscReal :: w
    PetscReal :: du_dt(3)
  end type stokes_auxvar_type

  !--------------------------------------------------------------------
  ! Stokes parameter container (mode-level settings)
  !--------------------------------------------------------------------
  type, public :: stokes_parameter_type
    PetscBool :: dummy_flag   ! extend with real parameters as needed
  end type stokes_parameter_type

  !--------------------------------------------------------------------
  ! Top-level Stokes aux container
  !--------------------------------------------------------------------
  type, public :: stokes_type
    PetscBool :: auxvars_up_to_date
    PetscInt  :: num_aux, num_aux_bc, num_aux_ss
    type(stokes_parameter_type), pointer :: param
    type(stokes_auxvar_type),    pointer :: auxvars(:)
    type(stokes_auxvar_type),    pointer :: auxvars_bc(:)
    type(stokes_auxvar_type),    pointer :: auxvars_ss(:)
    type(matrix_zeroing_type),   pointer :: matrix_zeroing
    ! Cached viscous diagonal for Rhie-Chow pressure stabilization.
    ! Allocated to grid%ngmax in StokesSetup; recomputed each
    ! residual/Jacobian evaluation via StokesComputeViscousDiag.
    PetscReal, pointer :: aP(:)
  end type stokes_type

  !--------------------------------------------------------------------
  ! Destroy interfaces
  !--------------------------------------------------------------------
  interface StokesAuxVarDestroy
    module procedure StokesAuxVarSingleDestroy
    module procedure StokesAuxVarArray1Destroy
    module procedure StokesAuxVarArray2Destroy
  end interface StokesAuxVarDestroy

  public :: StokesAuxCreate, &
            StokesAuxDestroy, &
            StokesAuxVarInit, &
            StokesAuxVarCopy, &
            StokesAuxVarCompute, &
            StokesAuxVarDestroy, &
            StokesAuxVarStrip

contains

! ************************************************************************** !

function StokesAuxCreate(option)
  !
  ! Allocate and initialize Stokes auxiliary object
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(stokes_type), pointer :: StokesAuxCreate
  type(stokes_type), pointer :: aux

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%num_aux     = 0
  aux%num_aux_bc  = 0
  aux%num_aux_ss  = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)
  nullify(aux%aP)

  allocate(aux%param)
  aux%param%dummy_flag = PETSC_FALSE

  StokesAuxCreate => aux

end function StokesAuxCreate

! ************************************************************************** !

subroutine StokesAuxVarInit(auxvar,option)
  !
  ! Initialize a per-cell Stokes auxiliary variable
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Option_module

  implicit none

  type(stokes_auxvar_type) :: auxvar
  type(option_type)        :: option

  auxvar%pressure = 0.d0
  auxvar%u = 0.d0
  auxvar%v = 0.d0
  auxvar%w = 0.d0
  auxvar%du_dt = 0.d0

end subroutine StokesAuxVarInit

! ************************************************************************** !

subroutine StokesAuxVarCopy(auxvar_src,auxvar_dst,option)
  !
  ! Copy a per-cell Stokes auxiliary variable
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Option_module

  implicit none

  type(stokes_auxvar_type) :: auxvar_src, auxvar_dst
  type(option_type)        :: option

  auxvar_dst%pressure = auxvar_src%pressure
  auxvar_dst%u = auxvar_src%u
  auxvar_dst%v = auxvar_src%v
  auxvar_dst%w = auxvar_src%w
  auxvar_dst%du_dt = auxvar_src%du_dt

end subroutine StokesAuxVarCopy

! ************************************************************************** !

subroutine StokesAuxVarCompute(x,auxvar,global_auxvar, &
                               option)
  !
  ! Compute auxiliary variables for each grid cell, given the primary DOFs.
  !
  ! Ucell      : primary unknowns in this cell, size STOKES_NUM_DOF = 4
  ! auxvar     : per-cell Stokes auxiliary variable
  ! global_aux : global auxvar for this cell
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Option_module
  use Global_Aux_module

  implicit none

  PetscReal :: x(STOKES_NUM_DOF)
  type(stokes_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option

  ! Pressure reference (for example)
  auxvar%pressure = x(STOKES_P_DOF)
  auxvar%u = x(STOKES_U_DOF)
  auxvar%v = x(STOKES_V_DOF)
  auxvar%w = x(STOKES_W_DOF)

  global_auxvar%temp = option%flow%reference_temperature
  global_auxvar%den_kg(option%liquid_phase) = &
    option%flow%reference_density(option%liquid_phase)
  global_auxvar%sat(option%liquid_phase) = option%flow%reference_saturation
  if (option%gas_phase > 0) then
    global_auxvar%sat(option%gas_phase) = 0.d0
  endif

end subroutine StokesAuxVarCompute

! ************************************************************************** !

subroutine StokesAuxVarSingleDestroy(auxvar)
  !
  ! Destroy single Stokes auxiliary variable
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  type(stokes_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call StokesAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine StokesAuxVarSingleDestroy

! ************************************************************************** !

subroutine StokesAuxVarArray1Destroy(auxvars)
  !
  ! Destroy array of Stokes auxiliary variables
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  type(stokes_auxvar_type), pointer :: auxvars(:)
  PetscInt :: i

  if (associated(auxvars)) then
    do i = 1, size(auxvars)
      call StokesAuxVarStrip(auxvars(i))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine StokesAuxVarArray1Destroy

! ************************************************************************** !

subroutine StokesAuxVarArray2Destroy(auxvars)
  !
  ! Destroy 2D array of Stokes auxiliary variables
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  type(stokes_auxvar_type), pointer :: auxvars(:,:)
  PetscInt :: i, j

  if (associated(auxvars)) then
    do j = 1, size(auxvars,2)
      do i = 1, size(auxvars,1)
        call StokesAuxVarStrip(auxvars(i,j))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine StokesAuxVarArray2Destroy

! ************************************************************************** !

subroutine StokesAuxVarStrip(auxvar)
  !
  ! If auxvar had dynamically allocated members, free them here.
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(stokes_auxvar_type) :: auxvar

  ! currently nothing to deallocate

end subroutine StokesAuxVarStrip

! ************************************************************************** !

subroutine StokesAuxDestroy(aux)
  !
  ! Destroy top-level aux container
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(stokes_type), pointer :: aux

  if (.not.associated(aux)) return

  call StokesAuxVarDestroy(aux%auxvars)
  call StokesAuxVarDestroy(aux%auxvars_bc)
  call StokesAuxVarDestroy(aux%auxvars_ss)
  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%aP)) then
    deallocate(aux%aP)
  endif
  nullify(aux%aP)

  if (associated(aux%param)) then
    ! nothing special to free inside param
  endif
  nullify(aux%param)

  deallocate(aux)
  nullify(aux)

end subroutine StokesAuxDestroy

end module Stokes_Aux_module
