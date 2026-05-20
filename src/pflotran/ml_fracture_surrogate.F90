module ML_Fracture_Surrogate_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use ML_Surrogate_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscReal, parameter :: SECONDS_PER_YEAR = DAYS_PER_YEAR * 24.d0 * 3600.d0

  ! Default r_norm offset / scale matching the 1D training deck.
  PetscReal, parameter :: DEFAULT_R_OFFSET = 2.55d0
  PetscReal, parameter :: DEFAULT_R_SCALE  = 8.3152d0

  ! Tortuosity exponent: D_eff = D_free * phi^3.05 = D_free * tortuosity * phi
  ! with tortuosity = phi^2.05, matching the QPAC formula in the design doc.
  PetscReal, parameter :: TORTUOSITY_EXPONENT = 2.05d0

  ! Floor for log10 inputs to avoid -Inf on zero velocity / aperture.
  PetscReal, parameter :: LOG10_FLOOR = 1.d-20

  type, public :: ml_fracture_surrogate_type
    PetscBool :: active
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: fracture_material_id
    PetscReal :: r_offset
    PetscReal :: r_scale
    type(ml_mlp_type) :: mlp
  end type ml_fracture_surrogate_type

  public :: MLFractureSurrogateInit, &
            MLFractureSurrogateReadInput, &
            MLFractureSurrogateSetup, &
            MLFractureSurrogateUpdate, &
            MLFractureSurrogateDestroy

contains

! ************************************************************************** !

subroutine MLFractureSurrogateInit(cfg)
  !
  ! Default-initializes the surrogate config to OFF.  Activation happens
  ! in MLFractureSurrogateReadInput when the FRACTURE_ML_SURROGATE block
  ! is encountered in the input deck.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-20
  !
  implicit none

  type(ml_fracture_surrogate_type) :: cfg

  cfg%active = PETSC_FALSE
  cfg%filename = ''
  cfg%fracture_material_id = UNINITIALIZED_INTEGER
  cfg%r_offset = DEFAULT_R_OFFSET
  cfg%r_scale = DEFAULT_R_SCALE
  call MLPCreate(cfg%mlp)

end subroutine MLFractureSurrogateInit

! ************************************************************************** !

subroutine MLFractureSurrogateReadInput(cfg, input, option)
  !
  ! Parses the FRACTURE_ML_SURROGATE nested block under SUBSURFACE_TRANSPORT:
  !
  !   FRACTURE_ML_SURROGATE
  !     FILENAME             <path>
  !     FRACTURE_MATERIAL_ID <int>
  !     R_OFFSET             <real>   ! optional, default 2.55
  !     R_SCALE              <real>   ! optional, default 8.3152
  !   /
  !
  ! Presence of the block sets cfg%active = PETSC_TRUE.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-20
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(ml_fracture_surrogate_type) :: cfg
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SUBSURFACE_TRANSPORT,FRACTURE_ML_SURROGATE'

  cfg%active = PETSC_TRUE

  call InputPushBlock(input, option)
  do
    call InputReadPflotranString(input, option)
    if (InputError(input)) exit
    if (InputCheckExit(input, option)) exit

    call InputReadCard(input, option, keyword)
    call InputErrorMsg(input, option, 'keyword', error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('FILENAME')
        call InputReadFilename(input, option, cfg%filename)
        call InputErrorMsg(input, option, keyword, error_string)
      case('FRACTURE_MATERIAL_ID')
        call InputReadInt(input, option, cfg%fracture_material_id)
        call InputErrorMsg(input, option, keyword, error_string)
      case('R_OFFSET')
        call InputReadDouble(input, option, cfg%r_offset)
        call InputErrorMsg(input, option, keyword, error_string)
      case('R_SCALE')
        call InputReadDouble(input, option, cfg%r_scale)
        call InputErrorMsg(input, option, keyword, error_string)
      case default
        call InputKeywordUnrecognized(input, keyword, error_string, option)
    end select
  enddo
  call InputPopBlock(input, option)

end subroutine MLFractureSurrogateReadInput

! ************************************************************************** !

subroutine MLFractureSurrogateSetup(cfg, option)
  !
  ! Validates the parsed config and loads MLP weights from HDF5.  Called
  ! once at startup from PMRTSetup when cfg%active.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-20
  !
  use Option_module

  implicit none

  type(ml_fracture_surrogate_type) :: cfg
  type(option_type) :: option

  if (len_trim(cfg%filename) == 0) then
    option%io_buffer = 'FRACTURE_ML_SURROGATE requires a FILENAME entry.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(cfg%fracture_material_id)) then
    option%io_buffer = 'FRACTURE_ML_SURROGATE requires FRACTURE_MATERIAL_ID.'
    call PrintErrMsg(option)
  endif
  if (cfg%r_scale == 0.d0) then
    option%io_buffer = 'FRACTURE_ML_SURROGATE R_SCALE must be non-zero.'
    call PrintErrMsg(option)
  endif

  call MLPLoad(cfg%mlp, cfg%filename, option)

end subroutine MLFractureSurrogateSetup

! ************************************************************************** !

subroutine MLFractureSurrogateUpdate(realization, cfg)
  !
  ! Per-timestep hook: for every local fracture cell, build the 4-input
  ! vector, run the MLP forward pass, and overlay the predicted porosity
  ! together with the derived permeability (cubic law) and tortuosity
  ! (phi^2.05).  Other cells are untouched; full RT chemistry still runs
  ! on every cell -- the surrogate just wins for fracture cells.
  !
  ! Velocity comes from PatchGetCellCenteredVelocities, which is local
  ! only.  Ghost cells are not updated here; they pick up new values via
  ! PFLOTRAN's normal property communication.
  !
  ! Author: Heeho Park
  ! Date: 2026-05-20
  !
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(ml_fracture_surrogate_type) :: cfg

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, allocatable :: velocities(:,:)
  PetscReal :: x_in(ML_N_IN)
  PetscReal :: phi_new, k_new, tortuosity_new
  PetscReal :: v_darcy, aperture, r, r_norm, time_y, k_cell
  PetscInt :: local_id, ghosted_id, i, iphase

  patch => realization%patch
  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars
  iphase = 1

  allocate(velocities(3, grid%nlmax))
  call PatchGetCellCenteredVelocities(patch, iphase, velocities)

  time_y = realization%option%time / SECONDS_PER_YEAR

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) /= cfg%fracture_material_id) cycle

    v_darcy = sqrt(velocities(1,local_id)**2 + &
                   velocities(2,local_id)**2 + &
                   velocities(3,local_id)**2)
    v_darcy = max(v_darcy, LOG10_FLOOR)

    k_cell = material_auxvars(ghosted_id)%permeability(perm_xx_index)
    aperture = sqrt(12.d0 * k_cell)
    aperture = max(aperture, LOG10_FLOOR)

    r = sqrt(grid%x(ghosted_id)**2 + grid%y(ghosted_id)**2)
    r_norm = (r - cfg%r_offset) / cfg%r_scale

    x_in(1) = log10(v_darcy)
    x_in(2) = log10(aperture)
    x_in(3) = time_y
    x_in(4) = r_norm

    call MLPForward(cfg%mlp, x_in, phi_new)

    k_new = phi_new**2 * aperture / 12.d0
    tortuosity_new = phi_new**TORTUOSITY_EXPONENT

    material_auxvars(ghosted_id)%porosity = phi_new
    do i = 1, size(material_auxvars(ghosted_id)%permeability)
      material_auxvars(ghosted_id)%permeability(i) = k_new
    enddo
    material_auxvars(ghosted_id)%tortuosity = tortuosity_new
  enddo

  deallocate(velocities)

end subroutine MLFractureSurrogateUpdate

! ************************************************************************** !

subroutine MLFractureSurrogateDestroy(cfg)
  !
  ! Author: Heeho Park
  ! Date: 2026-05-20
  !
  implicit none

  type(ml_fracture_surrogate_type) :: cfg

  call MLPDestroy(cfg%mlp)
  cfg%active = PETSC_FALSE

end subroutine MLFractureSurrogateDestroy

! ************************************************************************** !

end module ML_Fracture_Surrogate_module
