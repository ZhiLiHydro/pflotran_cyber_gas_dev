module ML_Fracture_Surrogate_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use ML_Surrogate_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscReal, parameter :: SECONDS_PER_YEAR = DAYS_PER_YEAR * 24.d0 * 3600.d0

  ! Tortuosity exponent: D_eff = D_free * phi^3.05 = D_free * tortuosity * phi
  ! with tortuosity = phi^2.05, matching the QPAC formula in the design doc.
  PetscReal, parameter :: TORTUOSITY_EXPONENT = 2.05d0

  ! Rate-integration contract (C_HUBER_SURROGATE_HANDOFF.md sections 0, 5, 8).
  PetscReal, parameter :: PHI_FLOOR   = 0.02d0    ! = MINIMUM_POROSITY
  PetscReal, parameter :: PHI_SEALED  = 0.0201d0  ! skip network below this
  PetscReal, parameter :: RATE_MAX    = 3.0d0     ! |dphi/dt| envelope [1/y]
  PetscReal, parameter :: START_TIME_Y = 1.0d0    ! surrogate active after 1 y

  type, public :: ml_fracture_surrogate_type
    PetscBool :: active
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: fracture_material_id   ! deck-supplied external ID
    PetscInt :: fracture_imat          ! internal index resolved at Setup
    PetscReal :: phi_initial           ! snapshot of deck POROSITY at Setup
    PetscReal :: k_initial             ! snapshot of deck PERM_XX  at Setup
    PetscReal :: aperture_m            ! sqrt(12*k_initial) [m], MLP input #2
    PetscInt :: id_h, id_ca, id_si, id_al, id_mg   ! primary species indices
    PetscReal, pointer :: phi_frac(:)  ! surrogate-owned porosity per local cell
    character(len=MAXSTRINGLENGTH) :: diag_file    ! optional 7-input CSV dump
    PetscInt :: diag_unit
    PetscInt :: diag_stride
    PetscInt :: diag_step
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
  ! Date: 2026-07-24
  !
  implicit none

  type(ml_fracture_surrogate_type) :: cfg

  cfg%active = PETSC_FALSE
  cfg%filename = ''
  cfg%fracture_material_id = UNINITIALIZED_INTEGER
  cfg%fracture_imat = UNINITIALIZED_INTEGER
  cfg%phi_initial = UNINITIALIZED_DOUBLE
  cfg%k_initial = UNINITIALIZED_DOUBLE
  cfg%aperture_m = UNINITIALIZED_DOUBLE
  cfg%id_h = UNINITIALIZED_INTEGER
  cfg%id_ca = UNINITIALIZED_INTEGER
  cfg%id_si = UNINITIALIZED_INTEGER
  cfg%id_al = UNINITIALIZED_INTEGER
  cfg%id_mg = UNINITIALIZED_INTEGER
  cfg%diag_file = ''
  cfg%diag_unit = -1
  cfg%diag_stride = 1
  cfg%diag_step = 0
  nullify(cfg%phi_frac)
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
  !     DIAGNOSTIC_FILE      <path>   ! optional: dump 7 inputs + raw rate
  !     DIAGNOSTIC_STRIDE    <int>    ! optional: dump every Nth step (def 1)
  !   /
  !
  ! Presence of the block sets cfg%active = PETSC_TRUE.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
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
      case('DIAGNOSTIC_FILE')
        call InputReadFilename(input, option, cfg%diag_file)
        call InputErrorMsg(input, option, keyword, error_string)
      case('DIAGNOSTIC_STRIDE')
        call InputReadInt(input, option, cfg%diag_stride)
        call InputErrorMsg(input, option, keyword, error_string)
      case default
        call InputKeywordUnrecognized(input, keyword, error_string, option)
    end select
  enddo
  call InputPopBlock(input, option)

end subroutine MLFractureSurrogateReadInput

! ************************************************************************** !

subroutine MLFractureSurrogateSetup(cfg, realization, option)
  !
  ! Validates the parsed config, resolves the deck-supplied external material
  ! ID to its internal index (what patch%imat stores per cell), resolves the
  ! five primary-species indices by name, allocates the surrogate-owned
  ! per-cell porosity, and loads MLP weights.  Called once from PMRTSetup.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Reaction_Aux_module
  use Option_module

  implicit none

  type(ml_fracture_surrogate_type) :: cfg
  class(realization_subsurface_type) :: realization
  type(option_type) :: option

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: name
  character(len=MAXSTRINGLENGTH) :: diag_name
  PetscInt :: i

  if (len_trim(cfg%filename) == 0) then
    option%io_buffer = 'FRACTURE_ML_SURROGATE requires a FILENAME entry.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(cfg%fracture_material_id)) then
    option%io_buffer = 'FRACTURE_ML_SURROGATE requires FRACTURE_MATERIAL_ID.'
    call PrintErrMsg(option)
  endif

  patch => realization%patch
  grid => patch%grid

  cfg%fracture_imat = UNINITIALIZED_INTEGER
  do i = 1, size(patch%material_property_array)
    if (patch%material_property_array(i)%ptr%external_id == &
        cfg%fracture_material_id) then
      cfg%fracture_imat = i
      exit
    endif
  enddo
  if (Uninitialized(cfg%fracture_imat)) then
    write(option%io_buffer,'(a,i0,a)') &
      'FRACTURE_ML_SURROGATE: no MATERIAL_PROPERTY with ID = ', &
      cfg%fracture_material_id, ' was found in the deck.'
    call PrintErrMsg(option)
  endif

  cfg%phi_initial = &
    patch%material_property_array(cfg%fracture_imat)%ptr%porosity
  cfg%k_initial = &
    patch%material_property_array(cfg%fracture_imat)%ptr%permeability(1,1)
  if (cfg%phi_initial <= 0.d0 .or. cfg%k_initial <= 0.d0) then
    option%io_buffer = 'FRACTURE_ML_SURROGATE: fracture MATERIAL_PROPERTY &
      &requires positive POROSITY and PERMEABILITY.'
    call PrintErrMsg(option)
  endif
  cfg%aperture_m = sqrt(12.d0 * cfg%k_initial)

  ! Primary-species indices into rt_auxvar%total, resolved by name.
  name = 'H+';         cfg%id_h  = ReactionAuxGetPriSpecIDFromName(name, &
                                     realization%reaction, option)
  name = 'Ca++';       cfg%id_ca = ReactionAuxGetPriSpecIDFromName(name, &
                                     realization%reaction, option)
  name = 'H4SiO4(aq)'; cfg%id_si = ReactionAuxGetPriSpecIDFromName(name, &
                                     realization%reaction, option)
  name = 'Al+++';      cfg%id_al = ReactionAuxGetPriSpecIDFromName(name, &
                                     realization%reaction, option)
  name = 'Mg++';       cfg%id_mg = ReactionAuxGetPriSpecIDFromName(name, &
                                     realization%reaction, option)

  ! Surrogate-owned porosity: fracture cells evolve from phi_initial.
  allocate(cfg%phi_frac(grid%nlmax))
  cfg%phi_frac = cfg%phi_initial

  ! Optional dump for external forward-pass verification: one file per rank,
  ! holding the raw MLP inputs and the pre-clamp network output.
  if (len_trim(cfg%diag_file) > 0) then
    write(diag_name,'(a,a,i0,a)') trim(cfg%diag_file), '.', option%myrank, '.csv'
    open(newunit=cfg%diag_unit, file=trim(diag_name), action='write', &
         status='replace')
    write(cfg%diag_unit,'(a)') &
      'time_y,cell_id,porosity,aperture,TotalH,Ca,Si,Al,Mg,rate_per_yr'
  endif

  call MLPLoad(cfg%mlp, cfg%filename, option)

end subroutine MLFractureSurrogateSetup

! ************************************************************************** !

subroutine MLFractureSurrogateUpdate(realization, cfg)
  !
  ! Per-timestep hook: for every local fracture cell, assemble the 7 B-huber
  ! features (current porosity, fixed aperture, TotalH, and the Ca/Si/Al/Mg
  ! local total concentrations), run the MLP to get dphi/dt, clamp to the
  ! clogging branch, Euler-step the surrogate porosity, and overlay porosity +
  ! derived permeability (cubic law) and tortuosity (phi^2.05) onto the material
  ! auxvar.  Full RT chemistry still runs on every cell; the surrogate
  ! overwrites fracture cells afterward.
  !
  ! The surrogate owns fracture porosity (cfg%phi_frac): the network consumes
  ! last step's value, so it must persist and not be recomputed from minerals.
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Material_Aux_module
  use Reactive_Transport_Aux_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(ml_fracture_surrogate_type) :: cfg

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  PetscReal :: x_in(ML_N_IN)
  PetscReal :: phi_old, phi_new, rate, dt_y, time_y
  PetscReal :: k_new, tortuosity_new, phi_ratio
  PetscInt :: local_id, ghosted_id, i, iphase

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  material_auxvars => patch%aux%Material%auxvars
  rt_auxvars => patch%aux%RT%auxvars
  iphase = 1                        ! aqueous phase

  time_y = option%time / SECONDS_PER_YEAR
  if (time_y < START_TIME_Y) return ! skip the sub-year equilibration transient
  dt_y = option%tran_dt / SECONDS_PER_YEAR
  cfg%diag_step = cfg%diag_step + 1

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) /= cfg%fracture_imat) cycle

    phi_old = cfg%phi_frac(local_id)
    if (phi_old <= PHI_SEALED) then
      phi_new = phi_old             ! cell is sealed; skip the network
    else
      x_in(1)  = phi_old
      x_in(2)  = cfg%aperture_m
      x_in(3)  = rt_auxvars(ghosted_id)%total(cfg%id_h,  iphase)  ! TotalH
      x_in(4)  = rt_auxvars(ghosted_id)%total(cfg%id_ca, iphase)
      x_in(5)  = rt_auxvars(ghosted_id)%total(cfg%id_si, iphase)
      x_in(6)  = rt_auxvars(ghosted_id)%total(cfg%id_al, iphase)
      x_in(7)  = rt_auxvars(ghosted_id)%total(cfg%id_mg, iphase)

      call MLPForward(cfg%mlp, x_in, rate)

      ! Dump the raw network output: the Python reference is unclamped, so the
      ! clamps below must not be applied to what gets compared.  (newunit gives
      ! a negative unit, so gate on diag_file, not diag_unit sign.)
      if (len_trim(cfg%diag_file) > 0) then
        if (mod(cfg%diag_step, cfg%diag_stride) == 0) then
          write(cfg%diag_unit,'(es17.9,a,i0,8(a,es24.16))') &
            time_y, ',', grid%nG2A(ghosted_id), &
            ',', x_in(1), ',', x_in(2), ',', x_in(3), ',', x_in(4), &
            ',', x_in(5), ',', x_in(6), ',', x_in(7), ',', rate
        endif
      endif

      rate = min(rate, 0.d0)          ! clogging-only branch (this contract)
      rate = max(rate, -RATE_MAX)     ! out-of-domain sinh envelope guard
      phi_new = max(phi_old + rate * dt_y, PHI_FLOOR)
    endif

    cfg%phi_frac(local_id) = phi_new

    phi_ratio = phi_new / cfg%phi_initial
    k_new = phi_ratio**3 * cfg%k_initial
    tortuosity_new = phi_new**TORTUOSITY_EXPONENT

    material_auxvars(ghosted_id)%porosity = phi_new
    do i = 1, size(material_auxvars(ghosted_id)%permeability)
      material_auxvars(ghosted_id)%permeability(i) = k_new
    enddo
    material_auxvars(ghosted_id)%tortuosity = tortuosity_new
  enddo

end subroutine MLFractureSurrogateUpdate

! ************************************************************************** !

subroutine MLFractureSurrogateDestroy(cfg)
  !
  ! Author: Heeho Park
  ! Date: 2026-07-24
  !
  implicit none

  type(ml_fracture_surrogate_type) :: cfg

  if (len_trim(cfg%diag_file) > 0) close(cfg%diag_unit)
  cfg%diag_unit = -1
  if (associated(cfg%phi_frac)) deallocate(cfg%phi_frac)
  nullify(cfg%phi_frac)
  call MLPDestroy(cfg%mlp)
  cfg%active = PETSC_FALSE

end subroutine MLFractureSurrogateDestroy

! ************************************************************************** !

end module ML_Fracture_Surrogate_module
