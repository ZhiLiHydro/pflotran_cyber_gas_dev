module Reaction_Sandbox_Simple_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_simple_type
    ! Aqueous species
    PetscInt :: species_Aaq_id
    PetscInt :: species_Baq_id
    PetscInt :: species_Caq_id
    PetscInt :: species_Daq_id
    PetscInt :: species_Eaq_id
    PetscInt :: species_Faq_id
    ! Immobile species (e.g. biomass)
    PetscInt :: species_Xim_id
    PetscInt :: species_Yim_id
  contains
    procedure, public :: Setup => SimpleSetup
    procedure, public :: Evaluate => SimpleEvaluate
  end type reaction_sandbox_simple_type

  public :: SimpleCreate

contains

! ************************************************************************** !

function SimpleCreate()
  !
  ! Allocates simple reaction object.
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  implicit none

  class(reaction_sandbox_simple_type), pointer :: SimpleCreate

  allocate(SimpleCreate)
  nullify(SimpleCreate%next)

end function SimpleCreate

! ************************************************************************** !

subroutine SimpleSetup(this,reaction,option)
  !
  ! Sets up the simple reaction with hardwired parameters
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/15

  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Option_module

  implicit none

  class(reaction_sandbox_simple_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  ! Aqueous species
  word = 'Aaq'
  this%species_Aaq_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Baq'
  this%species_Baq_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Caq'
  this%species_Caq_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Daq'
  this%species_Daq_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Eaq'
  this%species_Eaq_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  word = 'Faq'
  this%species_Faq_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)

  ! Immobile species
  word = 'Xim'
  this%species_Xim_id = &
    ReactionImGetSpeciesIDFromName(word,reaction%immobile,option)
  word = 'Yim'
  this%species_Yim_id = &
    ReactionImGetSpeciesIDFromName(word,reaction%immobile,option)


end subroutine SimpleSetup

! ************************************************************************** !

subroutine SimpleEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar,reaction, &
                          option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/15
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Reaction_Inhibition_Aux_module

  implicit none

  class(reaction_sandbox_simple_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume                 ! m^3 bulk volume
  PetscReal :: porosity               ! m^3 pore / m^3 bulk volume
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water

  PetscReal :: Aaq, Baq, Caq, Daq, Eaq, Faq  ! mol/L water
  PetscReal :: Xim, Yim  ! mol/m^3 bulk volume
  PetscReal :: Rate
  PetscReal :: RateA, RateB, RateC, RateD, RateE, RateF, RateX, RateY  ! mol/sec
  PetscReal :: stoichA, stoichB, stoichC, stoichD, stoichE, stoichF
  PetscReal :: stoichX, stoichY
  PetscReal :: k, kr        ! units are problem specific
  PetscReal :: K_Aaq, K_Baq ! [mol/L water]

  PetscReal :: Aaq_rate_constant
  PetscReal :: Aaq_half_saturation
  PetscReal :: Aaq_threshold
  PetscReal :: Aaq_inhibition
  PetscReal :: tempreal
  PetscReal :: dAaq_inhibition_dAaq
  PetscReal :: Aaq_monod
  PetscReal :: dAaq_monod_dAaq
  PetscReal :: Caq_rate_constant
  PetscReal :: Xim_rate_constant
  PetscReal :: Xim_decay_rate_constant
  PetscReal :: common_
  PetscReal :: kg_water
  PetscInt :: inhibition_method
  PetscReal :: kw_water_per_L_water
  PetscReal :: molarity_to_molality
  PetscReal :: immobile_to_aqueous
  PetscBool :: aqueous_biomass_units = PETSC_FALSE
  PetscBool :: aqueous_molarity_units = PETSC_TRUE
  PetscInt :: irow
  PetscInt :: icol
  PetscReal :: dcommon_dAaq
  PetscReal :: dcommon_dXim

  porosity = material_auxvar%porosity               ! [m^3 pore/m^3 bulk volume]
  liquid_saturation = global_auxvar%sat(iphase)     ! [m^3 water/m^3 pore]
  volume = material_auxvar%volume                   ! [m^3 bulk volume]
            ! multiplying by 1.d3 converts m^3 water -> L water
  kg_water = porosity*liquid_saturation*volume*global_auxvar%den_kg(iphase)
  L_water = porosity*liquid_saturation*volume*1.d3
  kw_water_per_L_water = global_auxvar%den_kg(iphase) * 1.d-3

  Aaq = rt_auxvar%total(this%species_Aaq_id,iphase) ! [mol/L water]
  Baq = rt_auxvar%total(this%species_Baq_id,iphase) ! [mol/L water]
  Caq = rt_auxvar%total(this%species_Caq_id,iphase) ! [mol/L water]
  Daq = rt_auxvar%total(this%species_Daq_id,iphase) ! [mol/L water]
  Eaq = rt_auxvar%total(this%species_Eaq_id,iphase) ! [mol/L water]
  Faq = rt_auxvar%total(this%species_Faq_id,iphase) ! [mol/L water]

  Xim = rt_auxvar%immobile(this%species_Xim_id)     ! [mol/m^3 bulk volume]
  Yim = rt_auxvar%immobile(this%species_Yim_id)     ! [mol/m^3 bulk volume]

  ! initialize all rates to zero
  Rate = 0.d0

  RateA = 0.d0
  RateB = 0.d0
  RateC = 0.d0
  RateD = 0.d0
  RateE = 0.d0
  RateF = 0.d0
  RateX = 0.d0
  RateY = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichA = 0.d0
  stoichB = 0.d0
  stoichC = 0.d0
  stoichD = 0.d0
  stoichE = 0.d0
  stoichF = 0.d0
  stoichX = 0.d0
  stoichY = 0.d0

  ! kinetic rate constants
  k = 0.d0
  kr = 0.d0

  ! Monod half-saturation constants
  K_Aaq = 0.d0
  K_Baq = 0.d0

  !----------------------------------------------------------------------------
  ! multiplicative Monod w/biomass (A + 2B -> C)
  Aaq_rate_constant = -2.66d-18     ! mol Aaq / mol Xim - sec
  Aaq_half_saturation = 2.93d-5     ! mol / L water
  Caq_rate_constant = 2.02d-18      ! mol Caq / mol Xim - sec
  Xim_rate_constant = 6.01d-6       ! 1 / sec
  Xim_decay_rate_constant = 1.d-12  ! 1 / sec

  tempreal = (Aaq_half_saturation + Aaq)
  Aaq_monod = Aaq / tempreal
  dAaq_monod_dAaq = Aaq_half_saturation /(tempreal*tempreal)

  Aaq_threshold = 1.d-9
  inhibition_method = INHIBITION_SMOOTHSTEP
!  inhibition_method = INHIBITION_THRESHOLD
  select case(inhibition_method)
    case(INHIBITION_SMOOTHSTEP)
      call ReactionInhibitionSmoothstep(Aaq,Aaq_threshold,0.1d0, &
                                        Aaq_inhibition,dAaq_inhibition_dAaq)
    case(INHIBITION_THRESHOLD)
      call ReactionInhibitionThreshold(Aaq,Aaq_threshold, &
                                       Aaq_inhibition,dAaq_inhibition_dAaq)
    case default
      Aaq_inhibition = 1.d0
      dAaq_inhibition_dAaq = 0.d0
  end select

  common_ = Xim * Aaq_monod * Aaq_inhibition * volume    ! mol Xim
  dcommon_dAaq = Xim * &
                 (dAaq_monod_dAaq * Aaq_inhibition + &
                  Aaq_monod * dAaq_inhibition_dAaq) * volume
  dcommon_dXim = Aaq_monod * Aaq_inhibition * volume

  RateA = Aaq_rate_constant * common_        ! mol Aaq / sec
  RateC = Caq_rate_constant * common_        ! mol Caq / sec
  RateX = Xim_rate_constant * common_  - &   ! mol Xim / sec
              Xim_decay_rate_constant * Xim * volume

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second
  Residual(this%species_Aaq_id) = Residual(this%species_Aaq_id) - RateA
  Residual(this%species_Baq_id) = Residual(this%species_Baq_id) - RateB
  Residual(this%species_Caq_id) = Residual(this%species_Caq_id) - RateC
  Residual(this%species_Daq_id) = Residual(this%species_Daq_id) - RateD
  Residual(this%species_Eaq_id) = Residual(this%species_Eaq_id) - RateE
  Residual(this%species_Faq_id) = Residual(this%species_Faq_id) - RateF
  Residual(this%species_Xim_id + reaction%offset_immobile) = &
    Residual(this%species_Xim_id + reaction%offset_immobile) - RateX
  Residual(this%species_Yim_id + reaction%offset_immobile) = &
    Residual(this%species_Yim_id + reaction%offset_immobile) - RateY

  if (compute_derivative) then
    ! dR_Aaq / dAaq
    irow = this%species_Aaq_id
    icol = irow
    Jacobian(irow,icol) = Jacobian(irow,icol) - & ! kg water / sec
      Aaq_rate_constant * dcommon_dAaq * kw_water_per_L_water
    ! dR_Aaq / dXim
    icol = this%species_Xim_id + reaction%offset_immobile
    Jacobian(irow,icol) = Jacobian(irow,icol) - &
      Aaq_rate_constant * dcommon_dXim
    ! dR_Caq / dAaq
    irow = this%species_Caq_id
    icol = this%species_Aaq_id
    Jacobian(irow,icol) = Jacobian(irow,icol) - & ! kg water / sec
      Caq_rate_constant * dcommon_dAaq * kw_water_per_L_water
    ! dR_Caq / dXim
    icol = this%species_Xim_id + reaction%offset_immobile
    Jacobian(irow,icol) = Jacobian(irow,icol) - &
      Caq_rate_constant * dcommon_dXim
    ! dR_Xim / dAaq
    irow = this%species_Xim_id + reaction%offset_immobile
    icol = this%species_Aaq_id
    Jacobian(irow,icol) = Jacobian(irow,icol) - & ! kg water / sec
      Xim_rate_constant * dcommon_dAaq * kw_water_per_L_water
    ! dR_Xim / dXim
    icol = irow
    Jacobian(irow,icol) = Jacobian(irow,icol) - &
      (Xim_rate_constant * dcommon_dXim - &
       Xim_decay_rate_constant * volume)
  endif

end subroutine SimpleEvaluate

end module Reaction_Sandbox_Simple_class
