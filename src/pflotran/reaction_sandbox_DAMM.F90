module Reaction_Sandbox_DAMM_class
  ! 
  ! Implementation of the Dual Arrhenius Michaelis Menten (DAMM) model for
  ! microbial greenhouse generation. Based on the model from the following studies:
  ! Davidson et al. (2012), The Dual Arrhenius and Michaelis–Menten kinetics model for decomposition
  ! of soil organic matter at hourly to seasonal time scales
  ! Sihi et al., (2020) Simultaneous numerical representation of soil microsite
  ! production and consumption of carbon dioxide, methane, and
  ! nitrous oxide using probability distribution functions
  !
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_DAMM_type

    PetscInt :: species_CH2O_id, species_O2_id, species_CO2_id
    PetscInt :: species_NO3_id, species_N2_id, species_HCO3_id
    PetscInt :: species_MNO2_id, species_MN_id
    PetscInt :: species_FEOH3_id, species_FE_id
    PetscInt :: species_SO4_id, species_H2S_id
    PetscInt :: species_CH4_id
    PetscInt :: species_D_id

    PetscInt :: respiration_id, methanogenesis_id
    PetscInt :: methane_oxidation_id
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
    PetscReal :: k_max
    PetscInt :: nrxn
    PetscReal :: K_CH2Oaq
    PetscReal :: K_CH2Oaq_methano
    PetscReal :: K_O2aq
    PetscReal :: I_CO2aq
    PetscReal :: yield
    PetscReal :: k_decay
    PetscReal :: k_CH2O
    PetscReal :: n
    PetscBool :: molarity_units
    PetscBool :: respiration
    PetscReal :: k_o2
    PetscReal :: k_c
    PetscReal :: r_o2_max
    PetscReal :: ae_o2
    PetscReal :: slope_o2
    PetscReal :: intercept_o2
    PetscReal :: alpha_o2
    PetscReal :: k_in_o2
    PetscBool :: methanogenesis
    PetscReal :: k_ch4
    PetscReal :: k_ch4ox
    PetscReal :: k_ch4ox_o2
    PetscReal :: r_ch4_max
    PetscReal :: r_ch4ox_max
    PetscReal :: ae_ch4
    PetscReal :: ae_ch4ox
    PetscReal :: slope_ch4
    PetscReal :: intercept_ch4
    PetscReal :: alpha_ch4
    PetscReal :: slope_ch4ox
    PetscReal :: intercept_ch4ox
    PetscReal :: alpha_ch4ox
    PetscReal :: k_in_ch4
    PetscBool :: methane_oxidation

    PetscReal, pointer :: stoich(:,:)
    PetscReal, pointer :: I_r(:)
  contains
    procedure, public :: ReadInput => DAMMRead
    procedure, public :: Setup => DAMMSetup
    procedure, public :: Evaluate => DAMMEvaluate
    procedure, public :: Destroy => DAMMDestroy
  end type reaction_sandbox_DAMM_type

  public :: DAMMCreate

contains

! ************************************************************************** !

function DAMMCreate()
  !
  ! Allocates example reaction object.
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !

  implicit none

  class(reaction_sandbox_DAMM_type), pointer :: DAMMCreate

  allocate(DAMMCreate)
  DAMMCreate%species_name = ''
  DAMMCreate%species_id = 0
  DAMMCreate%rate_constant = 0.d0
  DAMMCreate%nrxn = 0.d0
  DAMMCreate%k_max = UNINITIALIZED_DOUBLE
  DAMMCreate%K_CH2Oaq = UNINITIALIZED_DOUBLE
  DAMMCreate%K_O2aq = UNINITIALIZED_DOUBLE
  DAMMCreate%I_CO2aq = UNINITIALIZED_DOUBLE
  DAMMCreate%yield = UNINITIALIZED_DOUBLE
  DAMMCreate%k_decay = UNINITIALIZED_DOUBLE
  DAMMCreate%n = 1.d0
  DAMMCreate%molarity_units = PETSC_TRUE
  
  DAMMCreate%respiration = PETSC_FALSE
  DAMMCreate%r_o2_max = UNINITIALIZED_DOUBLE
  DAMMCreate%k_o2 = UNINITIALIZED_DOUBLE
  DAMMCreate%k_c = UNINITIALIZED_DOUBLE
  DAMMCreate%ae_o2 = UNINITIALIZED_DOUBLE
  DAMMCreate%slope_o2 = UNINITIALIZED_DOUBLE
  DAMMCreate%intercept_o2 = UNINITIALIZED_DOUBLE
  DAMMCreate%alpha_o2 = UNINITIALIZED_DOUBLE
  DAMMCreate%k_in_o2 = UNINITIALIZED_DOUBLE

  DAMMCreate%methanogenesis = PETSC_FALSE
  DAMMCreate%r_ch4_max = UNINITIALIZED_DOUBLE
  DAMMCreate%K_CH2Oaq_methano = UNINITIALIZED_DOUBLE
  DAMMCreate%ae_ch4 = UNINITIALIZED_DOUBLE
  DAMMCreate%slope_ch4 = UNINITIALIZED_DOUBLE
  DAMMCreate%intercept_ch4 = UNINITIALIZED_DOUBLE
  DAMMCreate%alpha_ch4 = UNINITIALIZED_DOUBLE
  DAMMCreate%k_in_ch4 = UNINITIALIZED_DOUBLE
  nullify(DAMMCreate%next)
  !DF: call MicrobialCreate

end function DAMMCreate

! ************************************************************************** !

subroutine DAMMRead(this,input,option)
  !
  ! Reads input deck for monod parameters
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_DAMM_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, error_string
  PetscReal :: K_CH2Oaq

  K_CH2Oaq = UNINITIALIZED_DOUBLE
  error_string = 'DAMM_GHG'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
           'CHEMISTRY,REACTION_SANDBOX,DAMM_GHG')
    call StringToUpper(word)

    select case(trim(word))
      case('REACTION')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'CONSTANT', &
                 'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
          call StringToUpper(word)

          select case(trim(word))
            case('RESPIRATION','AEROBIC_RESPIRATION')
              this%respiration = PETSC_TRUE
              this%nrxn = this%nrxn + 1
              call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'CONSTANT', &
                       'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
                call StringToUpper(word)
                select case(trim(word))
                  case('HALF_SATURATION_CONSTANT_O2','MONOD_CONSTANT_O2')
                    call InputReadDouble(input,option,this%k_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_o2, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                 &RESPIRATION,HALF_SATURATION_CONSTANT_O2',option)
                  case('HALF_SATURATION_CONSTANT_C','MONOD_CONSTANT_C')
                    call InputReadDouble(input,option,this%k_CH2Oaq)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_CH2Oaq, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &RESPIRATION,HALF_SATURATION_CONSTANT_C',option)
                  case('RMAX')
                    call InputReadDouble(input,option,this%r_o2_max)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%r_o2_max, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &RESPIRATION,ALPHA',option)
                  case('ACTIVATION_ENERGY')
                    call InputReadDouble(input,option,this%ae_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%ae_o2, &
                                       'J/mol','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &RESPIRATION,ACTIVATION_ENERGY',option)
                  case('MONOD_SLOPE')
                    call InputReadDouble(input,option,this%slope_o2)
                    call InputErrorMsg(input,option,word,error_string)
                  case('MONOD_INTERCEPT')
                    call InputReadDouble(input,option,this%intercept_o2)
                    call InputErrorMsg(input,option,word,error_string)
                  case('ALPHA')
                    call InputReadDouble(input,option,this%alpha_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%alpha_o2, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &RESPIRATION,ALPHA',option)
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'CHEMISTRY,REACTION_SANDBOX,DAMM,RESPIRATION',option)
                end select
              enddo
            case('METHANOGENESIS')
              this%nrxn = this%nrxn + 1
              this%methanogenesis = PETSC_TRUE
              call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'CONSTANT', &
                       'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
                call StringToUpper(word)
                select case(trim(word))
                  case('HALF_SATURATION_CONSTANT_C','MONOD_CONSTANT_C')
                    call InputReadDouble(input,option,this%k_CH2Oaq_methano)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_CH2Oaq_methano, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                 &METHANOGENESIS,HALF_SATURATION_CONSTANT_C',option)
                  case('HALF_SATURATION_CONSTANT','MONOD_CONSTANT')
                    call InputReadDouble(input,option,this%k_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_ch4, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANOGENESIS,HALF_SATURATION_CONSTANT',option)
                  case('RMAX')
                    call InputReadDouble(input,option,this%r_ch4_max)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%r_ch4_max, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANOGENESIS,ALPHA',option)
                  case('ACTIVATION_ENERGY')
                    call InputReadDouble(input,option,this%ae_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%ae_ch4, &
                                       'J/mol','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANOGENESIS,ACTIVATION_ENERGY',option)
                  case('MONOD_SLOPE')
                    call InputReadDouble(input,option,this%slope_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                  case('MONOD_INTERCEPT')
                    call InputReadDouble(input,option,this%intercept_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                  case('ALPHA')
                    call InputReadDouble(input,option,this%alpha_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%alpha_ch4, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANOGENESIS,ALPHA',option)
                  case('INHIBITION_CONSTANT')
                    call InputReadDouble(input,option,this%k_in_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_in_o2, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANOGENESIS,INHIBITION_CONSTANT',option)
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'CHEMISTRY,REACTION_SANDBOX,DAMM,METHANOGENESIS',option)
                end select
              enddo
            case('METHANE_OXIDATION')
              this%nrxn = this%nrxn + 1
              this%methane_oxidation = PETSC_TRUE
              call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'CONSTANT', &
                       'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
                call StringToUpper(word)
                select case(trim(word))
                  case('HALF_SATURATION_CONSTANT_CH4','MONOD_CONSTANT_CH4')
                    call InputReadDouble(input,option,this%k_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_ch4ox, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                 &METHANE_OXIDATION,MONOD_CONSTANT_CH4',option)
                  case('HALF_SATURATION_CONSTANT_O2','MONOD_CONSTANT_O2')
                    call InputReadDouble(input,option,this%k_ch4ox_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_ch4ox_o2, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                 &METHANE_OXIDATION,MONOD_CONSTANT_CH4',option)
                  case('RMAX')
                    call InputReadDouble(input,option,this%r_ch4ox_max)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%r_ch4ox_max, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANE_OXIDATION,RMAX',option)
                  case('ACTIVATION_ENERGY')
                    call InputReadDouble(input,option,this%ae_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%ae_ch4ox, &
                                       'J/mol','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANE_OXIDATION,ACTIVATION_ENERGY',option)
                  case('MONOD_SLOPE')
                    call InputReadDouble(input,option,this%slope_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                  case('MONOD_INTERCEPT')
                    call InputReadDouble(input,option,this%intercept_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                  case('ALPHA')
                    call InputReadDouble(input,option,this%alpha_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%alpha_ch4ox, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,DAMM,REACTION,&
                                                &METHANE_OXIDATION,ALPHA',option)
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'CHEMISTRY,REACTION_SANDBOX,DAMM,METHANE_OXIDATION',option)
                end select
              enddo
          end select
              call InputPopBlock(input,option)
        end do
        call InputPopBlock(input,option)

      case default
        call InputKeywordUnrecognized(input,word, &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
               &DECAY_AQUEOUS,CONSTANT',option)
        end select
  end do
  call InputPopBlock(input,option)


end subroutine DAMMRead

! ************************************************************************** !

subroutine DAMMSetup(this,reaction,option)
  !
  ! Sets up the DAMM reactions with parameters read from the
  ! input deck.
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !
  use Reaction_Aux_module, only : reaction_rt_type, ReactionAuxGetPriSpecIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_DAMM_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  ! Aqueous species
  this%respiration_id = 1
  this%methanogenesis_id = 2
  this%methane_oxidation_id = 3

  allocate(this%stoich(3,reaction%ncomp))

  this%stoich = 0.d0

  word = 'CH2O(aq)'
  this%species_CH2O_id = &
    ReactionAuxGetPriSpecIDFromName(word,reaction,option)
  if (this%respiration) then
    word = 'O2(aq)'
    this%species_O2_id = &
      ReactionAuxGetPriSpecIDFromName(word,reaction,option)
    word = 'CO2(aq)'
    this%species_CO2_id = &
      ReactionAuxGetPriSpecIDFromName(word,reaction,option)
    this%stoich(this%respiration_id,this%species_CH2O_id) = -1.d0
    this%stoich(this%respiration_id,this%species_O2_id) = -1.d0
    this%stoich(this%respiration_id,this%species_CO2_id) = 1.d0
  endif

  if (this%methanogenesis) then
    word = 'CH4(aq)'
    this%species_CH4_id = &
      ReactionAuxGetPriSpecIDFromName(word,reaction,option)
    if (.not. this%respiration) then
      word = 'CO2(aq)'
      this%species_CO2_id = &
        ReactionAuxGetPriSpecIDFromName(word,reaction,option)
    endif
    this%stoich(this%methanogenesis_id,this%species_CH2O_id) = -1.d0
    this%stoich(this%methanogenesis_id,this%species_CH4_id) = 0.5d0
    this%stoich(this%methanogenesis_id,this%species_CO2_id) = 0.5d0
  endif

  if (this%methane_oxidation) then
    if (.not. this%methanogenesis) then
      word = 'CH4(aq)'
      this%species_CH4_id = &
        ReactionAuxGetPriSpecIDFromName(word,reaction,option)
    endif
    if (.not. this%respiration) then
      word = 'CO2(aq)'
      this%species_CO2_id = &
        ReactionAuxGetPriSpecIDFromName(word,reaction,option)
    endif
    this%stoich(this%methane_oxidation_id,this%species_O2_id) = -1.d0
    this%stoich(this%methane_oxidation_id,this%species_CH4_id) = -1.d0
    this%stoich(this%methane_oxidation_id,this%species_CO2_id) = 1.d0
  endif

end subroutine DAMMSetup

! ************************************************************************** !

subroutine DAMMEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_DAMM_type) :: this
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
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: volume_L_water
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: CH2O,O2,CO2,CH4  ! [mole / L water] or [mole / kg water]
  PetscReal :: T, R

  allocate(this%I_r(3))

  T = global_auxvar%temp + 273.15d0 ! K
  R = IDEAL_GAS_CONSTANT            ! J/mol-K

  volume = material_auxvar%volume   ! m^3 bulk volume
  volume_L_water = volume * material_auxvar%porosity * 1.d3 * global_auxvar%sat(1)! L pore water
  if (this%molarity_units) then
    molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
 else
    
    molality_to_molarity = 1.d0
  endif

  CH2O = rt_auxvar%pri_molal(this%species_CH2O_id) * molality_to_molarity ! mol/L

  ! Respiration
  if (this%respiration) then
    O2   = rt_auxvar%pri_molal(this%species_O2_id) * molality_to_molarity
    CO2  = rt_auxvar%pri_molal(this%species_CO2_id) * molality_to_molarity
    if (Initialized(this%alpha_o2)) then
      this%r_o2_max = this%alpha_o2 * exp(-1.d0*this%ae_o2/(R*T)) 
      ! alpha_o2: [mol/L-s]
      ! ae_o2   : [J/mol]
      ! R       : [J/mol-K]
      ! T       : [K]
    endif
    if (Initialized(this%slope_o2)) then
      this%k_o2 = this%intercept_o2 + (T-273.15d0) * this%slope_o2
    endif
    ! I_r units: [mol/L-s] * [mol/L]/([mol/L]+[mol/L]) * [mol/L]/([mol/L]+[mol/L])
    this%I_r(this%respiration_id) = this%r_o2_max * CH2O / (this%k_CH2O + CH2O) * O2 / (this%k_o2 + O2)
    ! Residual units: [mol/s]
    Residual(this%species_O2_id) = Residual(this%species_O2_id) - this%I_r(this%respiration_id) * &
                                   volume_L_water * this%stoich(this%respiration_id,this%species_O2_id)
    Residual(this%species_CO2_id) = Residual(this%species_CO2_id) - this%I_r(this%respiration_id) * &
                                   volume_L_water * this%stoich(this%respiration_id,this%species_CO2_id)
    Residual(this%species_CH2O_id) = Residual(this%species_CH2O_id) - this%I_r(this%respiration_id) * &
                                   volume_L_water * this%stoich(this%respiration_id,this%species_CH2O_id)
  endif

  ! Methanogenesis
  if (this%methanogenesis) then
    CH4 = rt_auxvar%pri_molal(this%species_CH4_id) * molality_to_molarity
    if (.not. this%respiration) then
      O2   = rt_auxvar%pri_molal(this%species_O2_id) * molality_to_molarity
    endif
    if (Initialized(this%alpha_ch4)) then
      this%r_ch4_max = this%alpha_ch4 * exp(-1.d0*this%ae_ch4/(R*T))
    endif
    if (Initialized(this%slope_o2)) then
      this%k_ch4 = this%intercept_ch4 + (T-273.15d0) * this%slope_ch4
    endif
    this%I_r(this%methanogenesis_id) = this%r_ch4_max * CH2O / (this%k_CH2Oaq_methano + CH2O) * &
                                       this%k_in_o2 / (this%k_in_o2 + O2)
    Residual(this%species_CH4_id) = Residual(this%species_CH4_id) - this%I_r(this%methanogenesis_id) * &
                                   volume_L_water * this%stoich(this%methanogenesis_id,this%species_CH4_id)
    Residual(this%species_CO2_id) = Residual(this%species_CO2_id) - this%I_r(this%methanogenesis_id) * &
                                   volume_L_water * this%stoich(this%methanogenesis_id,this%species_CO2_id)
    Residual(this%species_CH2O_id) = Residual(this%species_CH2O_id) - this%I_r(this%methanogenesis_id) * &
                                   volume_L_water * this%stoich(this%methanogenesis_id,this%species_CH2O_id)
  endif

  ! Methane oxidation
  if (this%methane_oxidation) then
    if (.not. this%methanogenesis) then
      CH4 = rt_auxvar%pri_molal(this%species_CH4_id) * molality_to_molarity
    endif
    if (.not. this%respiration) then
      O2  = rt_auxvar%pri_molal(this%species_O2_id) * molality_to_molarity
      CO2 = rt_auxvar%pri_molal(this%species_CO2_id) * molality_to_molarity
    endif
    if (Initialized(this%alpha_ch4ox)) then
      this%r_ch4_max = this%alpha_ch4ox * exp(-1.d0*this%ae_ch4ox/(R*T))
    endif
    if (Initialized(this%slope_o2)) then
      this%k_ch4 = this%intercept_ch4ox + (T-273.15d0) * this%slope_ch4ox
    endif
    this%I_r(this%methane_oxidation_id) = this%r_ch4ox_max * CH4 / (this%k_ch4ox + CH4) * &
                                          O2 / (this%k_ch4ox_o2 + O2)
    Residual(this%species_CH4_id) = Residual(this%species_CH4_id) - this%I_r(this%methane_oxidation_id) * &
                                   volume_L_water * this%stoich(this%methane_oxidation_id,this%species_CH4_id)
    Residual(this%species_O2_id) = Residual(this%species_O2_id) - this%I_r(this%methane_oxidation_id) * &
                                   volume_L_water * this%stoich(this%methane_oxidation_id,this%species_O2_id)
    Residual(this%species_CO2_id) = Residual(this%species_CO2_id) - this%I_r(this%methane_oxidation_id) * &
                                   volume_L_water * this%stoich(this%methane_oxidation_id,this%species_CO2_id)
  endif

end subroutine DAMMEvaluate

! ************************************************************************** !

subroutine DAMMDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !
  implicit none

  class(reaction_sandbox_DAMM_type) :: this

end subroutine DAMMDestroy

end module Reaction_Sandbox_DAMM_class
