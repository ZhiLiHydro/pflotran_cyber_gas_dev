module Reaction_Sandbox_Particle_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_particle_type
    ! Implementation of 6 particle and emulation transport models listed 
    !   in table5.2 of 'Emulsion-based delivery of alkalinity for subsurface
    !   pH control' by Katherine Muller 2016

    ! Model index
    PetscInt :: n_model ! M-1~6 six models in Table5.2

    !ID number of particle species
    PetscInt :: species_aq_id ! Aqueous species
    PetscInt :: species_im_id ! Immobile species

    !Name of particle species
    character(len=MAXWORDLENGTH) :: name_aqueous
    character(len=MAXWORDLENGTH) :: name_immobile

    !Attachment rate (filtration model in Tufenkji & Elimelech, 2004)
    PetscReal :: diam_collector
    PetscReal :: diam_particle
    PetscReal :: hamaker_constant
    PetscReal :: BOLTZMANN_CONSTANT = 1.380649d-23 !J/K
    PetscReal :: density_particle
    PetscReal :: alpha_efficiency

    !Attachment/Detachment rates (Constant)
    PetscReal :: attachment_rate  ! 1/s
    PetscReal :: detachment_rate  ! 1/s

    PetscReal :: max_solid_conc  ! Smax [mol/m^3 bulk]
                                 ! was [mol/kg bulk] in the literature
    PetscReal :: beta            ! depth-dependent exponent [-]
    
  contains
    procedure, public :: ReadInput => ParticleRead
    procedure, public :: Setup => ParticleSetup
    procedure, public :: Evaluate => ParticleEvaluate
    procedure, public :: Destroy => ParticleDestroy
  end type reaction_sandbox_particle_type

  public :: ParticleCreate 

contains

! ************************************************************************** !

function ParticleCreate()
  !
  ! Allocates particle transport reaction object.
  !
  ! Author: Bryan He
  ! Date: 06/07/2024
  !

  implicit none

  class(reaction_sandbox_particle_type), pointer :: ParticleCreate

  allocate(ParticleCreate)
  ParticleCreate%n_model = 0
  ParticleCreate%name_aqueous = ''
  ParticleCreate%name_immobile = ''
  ParticleCreate%species_aq_id = 0
  ParticleCreate%species_im_id = 0
  ParticleCreate%max_solid_conc = 0.d0
  ParticleCreate%beta = 0.d0

  !Filtration Model
  ParticleCreate%diam_collector = 0.d0
  ParticleCreate%diam_particle = 0.d0
  ParticleCreate%hamaker_constant = 0.d0
  ParticleCreate%density_particle = 0.d0
  ParticleCreate%alpha_efficiency = 1.d0

  !Attachment rates
  ParticleCreate%attachment_rate = -1.d0
  !Detachment rates
  ParticleCreate%detachment_rate = 0.d0
  nullify(ParticleCreate%next)


end function ParticleCreate

! ************************************************************************** !

subroutine ParticleRead(this,input,option)
  !
  ! Reads input deck for particle transport reaction parameters
  !
  ! Author: Bryan He
  ! Date: 06/07/2024
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_particle_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT')
    call StringToUpper(word)

    select case(trim(word))
      case('PARTICLE_NAME_AQ')
        ! Mobile particle name 
        call InputReadWord(input,option,this%name_aqueous,PETSC_TRUE)
        call InputErrorMsg(input,option,'PARTICLE_NAME_AQ', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT')

      case('PARTICLE_NAME_IM')
        ! Immobile particle name 
        call InputReadWord(input,option,this%name_immobile,PETSC_TRUE)
        call InputErrorMsg(input,option,'PARTICLE_NAME_IM', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT')

      case('DEPOSTION_LANGMUIRIAN_BLOCKING')
        this%n_model = 1
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,word, &
                 'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT')
          call StringToUpper(word)  
          select case(trim(word))
            case('ATTACHMENT_RATE')
            ! Read kd in Eqn 5.13
              call ReadAttachmentRate(this,input,word,option)              
            case('MAX_SOLID_CONCENTRATION')
            ! Read Smax in Eqn 5.13 in mol/m^3 bulk
              call ReadMaxSolidConc(this,input,word,option)
            case default
              call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT,&
                     &DEPOSTION_LANGMUIRIAN_BLOCKING',option)
          end select    
        enddo
        call InputPopBlock(input,option)
              
!      case('DEPTH_DEPENDENT_RETENTION')
!         this%n_model = 2
!         do
!            select case(trim(word))
!      case('DEPTH_DEPENDENT_RETENTION_LANGMUIRIAN_BLOCKING')
!         this%n_model = 3
!      case('DEPOSTION_LANGMUIRIAN_BLOCKING_DEPTH_DEPENDENT_RETENTION')
!         this%n_model = 4 
!      case('DEPOSTION_WITHOUT_DISPERSION')
!         this%n_model = 5
!      case('DEPOSTION_WITH_SURFACE_CAPACITY')    
!         this%n_model = 6

      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ParticleRead

! ************************************************************************** !

subroutine ParticleSetup(this,reaction,option)
  !
  ! Sets up the particle transport reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: Bryan He
  ! Date: 06/07/2024
  !
  use Reaction_Aux_module
  use Option_module
  use Reaction_Immobile_Aux_module

  implicit none

  class(reaction_sandbox_particle_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%species_aq_id = &
    ReactionAuxGetPriSpecIDFromName(this%name_aqueous,reaction,option)
  this%species_im_id = &
    ReactionImGetSpeciesIDFromName(this%name_immobile,reaction%immobile,option)

end subroutine ParticleSetup

! ************************************************************************** !

subroutine ParticleEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: Bryan He
  ! Date: 06/07/2024
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_particle_type) :: this
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

  PetscReal :: volume                 ! m^3 bulk
  PetscReal :: porosity               ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water
  PetscReal :: temperature
  PetscReal :: bulk_density
  PetscReal :: rho_f
  PetscReal :: g
  PetscReal :: viscosity
  PetscReal :: psi_b,psi_d
  PetscReal :: dpsi_b_dc ![m^3_bulk/mol]
  PetscReal :: k_d ! deposition rate (attachment rate)
  PetscReal :: k_r ! remobilization rate (deachment rate)
  PetscReal :: smax
  PetscReal :: aq_conc !mol/kg_water
  PetscReal :: im_conc !mol/m^3_bulk
  PetscReal :: att_rate, det_rate  ! mol/sec
  PetscReal :: stoich_aq,stoich_im
  PetscReal :: exchange ! mol/(m^3_pore s) 
  PetscInt :: aq_id,im_id
  PetscReal :: tempreal

  ! Global stuff (Check global_aux.F90)
  volume = material_auxvar%volume
  porosity = material_auxvar%porosity
  bulk_density = material_auxvar%soil_particle_density
  liquid_saturation = global_auxvar%sat(iphase)
  L_water = porosity*liquid_saturation*volume*1.d3
  viscosity = 0.0008891
  ! Ns/m2 (This should be a field but not found in global_auxvar)

  temperature = global_auxvar%temp
  rho_f = global_auxvar%den_kg(iphase)
  g = EARTH_GRAVITY

  ! Assign concentrations of Vaq and Vim
  aq_conc = rt_auxvar%pri_molal(this%species_aq_id) ! mol/kg_water
  im_conc = rt_auxvar%immobile(this%species_im_id)  ! mol/m^3_bulk

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry 
  stoich_aq = -1.d0
  stoich_im = 1.d0

  if (this%n_model==1) then
  ! M-1 in table 5.2
    k_d = AttachmentRate(this,global_auxvar,iphase,temperature,viscosity, &
                        porosity,rho_f)
    smax = this%max_solid_conc
    call BlockingFunction(im_conc,smax,psi_b,dpsi_b_dc)
    ! eqn 5.13
    exchange = &    ![mol/(kg_bulk s)]
             liquid_saturation * k_d * psi_b * aq_conc * rho_f &
! unit     [m^3_water/m^3_pore] [1/s]  [-] [mol/kg_water] [kg_water/m^3_water]
             * porosity / bulk_density        
! unit   [m^3_pore/m^3_bulk] [kg_bulk/m^3_bulk]            

!    print *, 'exchange rate:',exchange 
!  elseif (this%n_model==2) then
!  elseif (this%n_model==3) then
!  elseif (this%n_model==4) then
!  elseif (this%n_model==5) then
!  elseif (this%n_model==6) then
  
  end if  

  att_rate = stoich_aq * exchange * volume * bulk_density  ! [mol/s]
! unit       [-]  [mol/(kg_bulk s)] [m^3_bulk] [kg_bulk/m^3_bulk]

  det_rate = stoich_im * exchange * volume * bulk_density  ! [mol/s]
! unit       [-]  [mol/(kg_bulk s)] [m^3_bulk] [kg_bulk/m^3_bulk]

  aq_id = this%species_aq_id
  im_id = this%species_im_id+reaction%offset_immobile
! residual
  Residual(aq_id) = Residual(aq_id) - att_rate
  Residual(im_id) = Residual(im_id) - det_rate   

! Jacobian
  if (compute_derivative) then    
    tempreal = & ! [kg_water/s]
             liquid_saturation * k_d * psi_b * &
!unit      [m^3_water/m^3_pore] [1/s]   [-]                 
             rho_f * volume * porosity 
!unit  [kg_water/m^3_water] [m^3_bulk] [m^3_pore/m^3_bulk]

    Jacobian(aq_id,aq_id) = Jacobian(aq_id,aq_id) - stoich_aq * tempreal
    Jacobian(im_id,aq_id) = Jacobian(im_id,aq_id) - stoich_im * tempreal   

    tempreal = & ! [m^3_bulk/s]
             liquid_saturation * k_d * dpsi_b_dc * & !bulk_density * & 
!unit      [m^3_water/m^3_pore] [1/s] [m^3_bulk/mol] [kg_bulk/m^3_bulk]
             aq_conc * rho_f * volume * porosity
!unit       [mol/L] [L/m^3_water] [m^3_bulk] [m^3_pore/m^3_bulk]
    Jacobian(aq_id,im_id) = Jacobian(aq_id,im_id) - stoich_aq * tempreal
    Jacobian(im_id,im_id) = Jacobian(im_id,im_id) - stoich_im * tempreal

  end if    

  ! Residual - 1D array storing Residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entries in units kg water/sec
end subroutine ParticleEvaluate

! ************************************************************************** !

subroutine ParticleDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Bryan He
  ! Date: 06/07/2024
  !
  implicit none

  class(reaction_sandbox_particle_type) :: this

end subroutine ParticleDestroy

subroutine ReadAttachmentRate(this,input,word,option)
  !   
  ! Read attachment rate from input file
  !   either a constant value or through 
  !   the filteration model (Tufenkji & Elimelech, 2004)
  ! Author: Bryan He
  ! Date: 06/08/2024
  !  
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_particle_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,word, &
         'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT')
  select case(trim(word))
    case('CONSTANT')
      call ReadAttachmentRateConstant(this,input,word,option)
    case('FILTRATION_MODEL')
      call ReadAttachmentRateModel(this,input,word,option)
    case default
      call InputKeywordUnrecognized(input,word, &
             'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT,&
             &ATTACHMENT_RATE',option)
  end select

end subroutine ReadAttachmentRate

subroutine ReadAttachmentRateConstant(this,input,word,option)
  !
  ! Read attachment rate from input file
  !   as a constant value
  ! Author: Bryan He
  ! Date: 06/08/2024
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_particle_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'CONSTANT', &
           'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT&
           &,ATTACHMENT_RATE')
    call StringToUpper(word)

    select case(trim(word))
      case('VALUE')
      ! Read the double precision rate constant
        call InputReadDouble(input,option,this%attachment_rate)
        call InputErrorMsg(input,option,'VALUE', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT&
               &,ATTACHMENT_RATE,CONSTANT')
      ! Read the units
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
        ! If units do not exist, assume default units of 1/s which are the
        ! standard internal PFLOTRAN units for this rate constant.
          input%err_buf = 'RATE CONSTANT UNITS assumed as 1/s'
          call InputDefaultMsg(input,option)
        else
        ! If units exist, convert to internal units of 1/s
          internal_units = 'unitless/sec'
          this%attachment_rate = this%attachment_rate * &
          UnitsConvertToInternal(word,internal_units, &
            'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT&
            &,ATTACHMENT_RATE,CONSTANT',option)
        endif
      case default
        call InputKeywordUnrecognized(input,word, &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT,&
               &ATTACHMENT_RATE',option)
    end select
  end do
  call InputPopBlock(input,option)

end subroutine ReadAttachmentRateConstant

subroutine ReadAttachmentRateModel(this,input,word,option)
  !
  ! Read attachment rate from input file
  !   as a filtration model
  ! Author: Bryan He
  ! Date: 06/08/2024
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_particle_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  call InputPushBlock(input,option)
  ! Gotta turn this on so Darcy velocity is stored in global
  ! and can be used in the reaction sandbox
  option%flow%store_darcy_vel = PETSC_TRUE
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'FILTRATION_MODEL', &
           'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT,&
           &ATTACHMENT_RATE')
    call StringToUpper(word)

    select case(trim(word))
      case('DIAMETER_COLLECTOR')
      ! Diameter of the collector, i.e., soil grain size [m]
        call InputReadDouble(input,option,this%diam_collector)
        call InputErrorMsg(input,option,'DIAMETER_COLLECTOR', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT, &
               &FILTRATION_MODEL')

      case('DIAMETER_PARTICLE')
      ! Diameter of the particle [m]
        call InputReadDouble(input,option,this%diam_particle)
        call InputErrorMsg(input,option,'DIAMETER_PARTICLE', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT, &
               &FILTRATION_MODEL')

      case('HAMAKER_CONSTANT')
      ! Hamaker constant particle-soil pair [Joules]
        call InputReadDouble(input,option,this%hamaker_constant)
        call InputErrorMsg(input,option,'HAMAKER_CONSTANT', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT, &
               &FILTRATION_MODEL')

      case('DENSITY_PARTICLE')
      ! Density of the particles [kg/m3]
        call InputReadDouble(input,option,this%density_particle)
        call InputErrorMsg(input,option,'DENSITY_PARTICLE', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT, &
               &FILTRATION_MODEL')

      case('ALPHA_EFFICIENCY')
      ! Collision/attachment efficiency [-]
        call InputReadDouble(input,option,this%alpha_efficiency)
        call InputErrorMsg(input,option,'ALPHA_EFFICIENCY', &
               'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT, &
               &FILTRATION_MODEL')

      ! Something else
      case default
        call InputKeywordUnrecognized(input,word, &
              'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT,&
              &ATTACHMENT_RATE,FILTRATION_MODEL',option)
    end select
  end do
  call InputPopBlock(input,option)

end subroutine ReadAttachmentRateModel

subroutine ReadMaxSolidConc(this,input,word,option)
  !
  ! Read max solid concentration from input file
  ! Author: Bryan He
  ! Date: 06/08/2024
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  class(reaction_sandbox_particle_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

! Read Smax in Eqn 5.13 in mol/m^3 bulk
  call InputReadDouble(input,option,this%max_solid_conc)
  call InputErrorMsg(input,option,word, &
         'CHEMISTRY,REACTION_SANDBOX,PARTICLE_TRANSPORT,&
         &MAX_SOLID_CONCENTRATION')

end subroutine ReadMaxSolidConc

function AttachmentRate(this,global_auxvar,iphase,temperature,viscosity, &
                       porosity,rho_f)
  !
  ! Compute attachment rate based on Tufenkji & Elimelech, 2004
  ! Following lines 730~830 in reaction_sandbox_bioTH.F90 
  ! Author: Bryan He
  ! Date: 06/07/2024
  !
  use Reaction_Aux_module

  implicit none
  class(reaction_sandbox_particle_type) :: this
  type(global_auxvar_type) :: global_auxvar

  ! Filtration model parameters for attachment
  PetscReal :: AttachmentRate 
  PetscReal :: dc,dp,Hamaker,rho_p, alpha, qMag, diffusionCoeff, kB
  PetscReal :: Gm, Gm5, Happel
  PetscReal :: NR, NPe, NvdW, Ngr
  PetscReal :: Eta_D, Eta_I, Eta_G, Eta_0
  PetscReal :: temperature,viscosity,g,porosity,rho_f 
  PetscInt :: iphase

  dc = 0.0
  dp = 0.0
  Hamaker = 0.0
  rho_p = 0.0
  alpha = 0.0
  g = EARTH_GRAVITY

  if (this%attachment_rate < 0.d0) then
    kB = this%BOLTZMANN_CONSTANT
    dc = this%diam_collector
    dp = this%diam_particle
    rho_p = this%density_particle
    alpha = this%alpha_efficiency
    Hamaker = this%hamaker_constant

    qMag = MAX(global_auxvar%darcy_vel(iphase),1.0d-20)

    diffusionCoeff = this%BOLTZMANN_CONSTANT*(temperature+T273K) / &
                   (3.0 * PI * viscosity * dp)

    ! Non-dimensional parameters
    !! Happel parameter As
    Gm = (1.0 - porosity)**(1./3.)
    Gm5 = Gm*Gm*Gm*Gm*Gm
    Happel = (2.0 * (1.0 - Gm5)) / (2.0 - (3.0*Gm) + (3.0*Gm5) - (2.0*Gm*Gm5))

    !! Aspect ratio
    NR = dp/dc

    !! Péclet number
    NPe = (qMag * dc) / (diffusionCoeff)

    !! van der Waals number
    NvdW = Hamaker/(kB*(temperature + T273K))

    !! Gravitational number
    NGr = PI/12.0 * (dp*dp*dp*dp) * (rho_p - rho_f) * g /&
          (kB*(temperature + T273K))
!    print *, 'Happel,NR,NPe,NvdW,NGr,temperature:',Happel,NR,NPe,NvdW,NGr,temperature
    ! Collector efficiencies
    ! ( see Tufenkji & Elimelech 2014
    !   DOI : 10.1021/es034049r )
    !! Transport by diffusion
    Eta_D = 2.4 &
            * Happel**(1./3.) &
            * NR**(-0.081) &
            * NPe**(-0.715) &
            * NvdW**(0.052)

    !! Transport by interception
    Eta_I = 0.55 &
            * Happel &
            * NR**(1.55) &
            * NPe**(-0.125) &
            * NvdW**(0.125)

    !! Transport due to gravity
    Eta_G = 0.475 &
            * NR**(-1.35) &
            * NPe**(-1.11) &
            * NvdW**(0.053) &
            * NGr**(1.11)

    !! Single collector efficiency
    Eta_0 = Eta_D + Eta_I + Eta_G
!    print *, 'Eta_D,Eta_I,Eta_G:',Eta_D,Eta_I,Eta_G
    ! Rate of attachment according to CFT
    AttachmentRate = 1.5 * (1 - porosity) * qMag * alpha * Eta_0 &
           / (dc * porosity)    
  else
    ! A constant rate of attachment
    AttachmentRate = this%attachment_rate
  end if

end function

function DetachmentRate(this)
  !
  ! Compute detachment rate
  ! Following lines 835~840 in reaction_sandbox_bioth.F90
  ! Author: Bryan He
  ! Date: 06/07/2024
  !
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_particle_type) :: this
  PetscReal :: DetachmentRate

  if (this%detachment_rate < 0.d0) then
    DetachmentRate = 0.0
  else
    DetachmentRate = this%detachment_rate
  end if

end function

subroutine BlockingFunction(im_conc,im_max_conc,psi_b,dpsi_b_dc)
  !
  ! Compute blocking function
  ! Eqn 5.7 in Muller 2015
  !
  ! Author: Bryan He
  ! Date: 06/07/2024 created
  !       06/11/2024 modified following Glenn's suggestion
  
  implicit none

  PetscReal :: psi_b, dpsi_b_dc 
  PetscReal :: im_conc,im_max_conc 

  dpsi_b_dc = min(0.d0,-1.d0/im_max_conc)

  psi_b = 1.d0 + im_conc*dpsi_b_dc

end subroutine BlockingFunction

function DepthDepRetenFunction(d50,dx,beta)
  !
  ! Compute depth dependent retention function
  ! Eqn 5.8 in Muller 2015 
  !
  ! Author: Bryan He
  ! Date: 06/07/2024
  !

  implicit none

  PetscReal :: DepthDepRetenFunction
  PetscReal :: d50,dx,beta

  DepthDepRetenFunction = ((d50+dx)/d50)**(-beta)

end function

end module Reaction_Sandbox_Particle_class
