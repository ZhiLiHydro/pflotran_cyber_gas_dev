module Reaction_Sandbox_Oregonator_class

#include "petsc/finclude/petscsys.h"
  use petscsys

! 1. Change all references to "Example" as desired to rename the module and
!    and subroutines within the module.

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.,
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_oregonator_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: species_X
    character(len=MAXWORDLENGTH) :: species_Y
    character(len=MAXWORDLENGTH) :: species_Z
    character(len=MAXWORDLENGTH) :: species_A
    character(len=MAXWORDLENGTH) :: species_B
    character(len=MAXWORDLENGTH) :: species_P
    PetscInt :: i_A
    PetscInt :: i_B
    PetscInt :: i_P
    PetscInt :: i_X
    PetscInt :: i_Y
    PetscInt :: i_Z

    PetscReal :: ref_temp
    PetscReal :: temp_slope
    
    PetscReal :: k_i,k_ii,k_iii,k_iv,k_v
  contains
    procedure, public :: ReadInput => OregonatorRead
    procedure, public :: Setup => OregonatorSetup
    procedure, public :: Evaluate => OregonatorEvaluate
    procedure, public :: Destroy => OregonatorDestroy
  end type reaction_sandbox_oregonator_type

  public :: OregonatorCreate

contains

! ************************************************************************** !

function OregonatorCreate()
  !
  ! Allocates oregonator reaction object.
  !
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  !

  implicit none

  class(reaction_sandbox_oregonator_type), pointer :: OregonatorCreate

! 4. Add code to allocate the object, initialize all variables to zero and
!    nullify all pointers. E.g.,
  allocate(OregonatorCreate)
  OregonatorCreate%species_x = ''
  OregonatorCreate%species_y = ''
  OregonatorCreate%species_z = ''
  OregonatorCreate%species_a = ''
  OregonatorCreate%species_b = ''
  OregonatorCreate%species_p = ''

  OregonatorCreate%i_x = UNINITIALIZED_INTEGER
  OregonatorCreate%i_y = UNINITIALIZED_INTEGER
  OregonatorCreate%i_z = UNINITIALIZED_INTEGER
  OregonatorCreate%i_a = UNINITIALIZED_INTEGER
  OregonatorCreate%i_b = UNINITIALIZED_INTEGER
  OregonatorCreate%i_p = UNINITIALIZED_INTEGER
  
  OregonatorCreate%k_i   = UNINITIALIZED_DOUBLE
  OregonatorCreate%k_ii  = UNINITIALIZED_DOUBLE
  OregonatorCreate%k_iii = UNINITIALIZED_DOUBLE
  OregonatorCreate%k_iv  = UNINITIALIZED_DOUBLE
  OregonatorCreate%k_v   = UNINITIALIZED_DOUBLE

  OregonatorCreate%ref_temp = UNINITIALIZED_DOUBLE
  Oregonatorcreate%temp_slope = UNINITIALIZED_DOUBLE
  nullify(OregonatorCreate%next)

end function OregonatorCreate

! ************************************************************************** !

subroutine OregonatorRead(this,input,option)
  !
  ! Reads input deck for oregonator reaction parameters (if any)
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_oregonator_type) :: this
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
                       'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
    call StringToUpper(word)

    select case(trim(word))
      case('SPECIES_A')
        call InputReadWord(input,option,this%species_a,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('SPECIES_B')
        call InputReadWord(input,option,this%species_b,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('SPECIES_P')
        call InputReadWord(input,option,this%species_p,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('SPECIES_X')
        call InputReadWord(input,option,this%species_x,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('SPECIES_Y')
        call InputReadWord(input,option,this%species_y,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('SPECIES_Z')
        call InputReadWord(input,option,this%species_z,PETSC_TRUE)
        call InputErrorMsg(input,option,'SPECIES_NAME', &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('K_I')
        call InputReadDouble(input,option,this%k_i)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%k_i, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &OREGONATOR,RATE_CONSTANT',option)
      case('K_II')
        call InputReadDouble(input,option,this%k_ii)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%k_ii, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &OREGONATOR,RATE_CONSTANT',option)
      case('K_III')
        call InputReadDouble(input,option,this%k_iii)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%k_iii, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &OREGONATOR,RATE_CONSTANT',option)
      case('K_IV')
        call InputReadDouble(input,option,this%k_iv)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%k_iv, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &OREGONATOR,RATE_CONSTANT',option)
      case('K_V')
        call InputReadDouble(input,option,this%k_v)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
        internal_units = 'unitless/sec'
        call InputReadAndConvertUnits(input,this%k_v, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &OREGONATOR,RATE_CONSTANT',option)
      case('REFERENCE_TEMP')
        call InputReadDouble(input,option,this%ref_temp)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')
      case('LINEAR_TEMPERATURE_DEPENDENCE')
        call InputReadDouble(input,option,this%temp_slope)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,OREGONATOR')

      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,OREGONATOR',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine OregonatorRead

! ************************************************************************** !

subroutine OregonatorSetup(this,reaction,option)
  !
  ! Sets up the oregonator reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Reaction_Aux_module
  use Option_module

  implicit none

  class(reaction_sandbox_oregonator_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize.
  this%i_a = &
    ReactionAuxGetPriSpecIDFromName(this%species_a,reaction,option)
  this%i_b = &
    ReactionAuxGetPriSpecIDFromName(this%species_b,reaction,option)
  this%i_p = &
    ReactionAuxGetPriSpecIDFromName(this%species_p,reaction,option)
  this%i_x = &
    ReactionAuxGetPriSpecIDFromName(this%species_x,reaction,option)
  this%i_y = &
    ReactionAuxGetPriSpecIDFromName(this%species_y,reaction,option)
  this%i_z = &
    ReactionAuxGetPriSpecIDFromName(this%species_z,reaction,option)
end subroutine OregonatorSetup

! ************************************************************************** !

subroutine OregonatorEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_oregonator_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscReal :: rx1_stoichA = -1.d0
  PetscReal :: rx1_stoichY = -1.d0
  PetscReal :: rx1_stoichX =  1.d0
  PetscReal :: rx1_stoichP =  1.d0

  PetscReal :: rx2_stoichX = -1.d0
  PetscReal :: rx2_stoichY = -1.d0
  PetscReal :: rx2_stoichP =  2.d0

  PetscReal :: rx3_stoichA = -1.d0
  PetscReal :: rx3_stoichX = -1.d0
  PetscReal :: rx3_stoichXp = -2.d0
  PetscReal :: rx3_stoichZ = -2.d0

  PetscReal :: rx4_stoichX = -2.d0
  PetscReal :: rx4_stoichA =  1.d0
  PetscReal :: rx4_stoichP =  1.d0

  PetscReal :: rx5_stoichB = -1.d0
  PetscReal :: rx5_stoichZ = -1.d0
  PetscReal :: rx5_stoichY =  0.5d0

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water

  PetscReal :: A, B, P, X, Y, Z
  PetscReal :: RateA, RateB, RateP, RateX, RateY, RateZ
  PetscReal :: molality_to_molarity

  PetscReal :: dT

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  molality_to_molarity = global_auxvar%den_kg(iphase) * 1.d-3

  A = rt_auxvar%pri_molal(this%i_A)*molality_to_molarity !mol/L
  B = rt_auxvar%pri_molal(this%i_B)*molality_to_molarity !mol/L
  P = rt_auxvar%pri_molal(this%i_P)*molality_to_molarity !mol/L
  X = rt_auxvar%pri_molal(this%i_X)*molality_to_molarity !mol/L
  Y = rt_auxvar%pri_molal(this%i_Y)*molality_to_molarity !mol/L
  Z = rt_auxvar%pri_molal(this%i_Z)*molality_to_molarity !mol/L

  if (Initialized(this%temp_slope) .and. Initialized(this%ref_temp)) then
    dT = (global_auxvar%temp - this%ref_temp)
    print *, 'dT, ref, slope', dT, this%ref_temp, this%temp_slope
    this%k_i   = this%k_i   + (dT * this%temp_slope)
    this%k_ii  = this%k_ii  + (dT * this%temp_slope)
    this%k_iii = this%k_iii + (dT * this%temp_slope)
    this%k_iv  = this%k_iv  + (dT * this%temp_slope)
    this%k_v   = this%k_v   + (dT * this%temp_slope)
  endif

  RateX = this%k_i   * A * Y      -   this%k_ii * X * Y + &
          this%k_iii * A * X      - 2*this%k_iv * X ** 2
  RateY = -1d0*this%k_i * A * Y   - this%k_ii * X * Y + &
          0.5d0*this%k_v * B * Z
  RateZ =  2d0*this%k_iii * A * X - this%k_v * B * Z
  RateA = -this%k_i * A * Y - this%k_iii * A * X + this%k_iv * X ** 2
  RateB = -this%k_v * B * Z
  RateP = this%k_i * A * Y + this%k_ii * X * Y + this%k_iv * X ** 2

  Residual(this%i_A) = Residual(this%i_A) - RateA * L_water
  Residual(this%i_B) = Residual(this%i_B) - RateB * L_water
  Residual(this%i_P) = Residual(this%i_P) - RateP * L_water
  Residual(this%i_X) = Residual(this%i_X) - RateX * L_water
  Residual(this%i_Y) = Residual(this%i_Y) - RateY * L_water
  Residual(this%i_Z) = Residual(this%i_Z) - RateZ * L_water

  ! Description of subroutine arguments:

  ! Residual - 1D array storing Residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entries in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analytical derivatives will
  !   be calculated.  The user must provide either the analytical derivatives
  !   or a numerical approximation unless always running with
  !   NUMERICAL_JACOBIAN defined within the NUMERICAL_METHODS TRANSPORT,
  !   NEWTON_SOLVER block of the input deck.  If the use of
  !   NUMERICAL_JACOBIAN is assumed, the user should provide an error
  !   message when compute_derivative is true.  E.g.,
  !
  !   if (compute_derivative) then
  !     option%io_buffer = 'NUMERICAL_JACOBIAN must be specified within &
  !       &the NEWTON_SOLVER block of NUMERICAL_METHODS TRANSPORT due to &
  !       &assumptions made in OregonatorEvaluate.'
  !     call PrintErrMsg(option)
  !   endif
  !
  ! rt_auxvar - Object holding chemistry information (e.g., concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g., saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - liquid density [mol/m^3]
  !     global_auxvar%den_kg(iphase) - liquid density [kg/m^3]
  !     global_auxvar%sat(iphase) - liquid saturation [m^3 water/m^3 pore]
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.,
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.

! 10. Add code for the Residual evaluation.

  ! Units of the Residual must be in moles/second.
  ! 1.d3 converts m^3 water -> L water

  ! Always "subtract" the contribution from the Residual.
!   Residual(this%species_id) = Residual(this%species_id) - &
!     (-1.d0) * & ! negative stoichiometry
!     this%rate_constant * &  ! 1/sec
!     L_water * & ! L water
!     rt_auxvar%total(this%species_id,iphase) ! mol/L water

!   if (compute_derivative) then

! ! 11. If using an analytical Jacobian, add code for the Jacobian evaluation.

!     ! Always "add" the contribution to the Jacobian.
!     ! Units = (mol/sec)*(kg water/mol) = kg water/sec
!     Jacobian(this%species_id,this%species_id) = &
!     Jacobian(this%species_id,this%species_id) - &
!       (-1.d0) * & ! negative stoichiometry
!       this%rate_constant * & ! 1/sec
!       L_water * & ! L water
!       ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) =
!       !   derivative of total component concentration with respect to the
!       !   free ion concentration of the same species.
!       ! kg water/L water
!       rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase)

!   endif

end subroutine OregonatorEvaluate

! ************************************************************************** !

subroutine OregonatorDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  implicit none

  class(reaction_sandbox_oregonator_type) :: this

! 12. Add code to deallocate dynamic members of reaction_sandbox_oregonator_type.

end subroutine OregonatorDestroy

end module Reaction_Sandbox_Oregonator_class
