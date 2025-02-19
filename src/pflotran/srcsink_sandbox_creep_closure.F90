!this is based on the LOPOS equations in GRS-204 by Noseck et al. (2005)
!as interpreted for the DECOVALEX 2023/2027 Salt Performance Assessment Tasks
!LaForce et al SAND2023-04005R and Jayne et al (2024) SAND2024XXXXX
module SrcSink_Sandbox_Creep_Closure_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class

  implicit none

  private

  PetscInt, parameter, public :: CC_LIQ_PRESSURE_INDEX   = 1
  PetscInt, parameter, public :: CC_GAS_PRESSURE_INDEX   = 2
  PetscInt, parameter, public :: CC_LIQ_SATURATION_INDEX = 3
  PetscInt, parameter, public :: CC_GAS_SATURATION_INDEX = 4
  PetscInt, parameter, public :: CC_POROSITY_INDEX       = 5
  PetscInt, parameter, public :: CC_TEMPERATURE_INDEX    = 6
  PetscInt, parameter, public :: CC_SALT_DEN_INDEX    = 6

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_creep_closure_type
    PetscReal :: p_0 !hydrostatic pressure at fully developed creep [Pa]
    PetscReal :: p_atm !atmospheric pressure [Pa]
    PetscReal :: pgzr !lithostatic rock pressure at reference (repository level) [Pa]
    PetscReal :: pgz  !lithostatic pressure at z (assumed to be equal to pgzr)
    PetscReal :: n    !stress exponent for convergence
    PetscReal :: x_F  !parameter for moisture-enhanced creep
    PetscReal :: phi_0 !initial porosity
    PetscReal :: phi_r !reference porosity
    PetscReal :: phi_g !final (fully converged) porosity
    PetscReal :: h0    !backfill support function
    PetscReal :: h1    !backfill support function
    PetscReal :: h2    !backfill support function
    PetscReal :: h3    !backfill support function
    PetscReal :: g0    !backfill support function
    PetscReal :: g1    !backfill support function
    PetscReal :: g2_wet!backfill support function for wet salt
    PetscReal :: g2_dry!backfill support function for dry salt
    PetscReal :: m     !stress exponent for backfill support
    PetscReal :: L     !lambda parameter for pre-history of convergence
    PetscReal :: Ls    !local variation of the salt convergence rate
    PetscReal :: K_ref !reference convergence rate [1/s]
    PetscReal :: K_0   !initial convergence rate at reference depth[1/s]
    PetscReal :: sum_volume
    PetscReal :: a ! temperature parameters, weighting between activation energies
    PetscReal :: Q_1 ! activation energy at lower temperature [units?] put in Q1=Q2=0 for isothermal creep
    PetscReal :: Q_2 ! activation energy at high temperature  [units?]
    PetscReal :: z_ref !reference depth [m]
    PetscReal :: T_ref !temperature (C) at reference depth 
    PetscReal ::  delT !geothermal gradient [C/m?]
    PetscReal ::  R    !gas constant 
    
  contains
    procedure, public :: ReadInput => CreepClosureRead
    procedure, public :: Setup => CreepClosureSetup
    procedure, public :: Evaluate => CreepClosureSrcSink
    procedure, public :: Destroy => CreepClosureDestroy
  end type srcsink_sandbox_creep_closure_type

  public :: CreepClosureCreate

contains

! ************************************************************************** !

function CreepClosureCreate()
  !
  ! Allocates a creep closure src/sink object.
  !
  ! Author: Tara LaForce, David Fukuyama
  ! Date: 01/22/25

  implicit none

  class(srcsink_sandbox_creep_closure_type), pointer :: CreepClosureCreate

  allocate(CreepClosureCreate)
  call SSSandboxBaseInit(CreepClosureCreate)
    CreepClosureCreate%p_0    = 1.0d6  !1 MPa initialize to values in python script
    CreepClosureCreate%p_atm  = 1.01d5 !1 atm in Pa
    CreepClosureCreate%pgzr   = 18.3d6 !Pa
    CreepClosureCreate%pgz    = 18.3d6 !Pa
    CreepClosureCreate%n      = 5.d0   ![-]
    CreepClosureCreate%x_F    = 4.d0   ![-]
    CreepClosureCreate%phi_0  = 0.1d0  ![-]
    CreepClosureCreate%phi_r  = 0.3d0  ![-]
    CreepClosureCreate%phi_g  = 4.d-6  ![-]
    CreepClosureCreate%h0     = 1.d0   ![-]
    CreepClosureCreate%h1     =-2.d0   ![-]
    CreepClosureCreate%h2     = -(3d0+2d0*CreepClosureCreate%h1)
    CreepClosureCreate%h3     = CreepClosureCreate%h1 + 2.d0
    CreepClosureCreate%g0     = 1.d0   ![-]
    CreepClosureCreate%g1     =-1.d0   ![-]
    CreepClosureCreate%g2_wet = 1000.d0 ![-]
    CreepClosureCreate%g2_dry = 100.d0 ![-]
    CreepClosureCreate%m      = 5.d0   ![-]
    CreepClosureCreate%L      = 0.2d0 ! lambda [-]
    CreepClosureCreate%Ls     = 0.2d0 ! L salt [-]
    CreepClosureCreate%K_ref  = 3.169d-16 !1.d-8/yr in 1/sec
    CreepClosureCreate%K_0    = 3.169d-9  !0.1d0 /yr in 1/sec
    CreepClosureCreate%sum_volume = 0.d0
    CreepClosureCreate%a = 0.d0 !0 to start using only Q1
    CreepClosureCreate%Q_1 = 5.2335d2 !kJ/mol at 29C from Le Comte, Paul. "Creep in rock salt." The Journal of Geology 73.3 (1965): 469-484.
    CreepClosureCreate%Q_2 = 1.2560d2 !ksame source, but at 300C
    CreepClosureCreate%z_ref = 8.5d1 !repository depth is 850
    CreepClosureCreate%T_ref = 3.625d1 !850 m deep at geothermal gradient given, surface temnp 15C
    CreepClosureCreate%delT = 2.5d-2 !25c/km is usual
    CreepClosureCreate%R = 8.3144598d0 !R in [J/mol/K]


end function CreepClosureCreate

! ************************************************************************** !

subroutine CreepClosureRead(this,input,option)
  !
  ! Reads input deck for creep closure src/sink parameters
  !
  ! Author: Tara LaForce, David Fukuyama
  ! Date: 01/22/25

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(srcsink_sandbox_creep_closure_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units, word2
  PetscBool :: found
  PetscBool :: inhibit_term_read

  inhibit_term_read = PETSC_FALSE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'SOURCE_SINK_SANDBOX,CREEP_CLOSURE')
    call StringToUpper(word)

    call SSSandboxBaseSelectCase(this,input,option,word,found)
    if (found) cycle

    select case(trim(word))

      case('P_0')
        call InputReadDouble(input,option,this%p_0)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
      case('P_ATM')
        call InputReadDouble(input,option,this%p_atm)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
      case('PGZR')
        call InputReadDouble(input,option,this%pgzr)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
      case('PGZ')
        call InputReadDouble(input,option,this%pgz)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
      case('N')
        call InputReadDouble(input,option,this%n)
        call InputErrorMsg(input,option,word, &
                           'SOURCE_SINK_SANDBOX,PRESSURE,'//trim(word))
      case('X_F')
        call InputReadDouble(input,option,this%x_f)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('PHI_0')
        call InputReadDouble(input,option,this%phi_0)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('PHI_R')
        call InputReadDouble(input,option,this%phi_r)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('PHI_G')
        call InputReadDouble(input,option,this%phi_g)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('H0')
        call InputReadDouble(input,option,this%h0)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('H1')
        call InputReadDouble(input,option,this%h1)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
        this%h2 = -(3d0+2d0*this%h1)
        this%h3 = this%h1 + 2.d0
      case('G0')
        call InputReadDouble(input,option,this%h2)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('G1')
        call InputReadDouble(input,option,this%g1)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('G2_WET')
        call InputReadDouble(input,option,this%g2_wet)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('G2_DRY')
        call InputReadDouble(input,option,this%g2_dry)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('M')
        call InputReadDouble(input,option,this%m)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('L')
        call InputReadDouble(input,option,this%l)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('LS')
        call InputReadDouble(input,option,this%ls)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('K_REF')
        call InputReadDouble(input,option,this%k_ref)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('K_0')
        call InputReadDouble(input,option,this%k_0)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('A')
        call InputReadDouble(input,option,this%a)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('Q_1')
        call InputReadDouble(input,option,this%q_1)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('Q_2')
        call InputReadDouble(input,option,this%q_2)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('Z_REF')
        call InputReadDouble(input,option,this%z_ref)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('T_REF')
        call InputReadDouble(input,option,this%t_ref)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('DELT')
        call InputReadDouble(input,option,this%delt)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SRCSINK_SANDBOX,PRESSURE',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine CreepClosureRead

! ************************************************************************** !

subroutine CreepClosureSetup(this,grid,region_list,material_auxvars,option)
  !
  ! Sets up the creep closure src/sink
  !
  ! Author: Tara LaForce, David Fukuyama
  ! Date: 01/22/25

  use Option_module
  use Grid_module
  use Material_Aux_module, only: material_auxvar_type
  use Region_module

  implicit none

  class(srcsink_sandbox_creep_closure_type) :: this
  type(grid_type) :: grid
  type(region_list_type) :: region_list
  type(material_auxvar_type) :: material_auxvars(:)
  type(option_type) :: option

  PetscInt :: icell
  PetscErrorCode :: ierr

  call SSSandboxBaseSetup(this,grid,region_list,material_auxvars,option)
  this%sum_volume = 0.d0
  if (associated(this%local_cell_ids)) then
    do icell = 1, size(this%local_cell_ids)
      this%sum_volume = this%sum_volume + &
        material_auxvars(grid%nL2G(this%local_cell_ids(icell)))%volume
    enddo
  endif
  call MPI_Allreduce(MPI_IN_PLACE,this%sum_volume,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                     ierr);CHKERRQ(ierr)

end subroutine CreepClosureSetup

! ************************************************************************** !

subroutine CreepClosureSrcSink(this,Residual,Jacobian,compute_derivative, &
                               material_auxvar,aux_real,option)
  !
  ! Evaluates src/sink storing residual and/or Jacobian
  !
  ! Author: Tara LaForce, David Fukuyama
  ! Date: 01/22/25

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module
  use String_module
  use Utility_module

  implicit none

  class(srcsink_sandbox_creep_closure_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  PetscReal :: p_l, p_g    ! [Pa]
  PetscReal :: s_l, s_g    ! [-]
  PetscReal :: phi         ! [-]
  PetscReal :: temp        ! [C]
  PetscReal :: cell_volume ! [m^3]
  PetscReal :: den_s       ! [kmol/m^3]
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)
  PetscReal :: p_avg
  PetscReal :: f_F
  PetscReal :: f_F0
  PetscReal :: f_R
  PetscReal :: f_P
  PetscReal :: f_T
  PetscReal :: f_Z
  PetscReal :: Kt
  PetscReal :: dphi
  PetscReal :: new_salt


  ! here's some variables that may be useful to you.
  p_l   = aux_real(CC_LIQ_PRESSURE_INDEX)
  p_g   = aux_real(CC_GAS_PRESSURE_INDEX)
  s_l   = aux_real(CC_LIQ_SATURATION_INDEX)
  s_g   = aux_real(CC_GAS_SATURATION_INDEX)
  phi   = aux_real(CC_POROSITY_INDEX)
  temp  = aux_real(CC_TEMPERATURE_INDEX)
  den_s = aux_real(CC_SALT_DEN_INDEX)

  cell_volume = material_auxvar%volume
  ! to get the total moles of salt per gridcell:
  ! cell_volume [m^3] * phi [-] * den_s [kmol/m^3] = [kmol]


  !These cutoffs are arbitrary to prevent high Pc from making p_avg blow up
  IF (s_l > 5.0D-2 .AND. s_g > 5.0D-2) THEN
     p_avg = p_l * s_l + p_g * s_g
  ELSE IF (s_l >= 5.0D-2) THEN
     p_avg = p_l
  ELSE
     p_avg = p_g
  ENDIF

  call calc_f_R(this,p_avg,f_R)
  call calc_f_P(this,p_avg,f_P)
  call calc_f_T(this,temp,f_T)
  call calc_f_F(this,phi,s_g,s_l,f_F)
  call calc_f_F(this,this%phi_0,s_g,s_l,f_F0)
  call calc_f_Z(this,phi,f_F0,f_Z)
  
  call calc_Kt(this,f_R,f_P,f_F,f_T,f_Z,Kt)

  call calc_dphi(this,phi,Kt,dphi)

  !print *, 'f_R', f_R
  !print *, 'f_P', f_P
  !print *, 'f_T', f_T
  !print *, 'f_F', f_F
  !print *, 'f_Z', f_Z
  
  new_salt = (cell_volume * den_s) * dphi

  Residual(THREE_INTEGER) =   Residual(THREE_INTEGER) - new_salt
  ! Residual(THREE_INTEGER) = ___amount of mass to add to residual_____ units: [kmol/s]
!  if (compute_derivative) then
!    Jacobian(this%iphase,this%iphase) =
!  endif

end subroutine CreepClosureSrcSink

! ************************************************************************** !

subroutine calc_f_R(this,p_avg,f_R)

  ! Calculates f_R (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: p_avg
  PetscReal, intent(out) :: f_R
  PetscReal :: arg

  !Note: p_l is used in the D2023 task specification  
  arg = (p_avg - this%p_atm) / (this%p_0)
  f_R = 1.0 + this%x_F * arg
  

end subroutine calc_f_R

! ************************************************************************** !

subroutine calc_f_P(this,p_avg,f_P)

  ! Calculates f_P (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: p_avg
  PetscReal, intent(out) :: f_P
  PetscReal :: arg

  arg = (this%pgz - p_avg) / (this%pgzr)
  f_P = (arg)**this%n
  

end subroutine calc_f_P

! ************************************************************************** !

subroutine calc_f_T(this,temp,f_T)

  ! Calculates f_T (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: temp
  PetscReal, intent(out) :: f_T
  PetscReal :: arg1
  PetscReal :: inner
  PetscReal :: arg2
  PetscReal :: arg3

  arg1 = 1.0 / (1.0 + this%a)
  !notice that inner would be different in the first and second term if z != z_ref
  inner = 1.0 / this%T_ref - 1.0 / temp
  arg2 = EXP((this%Q_1 / this%R) * inner) 
  arg3 =  1.0 + this%a * EXP(( (this%Q_2 - this%Q_1) / this%R) * inner)

  !print *, 'temp', temp
  !print *, 'inner', inner
  !print *, 'arg1', arg1
  !print *, 'arg2', arg2
  !print *, 'arg3', arg3

  f_T = arg1 * arg2 * arg3
  

end subroutine calc_f_T

! ************************************************************************** !

subroutine calc_f_Z(this,phi,f_F0,f_Z)

  ! Calculates f_Z (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: phi
  PetscReal, intent(in) :: f_F0
  PetscReal, intent(out) :: f_Z
  PetscReal :: arg
  PetscReal :: argA

  !Note: eventually this may depend on ln(V/Vo)
  argA = this%L * ( (this%K_0) / (this%K_ref * f_F0) - 1.0) 
  arg = (argA) / (this%L + 0.0) 
  f_Z = 1.0 + arg
  

end subroutine calc_f_Z

! ************************************************************************** !

subroutine calc_f_F(this,phi,s_g,s_l,f_F)

  ! Calculates f_F (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: phi
  PetscReal, intent(in) :: s_g
  PetscReal, intent(in) :: s_l
  PetscReal, intent(out) :: f_F
  PetscReal :: arg, h, g, g2_avg

  g2_avg = this%g2_wet * s_l + this%g2_dry * s_g

  arg = (phi - this%phi_g) / (this%phi_r - this%phi_g)
  h   = this%h0 + this%h1 * arg + this%h2 * arg ** 2 + this%h3 * arg**3
  g   = this%g0 + this%g1 * arg + g2_avg * arg**2
  f_F = (1.0 + h / (((phi - this%phi_g) * g) **(1.0/this%m)))**(-this%m)


end subroutine calc_f_F

! ************************************************************************** !
subroutine calc_Kt(this,f_R,f_P,f_F,f_T,f_Z,Kt)

  ! Calculates K(t) (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: f_R
  PetscReal, intent(in) :: f_P
  PetscReal, intent(in) :: f_F
  PetscReal, intent(in) :: f_T
  PetscReal, intent(in) :: f_Z
  PetscReal, intent(out) :: Kt

  Kt = this%Ls * this%K_ref * f_R * f_P * f_F * f_T * f_Z
  

end subroutine calc_Kt

! ************************************************************************** !
subroutine calc_dphi(this,phi,Kt,dphi)

  ! Calculates K(t) (Eq. 2-3)

  implicit none
  class(srcsink_sandbox_creep_closure_type) :: this
  PetscReal, intent(in) :: Kt
  PetscReal, intent(in) :: phi
  PetscReal, intent(out) :: dphi

  if (phi > this%phi_g) then
     dphi = -Kt
  else
     dphi = 0.0
  endif
  
end subroutine calc_dphi

! ************************************************************************** !

subroutine CreepClosureDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Tara LaForce, David Fukuyama
  ! Date: 01/22/25

  use Utility_module

  implicit none

  class(srcsink_sandbox_creep_closure_type) :: this

  call SSSandboxBaseDestroy(this)

end subroutine CreepClosureDestroy

end module SrcSink_Sandbox_Creep_Closure_class
