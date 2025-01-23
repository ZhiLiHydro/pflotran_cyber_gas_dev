module SrcSink_Sandbox_Creep_Closure_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class

  implicit none

  private

  PetscInt, parameter, public :: SS_PRES_LIQUID_PRESSURE_OFFSET = 1
  PetscInt, parameter, public :: SS_PRES_GAS_PRESSURE_OFFSET = 2
  PetscInt, parameter, public :: SS_PRES_LIQUID_DENSITY_OFFSET = 3
  PetscInt, parameter, public :: SS_PRES_GAS_DENSITY_OFFSET = 4

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_creep_closure_type
    PetscReal :: p_0
    PetscReal :: p_01
    PetscReal :: p_atm
    PetscReal :: pgzr
    PetscReal :: pgz
    PetscReal :: n    !DF: should this be an int?
    PetscReal :: x_F  !    change to PetscInt as needed
    PetscReal :: c_L
    PetscReal :: phi_0
    PetscReal :: phi_r
    PetscReal :: phi_g
    PetscReal :: h0
    PetscReal :: h1
    PetscReal :: h2
    PetscReal :: h3
    PetscReal :: g0
    PetscReal :: g1
    PetscReal :: g2
    PetscReal :: m
    PetscReal :: L
    PetscReal :: K_ref
    PetscReal :: K_0
    PetscReal :: sum_volume
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
    CreepClosureCreate%p_0    = 0.101d0 ! initialize to values in python script
    CreepClosureCreate%p_01   = 1.d0
    CreepClosureCreate%p_atm  = 0.101d0
    CreepClosureCreate%pgzr   = 18.3d0
    CreepClosureCreate%pgz    = 18.3d0
    CreepClosureCreate%n      = 5.d0
    CreepClosureCreate%x_F    = 4.d0
    CreepClosureCreate%c_L    = 1.d0
    CreepClosureCreate%phi_0  = 0.1d0
    CreepClosureCreate%phi_r  = 0.3d0
    CreepClosureCreate%phi_g  = 4.d-6
    CreepClosureCreate%h0     = 1.d0
    CreepClosureCreate%h1     =-2.d0
    CreepClosureCreate%h2     = -(3d0+2d0*CreepClosureCreate%h1)
    CreepClosureCreate%h3     = CreepClosureCreate%h1 + 2.d0
    CreepClosureCreate%g0     = 1.d0
    CreepClosureCreate%g1     =-1.d0
    CreepClosureCreate%g2     = 100.d0
    CreepClosureCreate%m      = 5.d0
    CreepClosureCreate%L      = 0.2d0
    CreepClosureCreate%K_ref  = 1.d-8
    CreepClosureCreate%K_0    = 0.1d0
    CreepClosureCreate%sum_volume = 0.d0

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
      case('P_01')
        call InputReadDouble(input,option,this%p_01)
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
      case('C_L')
        call InputReadDouble(input,option,this%c_L)
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
      case('G0')
        call InputReadDouble(input,option,this%h2)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('G1')
        call InputReadDouble(input,option,this%g1)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('G2')
        call InputReadDouble(input,option,this%g2)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('M')
        call InputReadDouble(input,option,this%m)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('L')
        call InputReadDouble(input,option,this%l)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('K_REF')
        call InputReadDouble(input,option,this%k_ref)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,PRESSURE,')
      case('K_0')
        call InputReadDouble(input,option,this%k_0)
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
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)

  print *, 'aux_real: ', aux_real(:)

  if (compute_derivative) then
!    Jacobian(this%iphase,this%iphase) =
  endif

end subroutine CreepClosureSrcSink

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
