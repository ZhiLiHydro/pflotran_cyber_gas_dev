module EOS_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: eos_state_type
    PetscReal :: pressure      = UNINITIALIZED_DOUBLE
    PetscReal :: temperature   = UNINITIALIZED_DOUBLE
    PetscReal :: concentration = UNINITIALIZED_DOUBLE
    PetscReal :: velocity      = UNINITIALIZED_DOUBLE
    PetscBool :: has_pressure      = PETSC_FALSE
    PetscBool :: has_temperature   = PETSC_FALSE
    PetscBool :: has_concentration = PETSC_FALSE
    PetscBool :: has_velocity      = PETSC_FALSE
  end type eos_state_type

  type, abstract, public :: eos_base_type
    type(eos_state_type) :: state
    PetscBool :: calculate_derivatives = PETSC_FALSE
    PetscReal :: density_dP   = UNINITIALIZED_DOUBLE
    PetscReal :: density_dT   = UNINITIALIZED_DOUBLE
    PetscReal :: viscosity_dP = UNINITIALIZED_DOUBLE
    PetscReal :: viscosity_dT = UNINITIALIZED_DOUBLE
  contains
    procedure, public :: LoadState    => EOSBaseLoadState
    procedure, public :: LoadPT       => EOSBaseLoadPT
    procedure, public :: LoadPTC      => EOSBaseLoadPTC
    procedure, public :: LoadPTCV     => EOSBaseLoadPTCV
    procedure, public :: LoadT        => EOSBaseLoadT
    procedure, public :: LoadVelocity => EOSBaseLoadVelocity
    procedure, public :: Reset        => EOSBaseReset
    procedure, public :: Density      => EOSBaseDensity
    procedure, public :: Viscosity    => EOSBaseViscosity
    procedure, public :: GetDwDp      => EOSBaseGetDwDp
    procedure, public :: GetDwDt      => EOSBaseGetDwDt
    procedure, public :: GetDvisDp    => EOSBaseGetDvisDp
    procedure, public :: GetDvisDt    => EOSBaseGetDvisDt
    procedure, public :: Init         => EOSBaseInit
    procedure, public :: Verify       => EOSBaseVerify
    procedure, public :: Strip        => EOSBaseStrip
    procedure, public :: ClearDerivatives => EOSBaseClearDerivatives
  end type eos_base_type

  public :: EOSBaseInit

contains

! ************************************************************************** !

subroutine EOSBaseInit(this)
  !
  ! Initializes an EOS base object
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_base_type) :: this

  this%state%pressure      = UNINITIALIZED_DOUBLE
  this%state%temperature   = UNINITIALIZED_DOUBLE
  this%state%concentration = UNINITIALIZED_DOUBLE
  this%state%velocity      = UNINITIALIZED_DOUBLE

  this%state%has_pressure      = PETSC_FALSE
  this%state%has_temperature   = PETSC_FALSE
  this%state%has_concentration = PETSC_FALSE
  this%state%has_velocity      = PETSC_FALSE

  this%calculate_derivatives = PETSC_FALSE
  call this%ClearDerivatives()

end subroutine EOSBaseInit

! ************************************************************************** !

subroutine EOSBaseClearDerivatives(this)
  !
  ! Resets cached property derivatives to UNINITIALIZED_DOUBLE
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_base_type) :: this

  this%density_dP   = UNINITIALIZED_DOUBLE
  this%density_dT   = UNINITIALIZED_DOUBLE
  this%viscosity_dP = UNINITIALIZED_DOUBLE
  this%viscosity_dT = UNINITIALIZED_DOUBLE

end subroutine EOSBaseClearDerivatives

! ************************************************************************** !

subroutine EOSBaseReset(this)
  implicit none

  class(eos_base_type) :: this

  call this%Init()

end subroutine EOSBaseReset

! ************************************************************************** !

subroutine EOSBaseLoadState(this, state)
  !
  ! Loads a full EOS state and updates validity flags
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_base_type) :: this
  type(eos_state_type) :: state

  this%state = state
  this%state%has_pressure = Initialized(state%pressure)
  this%state%has_temperature = Initialized(state%temperature)
  this%state%has_concentration = Initialized(state%concentration)
  this%state%has_velocity = Initialized(state%velocity)

end subroutine EOSBaseLoadState

! ************************************************************************** !

subroutine EOSBaseLoadPT(this, p, t)
  implicit none

  class(eos_base_type) :: this
  PetscReal, intent(in) :: p, t

  this%state%pressure = p
  this%state%temperature = t
  this%state%has_pressure = PETSC_TRUE
  this%state%has_temperature = PETSC_TRUE

end subroutine EOSBaseLoadPT

! ************************************************************************** !

subroutine EOSBaseLoadPTC(this, p, t, c)
  implicit none

  class(eos_base_type) :: this
  PetscReal, intent(in) :: p, t, c

  call this%LoadPT(p, t)
  this%state%concentration = c
  this%state%has_concentration = PETSC_TRUE

end subroutine EOSBaseLoadPTC

! ************************************************************************** !

subroutine EOSBaseLoadPTCV(this, p, t, c, v)
  implicit none

  class(eos_base_type) :: this
  PetscReal, intent(in) :: p, t, c, v

  call this%LoadPTC(p, t, c)
  this%state%velocity = v
  this%state%has_velocity = PETSC_TRUE

end subroutine EOSBaseLoadPTCV

! ************************************************************************** !

subroutine EOSBaseLoadT(this, t)
  implicit none

  class(eos_base_type) :: this
  PetscReal, intent(in) :: t

  this%state%temperature = t
  this%state%has_temperature = PETSC_TRUE

end subroutine EOSBaseLoadT

! ************************************************************************** !

subroutine EOSBaseLoadVelocity(this, v)
  implicit none

  class(eos_base_type) :: this
  PetscReal, intent(in) :: v

  this%state%velocity = v
  this%state%has_velocity = PETSC_TRUE

end subroutine EOSBaseLoadVelocity

! ************************************************************************** !

function EOSBaseDensity(this) result(dw)
  implicit none

  class(eos_base_type) :: this
  PetscReal :: dw

  dw = UNINITIALIZED_DOUBLE

end function EOSBaseDensity

! ************************************************************************** !

function EOSBaseViscosity(this) result(vis)
  implicit none

  class(eos_base_type) :: this
  PetscReal :: vis

  vis = UNINITIALIZED_DOUBLE

end function EOSBaseViscosity

! ************************************************************************** !

function EOSBaseGetDwDp(this) result(dwp)
  implicit none

  class(eos_base_type) :: this
  PetscReal :: dwp

  dwp = this%density_dP

end function EOSBaseGetDwDp

! ************************************************************************** !

function EOSBaseGetDwDt(this) result(dwt)
  implicit none

  class(eos_base_type) :: this
  PetscReal :: dwt

  dwt = this%density_dT

end function EOSBaseGetDwDt

! ************************************************************************** !

function EOSBaseGetDvisDp(this) result(dvisp)
  implicit none

  class(eos_base_type) :: this
  PetscReal :: dvisp

  dvisp = this%viscosity_dP

end function EOSBaseGetDvisDp

! ************************************************************************** !

function EOSBaseGetDvisDt(this) result(dvist)
  implicit none

  class(eos_base_type) :: this
  PetscReal :: dvist

  dvist = this%viscosity_dT

end function EOSBaseGetDvisDt

! ************************************************************************** !

subroutine EOSBaseVerify(this, name, option)
  use Option_module

  implicit none

  class(eos_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

end subroutine EOSBaseVerify

! ************************************************************************** !

subroutine EOSBaseStrip(this)
  implicit none

  class(eos_base_type) :: this

end subroutine EOSBaseStrip

end module EOS_Base_class