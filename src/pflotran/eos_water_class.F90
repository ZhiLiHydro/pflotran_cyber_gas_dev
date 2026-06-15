module EOS_Water_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use EOS_Base_class
  use EOS_Water_module

  implicit none

  private

  PetscInt, parameter, public :: EOS_WATER_DENSITY_CONSTANT = 1
  PetscInt, parameter, public :: EOS_WATER_DENSITY_IFC67    = 2
  PetscInt, parameter, public :: EOS_WATER_DENSITY_IF97     = 3

  PetscInt, parameter, public :: EOS_WATER_VISCOSITY_CONSTANT = 1
  PetscInt, parameter, public :: EOS_WATER_VISCOSITY_DEFAULT  = 2

  character(len=MAXWORDLENGTH), save :: pending_density_keyword = ''
  PetscReal, save :: pending_density_aux(10) = 0.d0
  PetscBool, save :: pending_density_configured = PETSC_FALSE

  character(len=MAXWORDLENGTH), save :: pending_viscosity_keyword = ''
  PetscReal, save :: pending_viscosity_aux(10) = 0.d0
  PetscBool, save :: pending_viscosity_configured = PETSC_FALSE

  ! Match EOS_Water_module kernel signatures for procedure pointers.
  interface
    subroutine EOSWaterClassDensityKernel(t, p, calculate_derivatives, &
                                          dw, dwmol, dwp, dwt, ierr, &
                                          table_idxs)
      use petscsys
      implicit none
      PetscReal, intent(in) :: t, p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw, dwmol, dwp, dwt
      PetscErrorCode, intent(inout) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSWaterClassDensityKernel
    subroutine EOSWaterClassViscosityKernel(t, p, ps, dps_dt, &
                                            calculate_derivatives, vw, &
                                            dvw_dt, dvw_dp, ierr, &
                                            table_idxs)
      use petscsys
      implicit none
      PetscReal, intent(in) :: t, p, ps, dps_dt
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: vw
      PetscReal, intent(out) :: dvw_dt, dvw_dp
      PetscErrorCode, intent(inout) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSWaterClassViscosityKernel
    subroutine EOSWaterClassSatPresKernel(t, calculate_derivatives, &
                                          ps, dps_dt, ierr)
      use petscsys
      implicit none
      PetscReal, intent(in) :: t
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: ps, dps_dt
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterClassSatPresKernel
  end interface

  type, public, extends(eos_base_type) :: eos_water_type
    PetscInt  :: density_model = EOS_WATER_DENSITY_IFC67
    PetscInt  :: viscosity_model = EOS_WATER_VISCOSITY_DEFAULT
    PetscReal :: constant_density = 1000.d0
    PetscReal :: constant_viscosity = 1.d-3
    procedure(EOSWaterClassDensityKernel), pointer, nopass :: &
      DensityEval => null()
    procedure(EOSWaterClassViscosityKernel), pointer, nopass :: &
      ViscosityEval => null()
    procedure(EOSWaterClassSatPresKernel), pointer, nopass :: &
      SatPresEval => null()
  contains
    procedure, public :: Init            => EOSWaterClassInit
    procedure, public :: Density        => EOSWaterClassDensity
    procedure, public :: Viscosity      => EOSWaterClassViscosity
    procedure, public :: SetDensityModel => EOSWaterSetDensityModel
    procedure, public :: SetViscosityModel => EOSWaterSetViscosityModel
  end type eos_water_type

  public :: EOSWaterCreate, &
            EOSWaterClassInit, &
            EOSWaterClassSetDensity, &
            EOSWaterClassSetViscosity, &
            EOSWaterApplyPendingConfig, &
            EOSWaterSetDensityModel, &
            EOSWaterSetViscosityModel

contains

! ************************************************************************** !

subroutine EOSWaterClassInit(this)
  !
  ! Initializes a water EOS object
  !
  ! Author: Glenn Hammond
  ! Date: 06/15/26
  !
  use EOS_Base_class, only : EOSBaseInit

  implicit none

  class(eos_water_type) :: this

  call EOSBaseInit(this)

  this%density_model = EOS_WATER_DENSITY_IFC67
  this%viscosity_model = EOS_WATER_VISCOSITY_DEFAULT
  this%constant_density = 1000.d0
  this%constant_viscosity = 1.d-3
  this%DensityEval => EOSWaterDensityIFC67
  this%ViscosityEval => EOSWaterViscosity1
  this%SatPresEval => EOSWaterSaturationPressureIFC67

end subroutine EOSWaterClassInit

! ************************************************************************** !

function EOSWaterCreate() result(water)
  !
  ! Creates a water EOS object
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_water_type), pointer :: water

  allocate(water)
  call water%Init()
  call EOSWaterApplyPendingConfig(water)

end function EOSWaterCreate

! ************************************************************************** !

subroutine EOSWaterApplyPendingConfig(water)
  !
  ! Applies EOS input read before the instance was allocated
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_water_type) :: water

  if (pending_density_configured) then
    call water%SetDensityModel(pending_density_keyword, &
                               pending_density_aux)
    pending_density_configured = PETSC_FALSE
    pending_density_keyword = ''
  endif
  if (pending_viscosity_configured) then
    call water%SetViscosityModel(pending_viscosity_keyword, &
                                 pending_viscosity_aux)
    pending_viscosity_configured = PETSC_FALSE
    pending_viscosity_keyword = ''
  endif

end subroutine EOSWaterApplyPendingConfig

! ************************************************************************** !

subroutine EOSWaterClassSetDensity(keyword, aux)
  !
  ! Stores water density model input for later instance configuration
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  character(len=*) :: keyword
  PetscReal :: aux(10)

  pending_density_keyword = keyword
  pending_density_aux(1:10) = aux(1:10)
  pending_density_configured = PETSC_TRUE

end subroutine EOSWaterClassSetDensity

! ************************************************************************** !

subroutine EOSWaterClassSetViscosity(keyword, aux)
  !
  ! Stores water viscosity model input for later instance configuration
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  character(len=*) :: keyword
  PetscReal :: aux(10)

  pending_viscosity_keyword = keyword
  pending_viscosity_aux(1:10) = aux(1:10)
  pending_viscosity_configured = PETSC_TRUE

end subroutine EOSWaterClassSetViscosity

! ************************************************************************** !

subroutine EOSWaterSetDensityModel(this, keyword, aux)
  implicit none

  class(eos_water_type) :: this
  character(len=*) :: keyword
  PetscReal :: aux(10)

  character(len=MAXWORDLENGTH) :: word

  word = trim(keyword)

  select case(trim(word))
    case('CONSTANT')
      this%density_model = EOS_WATER_DENSITY_CONSTANT
      this%constant_density = aux(1)
      constant_density = aux(1)
      this%DensityEval => EOSWaterDensityConstant
    case('IF97')
      this%density_model = EOS_WATER_DENSITY_IF97
      this%DensityEval => EOSWaterDensityIF97
    case('DEFAULT','IFC67')
      this%density_model = EOS_WATER_DENSITY_IFC67
      this%DensityEval => EOSWaterDensityIFC67
    case default
      this%density_model = EOS_WATER_DENSITY_IFC67
      this%DensityEval => EOSWaterDensityIFC67
  end select

end subroutine EOSWaterSetDensityModel

! ************************************************************************** !

subroutine EOSWaterSetViscosityModel(this, keyword, aux)
  !
  ! Configures water viscosity model on an EOS instance
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_water_type) :: this
  character(len=*) :: keyword
  PetscReal :: aux(10)

  character(len=MAXWORDLENGTH) :: word

  word = trim(keyword)

  select case(trim(word))
    case('CONSTANT')
      this%viscosity_model = EOS_WATER_VISCOSITY_CONSTANT
      this%constant_viscosity = aux(1)
      constant_viscosity = aux(1)
      this%ViscosityEval => EOSWaterViscosityConstant
    case('DEFAULT')
      this%viscosity_model = EOS_WATER_VISCOSITY_DEFAULT
      this%ViscosityEval => EOSWaterViscosity1
    case default
      this%viscosity_model = EOS_WATER_VISCOSITY_DEFAULT
      this%ViscosityEval => EOSWaterViscosity1
  end select

end subroutine EOSWaterSetViscosityModel

! ************************************************************************** !

function EOSWaterClassDensity(this) result(dw)
  !
  ! Computes water density via the configured EOS_Water_module kernel
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_water_type) :: this
  PetscReal :: dw

  PetscReal :: t, p, dwmol, dwp, dwt
  PetscErrorCode :: ierr

  t = this%state%temperature
  p = this%state%pressure
  ierr = 0

  if (.not.this%calculate_derivatives) then
    call this%ClearDerivatives()
    call this%DensityEval(t, p, PETSC_FALSE, dw, dwmol, dwp, dwt, ierr)
  else
    call this%DensityEval(t, p, PETSC_TRUE, dw, dwmol, &
                           this%density_dP, this%density_dT, ierr)
  endif

end function EOSWaterClassDensity

! ************************************************************************** !

function EOSWaterClassViscosity(this) result(vis)
  !
  ! Computes water viscosity via the configured EOS_Water_module kernel
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/26
  !
  implicit none

  class(eos_water_type) :: this
  PetscReal :: vis

  PetscReal :: t, p, ps, dps_dt, dvis_dt, dvis_dp
  PetscErrorCode :: ierr

  t = this%state%temperature
  p = this%state%pressure
  ierr = 0

  if (.not.this%calculate_derivatives) then
    call this%ClearDerivatives()
    call this%SatPresEval(t, PETSC_FALSE, ps, dps_dt, ierr)
    call this%ViscosityEval(t, p, ps, dps_dt, PETSC_FALSE, vis, &
                             dvis_dt, dvis_dp, ierr)
  else
    call this%SatPresEval(t, PETSC_TRUE, ps, dps_dt, ierr)
    call this%ViscosityEval(t, p, ps, dps_dt, PETSC_TRUE, vis, &
                             this%viscosity_dT, this%viscosity_dP, ierr)
  endif

end function EOSWaterClassViscosity

end module EOS_Water_class