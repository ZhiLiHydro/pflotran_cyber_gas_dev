module PM_Stokes_TS_class

#include "petsc/finclude/petscts.h"
  use petscts
  use Stokes_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_Stokes_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_stokes_type) :: pm_stokes_ts_type
  contains
    procedure, public :: IFunction => PMStokesTSIFunction
    procedure, public :: IJacobian => PMStokesTSIJacobian
    procedure, public :: RHSFunction => PMStokesTSRHSFunction
    procedure, public :: InitializeTimestep => PMStokesTSInitializeTimestep
  end type pm_stokes_ts_type

  public :: PMStokesTSCreate

contains

! ************************************************************************** !

function PMStokesTSCreate()
  !
  ! Creates Stokes TS process model shell
  !
  ! Author: Satish Karra
  ! Date: 02/06/26
  !

  implicit none

  class(pm_stokes_ts_type), pointer :: PMStokesTSCreate
  class(pm_stokes_ts_type), pointer :: stokes_ts_pm

  allocate(stokes_ts_pm)
  call PMStokesInit(stokes_ts_pm)
  stokes_ts_pm%name = 'Stokes TS Flow'
  stokes_ts_pm%header = 'STOKES TS FLOW'

  PMStokesTSCreate => stokes_ts_pm

end function PMStokesTSCreate

! ************************************************************************** !

subroutine PMStokesTSIFunction(this,ts,time,U,Udot,F,ierr)
  !
  ! TS IFunction for Stokes IMEX split
  !
  ! Author: Satish Karra
  ! Date: 02/06/26
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none

  class(pm_stokes_ts_type) :: this

  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  Vec :: F
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  call VecZeroEntries(F,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,U,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationGlobalToLocal(discretization,Udot,field%flow_xxdot_loc, &
                                   NFLOWDOF)

  call StokesUpdateAuxVars(realization)

  call StokesResidualInternalConn(F,realization,this%rho,this%mu, &
                                  PETSC_FALSE,ierr)
  call StokesResidualBoundaryConn(F,realization,this%rho,this%mu,ierr)
  call StokesResidualSourceSink(F,realization,this%rho,ierr)
  call StokesIFunctionAccumulation(F,realization,this%rho,ierr)
  call StokesResidualPinPressure(F,U,realization,0.d0,ierr)
  call StokesResidualBoundaryFinal(F,U,realization,ierr)

end subroutine PMStokesTSIFunction

! ************************************************************************** !

subroutine PMStokesTSRHSFunction(this,ts,time,xx,ff,ierr)
  !
  ! TS RHSFunction for explicit convection
  !
  ! Author: Satish Karra
  ! Date: 02/06/26
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none

  class(pm_stokes_ts_type) :: this

  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  call VecZeroEntries(ff,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)

  if (this%add_convection) then
    call StokesResidualConvectionInternalConn(ff,realization,this%rho,ierr)
  endif

  call VecScale(ff,-1.d0,ierr);CHKERRQ(ierr)

end subroutine PMStokesTSRHSFunction

! ************************************************************************** !

subroutine PMStokesTSIJacobian(this,ts,time,U,Udot,shift,A,B,ierr)
  !
  ! TS IJacobian for Stokes IMEX split
  !
  ! Author: Satish Karra
  ! Date: 02/06/26
  !
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none

  class(pm_stokes_ts_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  PetscReal :: shift
  Mat :: A, B
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization
  Mat :: J

  realization => this%realization
  field => realization%field
  discretization => realization%discretization

  J = B

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call StokesJacobianInternalConn(J,realization,this%mu,this%rho, &
                                      PETSC_FALSE,ierr)
  call StokesIJacobianAccumulation(J,shift,realization,this%rho,ierr)
  call StokesJacobianBoundaryFacesConn(J,realization,this%mu,ierr)
  call StokesJacobianPinPressure(J,realization,ierr)
  call StokesJacobianBoundaryConn(J,realization,this%rho,this%mu,ierr)

  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (A /= B) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  endif

end subroutine PMStokesTSIJacobian

! ************************************************************************** !

subroutine PMStokesTSInitializeTimestep(this)
  !
  ! TS timestep initialization
  !
  ! Author: Satish Karra
  ! Date: 02/06/26
  !
  implicit none

  class(pm_stokes_ts_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMStokesTSInitializeTimestep

end module PM_Stokes_TS_class
