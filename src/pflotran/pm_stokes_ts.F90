module PM_Stokes_TS_class

#include "petsc/finclude/petscts.h"
  use petscts
  use Stokes_module
  use Stokes_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_Stokes_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_stokes_type) :: pm_stokes_ts_type
  contains
    procedure, public :: UpdateAuxVars => PMStokesTSUpdateAuxVars
    procedure, public :: IFunction => PMStokesTSIFunction
    procedure, public :: IJacobian => PMStokesTSIJacobian
    procedure, public :: RHSFunction => PMStokesTSRHSFunction
    procedure, public :: InitializeTimestep => PMStokesTSInitializeTimestep
  end type pm_stokes_ts_type

  public :: PMStokesTSCreate, &
            PMStokesTSUpdateAuxVarsPatch

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
  stokes_ts_pm%residual_abs_inf_tol    = 1.0d-6
  stokes_ts_pm%residual_scaled_inf_tol = 1.0d-6
  stokes_ts_pm%abs_update_inf_tol      = 1.0d-2
  stokes_ts_pm%rel_update_inf_tol      = 1.0d-4
  stokes_ts_pm%rho = 998.32d0
  stokes_ts_pm%mu  = 1.002d-3
  stokes_ts_pm%add_convection = PETSC_FALSE

  PMStokesTSCreate => stokes_ts_pm

end function PMStokesTSCreate

! ************************************************************************** !

subroutine PMStokesTSUpdateAuxVars(this)
  !
  ! Author: Satish Karra
  ! Date: 02/06/26

  implicit none

  class(pm_stokes_ts_type) :: this

  call PMStokesTSUpdateAuxVarsPatch(this%realization)

end subroutine PMStokesTSUpdateAuxVars

! ************************************************************************** !

subroutine PMStokesTSUpdateAuxVarsPatch(realization)
  !
  ! Author: Satish Karra
  ! Date: 02/06/26

  use Realization_Subsurface_class
  use Field_module
  use Grid_module
  use Patch_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(stokes_auxvar_type), pointer :: stokes_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt :: ghosted_id
  PetscReal, pointer :: xx_loc_p(:),xxdot_loc_p(:)
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  stokes_auxvars => patch%aux%Stokes%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  ! 1. Update auxvars based on new values of pressure
  call StokesUpdateAuxVars(realization)

  ! 2. Update auxvars based on new values
  call VecGetArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_xxdot_loc,xxdot_loc_p, &
                          ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    stokes_auxvars(ghosted_id)%du_dt(1) = xxdot_loc_p((ghosted_id-1)*option%nflowdof+1)
    stokes_auxvars(ghosted_id)%du_dt(2) = xxdot_loc_p((ghosted_id-1)*option%nflowdof+2)
    stokes_auxvars(ghosted_id)%du_dt(3) = xxdot_loc_p((ghosted_id-1)*option%nflowdof+3)
  enddo

  call VecRestoreArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%flow_xxdot_loc,xxdot_loc_p, &
                              ierr);CHKERRQ(ierr)

end subroutine PMStokesTSUpdateAuxVarsPatch

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
  use Matrix_Zeroing_module

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

  call StokesApplyPrescribedConditions(realization)
  call PMStokesTSUpdateAuxVarsPatch(realization)
  call StokesUpdateTransportVelocities(realization)

  call StokesResidualInternalConn(F,realization,this%rho,this%mu, &
                                  this%add_convection,ierr)
  call StokesResidualBoundaryConn(F,realization,this%rho,this%mu,ierr)
  call StokesResidualSourceSink(F,realization,this%rho,ierr)
  call StokesIFunctionAccumulation(F,realization,this%rho,ierr)
  call StokesResidualPinPressure(F,U,realization,0.d0,ierr)
  call StokesResidualBoundaryFinal(F,U,realization,ierr)
  if (associated(realization%patch%aux%Stokes)) then
    call MatrixZeroingZeroVecEntries(realization%patch%aux%Stokes%matrix_zeroing,F)
  endif

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

  ! Fully implicit STOKES_TS (BEULER): convection is assembled in IFunction/IJacobian.
  ! Keep RHS identically zero to avoid double-counting explicit terms.

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
  use Matrix_Zeroing_module, only : MatrixZeroingZeroMatEntries

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

  ! Follow Richards/TH TS convention: assemble Jacobian in B
  J = B

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,U,field%flow_xx_loc,NFLOWDOF)
  call DiscretizationGlobalToLocal(discretization,Udot,field%flow_xxdot_loc,&
                                   NFLOWDOF)
  call StokesApplyPrescribedConditions(realization)
  call PMStokesTSUpdateAuxVarsPatch(realization)
  call StokesUpdateTransportVelocities(realization)

  call StokesJacobianInternalConn(J,realization,this%mu,this%rho,&
                                     this%add_convection,ierr)
  call StokesJacobianBoundaryFacesConn(J,realization,this%mu,ierr)
  call StokesJacobianSourceSink(J,realization,this%rho,ierr)
  call StokesIJacobianAccumulation(J,shift,realization,this%rho,ierr)

  ! Match steady STOKES assembly sequence: complete ADD_VALUES exchange
  ! before any row zeroing operations.
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (associated(realization%patch%aux%Stokes)) then
    call MatrixZeroingZeroMatEntries(realization%patch%aux%Stokes%matrix_zeroing,J)
  endif

  call StokesJacobianPinPressure(J,realization,ierr)
  call StokesJacobianBoundaryConn(J,realization,this%rho,this%mu,ierr)

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
  use Stokes_module, only : StokesApplyPrescribedConditions
  implicit none

  class(pm_stokes_ts_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  if (associated(this%solver)) then
    ! TS nonlinear solve can legitimately converge at iter 0.
    this%solver%newton_min_iterations = 0
  endif
  call StokesApplyPrescribedConditions(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMStokesTSInitializeTimestep

end module PM_Stokes_TS_class
