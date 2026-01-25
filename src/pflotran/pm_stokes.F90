module PM_Stokes_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PFLOTRAN_Constants_module
  use Stokes_module

  implicit none

  private

  PetscInt, parameter :: ABS_UPDATE_INDEX      = 1
  PetscInt, parameter :: REL_UPDATE_INDEX      = 2
  PetscInt, parameter :: RESIDUAL_INDEX        = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  PetscInt, parameter :: MAX_INDEX             = SCALED_RESIDUAL_INDEX

  type, public, extends(pm_subsurface_flow_type) :: pm_stokes_type
    PetscBool :: converged_flag(MAX_INDEX)
    PetscInt  :: converged_cell(MAX_INDEX)
    PetscReal :: converged_real(MAX_INDEX)
    PetscReal :: residual_abs_inf_tol
    PetscReal :: residual_scaled_inf_tol
    PetscReal :: abs_update_inf_tol
    PetscReal :: rel_update_inf_tol
    PetscReal :: rho
    PetscReal :: mu
  contains
    procedure, public :: ReadSimulationOptionsBlock => PMStokesReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock            => PMStokesReadNewtonSelectCase
    procedure, public :: Setup                      => PMStokesSetup
    procedure, public :: InitializeRun              => PMStokesInitializeRun
    procedure, public :: InitializeTimestep         => PMStokesInitializeTimestep
    procedure, public :: Residual                   => PMStokesResidual
    procedure, public :: Jacobian                   => PMStokesJacobian
    procedure, public :: UpdateTimestep             => PMStokesUpdateTimestep
    procedure, public :: PreSolve                   => PMStokesPreSolve
    procedure, public :: PostSolve                  => PMStokesPostSolve
    procedure, public :: CheckUpdatePre             => PMStokesCheckUpdatePre
    procedure, public :: CheckUpdatePost            => PMStokesCheckUpdatePost
    procedure, public :: CheckConvergence           => PMStokesCheckConvergence
    procedure, public :: TimeCut                    => PMStokesTimeCut
    procedure, public :: UpdateSolution             => PMStokesUpdateSolution
    procedure, public :: UpdateAuxVars              => PMStokesUpdateAuxVars
    procedure, public :: MaxChange                  => PMStokesMaxChange
    procedure, public :: InputRecord                => PMStokesInputRecord
    procedure, public :: Destroy                    => PMStokesDestroy
  end type pm_stokes_type

  public :: PMStokesCreate, PMStokesInit, PMStokesDestroy, PMStokesCheckConvergence

contains

! ************************************************************************** !

function PMStokesCreate()
  !
  ! Creates a new Stokes process model
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  class(pm_stokes_type), pointer :: PMStokesCreate
  class(pm_stokes_type), pointer :: this

  allocate(this)
  call PMStokesInit(this)

  this%name   = 'Stokes Flow'
  this%header = 'STOKES FLOW'
  this%residual_abs_inf_tol    = 1.0d-6
  this%residual_scaled_inf_tol = 1.0d-6
  this%abs_update_inf_tol      = 1.0d-2
  this%rel_update_inf_tol      = 1.0d-4
  this%rho = 1000.0d0
  this%mu  = 1.0d-3

  PMStokesCreate => this

end function PMStokesCreate

! ************************************************************************** !

subroutine PMStokesInit(this)
  !
  ! Initializes Stokes process model
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  class(pm_stokes_type) :: this

  call PMSubsurfaceFlowInit(this)

end subroutine PMStokesInit

! ************************************************************************** !

subroutine PMStokesReadSimOptionsBlock(this,input)
  !
  ! Reads simulation options block for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_stokes_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer     :: option
  PetscBool :: found
  PetscReal :: tempreal

  option => this%option
  error_string = 'Stokes Options'
  input%ierr = INPUT_ERROR_NONE

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('DENSITY')
        call InputReadDouble(input,option,tempreal); this%rho = tempreal
      case('VISCOSITY')
        call InputReadDouble(input,option,tempreal); this%mu  = tempreal
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMStokesReadSimOptionsBlock

! ************************************************************************** !

subroutine PMStokesReadNewtonSelectCase(this,input,keyword,found, &
                                        error_string,option)
  !
  ! Reads Newton solver options for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_stokes_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer     :: option
  PetscBool :: found
  PetscReal :: tempreal

  option      => this%option
  error_string = 'Stokes Newton Solver'
  found = PETSC_FALSE

  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('RESIDUAL_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      this%residual_abs_inf_tol    = tempreal
      this%residual_scaled_inf_tol = tempreal
    case('RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      this%residual_abs_inf_tol    = tempreal
    case('RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      this%residual_scaled_inf_tol = tempreal
    case('UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      this%abs_update_inf_tol      = tempreal
      this%rel_update_inf_tol      = tempreal
    case('ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      this%abs_update_inf_tol      = tempreal
    case('REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      this%rel_update_inf_tol      = tempreal
    case default
      found = PETSC_FALSE
  end select

end subroutine PMStokesReadNewtonSelectCase

! ************************************************************************** !

subroutine PMStokesSetup(this)
  !
  ! Sets up Stokes process model
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  class(pm_stokes_type) :: this

  call this%SetRealization()
  call StokesSetup(this%realization,this%rho,this%mu)
  call PMSubsurfaceFlowSetup(this)

end subroutine PMStokesSetup

! ************************************************************************** !

recursive subroutine PMStokesInitializeRun(this)
  !
  ! Initializes Stokes run
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  class(pm_stokes_type) :: this

  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMStokesInitializeRun

! ************************************************************************** !

subroutine PMStokesInitializeTimestep(this)
  !
  ! Initializes Stokes timestep
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Stokes_module, only : StokesInitializeTimestep

  implicit none

  class(pm_stokes_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call StokesInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMStokesInitializeTimestep

! ************************************************************************** !

subroutine PMStokesPreSolve(this)
  !
  ! Pre-solve operations for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  class(pm_stokes_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMStokesPreSolve

! ************************************************************************** !

subroutine PMStokesPostSolve(this)
  !
  ! Post-solve operations for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  class(pm_stokes_type) :: this

end subroutine PMStokesPostSolve

! ************************************************************************** !

subroutine PMStokesUpdateTimestep(this,update_dt, &
                                  dt,dt_min,dt_max,iacceleration, &
                                  num_newton_iterations,tfac, &
                                  time_step_max_growth_factor)
  !
  ! Updates timestep for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none

  class(pm_stokes_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt, dt_min, dt_max
  PetscInt  :: iacceleration, num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)

end subroutine PMStokesUpdateTimestep

! ************************************************************************** !

subroutine PMStokesResidual(this,snes,xx,r,ierr)
  !
  ! Computes Stokes residual
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Debug_module

  implicit none

  class(pm_stokes_type) :: this
  SNES :: snes
  Vec  :: xx, r
  PetscErrorCode :: ierr
  type(debug_type) :: dbg

  dbg = this%debug
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call StokesResidual(snes,xx,r,this%realization,this%rho,this%mu,dbg,ierr)

end subroutine PMStokesResidual

! ************************************************************************** !

subroutine PMStokesJacobian(this,snes,xx,A,B,ierr)
  !
  ! Computes Stokes Jacobian
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Debug_module

  implicit none

  class(pm_stokes_type) :: this
  SNES :: snes
  Vec  :: xx
  Mat  :: A, B
  PetscErrorCode :: ierr
  type(debug_type) :: dbg

  dbg = this%debug
  call StokesJacobian(snes,xx,A,B,this%realization,this%rho,this%mu,dbg,ierr)

end subroutine PMStokesJacobian

! ************************************************************************** !

subroutine PMStokesCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Checks update before applying for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  class(pm_stokes_type) :: this
  SNES :: snes
  Vec  :: X, dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  changed = PETSC_FALSE

end subroutine PMStokesCheckUpdatePre

! ************************************************************************** !

subroutine PMStokesCheckUpdatePost(this,snes,X0,dX,X1,dX_changed,X1_changed,ierr)
  !
  ! Checks update after applying for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Grid_module
  use Patch_module
  use Option_module
  use Field_module

  implicit none

  class(pm_stokes_type) :: this
  SNES :: snes
  Vec  :: X0, dX, X1
  PetscBool :: dX_changed, X1_changed
  PetscErrorCode :: ierr
  PetscReal, pointer :: X0_p(:), dX_p(:)
  type(grid_type),   pointer :: grid
  type(patch_type),  pointer :: patch
  type(option_type), pointer :: option
  type(field_type),  pointer :: field
  PetscInt :: lid, gid, nat_id
  PetscReal :: dX_, X0_, dX_X0
  PetscBool :: conv_abs_flag, conv_rel_flag
  PetscReal :: conv_abs_real, conv_rel_real
  PetscInt  :: conv_abs_cell, conv_rel_cell
  PetscBool :: conv_abs, conv_rel

  patch  => this%realization%patch
  grid   => patch%grid
  option => this%option
  field  => this%realization%field

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayRead(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)

  conv_abs_flag = PETSC_TRUE
  conv_rel_flag = PETSC_TRUE
  conv_abs_real = 0.d0
  conv_rel_real = 0.d0
  conv_abs_cell = 0
  conv_rel_cell = 0

  do lid = 1, grid%nlmax
    gid    = grid%nL2G(lid)
    nat_id = grid%nG2A(gid)
    if (patch%imat(gid) <= 0) cycle

    dX_ = dabs(dX_p(lid))
    X0_ = dabs(X0_p(lid))
    if (X0_ < 1.d-40) X0_ = 1.d-40
    dX_X0 = dX_/X0_

    conv_abs = (dX_   <= this%abs_update_inf_tol)
    conv_rel = (dX_X0 <= this%rel_update_inf_tol)

    if (.not.conv_abs) then
      conv_abs_flag = PETSC_FALSE
      if (dX_ > conv_abs_real) then
        conv_abs_real = dX_
        conv_abs_cell = nat_id
      endif
    endif
    if (.not.conv_rel) then
      conv_rel_flag = PETSC_FALSE
      if (dX_X0 > conv_rel_real) then
        conv_rel_real = dX_X0
        conv_rel_cell = nat_id
      endif
    endif
  enddo

  call VecRestoreArrayRead(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(ABS_UPDATE_INDEX) = conv_abs_flag
  this%converged_flag(REL_UPDATE_INDEX) = conv_rel_flag
  this%converged_real(ABS_UPDATE_INDEX) = conv_abs_real
  this%converged_real(REL_UPDATE_INDEX) = conv_rel_real
  this%converged_cell(ABS_UPDATE_INDEX) = conv_abs_cell
  this%converged_cell(REL_UPDATE_INDEX) = conv_rel_cell

end subroutine PMStokesCheckUpdatePost

! ************************************************************************** !

subroutine PMStokesCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
  !
  ! Checks convergence for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use String_module

  implicit none

  class(pm_stokes_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm, unorm, fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr
  type(option_type), pointer :: option
  type(field_type),  pointer :: field
  type(grid_type),   pointer :: grid
  type(patch_type),  pointer :: patch
  PetscReal, pointer :: r_p(:), accum2_p(:)
  PetscInt :: lid, gid, nat_id
  PetscReal :: R, A, R_A
  PetscReal, parameter :: A_zero = 1.d-15
  PetscBool :: conv_abs_res_flag, conv_scaled_res_flag
  PetscReal :: conv_abs_res_real, conv_scaled_res_real
  PetscInt  :: conv_abs_res_cell, conv_scaled_res_cell
  PetscMPIInt :: mpi_int
  character(len=MAXSTRINGLENGTH) :: string
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
       ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  option => this%option
  patch  => this%realization%patch
  field  => this%realization%field
  grid   => patch%grid

  call VecGetArrayRead(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)

  conv_abs_res_flag    = PETSC_TRUE
  conv_scaled_res_flag = PETSC_TRUE
  conv_abs_res_real    = 0.d0
  conv_scaled_res_real = 0.d0
  conv_abs_res_cell    = 0
  conv_scaled_res_cell = 0

  do lid = 1, grid%nlmax
    gid    = grid%nL2G(lid)
    nat_id = grid%nG2A(gid)
    if (patch%imat(gid) <= 0) cycle

    R   = dabs(r_p(lid))
    A   = dabs(accum2_p(lid))
    R_A = R / max(A,A_zero)

    if (R > this%residual_abs_inf_tol) then
      conv_abs_res_flag = PETSC_FALSE
      if (R > conv_abs_res_real) then
        conv_abs_res_real = R
        conv_abs_res_cell = nat_id
      endif
    endif

    if (R_A > this%residual_scaled_inf_tol) then
      conv_scaled_res_flag = PETSC_FALSE
      if (R_A > conv_scaled_res_real) then
        conv_scaled_res_real = R_A
        conv_scaled_res_cell = nat_id
      endif
    endif
  enddo

  call VecRestoreArrayRead(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)

  this%converged_flag(RESIDUAL_INDEX)        = conv_abs_res_flag
  this%converged_real(RESIDUAL_INDEX)        = conv_abs_res_real
  this%converged_cell(RESIDUAL_INDEX)        = conv_abs_res_cell
  this%converged_flag(SCALED_RESIDUAL_INDEX) = conv_scaled_res_flag
  this%converged_real(SCALED_RESIDUAL_INDEX) = conv_scaled_res_real
  this%converged_cell(SCALED_RESIDUAL_INDEX) = conv_scaled_res_cell

  mpi_int = MAX_INDEX
  call MPI_Allreduce(MPI_IN_PLACE,this%converged_flag,mpi_int,MPI_INTEGER, &
                     MPI_LAND,option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)

  option%convergence = CONVERGENCE_CONVERGED
  if (.not.all(this%converged_flag)) option%convergence = CONVERGENCE_KEEP_ITERATING

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)

end subroutine PMStokesCheckConvergence

! ************************************************************************** !

subroutine PMStokesTimeCut(this)
  !
  ! Handles time step cut for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  class(pm_stokes_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call PMSubsurfaceFlowTimeCutPostInit(this)

end subroutine PMStokesTimeCut

! ************************************************************************** !

subroutine PMStokesUpdateSolution(this)
  !
  ! Updates solution for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Stokes_module, only : StokesUpdateSolution

  implicit none

  class(pm_stokes_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call StokesUpdateSolution(this%realization)

end subroutine PMStokesUpdateSolution

! ************************************************************************** !

subroutine PMStokesUpdateAuxVars(this)
  !
  ! Updates auxiliary variables for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Stokes_module, only : StokesUpdateAuxVars

  implicit none

  class(pm_stokes_type) :: this

  call StokesUpdateAuxVars(this%realization)

end subroutine PMStokesUpdateAuxVars

! ************************************************************************** !

subroutine PMStokesMaxChange(this)
  !
  ! Computes maximum change for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Stokes_module, only : StokesMaxChange
  use Option_module

  implicit none

  class(pm_stokes_type) :: this

  call StokesMaxChange(this%realization,this%max_pressure_change)
  write(this%option%io_buffer,'("  --> max change: d(max DOF)= ",1pe12.4)') &
        this%max_pressure_change
  call PrintMsg(this%option)

end subroutine PMStokesMaxChange

! ************************************************************************** !

subroutine PMStokesInputRecord(this)
  !
  ! Writes input record for Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  implicit none

  class(pm_stokes_type) :: this
  PetscInt :: id

  id = INPUT_RECORD_UNIT
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'stokes'

end subroutine PMStokesInputRecord

! ************************************************************************** !

subroutine PMStokesDestroy(this)
  !
  ! Destroys Stokes process model
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Stokes_module, only : StokesDestroy

  implicit none

  class(pm_stokes_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif


  call StokesDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMStokesDestroy

end module PM_Stokes_class