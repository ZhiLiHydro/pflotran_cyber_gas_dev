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
    PetscBool :: add_convection
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
  this%rho = 998.32d0 ! at 20 deg C
  this%mu  = 1.002d-3 ! at 20 deg C
  this%add_convection = PETSC_FALSE

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
        option%flow%reference_density(option%liquid_phase) = tempreal
      case('VISCOSITY')
        call InputReadDouble(input,option,tempreal); this%mu  = tempreal
      case('ADD_CONVECTION')
        this%add_convection = PETSC_TRUE
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
  call PMStokesRegisterFieldSplit(this)

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
  call StokesInitializeTimestep(this%realization,this%rho)
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
  call StokesApplyPrescribedConditions(this%realization)

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
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac

  if (update_dt .and. iacceleration /= 0) then
    if (iacceleration > 0) then
      fac = 0.5d0
      if (num_newton_iterations >= iacceleration) then
        fac = 0.33d0
        ut = 0.d0
      else
        up = this%pressure_change_governor/(this%max_pressure_change+0.1)
        ut = up
      endif
      dtt = fac * dt * (1.d0 + ut)
    else
      ifac = max(min(num_newton_iterations,size(tfac)),1)
      dt_tfac = tfac(ifac) * dt

      fac = 0.5d0
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      dt_p = fac * dt * (1.d0 + up)

      dtt = min(dt_tfac,dt_p)
    endif

    dtt = min(time_step_max_growth_factor*dt,dtt)
    if (dtt > dt_max) dtt = dt_max

    dtt = max(dtt,dt_min)
    dt = dtt
  endif

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
  call StokesResidual(snes,xx,r,this%realization,this%rho,this%mu, &
                      this%add_convection,dbg,ierr)

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
  call StokesJacobian(snes,xx,A,B,this%realization,this%rho,this%mu, &
                      this%add_convection,dbg,ierr)

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
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: dX_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  PetscInt :: local_id, ghosted_id, natural_id, ival, idof
  PetscReal :: dX_, X0_, dX_X0
  PetscBool :: converged_abs_update_flag
  PetscBool :: converged_rel_update_flag
  PetscInt :: converged_abs_update_cell
  PetscInt :: converged_rel_update_cell
  PetscReal :: converged_abs_update_real
  PetscReal :: converged_rel_update_real
  PetscBool :: converged_absolute
  PetscBool :: converged_relative

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayRead(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)
  converged_abs_update_flag = PETSC_TRUE
  converged_rel_update_flag = PETSC_TRUE
  converged_abs_update_cell = ZERO_INTEGER
  converged_rel_update_cell = ZERO_INTEGER
  converged_abs_update_real = 0.d0
  converged_rel_update_real = 0.d0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    do idof = 1, option%nflowdof
      ival = (local_id-1)*option%nflowdof + idof
      ! infinity norms on update
      converged_absolute = PETSC_TRUE
      converged_relative = PETSC_TRUE
      dX_ = dabs(dX_p(ival))
      X0_ = dabs(X0_p(ival))
      X0_ = max(X0_,1.d-40)
      dX_X0 = dX_/X0_
      if (dX_ > this%abs_update_inf_tol) then
        converged_absolute = PETSC_FALSE
      endif
      if (converged_abs_update_real < dX_) then
        converged_abs_update_real = dX_
        converged_abs_update_cell = natural_id
      endif
      if (dX_X0 > this%rel_update_inf_tol) then
        converged_relative = PETSC_FALSE
      endif
      if (converged_rel_update_real < dX_X0) then
        converged_rel_update_real = dX_X0
        converged_rel_update_cell = natural_id
      endif
      ! only enter this condition if both are not converged
      if (.not.(converged_absolute .or. converged_relative)) then
        if (.not.converged_absolute) then
          converged_abs_update_flag = PETSC_FALSE
        endif
        if (.not.converged_relative) then
          converged_rel_update_flag = PETSC_FALSE
        endif
      endif
    enddo
  enddo
  call VecRestoreArrayRead(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(ABS_UPDATE_INDEX) = converged_abs_update_flag
  this%converged_flag(REL_UPDATE_INDEX) = converged_rel_update_flag
  this%converged_real(ABS_UPDATE_INDEX) = converged_abs_update_real
  this%converged_real(REL_UPDATE_INDEX) = converged_rel_update_real
  this%converged_cell(ABS_UPDATE_INDEX) = converged_abs_update_cell
  this%converged_cell(REL_UPDATE_INDEX) = converged_rel_update_cell

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
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscInt :: local_id, ghosted_id, natural_id, idof, ival
  PetscInt :: itol
  PetscReal :: R, A, R_A
  PetscReal, parameter :: A_zero = 1.d-15
  PetscBool :: converged_abs_residual_flag
  PetscReal :: converged_abs_residual_real
  PetscInt :: converged_abs_residual_cell
  PetscBool :: converged_scaled_residual_flag
  PetscReal :: converged_scaled_residual_real
  PetscInt :: converged_scaled_residual_cell
  PetscBool :: converged_absolute
  PetscBool :: converged_scaled
  PetscMPIInt :: mpi_int
  character(len=MAXSTRINGLENGTH) :: string
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid

  if (this%check_post_convergence) then
    call VecGetArrayRead(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayRead(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
    converged_abs_residual_flag = PETSC_TRUE
    converged_abs_residual_real = 0.d0
    converged_abs_residual_cell = ZERO_INTEGER
    converged_scaled_residual_flag = PETSC_TRUE
    converged_scaled_residual_real = 0.d0
    converged_scaled_residual_cell = ZERO_INTEGER
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      do idof = 1, option%nflowdof
        ival = (local_id-1)*option%nflowdof + idof
        converged_absolute = PETSC_TRUE
        converged_scaled = PETSC_TRUE
        ! infinity norms on residual
        R = dabs(r_p(ival))
        A = dabs(accum2_p(ival))
        R_A = R/A
        if (R > this%residual_abs_inf_tol) then
          converged_absolute = PETSC_FALSE
        endif
        ! find max value regardless of convergence
        if (converged_abs_residual_real < R) then
          converged_abs_residual_real = R
          converged_abs_residual_cell = natural_id
        endif
        if (A > A_zero) then
          if (R_A > this%residual_scaled_inf_tol) then
            converged_scaled = PETSC_FALSE
          endif
          ! find max value regardless of convergence
          if (converged_scaled_residual_real < R_A) then
            converged_scaled_residual_real = R_A
            converged_scaled_residual_cell = natural_id
          endif
        endif
        ! only enter this condition if both are not converged
        if (.not.(converged_absolute .or. converged_scaled)) then
          if (.not.converged_absolute) then
            converged_abs_residual_flag = PETSC_FALSE
          endif
          if (.not.converged_scaled) then
            converged_scaled_residual_flag = PETSC_FALSE
          endif
        endif
      enddo
    enddo
    call VecRestoreArrayRead(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayRead(field%flow_accum2,accum2_p, &
                                ierr);CHKERRQ(ierr)

    this%converged_flag(RESIDUAL_INDEX) = converged_abs_residual_flag
    this%converged_real(RESIDUAL_INDEX) = converged_abs_residual_real
    this%converged_cell(RESIDUAL_INDEX) = converged_abs_residual_cell
    this%converged_flag(SCALED_RESIDUAL_INDEX) = converged_scaled_residual_flag
    this%converged_real(SCALED_RESIDUAL_INDEX) = converged_scaled_residual_real
    this%converged_cell(SCALED_RESIDUAL_INDEX) = converged_scaled_residual_cell
    mpi_int = MAX_INDEX
    ! do not perform an all reduce on cell id as this info is not printed
    ! in parallel
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_flag,mpi_int,MPI_INTEGER, &
                       MPI_LAND,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                       ierr);CHKERRQ(ierr)

    ! Only declare convergence after at least one Newton iteration.
    ! At it=0 (initial residual evaluation before any solve), the residuals
    ! may be trivially zero due to prescribed BCs / zero IC, leading to
    ! spurious convergence before the system is actually solved.
    if (it > 0) then
      option%convergence = CONVERGENCE_CONVERGED
    endif

    do itol = 1, MAX_INDEX
      if (.not.this%converged_flag(itol)) then
        option%convergence = CONVERGENCE_KEEP_ITERATING
        if (this%logging_verbosity > 0) then
          string = '   ' // trim(tol_string(itol)) // ', Liquid Pressure'
          if (option%comm%size == 1) then
            string = trim(string) // ' (' // &
              trim(StringFormatInt(this%converged_cell(itol))) &
              // ')'
          endif
          string = trim(string) // ' : ' // &
            StringFormatDouble(this%converged_real(itol))
          call PrintMsg(option,string)
        endif
      endif
    enddo
    if (this%logging_verbosity > 0 .and. it > 0 .and. &
        option%convergence == CONVERGENCE_CONVERGED) then
      string = '   Converged'
      call PrintMsg(option,string)
      write(string,'(4x," R:",es8.1)') this%converged_real(RESIDUAL_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"SR:",es8.1)') &
        this%converged_real(SCALED_RESIDUAL_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"AU:",es8.1)') this%converged_real(ABS_UPDATE_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"RU:",es8.1)') this%converged_real(REL_UPDATE_INDEX)
      call PrintMsg(option,string)
    endif
  endif

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

  use Stokes_module, only : StokesUpdateAuxVars, &
                            StokesUpdateTransportVelocities

  implicit none

  class(pm_stokes_type) :: this

  call StokesUpdateAuxVars(this%realization)
  call StokesUpdateTransportVelocities(this%realization)

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
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_UPDATE: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMStokesInputRecord

! ************************************************************************** !

subroutine PMStokesRegisterFieldSplit(this)
  !
  ! Register velocity (u) and pressure (p) fields for PCFIELDSPLIT.
  ! Equivalent to the C code:
  !   const PetscInt ufields[] = {0,1,2}, pfields[] = {3};
  !   KSPGetPC(ksp_S, &pc);
  !   PCFieldSplitSetBlockSize(pc, 4);
  !   PCFieldSplitSetFields(pc, "u", 3, ufields, ufields);
  !   PCFieldSplitSetFields(pc, "p", 1, pfields, pfields);
  !
  ! Author: Satish Karra
  ! Date: 01/28/26
  !

  implicit none
  class(pm_stokes_type)  :: this
  KSP        :: ksp_flow
  PC         :: pc
  PetscErrorCode :: ierr
  PetscInt   :: ufields(3), pfields(1)
  character(len=*), parameter :: field_u = 'u'
  character(len=*), parameter :: field_p = 'p'

  ksp_flow = this%solver%ksp
  pc = this%solver%pc

  ! Define DOF positions inside the 4-wide block:
  ! [0] -> u, [1] -> v, [2] -> w, [3] -> p
  ufields(1) = 0
  ufields(2) = 1
  ufields(3) = 2
  pfields(1) = 3

  ! Get the preconditioner context
  call KSPGetPC(ksp_flow, pc, ierr);CHKERRQ(ierr)

  ! Set the block size (4 DOFs per cell: u,v,w,p)
  call PCFieldSplitSetBlockSize(pc, 4, ierr);CHKERRQ(ierr)

  ! Define velocity field "u" as DOFs 0,1,2 within each block
  call PCFieldSplitSetFields(pc, field_u, 3, ufields, ufields, ierr);CHKERRQ(ierr)

  ! Define pressure field "p" as DOF 3 within each block
  call PCFieldSplitSetFields(pc, field_p, 1, pfields, pfields, ierr);CHKERRQ(ierr)

end subroutine PMStokesRegisterFieldSplit

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
