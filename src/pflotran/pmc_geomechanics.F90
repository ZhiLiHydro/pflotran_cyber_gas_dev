module PMC_Geomechanics_class

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PMC_Base_class
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pmc_base_type) :: pmc_geomechanics_type
    class(realization_subsurface_type), pointer :: subsurf_realization
    class(realization_geomech_type), pointer :: geomech_realization
  contains
    procedure, public :: Init => PMCGeomechanicsInit
    procedure, public :: SetupSolvers => PMCGeomechanicsSetupSolvers
    procedure, public :: InitializeRun => PMCGeomechanicsInitializeRun
    procedure, public :: RunToTime => PMCGeomechanicsRunToTime
    procedure, public :: GetAuxData => PMCGeomechanicsGetAuxData
    procedure, public :: SetAuxData => PMCGeomechanicsSetAuxData
    procedure, public :: Destroy => PMCGeomechanicsDestroy
  end type pmc_geomechanics_type

  public :: PMCGeomechanicsCreate

contains

! ************************************************************************** !

function PMCGeomechanicsCreate()
  !
  ! This routine allocates and initializes a new object.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  implicit none

  class(pmc_geomechanics_type), pointer :: PMCGeomechanicsCreate

  class(pmc_geomechanics_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCGeomechanicsCreate%Create()'
#endif

  allocate(pmc)
  call pmc%Init()

  PMCGeomechanicsCreate => pmc

end function PMCGeomechanicsCreate

! ************************************************************************** !

subroutine PMCGeomechanicsInit(this)
  !
  ! This routine initializes a new process model coupler object.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  implicit none

  class(pmc_geomechanics_type) :: this

#ifdef DEBUG
  print *, 'PMCGeomechanics%Init()'
#endif

  call PMCBaseInit(this)
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsInit

! ************************************************************************** !

subroutine PMCGeomechanicsSetupSolvers(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !
  use Geomechanics_Discretization_module
  use Timestepper_KSP_class
  use PM_Base_class
  use PM_Base_Pointer_module
  use Option_module
  use Solver_module

  implicit none

  class(pmc_geomechanics_type) :: this

  class(realization_geomech_type), pointer :: geomech_realization
  class(geomech_discretization_type), pointer :: geomech_discretization
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeomechanicsSetupSolvers')
#endif

  option => this%option
  solver => this%timestepper%solver
  geomech_realization => this%geomech_realization
  geomech_discretization => geomech_realization%geomech_discretization

  select type(ts=>this%timestepper)
    class is(timestepper_KSP_type)
    class default
      option%io_buffer = 'A KSP timestepper must be used for geomechanics.'
      call PrintErrMsg(option)
  end select

  call SolverCreateKSP(solver,option%mycomm)

  call PrintMsg(option,"  Beginning setup of GEOMECH KSP")
  call KSPSetOptionsPrefix(solver%ksp,"geomech_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  if (solver%M_mat_type == MATAIJ) then
    option%io_buffer = 'AIJ matrix not supported for geomechanics.'
    call PrintErrMsg(option)
  endif


  if (Uninitialized(solver%Mpre_mat_type) .and. &
      Uninitialized(solver%M_mat_type)) then
    ! Matrix types not specified, so set to default.
    solver%Mpre_mat_type = MATBAIJ
    solver%M_mat_type = solver%Mpre_mat_type
  else if (Uninitialized(solver%Mpre_mat_type)) then
    if (solver%M_mat_type == MATMFFD) then
      solver%Mpre_mat_type = MATBAIJ
    else
      solver%Mpre_mat_type = solver%M_mat_type
    endif
  else if (Uninitialized(solver%M_mat_type)) then
    solver%M_mat_type = solver%Mpre_mat_type
  endif

  call GeomechDiscretizationCreateMatrix(geomech_realization% &
                                         geomech_discretization,NGEODOF, &
                                         solver%Mpre_mat_type, &
                                         solver%Mpre,option)

  call MatSetOptionsPrefix(solver%Mpre,"geomech_",ierr);CHKERRQ(ierr)
  solver%M = solver%Mpre
  geomech_realization%geomech_field%A = solver%M


  ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
  if (option%verbosity >= 2) then
    string = '-geomech_ksp_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                  ierr);CHKERRQ(ierr)
    string = '-geomech_ksp_monitor'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                  ierr);CHKERRQ(ierr)
  endif


  ! call KSPSetOperators(solver%ksp,solver%M,solver%Mpre,ierr);CHKERRQ(ierr)

  call SolverSetKSPOptions(solver,option)

  call PrintMsg(option,"  Finished setting up GEOMECH KSP ")

  ! jaa - account for initialization
  if (this%option%geomechanics%split_scheme == &
      GEOMECH_FIXED_STRAIN_SPLIT .or. &
      this%option%geomechanics%split_scheme == &
      GEOMECH_FIXED_STRESS_SPLIT .or. &
      this%option%geomechanics%split_scheme == &
      GEOMECH_DRAINED_SPLIT .or. &
      this%option%geomechanics%split_scheme == &
      GEOMECH_UNDRAINED_SPLIT) &
    this%timestepper%steps = -1

end subroutine PMCGeomechanicsSetupSolvers

! ************************************************************************** !


recursive subroutine PMCGeomechanicsInitializeRun(this)
  !
  ! Initializes the geomechanics process model coupler
  !
  ! Author: Glenn Hammond
  ! Date: 03/31/25
  use Timestepper_Base_class
  implicit none
  class(pmc_geomechanics_type) :: this
  PetscReal :: target_time
  PetscInt :: local_stop_flag
  call PMCBaseInitializeRun(this)
  target_time = 0.d0
  local_stop_flag = TS_CONTINUE
  call this%RunToTime(target_time,local_stop_flag)

end subroutine PMCGeomechanicsInitializeRun

! ************************************************************************** !


recursive subroutine PMCGeomechanicsRunToTime(this,sync_time,stop_flag)
  !
  ! This routine runs the geomechanics simulation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  use Timestepper_Base_class
  use Option_module
  use PM_Base_class
  use Output_Geomechanics_module

  implicit none

  class(pmc_geomechanics_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  PetscInt :: local_stop_flag
  PetscBool :: sync_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_flag
  PetscBool :: peer_already_run_to_time

  class(pm_base_type), pointer :: cur_pm

  if (stop_flag == TS_STOP_FAILURE) return

  ! jaa:  want to sync geomech dt with flow_dt (master)
  this%timestepper%dt = this%option%flow_dt
  if (this%option%geomechanics%split_scheme == GEOMECH_FIXED_STRAIN_SPLIT)&
    this%timestepper%dt = this%option%flow_dt

  call this%PrintHeader()
  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call PrintVerboseMsg(this%option)

  ! Get data of other process-model
  call this%GetAuxData()

  local_stop_flag = 0

  call SetOutputFlags(this)
  sync_flag = PETSC_FALSE
  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE

  call this%timestepper%SetTargetTime(sync_time,this%option,local_stop_flag, &
                                      sync_flag, &
                                      snapshot_plot_flag, &
                                      observation_plot_flag, &
                                      massbal_plot_flag,checkpoint_flag)
  call this%timestepper%StepDT(this%pm_list,local_stop_flag)

  ! Check if it is initial solve
  ! jaa - do we need this?
  !if (this%timestepper%steps == 1) then
  !  this%option%geomechanics%initial_flag = PETSC_TRUE
  !endif

  if (this%option%geomechanics%skip_stepdt_flag) then
    ! reset
    this%option%geomechanics%skip_stepdt_flag = PETSC_FALSE
    print *, 'skipping PMCGeomechanicsRunToTime ...'
    return
  endif

  ! Have to loop over all process models coupled in this object and update
  ! the time step size.  Still need code to force all process models to
  ! use the same time step size if tightly or iteratively coupled.
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    ! Geomechanics PM does not have an associate time
    !call this%timestepper%UpdateDT(cur_pm)
    cur_pm => cur_pm%next
  enddo

  ! Run underlying process model couplers
  if (associated(this%child)) then
    ! Set data needed by process-model
    call this%SetAuxData()
    call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
    call this%GetAuxData()
  endif

  peer_already_run_to_time = PETSC_FALSE
  if (sync_flag .and. associated(this%peer)) then
    ! synchronize peers
    call this%SetAuxData()
    ! Run neighboring process model couplers
    call this%peer%RunToTime(this%timestepper%target_time, &
                             local_stop_flag)
    peer_already_run_to_time = PETSC_TRUE
    call this%GetAuxData()
  endif

  if (this%timestepper%time_step_cut_flag) then
    snapshot_plot_flag = PETSC_FALSE
  endif
  ! however, if we are using the modulus of the output_option%imod, we may
  ! still print
  if (mod(this%timestepper%steps,this%pm_list% &
          output_option%periodic_snap_output_ts_imod) == 0) then
    snapshot_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_obs_output_ts_imod) == 0) then
    observation_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_msbl_output_ts_imod) == 0) then
    massbal_plot_flag = PETSC_TRUE
  endif

  call OutputGeomechanics(this%geomech_realization,snapshot_plot_flag, &
                          observation_plot_flag,massbal_plot_flag)
  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer) .and. .not.peer_already_run_to_time) then
    call this%SetAuxData()
    call this%peer%RunToTime(sync_time,local_stop_flag)
    call this%GetAuxData()
  endif

  stop_flag = max(stop_flag,local_stop_flag)

end subroutine PMCGeomechanicsRunToTime

! ************************************************************************** !

subroutine PMCGeomechanicsSetAuxData(this)
  !
  ! This routine updates data in simulation_aux that is required by other
  ! process models.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !
  use Option_module
  use Grid_module
  use Discretization_module
  use Geomechanics_Subsurface_Properties_module
  use Parameter_module
  use Global_Aux_module

  implicit none

  class(pmc_geomechanics_type) :: this

  type(grid_type), pointer :: grid
  PetscInt :: local_id, ghosted_id
  PetscScalar, pointer :: por0_p(:)
  PetscScalar, pointer :: por_p(:)
  PetscScalar, pointer :: perm0_p(:)
  PetscScalar, pointer :: perm_p(:)
  PetscScalar, pointer :: strain_p(:)
  PetscScalar, pointer :: stress_p(:)
  PetscScalar, pointer :: press_p(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: local_stress(6), local_strain(6), local_pressure
  PetscErrorCode :: ierr
  PetscReal :: por_new
  PetscReal :: perm_new
  PetscInt :: i
  PetscInt :: parameter_index
  ! jaa
  PetscInt :: strainv_0_id, stressv_0_id, press_0_id
  PetscInt :: stressv_id, stressv_old_id
  PetscInt :: strainv_id, strainv_old_id
  PetscInt :: press_old_id, porosity_id
  PetscInt :: temp_porosity_id, porosity_old_id
  type(global_auxvar_type), pointer :: global_auxvars(:)

  Vec :: geomech_vec
  Vec :: subsurf_vec
  PetscReal :: biot_coeff, youngs_mod, poissons_ratio, dr_bulk_modulus
  PetscReal :: strain_vol, stress_vol
  PetscReal :: press_init, porosity_0, strain_vol_0
  PetscInt :: imat

#if GEOMECH_DEBUG
  PetscViewer :: viewer
  print *, 'PMCGeomechSetAuxData'
#endif


  ! If at initialization stage, do nothing
!  if (this%timestepper%steps == 0) return

  select type(pmc => this)
    class is(pmc_geomechanics_type)
      if (this%option%geomechanics%flow_coupling == &
            GEOMECH_TWO_WAY_COUPLED .or. &
          this%option%geomechanics%geophysics_coupling == &
            GEOMECH_ERT_COUPLING) then

        grid => pmc%subsurf_realization%patch%grid
        global_auxvars => pmc%subsurf_realization% &
                          patch%aux%Global%auxvars

        ! Find the number of geomech grid nodes for each flow cell
        call VecDuplicate(pmc%geomech_realization%geomech_field%strain, &
                          geomech_vec,ierr);CHKERRQ(ierr)
        call VecSet(geomech_vec,1.d0,ierr);CHKERRQ(ierr)

        call VecDuplicate(pmc%sim_aux%subsurf_strain,subsurf_vec, &
                          ierr);CHKERRQ(ierr)
        call VecSet(subsurf_vec,0.d0,ierr);CHKERRQ(ierr)

        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf,geomech_vec, &
                             subsurf_vec,ADD_VALUES,SCATTER_FORWARD, &
                             ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf,geomech_vec, &
                           subsurf_vec,ADD_VALUES,SCATTER_FORWARD, &
                           ierr);CHKERRQ(ierr)


#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(this%option%mycomm, &
                            'subsurf_vec_adjacency_count.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(subsurf_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

       ! Save strain dataset in sim_aux%subsurf_strain
        call VecSet(pmc%sim_aux%subsurf_strain,0.d0,ierr);CHKERRQ(ierr)
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%strain, &
                             pmc%sim_aux%subsurf_strain,ADD_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%strain, &
                           pmc%sim_aux%subsurf_strain,ADD_VALUES, &
                           SCATTER_FORWARD,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(this%option%mycomm, &
                            'subsurf_strain_vector_before_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_strain,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        ! Save stress dataset in sim_aux%subsurf_stress
        call VecSet(pmc%sim_aux%subsurf_stress,0.d0,ierr);CHKERRQ(ierr)
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%stress, &
                             pmc%sim_aux%subsurf_stress,ADD_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%stress, &
                           pmc%sim_aux%subsurf_stress,ADD_VALUES, &
                           SCATTER_FORWARD,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(this%option%mycomm, &
                            'subsurf_stress_vector_before_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_stress,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        ! Calculate the average stress and strain
        call VecPointwiseDivide(pmc%sim_aux%subsurf_strain, &
                                pmc%sim_aux%subsurf_strain,subsurf_vec, &
                                ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(this%option%mycomm, &
                            'subsurf_strain_vector_after_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_strain,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        call VecPointwiseDivide(pmc%sim_aux%subsurf_stress, &
                                pmc%sim_aux%subsurf_stress,subsurf_vec, &
                                ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(this%option%mycomm, &
                            'subsurf_stress_vector_after_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_stress,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        ! Strain
        call VecGetArray(pmc%sim_aux%subsurf_strain,strain_p, &
                            ierr);CHKERRQ(ierr)
        ! Stress
        call VecGetArray(pmc%sim_aux%subsurf_stress,stress_p, &
                            ierr);CHKERRQ(ierr)

        if (this%option%geomechanics%geophysics_coupling == &
                                GEOMECH_ERT_COUPLING) then
          ! stress
          call VecGetArray(pmc%subsurf_realization%field%work, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            vec_ptr(local_id) = 0.d0
            do i = 1, 3
              vec_ptr(local_id) = vec_ptr(local_id) + &
                              stress_p((local_id - 1)*SIX_INTEGER + i)
            enddo
            vec_ptr(local_id) = vec_ptr(local_id) / 3.d0
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          call DiscretizationGlobalToLocal( &
                 pmc%subsurf_realization%discretization, &
                 pmc%subsurf_realization%field%work, &
                 pmc%subsurf_realization%field%work_loc,ONEDOF)
          call VecGetArray(pmc%subsurf_realization%field%work_loc, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          parameter_index = ParameterGetIDFromName('geomechanics_stress', &
                                             pmc%subsurf_realization%option)
          do ghosted_id = 1, grid%ngmax
            pmc%subsurf_realization%patch%aux% &
              Global%auxvars(ghosted_id)%parameters(parameter_index) = &
                vec_ptr(ghosted_id)
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work_loc, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          ! strain
          call VecGetArray(pmc%subsurf_realization%field%work, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            vec_ptr(local_id) = 0.d0
            do i = 1, 3
              vec_ptr(local_id) = vec_ptr(local_id) + &
                              strain_p((local_id - 1)*SIX_INTEGER + i)
            enddo
            !vec_ptr(local_id) = vec_ptr(local_id) / 3.d0
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          call DiscretizationGlobalToLocal( &
                 pmc%subsurf_realization%discretization, &
                 pmc%subsurf_realization%field%work, &
                 pmc%subsurf_realization%field%work_loc,ONEDOF)
          call VecGetArray(pmc%subsurf_realization%field%work_loc, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          parameter_index = ParameterGetIDFromName('geomechanics_strain', &
                                             pmc%subsurf_realization%option)
          do ghosted_id = 1, grid%ngmax
            pmc%subsurf_realization%patch%aux% &
              Global%auxvars(ghosted_id)%parameters(parameter_index) = &
                vec_ptr(ghosted_id)
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work_loc, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
        endif

        if (this%option%geomechanics%flow_coupling == &
            GEOMECH_TWO_WAY_COUPLED) then

          ! Update porosity dataset in sim_aux%subsurf_por
          call VecGetArray(pmc%sim_aux%subsurf_por0,por0_p, &
                              ierr);CHKERRQ(ierr)
          call VecGetArray(pmc%sim_aux%subsurf_por,por_p,ierr);CHKERRQ(ierr)
          ! Perm
          call VecGetArray(pmc%sim_aux%subsurf_perm0,perm0_p, &
                              ierr);CHKERRQ(ierr)
          call VecGetArray(pmc%sim_aux%subsurf_perm,perm_p, &
                              ierr);CHKERRQ(ierr)
          ! Flow
          call VecGetArray(pmc%subsurf_realization%field%flow_xx,press_p, &
                              ierr);CHKERRQ(ierr)


          do local_id = 1, grid%nlmax
            do i = 1, SIX_INTEGER
              local_stress(i) = stress_p((local_id - 1)*SIX_INTEGER + i)
              local_strain(i) = strain_p((local_id - 1)*SIX_INTEGER + i)
            enddo
            if (this%option%iflowmode == TH_MODE) then
              local_pressure = press_p(local_id*this%option%nflowdof-1)
            else ! richards mode
              local_pressure = press_p(local_id)
            endif
            ! Update porosity based on stress/strain
            call GeomechanicsSubsurfacePropsPoroEvaluate( &
                  grid, &
                  pmc%subsurf_realization%patch%aux%Material%auxvars(local_id), &
                  por0_p(local_id),local_stress,local_strain,local_pressure, &
                  por_new,this%option)
            por_p(local_id) = por_new
            ! Update permeability based on stress/strain
            call GeomechanicsSubsurfacePropsPermEvaluate( &
                  grid, &
                  pmc%subsurf_realization%patch%aux%Material%auxvars(local_id), &
                  perm0_p(local_id),local_stress,local_strain,local_pressure, &
                  perm_new,this%option)
            perm_p(local_id) = perm_new
            ! jaa begin -------------------------------------------------------
            if (this%option%geomechanics%split_scheme == &
              GEOMECH_FIXED_STRESS_SPLIT) then
              ! ID's are registered in factory_subsurface.F90
              ! get id's for global_aux%parameters
              strainv_0_id = ParameterGetIDFromName( &
                                 'geomech_vol_strain_0', &
                                 pmc%subsurf_realization%option)
              stressv_0_id = ParameterGetIDFromName( &
                                 'geomech_vol_stress_0', &
                                 pmc%subsurf_realization%option)
              strainv_id = ParameterGetIDFromName( &
                                 'geomech_vol_strain', &
                                 pmc%subsurf_realization%option)
              strainv_old_id = ParameterGetIDFromName( &
                                 'geomech_vol_strain_old', &
                                 pmc%subsurf_realization%option)
              stressv_id = ParameterGetIDFromName( &
                                 'geomech_vol_stress', &
                                 pmc%subsurf_realization%option)
              stressv_old_id = ParameterGetIDFromName( &
                                 'geomech_vol_stress_old', &
                                 pmc%subsurf_realization%option)
              press_old_id = ParameterGetIDFromName( &
                                 'geomech_old_press', &
                                 pmc%subsurf_realization%option)
              press_0_id = ParameterGetIDFromName( &
                                 'geomech_press_0', &
                                 pmc%subsurf_realization%option)
              porosity_id = ParameterGetIDFromName( &
                                 'geomech_porosity', &
                                 pmc%subsurf_realization%option)
              porosity_old_id = ParameterGetIDFromName( &
                                  'geomech_old_porosity', &
                                  pmc%subsurf_realization%option)
              temp_porosity_id = ParameterGetIDFromName( &
                                   'geomech_temp_porosity', &
                                   pmc%subsurf_realization%option)
              ghosted_id = grid%nL2G(local_id)
              imat = pmc%subsurf_realization%patch%aux%Material% &
                     auxvars(local_id)%id
              ! calc volumetric strain
              strain_vol = local_strain(1) + local_strain(2) + local_strain(3)
              ! calc volumetric stress
              stress_vol = (local_stress(1) + local_stress(2) + &
                            local_stress(3))/3.d0
              youngs_mod = this%option%geomechanics%youngs_modulus(imat)
              poissons_ratio = this%option%geomechanics%poissons_ratio(imat)
              biot_coeff = this%option%geomechanics%biot_coefficient(imat)
              dr_bulk_modulus = youngs_mod / &
                         (3.d0 * (1.d0 - (2.d0 * poissons_ratio)))
              print *, 'timestepper steps: ', this%timestepper%steps
              if (this%timestepper%steps < 0) then ! before init
                global_auxvars(ghosted_id)%parameters(stressv_id) = 0.d0
                global_auxvars(ghosted_id)%parameters(stressv_old_id) = 0.d0
                global_auxvars(ghosted_id)%parameters(press_0_id) = local_pressure
              elseif (this%timestepper%steps < 1) then ! at init
                press_init = global_auxvars(ghosted_id)%parameters(press_0_id)
                global_auxvars(ghosted_id)%parameters(strainv_0_id) = strain_vol
                global_auxvars(ghosted_id)%parameters(stressv_0_id) = stress_vol
                global_auxvars(ghosted_id)%parameters(stressv_old_id) = 0.d0
                global_auxvars(ghosted_id)%parameters(stressv_id) = &
                  ! with geomech
                  stress_vol
                  ! without geomech
                  !(biot_coeff * (local_pressure-press_init)/3.d0)
                global_auxvars(ghosted_id)%parameters(porosity_old_id) = &
                  pmc%subsurf_realization%patch%aux%Material% &
                  auxvars(local_id)%porosity_0
              else ! post init
                press_init = global_auxvars(ghosted_id)%parameters(press_0_id)
                global_auxvars(ghosted_id)%parameters(stressv_old_id) = &
                  global_auxvars(ghosted_id)%parameters(stressv_id)
                global_auxvars(ghosted_id)%parameters(stressv_id) = &
                  ! with geomech
                  stress_vol
                  ! without geomech
                  !(biot_coeff*(local_pressure-press_init)/3.d0)
                ! use temp porosity for flow sol. at fixed
                global_auxvars(ghosted_id)%parameters(porosity_old_id) = &
                  global_auxvars(ghosted_id)%parameters(temp_porosity_id)
              endif
              ! jaa - using Jeff's paper
              porosity_0 = pmc%subsurf_realization%patch%aux%Material% &
                           auxvars(local_id)%porosity_0
              strain_vol_0 = global_auxvars(ghosted_id)%parameters(strainv_0_id)
              por_p(local_id) = porosity_0 +&
                                (biot_coeff*(strain_vol-strain_vol_0)) + &
                                (((biot_coeff-porosity_0)* &
                                (1.d0-biot_coeff)/dr_bulk_modulus)* &
                                (local_pressure-press_init))
              global_auxvars(ghosted_id)%parameters(porosity_id) = por_p(local_id)
              !print *, 'id, por_p, del_volstrain, press ', local_id, por_p(local_id), &
              !       strain_vol - strain_vol_0, local_pressure
              print *, 'local id, press, por_p', local_id, local_pressure, por_p(local_id)
              global_auxvars(ghosted_id)%parameters(press_old_id) = local_pressure
            ! below changes will need to be part of undrained, drained, fixed strain
            ! bc biot coef should be part of porosity update
            elseif (this%option%geomechanics%split_scheme == &
                    GEOMECH_UNDRAINED_SPLIT &
                    .or. this%option%geomechanics%split_scheme == &
                    GEOMECH_DRAINED_SPLIT &
                    .or. this%option%geomechanics%split_scheme == &
                    GEOMECH_FIXED_STRAIN_SPLIT &
                   ) then
              strainv_0_id = ParameterGetIDFromName( &
                               'geomech_vol_strain_0', &
                               pmc%subsurf_realization%option)
              ghosted_id = grid%nL2G(local_id)
              strain_vol = local_strain(1) + local_strain(2) + local_strain(3)
              imat = pmc%subsurf_realization%patch%aux%Material% &
                     auxvars(local_id)%id
              biot_coeff = this%option%geomechanics%biot_coefficient(imat)
              if (this%timestepper%steps == 0) then ! at init
                global_auxvars(ghosted_id)%parameters(strainv_0_id) = strain_vol
                strain_vol_0 = strain_vol
              elseif (this%timestepper%steps > 0) then
                strain_vol_0 = global_auxvars(ghosted_id)%parameters(strainv_0_id)
              else ! before init
                cycle
              endif
              porosity_0 = pmc%subsurf_realization%patch%aux%Material% &
                           auxvars(local_id)%porosity_0
              !print *, 'id, porosity_0, biot_coeff, strain_vol, strain_vol_0', &
              !  local_id, porosity_0, biot_coeff, strain_vol, strain_vol_0
              por_p(local_id) = porosity_0 +&
                                (biot_coeff*(strain_vol-strain_vol_0))
              print *, 'local id, press, por_p', local_id, local_pressure, por_p(local_id)
            !else ! print pressure for other split schemes
            !  print *, 'local id, press', local_id, local_pressure
            endif
              ! jaa end --------------------------------------------------------
          enddo ! local_id loop

          call VecRestoreArray(pmc%sim_aux%subsurf_por0,por0_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArray(pmc%sim_aux%subsurf_por,por_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArray(pmc%subsurf_realization%field%flow_xx,press_p, &
                                  ierr);CHKERRQ(ierr)

          call VecRestoreArray(pmc%sim_aux%subsurf_perm0,perm0_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArray(pmc%sim_aux%subsurf_perm,perm_p, &
                                  ierr);CHKERRQ(ierr)
        endif

        call VecRestoreArray(pmc%sim_aux%subsurf_stress,stress_p, &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArray(pmc%sim_aux%subsurf_strain,strain_p, &
                                ierr);CHKERRQ(ierr)

        call VecDestroy(geomech_vec,ierr);CHKERRQ(ierr)
        call VecDestroy(subsurf_vec,ierr);CHKERRQ(ierr)

      endif

  end select

end subroutine PMCGeomechanicsSetAuxData

! ************************************************************************** !

subroutine PMCGeomechanicsGetAuxData(this)
  !
  ! This routine updates data for geomechanics simulation from other process
  ! models.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  use Option_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module

  implicit none

  class(pmc_geomechanics_type) :: this

  PetscErrorCode :: ierr

#if GEOMECH_DEBUG
print *, 'PMCGeomechanicsGetAuxData'
#endif

  select type(pmc => this)
    class is(pmc_geomechanics_type)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_pres, &
                           pmc%geomech_realization%geomech_field%press, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                         pmc%sim_aux%subsurf_pres, &
                         pmc%geomech_realization%geomech_field%press, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_temp, &
                           pmc%geomech_realization%geomech_field%temp, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                         pmc%sim_aux%subsurf_temp, &
                         pmc%geomech_realization%geomech_field%temp, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_fluid_den, &
                      pmc%geomech_realization%geomech_field%fluid_density, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_fluid_den, &
                      pmc%geomech_realization%geomech_field%fluid_density, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_por, &
                      pmc%geomech_realization%geomech_field%porosity, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_por, &
                      pmc%geomech_realization%geomech_field%porosity, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%press, &
                            pmc%geomech_realization%geomech_field%press_loc, &
                            ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%temp, &
                            pmc%geomech_realization%geomech_field%temp_loc, &
                            ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                    pmc%geomech_realization%geomech_discretization, &
                    pmc%geomech_realization%geomech_field%fluid_density, &
                    pmc%geomech_realization%geomech_field%fluid_density_loc, &
                    ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                    pmc%geomech_realization%geomech_discretization, &
                    pmc%geomech_realization%geomech_field%porosity, &
                    pmc%geomech_realization%geomech_field%porosity_loc, &
                    ONEDOF)

  end select

end subroutine PMCGeomechanicsGetAuxData

! ************************************************************************** !

subroutine PMCGeomechanicsStrip(this)
  !
  ! Deallocates members of PMC Geomechanics.
  !
  ! Author: Satish Karra
  ! Date: 06/01/16

  implicit none

  class(pmc_geomechanics_type) :: this

  call PMCBaseStrip(this)
  ! realizations destroyed elsewhere
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsStrip

! ************************************************************************** !

recursive subroutine PMCGeomechanicsDestroy(this)
  !
  ! Author: Satish Karra
  ! Date: 06/01/16
  !
  use Option_module

  implicit none

  class(pmc_geomechanics_type) :: this

#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeomechanics%Destroy()')
#endif

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif

  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif

  call PMCGeomechanicsStrip(this)

end subroutine PMCGeomechanicsDestroy

end module PMC_Geomechanics_class
