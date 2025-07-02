module Geothermal_Well_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use TH_Aux_module
  use PM_Well_class
  use PM_Base_class
  use Well_Grid_module
  use Option_module
  use Input_Aux_module
  use Condition_module
  use Geometry_module
  use Strata_module
  use String_module
  use Grid_module
  use Patch_module
  use Material_Aux_module
  use Global_Aux_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !    TOP
  !   -------- ql_bc(2)
  !   | i=n|
  !   -------- liq%ql(k=n-1)
  !   |    |
  !   ------
  !   |    |
  !   -------- liq%ql(k=3)
  !   | i=3|
  !   -------- liq%ql(k=2)
  !   | i=2|
  !   -------- liq%ql(k=1)
  !   | i=1|
  !   -------- ql_bc(1)
  !    BOT
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Geothermal Quasi-Implicit Well Model
  type, public, extends(pm_well_qi_type) :: pm_well_geoth_qi_type
    contains

    procedure, public :: ReadPMBlock => PMWellReadPMBlockGeoThQI
    procedure, public :: ReadSimulationOptionsBlock => &
                         PMWellReadSimOptionsBlockGeoThQI
    procedure, public :: Setup => PMWellSetupGeoThQI
    procedure, public :: InitializeTimestep => PMWellInitializeTimestepGeoThQI
    procedure, public :: Solve => PMWellSolveGeoThQI
    procedure, public :: SolveFlow => PMWellSolveFlowGeoThQI
    procedure, public :: ModifyFlowResidual => PMWellModifyFlowResGeoThQI
    procedure, public :: ModifyFlowJacobian => PMWellModifyFlowJacGeoThQI
    procedure, public :: UpdateFlowRates => PMWellUpdateFlowRatesGeoThQI
    procedure, public :: UpdateFlowProperties => UpdateFlowPropertiesGeoThQI
  end type pm_well_geoth_qi_type

  public :: PMWellGeoThQICreate

  contains

! ************************************************************************** !

function PMWellGeoThQICreate()
  !
  ! Creates the geothermal quasi-implicit well process model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_geoth_qi_type), pointer :: PMWellGeoThQICreate
  class(pm_well_geoth_qi_type), pointer :: pm_well

  allocate(pm_well)
  call PMWellQIInit(pm_well)

  pm_well%well%well_model_type = WELL_MODEL_TH_GEOTHERMAL
  pm_well%flow_coupling = QUASI_IMPLICIT_WELL

  nullify(pm_well%next_well)

  PMWellGeoThQICreate => pm_well

end function PMWellGeoThQICreate

! ************************************************************************** !

subroutine PMWellSetupGeoTh(pm_well)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  use Utility_module
  use Coupler_module
  use Connection_module
  use Input_Aux_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: res_grid
  type(well_grid_type), pointer :: well_grid
  type(coupler_type), pointer :: source_sink
  type(input_type) :: input_dummy
  class(realization_subsurface_type), pointer :: realization
  type(flow_sub_condition_type), pointer :: rate,temperature,pressure
  type(TH_auxvar_type), pointer :: TH_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: sum_connection,num_new_source_sinks, offset
  PetscInt :: k

  call PMWellSetupBase(pm_well)
  option => pm_well%option
  realization => pm_well%realization
  patch => realization%patch
  res_grid => patch%grid
  well_grid => pm_well%well_grid

  pm_well%flow_soln%ndof = pm_well%nphase

  if (option%itranmode /= NULL_MODE) then
    option%io_buffer ="The WELLBORE_MODEL, type WELL_MODEL_TH_GEOTHERMAL &
      &doesn't yet work with a transport mode. You need to remove the &
      &transport process model from your simulation for the WELLBORE_MODEL, &
      &type WELL_MODEL_TH_GEOTHERMAL to work."
    call PrintErrMsg(option)
  endif
  if (option%iflowmode /= TH_MODE) then
    option%io_buffer ="The WELLBORE_MODEL, type WELL_MODEL_TH_GEOTHERMAL &
      &can only be used with TH flow mode."
    call PrintErrMsg(option)
  endif

  ! add a reservoir src/sink coupler for each well segment
  num_new_source_sinks = 0
  if (associated(patch%aux%Global%auxvars_ss)) then
    offset = patch%aux%Global%num_aux_ss
  else
    offset = 0
  endif
  do k = 1,well_grid%nsegments

    source_sink => CouplerCreate()
    source_sink%itype = SRC_SINK_COUPLER_TYPE
    source_sink%name = trim(pm_well%name) // '_wellseg_' // StringWrite(k)
    source_sink%well_name = trim(pm_well%name)
    source_sink%flow_condition_name = trim(pm_well%name) // &
                          '_wellseg_' // StringWrite(k) // '_flow_ss'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
    allocate(source_sink%flow_condition%itype(TWO_INTEGER))
    source_sink%flow_condition%num_sub_conditions = TWO_INTEGER
    allocate(source_sink%flow_condition%sub_condition_ptr(THREE_INTEGER))

    if (well_grid%h_rank_id(k) /= option%myrank) cycle
    source_sink%flow_condition%well => FlowSubConditionCreate(ONE_INTEGER)

    ! ----- flow: begin ------------------
    source_sink%flow_condition%itype(TH_PRESSURE_DOF) = MASS_RATE_SS
    nullify(source_sink%flow_condition%sub_condition_ptr(1)%ptr)
    rate => FlowSubConditionCreate(ONE_INTEGER)
    rate%name = trim(pm_well%name) // '_wellseg_' // StringWrite(k) &
                       // '_rate_ss'
    rate%units = 'kg/sec'
    rate%itype = MASS_RATE_SS
    allocate(rate%dataset%rarray(1))
    rate%dataset%rarray(:) = 0.d0
    source_sink%flow_condition%rate => rate
    source_sink%flow_condition%sub_condition_ptr(1)%ptr => rate
    nullify(source_sink%flow_condition%sub_condition_ptr(2)%ptr)
    pressure => FlowSubConditionCreate(ONE_INTEGER)
    pressure%name = trim(pm_well%name) // '_wellseg_' // StringWrite(k) &
                       // '_pressure_ss'
    pressure%units = 'Pa'
    pressure%itype = DIRICHLET_BC
    allocate(pressure%dataset%rarray(1))
    pressure%dataset%rarray(:) = 0.d0
    source_sink%flow_condition%pressure => pressure
    source_sink%flow_condition%sub_condition_ptr(2)%ptr => pressure
    ! ----- flow: end --------------------

    ! ----- temperature: begin ------------------
    source_sink%flow_condition%itype(TH_TEMPERATURE_DOF) = DIRICHLET_BC ! or ZERO_GRADIENT_BC?
    nullify(source_sink%flow_condition%sub_condition_ptr(3)%ptr)
    temperature => FlowSubConditionCreate(ONE_INTEGER)
    temperature%name = trim(pm_well%name) // '_wellseg_' // StringWrite(k) &
                       // '_temp_ss'
    temperature%units = 'C'
    temperature%itype = DIRICHLET_BC ! or ZERO_GRADIENT_BC?
    allocate(temperature%dataset%rarray(1))
    temperature%dataset%rarray(:) = 25.d0
    source_sink%flow_condition%temperature => temperature
    source_sink%flow_condition%sub_condition_ptr(3)%ptr => temperature
    ! ----- temperature: end --------------------

    source_sink%connection_set => &
      ConnectionCreate(1,GENERIC_CONNECTION_TYPE,res_grid%itype)
    source_sink%connection_set%id_dn = well_grid%h_local_id(k)
    if (well_grid%h_local_id(k) < 0) then
      source_sink%connection_set%num_connections = 0
    else
      num_new_source_sinks = num_new_source_sinks + 1
      source_sink%connection_set%offset = num_new_source_sinks + offset - 1
    endif

    call CouplerAddToList(source_sink,patch%source_sink_list)
    nullify(source_sink)
  enddo

  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > num_new_source_sinks) then
    ! there were some non-well ss's in TH flow conditions
    ! figure out how to add to them, or copy what's there into a new one
    ! that also includes the well ss's
    option%io_buffer = 'There are non-well source_sink objects. Stopping!'
    call PrintErrMsg(option)
  else
    ! all ss's are well ss's, carry out line 375 in th.F90
    allocate(TH_auxvars_ss(sum_connection))
    allocate(global_auxvars_ss(sum_connection))
    do k = 1, sum_connection
      call THAuxVarInit(TH_auxvars_ss(k),option)
      call GlobalAuxVarInit(global_auxvars_ss(k),option)
    enddo
    patch%aux%TH%auxvars_ss => TH_auxvars_ss
    patch%aux%Global%auxvars_ss => global_auxvars_ss
    patch%aux%TH%num_aux_ss = sum_connection
    patch%aux%Global%num_aux_ss = sum_connection
  endif
  !call THUpdateAuxVars(realization) ! can't call, or there is a circular dependency

  ! Next, re-set up patch%aux%TH%auxvars_ss => TH_auxvars_ss. This was already
  ! done for TH mode in th.F90, line 372. You will need to add the new set of
  ! auxvars_ss for all well segments. Then call AuxvarCompute again.

end subroutine PMWellSetupGeoTh

! ************************************************************************** !

subroutine PMWellSetupGeoThQI(this)

  implicit none

  class(pm_well_geoth_qi_type) :: this

  call PMWellSetupGeoTh(this)

end subroutine PMWellSetupGeoThQI

! ************************************************************************** !

subroutine PMWellReadSimOptionsBlockGeoThQI(this,input)

  implicit none

  type(input_type), pointer :: input
  class(pm_well_geoth_qi_type) :: this

  call PMWellReadSimOptionsBlockBase(this,input)

end subroutine PMWellReadSimOptionsBlockGeoThQI

! ************************************************************************** !

subroutine PMWellReadPMBlockGeoTh(pm_well,input)
  !
  ! Reads input file parameters associated with the well process model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(input_type), pointer :: input

  ! PLACEHOLDER (this is called by pm_well%ReadPMBlock => PMWellReadPMBlockGeoThQI)

  call PMWellReadPMBlockBase(pm_well,input)

end subroutine PMWellReadPMBlockGeoTh

! ************************************************************************** !

subroutine PMWellReadPMBlockGeoThQI(this,input)

  implicit none

  class(pm_well_geoth_qi_type) :: this
  type(input_type), pointer :: input

  call PMWellReadPMBlockGeoTh(this,input)

end subroutine PMWellReadPMBlockGeoThQI

! ************************************************************************** !

subroutine PMWellInitializeTimestepGeoTh(pm_well)
  !
  ! Initializes and takes the time step for the well process model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  PetscReal :: curr_time

  curr_time = pm_well%option%time - pm_well%option%flow_dt
  pm_well%dt_flow = pm_well%realization%option%flow_dt

  if (any(pm_well%well_grid%h_rank_id == pm_well%option%myrank)) then
    call PMWellInitializeTimestepFlow(pm_well,curr_time)
  endif

end subroutine PMWellInitializeTimestepGeoTh

! ************************************************************************** !

subroutine PMWellInitializeTimestepGeoThQI(this)

  implicit none

  class(pm_well_geoth_qi_type) :: this

  call PMWellInitializeTimestepGeoTh(this)

end subroutine PMWellInitializeTimestepGeoThQI

! ************************************************************************** !

subroutine PMWellInitializeTimestepFlow(pm_well,curr_time)
  !
  ! Initializes and takes the time step for the well process model - flow.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  use Option_module

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: curr_time

  PetscInt :: i

  write(*,*) 'PMWellInitializeTimestepFlow is called within geothermal_well.F90'

  ! update the reservoir object with current reservoir properties
  call PMWellCopyReservoir(pm_well%well%reservoir,pm_well%well%reservoir_save, &
                           pm_well%transport)

  call PMWellUpdateReservoirGeoTh(pm_well,UNINITIALIZED_INTEGER)

  call PMWellComputeWellIndex(pm_well)

  call PMWellUpdateStrata(pm_well,curr_time)

  if (initialize_well_flow) then
      ! enter here if its the very first timestep
      call PMWellInitializeWellFlow(pm_well)
  endif

  do i = 1,pm_well%option%nflowdof
    call PMWellCopyWell(pm_well%well,pm_well%well_pert(i),pm_well%transport)
  enddo

  ! TODO: Get rid of bh_p_set_by_reservoir, bh_sg_set_by_reservoir variables
  if (pm_well%well%bh_p_set_by_reservoir) then
      pm_well%well%bh_p = pm_well%well%reservoir%p_l(1)
  endif
  if (pm_well%well%bh_sg_set_by_reservoir) then
      pm_well%well%bh_sg = pm_well%well%reservoir%s_g(1)
  endif

  call pm_well%UpdateFlowProperties(PETSC_FALSE,ZERO_INTEGER)

end subroutine PMWellInitializeTimestepFlow

! ************************************************************************** !

subroutine PMWellUpdateReservoirGeoTh(pm_well,TH_update_index)
  !
  ! Updates the reservoir properties for the well process model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: TH_update_index

  type(well_reservoir_type), pointer :: reservoir
  type(grid_type), pointer :: res_grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(well_comm_type), pointer :: well_comm
  type(material_auxvar_type), pointer :: material_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(TH_auxvar_type), pointer :: TH_auxvar
  PetscInt :: ghosted_id
  PetscInt :: k,vec_size
  PetscErrorCode :: ierr

  option => pm_well%option
  well_comm => pm_well%well_comm
  reservoir => pm_well%well%reservoir
  patch => pm_well%realization%patch
  res_grid => patch%grid

  reservoir%tmp_flow = -MAX_DOUBLE

  do k = 1,pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_rank_id(k) /= option%myrank) cycle

    ghosted_id = pm_well%well_grid%h_ghosted_id(k)

    material_auxvar => patch%aux%material%auxvars(ghosted_id)
    global_auxvar => patch%aux%Global%auxvars(ghosted_id)
    TH_auxvar => patch%aux%TH%auxvars(ghosted_id)

    ! TODO: Do we need to extend the well_reservoir_type to include TH mode
    ! reservoir variables?

    ! reservoir%tmp_flow(k,?) = TH_auxvar%pc     ! liquid capillary pressure
    ! reservoir%tmp_flow(k,?) = TH_auxvar%avgmw  ! avg. mw
    ! reservoir%tmp_flow(k,?) = TH_auxvar%u      ! liquid energy

    reservoir%tmp_flow(k,1) = global_auxvar%pres(1)     ! liquid pressure
    reservoir%tmp_flow(k,2) = global_auxvar%sat(1)      ! liquid saturation
    reservoir%tmp_flow(k,3) = 1.d0-global_auxvar%sat(1) ! gas saturation
    reservoir%tmp_flow(k,4) = TH_auxvar%kvr              ! liquid mobility
    reservoir%tmp_flow(k,5) = global_auxvar%den_kg(1)   ! liquid density
    reservoir%tmp_flow(k,6) = TH_auxvar%vis              ! liquid viscosity
    reservoir%tmp_flow(k,7) = global_auxvar%temp        ! liq temperature
    reservoir%tmp_flow(k,8) = TH_auxvar%h                ! liquid enthalpy

    reservoir%tmp_flow(k,14) = material_auxvar%permeability(1) ! kx
    reservoir%tmp_flow(k,15) = material_auxvar%permeability(2) ! ky
    reservoir%tmp_flow(k,16) = material_auxvar%permeability(3) ! kz
    reservoir%tmp_flow(k,17) = material_auxvar%volume

    if (res_grid%itype == STRUCTURED_GRID) then
      reservoir%tmp_flow(k,18) = res_grid%structured_grid%dx(ghosted_id) ! dx
      reservoir%tmp_flow(k,19) = res_grid%structured_grid%dy(ghosted_id) ! dy
      reservoir%tmp_flow(k,20) = res_grid%structured_grid%dz(ghosted_id) ! dz
    else
      reservoir%tmp_flow(k,20) = pm_well%well_grid%res_dz(k) ! dz
      reservoir%tmp_flow(k,18) = sqrt(material_auxvar%volume/ &
                                  reservoir%tmp_flow(k,20)) ! dx
      reservoir%tmp_flow(k,19) = reservoir%tmp_flow(k,18) ! dy
    endif

  enddo

  vec_size = pm_well%well_grid%nsegments
  if (well_comm%commsize > 1) then
    ! Updates reservoir property vector in place using the maximum value.
    ! The rank-updated value will be larger then the initialized value (-MAX_DOUBLE)
    call MPI_Allreduce(MPI_IN_PLACE,reservoir%tmp_flow,vec_size*20,&
                       MPI_DOUBLE_PRECISION,MPI_MAX,pm_well%well_comm%comm,ierr)
  endif

  do k = 1,pm_well%well_grid%nsegments
    reservoir%p_l(k) = reservoir%tmp_flow(k,1)
    reservoir%s_l(k) = reservoir%tmp_flow(k,2)
    reservoir%s_g(k) = reservoir%tmp_flow(k,3)
    reservoir%mobility_l(k) = reservoir%tmp_flow(k,4)
    reservoir%den_l(k) = reservoir%tmp_flow(k,5)
    reservoir%visc_l(k) = reservoir%tmp_flow(k,6)
    reservoir%temp(k) = reservoir%tmp_flow(k,7)
    reservoir%H_l(k) = reservoir%tmp_flow(k,8)

    reservoir%kx(k) = reservoir%tmp_flow(k,14)
    reservoir%ky(k) = reservoir%tmp_flow(k,15)
    reservoir%kz(k) = reservoir%tmp_flow(k,16)
    reservoir%volume(k) = reservoir%tmp_flow(k,17)
    reservoir%dx(k) = reservoir%tmp_flow(k,18)
    reservoir%dy(k) = reservoir%tmp_flow(k,19)
    reservoir%dz(k) = reservoir%tmp_flow(k,20)
  enddo

end subroutine PMWellUpdateReservoirGeoTh

! ************************************************************************** !

subroutine UpdateFlowPropertiesGeoTh(pm_well,pert,index)
  !
  ! Updates flow well object properties, when TH is the flow mode.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscBool :: pert
  PetscInt :: index

  ! PLACEHOLDER (this is pm_well%UpdateFlowProperties => UpdateFlowPropertiesGeoThQI)
  ! see UpdateFlowPropertiesWIPP() in wipp_well.F90

end subroutine UpdateFlowPropertiesGeoTh

! ************************************************************************** !

subroutine UpdateFlowPropertiesGeoThQI(this,pert,index)

  implicit none

  class(pm_well_geoth_qi_type) :: this
  PetscBool :: pert
  PetscInt :: index

  call UpdateFlowPropertiesGeoTh(this,pert,index)

end subroutine UpdateFlowPropertiesGeoThQI

! ************************************************************************** !

subroutine PMWellUpdateFlowRatesGeoThQI(this,well_pert,res_pert,segment_index,&
                                        ierr)
  !
  ! This subroutine performs the well rate computation when called from the
  ! fully- or quasi-coupled source/sink update.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_geoth_qi_type) :: this
  PetscInt :: well_pert
  PetscInt :: res_pert
  PetscInt :: segment_index
  PetscErrorCode :: ierr

  ! PLACEHOLDER (this is pm_well%UpdateFlowRates)

end subroutine PMWellUpdateFlowRatesGeoThQI

! ************************************************************************** !

subroutine PMWellNewtonFlow(pm_well)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  use Utility_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! PLACEHOLDER

end subroutine PMWellNewtonFlow

! ************************************************************************** !

subroutine PMWellResidualFlow(pm_well)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! PLACEHOLDER

end subroutine PMWellResidualFlow

! ************************************************************************** !

subroutine PMWellJacobianFlow(pm_well)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! PLACEHOLDER

end subroutine PMWellJacobianFlow

! ************************************************************************** !

subroutine PMWellFillJacFlow(pm_well,Jac,Jtmp,id1,id2)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jtmp(pm_well%nphase,pm_well%nphase), &
               Jac(pm_well%nphase*pm_well%well_grid%nsegments, &
                   pm_well%nphase*pm_well%well_grid%nsegments)
  PetscInt :: id1,id2

  PetscInt :: i,j

  do i = 1,pm_well%nphase
    do j = 1,pm_well%nphase
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) = &
      Jac((id1-1)*pm_well%nphase+i,(id2-1)*pm_well%nphase+j) + Jtmp(i,j)
    enddo
  enddo

end subroutine PMWellFillJacFlow

! ************************************************************************** !

subroutine PMWellAccumulationFlow(pm_well,well,id,Res)
  !
  ! Computes the accumulation term for the flow residual based on the
  ! chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0

  ! PLACEHOLDER

end subroutine PMWellAccumulationFlow

! ************************************************************************** !

subroutine PMWellAccumDerivative(pm_well,local_id,Jac)
  !
  ! Computes the derivative of the accumulation term for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: local_id
  PetscReal :: Jac(pm_well%nphase,pm_well%nphase)

  PetscInt :: idof, irow
  PetscReal :: res(pm_well%nphase),res_pert(pm_well%nphase)

  call PMWellAccumulationFlow(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellAccumulationFlow(pm_well,pm_well%well_pert(idof),local_id, &
                                res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellAccumDerivative

! ************************************************************************** !

subroutine PMWellSrcSink(pm_well,well,id,Res)
  !
  ! Computes the source/sink term for the flow residual based on the
  ! chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well
  PetscInt :: id
  PetscReal :: Res(pm_well%nphase)

  Res = 0.d0

  ! PLACEHOLDER

end subroutine PMWellSrcSink

! ************************************************************************** !

subroutine PMWellSrcSinkDerivative(pm_well,local_id,Jac)
  !
  ! Computes the derivative of the source/sink term for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: local_id
  PetscReal :: Jac(pm_well%nphase,pm_well%nphase)

  PetscInt :: idof, irow
  PetscReal :: res(pm_well%nphase),res_pert(pm_well%nphase)

  call PMWellSrcSink(pm_well,pm_well%well,local_id,res)

  do idof = 1, pm_well%nphase
    call PMWellSrcSink(pm_well,pm_well%well_pert(idof),local_id,res_pert)
    do irow = 1, pm_well%nphase
      Jac(irow,idof) = (res_pert(irow)-res(irow))/pm_well%pert(local_id,idof)
    enddo !irow
  enddo ! idof

end subroutine PMWellSrcSinkDerivative

! ************************************************************************** !

subroutine PMWellFlux(pm_well,well_up,well_dn,iup,idn,Res,save_flux)
  !
  ! Computes the internal flux terms for the residual based on
  ! the chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well_up, well_dn
  PetscInt :: iup, idn
  PetscReal :: Res(pm_well%nphase)
  PetscBool :: save_flux

  ! PLACEHOLDER

  ! TODO: Consider changing the variables iup and idn to be consistent with
  !       the well model segment ordering, and consistent with up the well
  !       and down the well, i.e. iup is the cell towards the top of the well
  !       and idn is the cell towards the bottom of the well.

end subroutine PMWellFlux

! ************************************************************************** !

subroutine PMWellFluxDerivative(pm_well,iup,idn,Jup,Jdn)
  !
  ! Computes the derivative of the internal flux terms for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: iup, idn, iphase, irow
  PetscReal :: Jup(pm_well%nphase,pm_well%nphase), &
               Jdn(pm_well%nphase,pm_well%nphase)

  PetscReal :: res_up(pm_well%nphase),res_dn(pm_well%nphase), &
               res_pert(pm_well%nphase)

  call PMWellFlux(pm_well,pm_well%well,pm_well%well,iup,idn,res_up,PETSC_FALSE)

  res_dn = res_up

  ! upgradient derivatives
  do iphase = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well_pert(iphase),pm_well%well,iup,idn, &
                    res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jup(irow,iphase) = (res_pert(irow)-res_up(irow)) / &
                         pm_well%pert(iup,iphase)
    enddo !irow
  enddo

  ! downgradient derivatives
  do iphase = 1,pm_well%nphase
    call PMWellFlux(pm_well,pm_well%well,pm_well%well_pert(iphase),iup,idn, &
                    res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jdn(irow,iphase) = (res_pert(irow)-res_dn(irow)) / &
                         pm_well%pert(idn,iphase)
    enddo !irow
  enddo

end subroutine PMWellFluxDerivative

! ************************************************************************** !

subroutine PMWellBCFlux(pm_well,well,Res,save_flux)
  !
  ! Computes the boundary flux terms for the residual based on the
  ! chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  type(well_type) :: well
  PetscReal :: Res(2*pm_well%nphase)
  PetscBool :: save_flux

  ! PLACEHOLDER

end subroutine PMWellBCFlux

! ************************************************************************** !

subroutine PMWellBCFluxDerivative(pm_well,Jtop,Jbtm)
  !
  ! Computes the derivative of the boundary flux terms for the Jacobian,
  ! based on the chosen well model.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: Jtop(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof), &
               Jbtm(pm_well%flow_soln%ndof,pm_well%flow_soln%ndof)

  PetscInt :: idof, irow
  PetscReal :: res(2*pm_well%flow_soln%ndof),res_pert(2*pm_well%flow_soln%ndof)

  call PMWellBCFlux(pm_well,pm_well%well,res,PETSC_FALSE)

  Jtop = 0.d0
  Jbtm = 0.d0

  ! downgradient derivatives
  do idof = 1,pm_well%nphase
    call PMWellBCFlux(pm_well,pm_well%well_pert(idof),res_pert,PETSC_FALSE)
    do irow = 1, pm_well%nphase
      Jbtm(irow,idof) = (res_pert(irow)-res(irow)) / &
                        pm_well%pert(1,idof)
    enddo !irow
    do irow = 1, pm_well%nphase
      Jtop(irow,idof) = (res_pert(irow + pm_well%flow_soln%ndof)- &
                         res(irow + pm_well%flow_soln%ndof)) / &
                         pm_well%pert(pm_well%well_grid%nsegments,idof)
    enddo !irow
  enddo

end subroutine PMWellBCFluxDerivative

! ************************************************************************** !

subroutine PMWellPerturb(pm_well)
  !
  ! Calculates the state variables for the perturbed well system.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! PLACEHOLDER

end subroutine PMWellPerturb

! ************************************************************************** !

subroutine PMWellModifyFlowResGeoThQI(this,residual)
  !
  ! This subroutine computes the well contribution to the reservoir residual
  ! when called from the quasi-coupled source/sink update in TH flow mode.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_geoth_qi_type) :: this
  PetscReal, pointer :: residual(:)

  PetscReal :: Res(this%flow_soln%ndof)

  ! PLACEHOLDER (this is pm_well%ModifyFlowResidual)

  Res(:) = 0.d0

  ! PLACEHOLDER (this is pm_well%ModifyFlowResidual)

end subroutine PMWellModifyFlowResGeoThQI

! ************************************************************************** !

subroutine PMWellModifyFlowJacGeoThQI(this,Jac,ierr)
  !
  ! This subroutine computes the well contribution to the reservoir
  ! Jacobian when called from the fully- or quasi-coupled source/sink update.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_geoth_qi_type) :: this
  Mat :: Jac
  PetscErrorCode :: ierr

  ! PLACEHOLDER (this is pm_well%ModifyFlowJacobian)

end subroutine PMWellModifyFlowJacGeoThQI

! ************************************************************************** !

subroutine PMWellPreSolveFlow(pm_well)
    !
    ! Author: Jennifer M. Frederick
    ! Date: 08/04/2021
    !

    implicit none

    class(pm_well_sequential_type) :: pm_well

    character(len=MAXSTRINGLENGTH) :: out_string
    PetscReal :: cur_time, cur_time_converted
    PetscReal :: dt_converted

    pm_well%flow_soln%not_converged = PETSC_TRUE
    pm_well%flow_soln%converged = PETSC_FALSE

    cur_time = pm_well%option%time + pm_well%cumulative_dt_flow
    cur_time_converted = cur_time/pm_well%output_option%tconv
    dt_converted = pm_well%dt_flow/pm_well%output_option%tconv

    if (pm_well%print_output) then
      write(out_string,'(" FLOW Step ",i6,"   Time =",1pe12.5,"   Dt =", &
                         1pe12.5," ",a4)') &
                       (pm_well%flow_soln%n_steps+1),cur_time_converted, &
                       dt_converted,pm_well%output_option%tunit
      call PMWellPrint(pm_well,out_string)
    endif

  end subroutine PMWellPreSolveFlow

! ************************************************************************** !

subroutine PMWellSolveGeoTh(pm_well,time,qi_coupling,ierr)

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscReal :: time
  PetscBool :: qi_coupling
  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: curr_time, curr_time_converted

  curr_time = pm_well%option%time
  curr_time_converted = curr_time/pm_well%output_option%tconv

  ierr = 0 ! If pm_well is not set to zero, TS_STOP_FAILURE occurs if the solve
           ! routines are not entered, either due to an inactive well or due
           ! to being on a process that doesn't contain a well segment.

  write(out_string,'(" FLOW Step          Quasi-implicit wellbore flow &
                      &coupling is being used.")')
  call PrintMsg(pm_well%option,out_string)

  call PMWellCalcCumulativeQFlux(pm_well)

end subroutine PMWellSolveGeoTh

! ************************************************************************** !

subroutine PMWellSolveGeoThQI(this,time,ierr)

  implicit none

  class(pm_well_geoth_qi_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ! This is called when the well model is executed, after TH flow and
  ! after the fracture process model. It doesn't do much, since all the
  ! action should be happening with the TH mode solve calling the well model.

  call PMWellSolveGeoTh(this,time,this%update_for_flow_qi_coupling,ierr)

end subroutine PMWellSolveGeoThQI

! ************************************************************************** !

subroutine GeoThWellSolveFlowSequential(pm_well,perturbation_index,ierr)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  PetscReal :: res_fixed(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)
  PetscReal :: res(pm_well%flow_soln%ndof)
  PetscInt :: n_iter
  PetscInt :: i,istart, iend

  ! PLACEHOLDER (this is pm_well%SolveFlow => PMWellSolveFlowGeoThQI)
  ! see WIPPWellSolveFlowSequential in wipp_well.F90

  ! update well index
  call PMWellComputeWellIndex(pm_well)

  do while (pm_well%cumulative_dt_flow < pm_well%realization%option%flow_dt)

    ! update the well src/sink Q vector at start of time step
    call PMWellUpdateWellQ(pm_well%well,pm_well%well%reservoir)

    call PMWellPreSolveFlow(pm_well)

    ! Fixed accumulation term
    res_fixed = 0.d0
    res = 0.d0
    do i = 1,pm_well%well_grid%nsegments
      call PMWellAccumulationFlow(pm_well,pm_well%well,i,res)
      istart = pm_well%flow_soln%ndof*(i-1)+1
      iend = pm_well%flow_soln%ndof*i
      res_fixed(istart:iend) = -1.d0 * res * pm_well%dt_flow
    enddo

    n_iter = 0

    do while (pm_well%flow_soln%not_converged)

      call PMWellNewtonFlow(pm_well)

      call PMWellCheckConvergenceFlow(pm_well,n_iter,res_fixed)

    enddo

    ! hack so we get out of while loop
    pm_well%cumulative_dt_flow = pm_well%realization%option%flow_dt

    call PMWellPostSolveFlow(pm_well)

  enddo

end subroutine GeoThWellSolveFlowSequential

! ************************************************************************** !

subroutine PMWellSolveFlowGeoThQI(this,perturbation_index,ierr)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_geoth_qi_type) :: this
  PetscInt :: perturbation_index
  PetscErrorCode :: ierr

  ! This is called by TH mode so that the well model can solve for pressure
  ! and the exchange fluxes.

  call GeoThWellSolveFlowSequential(this,perturbation_index,ierr)

end subroutine PMWellSolveFlowGeoThQI

! ************************************************************************** !

subroutine PMWellPostSolveFlow(pm_well)
  !
  ! Author: Jennifer M. Frederick, SNL
  ! Date: 06/25/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  character(len=MAXSTRINGLENGTH) :: out_string
  PetscReal :: cur_time_converted

  cur_time_converted = pm_well%option%time/pm_well%output_option%tconv

  WRITE(out_string,'(" PM Well FLOW Step Complete!    Time=",1pe12.5," &
                    &",a4,"Total Newton Its =",i8)') &
                    cur_time_converted,pm_well%output_option%tunit, &
                    pm_well%flow_soln%n_newton
  call PMWellPrint(pm_well,out_string)

end subroutine PMWellPostSolveFlow

! ************************************************************************** !

subroutine PMWellUpdateSolutionFlow(pm_well)
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! PLACEHOLDER

end subroutine PMWellUpdateSolutionFlow

! ************************************************************************** !

subroutine PMWellCheckConvergenceFlow(pm_well,n_iter,fixed_accum)
  !
  ! Checks flow solution convergence against prescribed tolerances.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well
  PetscInt :: n_iter
  PetscReal :: fixed_accum(pm_well%flow_soln%ndof*pm_well%well_grid%nsegments)

  type(well_soln_flow_type), pointer :: flow_soln
  character(len=MAXSTRINGLENGTH) :: out_string
  character(len=MAXSTRINGLENGTH) :: rsn_string

  flow_soln => pm_well%flow_soln

  n_iter = n_iter + 1
  flow_soln%n_newton = flow_soln%n_newton + 1

  ! PLACEHOLDER

  ! pretend the soln converged
  flow_soln%converged = PETSC_TRUE
  flow_soln%not_converged = PETSC_FALSE
  rsn_string = 'dummy!'

  if (pm_well%print_output) then
    if (flow_soln%converged) then
      out_string = ' WELL FLOW Solution converged!  ---> ' // trim(rsn_string)
      call PMWellPrint(pm_well,out_string)
    endif
  endif
  
end subroutine PMWellCheckConvergenceFlow

! ************************************************************************** !

subroutine PMWellUpdateWellQ(well,reservoir)
  !
  ! Updates the src/sink vector for the fluid object.
  !
  ! Author: Fill Me In
  ! Date: 01/01/2025
  !

  implicit none

  type(well_type) :: well
  type(well_reservoir_type), pointer :: reservoir

  ! PLACEHOLDER

end subroutine PMWellUpdateWellQ

! ************************************************************************** !

subroutine PMWellCutTimestepFlow(pm_well)
  !
  ! Author: Michael Nole
  ! Date: 01/24/2022
  !

  implicit none

  class(pm_well_sequential_type) :: pm_well

  ! TODO: could make this smarter or call smarter timestepping routines
  pm_well%dt_flow = pm_well%dt_flow / pm_well%flow_soln%ts_cut_factor
  pm_well%dt_flow = max(pm_well%dt_flow,pm_well%min_dt_flow)
  pm_well%well%pl = pm_well%flow_soln%prev_soln%pl
  pm_well%well%gas%s = pm_well%flow_soln%prev_soln%sg
  call pm_well%UpdateFlowProperties(PETSC_FALSE,ZERO_INTEGER)

end subroutine PMWellCutTimestepFlow

! ************************************************************************** !

subroutine PMWellPrint(pm_well,string)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 03/06/2025

  implicit none

  class(pm_well_sequential_type) :: pm_well
  character(len=MAXSTRINGLENGTH) :: string

  if (MINVAL(pm_well%well_comm%petsc_rank_list) == pm_well%option%myrank) then
    write(STDOUT_UNIT,'(a)') trim(string)
  endif

end subroutine PMWellPrint

! ************************************************************************** !

end module Geothermal_Well_class

! TODO:
!
! How to initialize the well variables?
!
! Add TH to factory_subsurface.F90(?)
!
! Look at Julia Camargo commits (2303b58) - she is adding an implicit
! hydrostatic well.