module Geothermal_Well_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
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
  use Material_Aux_module
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

  use Grid_module
  use Utility_module
  use Coupler_module
  use Connection_module
  use Input_Aux_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Global_Aux_module

  implicit none

  class(pm_well_sequential_type) :: pm_well

  type(option_type), pointer :: option
  type(grid_type), pointer :: res_grid
  type(well_grid_type), pointer :: well_grid
  type(coupler_type), pointer :: source_sink
  type(input_type) :: input_dummy
  class(realization_subsurface_type), pointer :: realization
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: num_new_source_sinks, offset

  call PMWellSetupBase(pm_well)
  option => pm_well%option
  realization => pm_well%realization
  res_grid => realization%patch%grid
  well_grid => pm_well%well_grid

  pm_well%flow_soln%ndof = pm_well%nphase

  if (option%itranmode /= NULL_MODE) then
    option%io_buffer ="The WELLBORE_MODEL, type WELL_MODEL_TH_GEOTHERMAL &
      &doesn't yet work with a transport mode. You need to remove the &
      &transport process model from your simulation for the WELLBORE_MODEL, &
      &type WELL_MODEL_TH_GEOTHERMAL to work."
    call PrintErrMsg(option)
  endif

  ! add a reservoir src/sink coupler for each well segment
  num_new_source_sinks = 0
  if (associated(realization%patch%aux%Global%auxvars_ss)) then
    offset = realization%patch%aux%Global%num_aux_ss
  else
    offset = 0
  endif
  do k = 1,well_grid%nsegments

    source_sink => CouplerCreate()
    source_sink%itype = SRC_SINK_COUPLER_TYPE
    source_sink%name = trim(pm_well%name) // '_well_segment_' // StringWrite(k)

    ! ----- flow: begin ------------------
    source_sink%flow_condition_name = trim(pm_well%name) // &
                          '_well_segment_' // StringWrite(k) // '_flow_srcsink'
    source_sink%flow_condition => FlowConditionCreate(option)
    source_sink%flow_condition%name = source_sink%flow_condition_name
    if (well_grid%h_rank_id(k) /= option%myrank) cycle

    source_sink%flow_condition%general => FlowGeneralConditionCreate(option)
    string = 'RATE'
    source_sink%flow_condition%general%rate => FlowGeneralSubConditionPtr( &
      input_dummy,string,source_sink%flow_condition%general,option)
    source_sink%flow_condition%general%rate%itype = SCALED_MASS_RATE_SS
    source_sink%flow_condition%general%liquid_pressure => &
          FlowGeneralSubConditionPtr(input_dummy,string,source_sink% &
                                      flow_condition%general,option)
    source_sink%flow_condition%general%gas_pressure => &
          FlowGeneralSubConditionPtr(input_dummy,string,source_sink% &
                                    flow_condition%general,option)
    allocate(source_sink%flow_condition%general%rate%dataset%rarray(4))
    source_sink%flow_condition%general%rate%dataset%rarray(:) = 0.d0

    source_sink%flow_condition%well => FlowSubConditionCreate(ONE_INTEGER)
    ! ----- flow: end --------------------

    source_sink%connection_set => &
      ConnectionCreate(1,GENERIC_CONNECTION_TYPE,res_grid%itype)
    source_sink%connection_set%id_dn = well_grid%h_local_id(k)
    if (well_grid%h_local_id(k) < 0) then
      source_sink%connection_set%num_connections = 0
    else
      num_new_source_sinks = num_new_source_sinks + 1
      source_sink%connection_set%offset = num_new_source_sinks + offset - 1
    endif

    call CouplerAddToList(source_sink, &
                          pm_well%realization%patch%source_sink_list)
    nullify(source_sink)
  enddo

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

  ! PLACEHOLDER

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

  ! PLACEHOLDER (called by PMWellSolveGeoThQI)

end subroutine PMWellSolveGeoTh

! ************************************************************************** !

subroutine PMWellSolveGeoThQI(this,time,ierr)

  implicit none

  class(pm_well_geoth_qi_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

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

  ! PLACEHOLDER (this is pm_well%SolveFlow => PMWellSolveFlowGeoThQI)

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

  ! PLACEHOLDER
  
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
! Add case(TH_MODE) to pm_well.F90,PMWellUpdateReservoirSrcSinkFlow(). Note, a
! placeholder has already been inserted into the right spot, you just need to
! put in the content for updating the source_sink with Q values and pm_well 
! with th auxvar values.
!
! Add TH to patch.F90, factory_subsurface.F90(?)
! 