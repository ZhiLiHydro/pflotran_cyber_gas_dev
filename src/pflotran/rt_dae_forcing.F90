module RT_DAE_Forcing_module
  !
  ! Forcing of the GIRT_DAE model (RT_DAE_Model_class): what an evaluation
  ! at a time reads besides the state. It is captured from the realization
  ! one waypoint interval (segment) at a time: the cells' flow state and
  ! porosity, face and source fluxes, and boundary and source chemistry.
  ! The records of the segments serve an adjoint's replay (TS_Model_class).
  ! Coupled to MODE RICHARDS_DAE, the flow state comes from the flow model
  ! (RTDAEForcingFromFlowState) and only the chemistry is forcing. The only
  ! other case has no flow mode (the factory rejects classic flow): the
  ! flow state is the initial one, constant.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Realization_Subsurface_class
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use RT_DAE_Chem_module
  use TS_Flow_State_module
  use TS_Connections_module

  implicit none

  private

  ! The forcing of one waypoint interval [t_start, t_end], immutable once
  ! RTDAEForcingBegin returns. Its flow state, constant, is held by the
  ! forcing (gaux, gaux_bc).
  type, public :: rt_dae_segment_type
    PetscInt :: id = 0
    PetscReal :: t_start = 0.d0
    PetscReal :: t_end = 0.d0
    PetscReal, allocatable :: por(:)        ! (ngmax) porosity
    ! face fluxes at t_start (vel_*) and t_end (vel_*_end), linear in
    ! between: exact for a LINEAR uniform velocity dataset, whose knots are
    ! waypoints
    PetscBool :: transient_velocity = PETSC_FALSE
    ! Darcy fluxes [m/s], (nphase,nconn_int) and (nphase,nconn_bc)
    PetscReal, allocatable :: vel_int(:,:)
    PetscReal, allocatable :: vel_bc(:,:)
    PetscReal, allocatable :: vel_int_end(:,:)
    PetscReal, allocatable :: vel_bc_end(:,:)
    PetscReal, allocatable :: total_bc(:,:) ! (naqcomp,nconn_bc) [mol/L]
    PetscReal, allocatable :: total_ss(:,:) ! (naqcomp,nconn_ss)
    ! ghosted cell-centered Darcy velocities, for transverse dispersion only
    PetscReal, allocatable :: ccv(:,:,:)      ! (3,nphase,ngmax)
    PetscReal, allocatable :: ccv_end(:,:,:)
  end type rt_dae_segment_type

  type, public :: rt_dae_forcing_type
    class(realization_subsurface_type), pointer :: realization => null()
    type(option_type), pointer :: option => null()
    class(reaction_rt_type), pointer :: reaction => null()
    ! the model's connection table
    type(ts_connections_type), pointer :: conn => null()
    PetscInt :: ngmax = 0
    PetscInt :: naqcomp = 0
    PetscInt :: nphase = 0
    ! the segment (during a replay, the record in use) and the records
    type(rt_dae_segment_type) :: seg
    type(rt_dae_segment_type), allocatable :: records(:)
    ! flow state of the cells and boundaries
    type(global_auxvar_type), pointer :: gaux(:) => null()
    type(global_auxvar_type), pointer :: gaux_bc(:) => null()
    ! face fluxes and cell-centered velocities (null without transverse
    ! dispersion) at forcing_time; pointers, as RTComputeTransportCoefs
    ! takes them
    PetscReal, pointer :: vel_int(:,:) => null()
    PetscReal, pointer :: vel_bc(:,:) => null()
    PetscReal, pointer :: ccv(:,:,:) => null()
    ! the values above are at forcing_time of segment forcing_segment or,
    ! coupled, from the flow state of id forcing_flow_id
    PetscReal :: forcing_time = -MAX_DOUBLE
    PetscInt :: forcing_segment = -1
    PetscInt :: forcing_flow_id = -1
    ! chemistry scratch of the boundary totals
    type(reactive_transport_auxvar_type), pointer :: scratch => null()
    type(rt_dae_gamma_type) :: gam
  end type rt_dae_forcing_type

  public :: RTDAEForcingInit, &
            RTDAEForcingBegin, &
            RTDAEForcingStore, &
            RTDAEForcingRestore, &
            RTDAEForcingInterpolate, &
            RTDAEForcingFromFlowState, &
            RTDAEForcingDestroy

contains

! ************************************************************************** !

subroutine RTDAEForcingInit(this,realization,conn)
  !
  ! Allocates the flow state, a copy of the realization's, and the face
  ! fluxes over the connections of conn.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
  use Patch_module
  use Reactive_Transport_Aux_module, only : RTAuxVarInit

  implicit none

  type(rt_dae_forcing_type) :: this
  class(realization_subsurface_type), pointer :: realization
  type(ts_connections_type), target :: conn

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscInt :: ghosted_id, iconn

  this%realization => realization
  option => realization%option
  patch => realization%patch
  this%option => option
  this%reaction => realization%reaction
  this%conn => conn
  this%ngmax = patch%grid%ngmax
  this%naqcomp = this%reaction%naqcomp
  this%nphase = patch%aux%RT%rt_parameter%nphase

  allocate(this%gaux(this%ngmax))
  do ghosted_id = 1, this%ngmax
    call GlobalAuxVarInit(this%gaux(ghosted_id),option)
    call GlobalAuxVarCopy(patch%aux%Global%auxvars(ghosted_id), &
                          this%gaux(ghosted_id),option)
  enddo
  allocate(this%gaux_bc(max(conn%nconn_bc,1)))
  do iconn = 1, conn%nconn_bc
    call GlobalAuxVarInit(this%gaux_bc(iconn),option)
    call GlobalAuxVarCopy(patch%aux%Global%auxvars_bc(iconn), &
                          this%gaux_bc(iconn),option)
  enddo

  allocate(this%vel_int(this%nphase,max(conn%nconn_int,1)))
  allocate(this%vel_bc(this%nphase,max(conn%nconn_bc,1)))
  this%vel_int = 0.d0
  this%vel_bc = 0.d0
  if (patch%aux%RT%rt_parameter%calculate_transverse_dispersion) then
    allocate(this%ccv(3,this%nphase,this%ngmax))
    this%ccv = 0.d0
  endif

  allocate(this%scratch)
  call RTAuxVarInit(this%scratch,this%reaction,PETSC_FALSE,option)
  call RTDAEGammaInit(this%gam,this%reaction)

end subroutine RTDAEForcingInit

! ************************************************************************** !

subroutine RTDAEForcingBegin(this,coupled,t_start,t_end,discontinuous)
  !
  ! Captures the forcing on [t_start, t_end] into a new segment;
  ! discontinuous: it jumps from the previous one (the integrator history
  ! must restart). Reads the realization at t_start, then at t_end for a
  ! LINEAR uniform velocity dataset (TimeStorageUpdate cannot go back).
  ! Restores option%time; leaves the realization's transport conditions and
  ! boundary chemistry at t_start and its uniform velocities at the last
  ! time read. coupled: the flow state is not forcing (see
  ! RTDAEForcingFromFlowState).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Condition_module, only : TranConditionUpdate

  implicit none

  type(rt_dae_forcing_type) :: this
  PetscBool :: coupled
  PetscReal :: t_start
  PetscReal :: t_end
  PetscBool :: discontinuous

  type(rt_dae_segment_type) :: new
  PetscReal :: saved_time

  new%id = this%seg%id + 1
  new%t_start = t_start
  new%t_end = t_end

  saved_time = this%option%time
  this%option%time = t_start
  ! conditions and uniform velocities at t_start; their knots are waypoints
  call TranConditionUpdate(this%realization%transport_conditions,this%option)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif
  call RTDAECaptureFlow(this,new)
  call RTDAECaptureBoundaries(this,new)
  call RTDAECaptureSources(this,new)
  ! last: a LINEAR uniform velocity dataset moves on to t_end
  call RTDAECaptureFluxes(this,new)
  this%option%time = saved_time
  if (.not.coupled) call RTDAECheckForcing(this)

  discontinuous = PETSC_TRUE
  if (this%seg%id > 0) discontinuous = &
    RTDAESegmentJumps(this%seg,new,coupled,this%option)
  this%seg = new

  ! force the forcing to be rebuilt
  this%forcing_segment = -1
  this%forcing_time = -MAX_DOUBLE

end subroutine RTDAEForcingBegin

! ************************************************************************** !

subroutine RTDAECaptureFlow(this,seg)
  !
  ! Copies the cells' flow state into gaux and their porosity into seg.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module

  implicit none

  type(rt_dae_forcing_type) :: this
  type(rt_dae_segment_type) :: seg

  type(patch_type), pointer :: patch
  PetscInt :: ghosted_id

  patch => this%realization%patch
  allocate(seg%por(this%ngmax))
  do ghosted_id = 1, this%ngmax
    call GlobalAuxVarCopy(patch%aux%Global%auxvars(ghosted_id), &
                          this%gaux(ghosted_id),this%option)
    seg%por(ghosted_id) = patch%aux%Material%auxvars(ghosted_id)%porosity
  enddo

end subroutine RTDAECaptureFlow

! ************************************************************************** !

subroutine RTDAECaptureFluxes(this,seg)
  !
  ! Face fluxes at the segment start and, for a LINEAR uniform velocity
  ! dataset, at its end (see RTDAEForcingBegin);
  ! ghosted cell-centered velocities for transverse dispersion.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Reactive_Transport_module, only : RTGhostedCellCenteredVelocities

  implicit none

  type(rt_dae_forcing_type) :: this
  type(rt_dae_segment_type) :: seg

  type(patch_type), pointer :: patch
  PetscBool :: transverse
  PetscReal, pointer :: ccv(:,:,:)

  patch => this%realization%patch
  transverse = patch%aux%RT%rt_parameter%calculate_transverse_dispersion
  nullify(ccv)

  allocate(seg%vel_int(this%nphase,max(this%conn%nconn_int,1)))
  allocate(seg%vel_bc(this%nphase,max(this%conn%nconn_bc,1)))
  seg%vel_int = 0.d0
  seg%vel_bc = 0.d0
  call CaptureVelocities(seg%vel_int,seg%vel_bc)
  seg%vel_int_end = seg%vel_int
  seg%vel_bc_end = seg%vel_bc
  if (transverse) then
    call RTGhostedCellCenteredVelocities(this%realization,ccv)
    seg%ccv = ccv
    deallocate(ccv)
    seg%ccv_end = seg%ccv
  endif
  ! a STEP dataset is constant inside a segment (its knots are waypoints);
  ! evaluated at t_end it would already give the next piece
  if (RTDAELinearVelocityDataset(this) .and. seg%t_end > seg%t_start) then
    this%option%time = seg%t_end
    call RealizUpdateUniformVelocity(this%realization)
    call CaptureVelocities(seg%vel_int_end,seg%vel_bc_end)
    if (transverse) then
      call RTGhostedCellCenteredVelocities(this%realization,ccv)
      seg%ccv_end = ccv
      deallocate(ccv)
    endif
  endif
  seg%transient_velocity = any(seg%vel_int_end /= seg%vel_int) .or. &
                           any(seg%vel_bc_end /= seg%vel_bc)

contains

  subroutine CaptureVelocities(vel_int,vel_bc)
    ! the patch's face fluxes of the transport phases
    PetscReal :: vel_int(:,:), vel_bc(:,:)
    if (this%conn%nconn_int > 0) &
      vel_int(:,1:this%conn%nconn_int) = &
        patch%internal_velocities(1:this%nphase,1:this%conn%nconn_int)
    if (this%conn%nconn_bc > 0) &
      vel_bc(:,1:this%conn%nconn_bc) = &
        patch%boundary_velocities(1:this%nphase,1:this%conn%nconn_bc)
  end subroutine CaptureVelocities

end subroutine RTDAECaptureFluxes

! ************************************************************************** !

function RTDAELinearVelocityDataset(this)
  !
  ! Whether the uniform velocity varies linearly in time between its knots.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(rt_dae_forcing_type) :: this
  PetscBool :: RTDAELinearVelocityDataset

  RTDAELinearVelocityDataset = PETSC_FALSE
  if (.not.associated(this%realization%uniform_velocity_dataset)) return
  if (.not.associated(this%realization%uniform_velocity_dataset% &
                      time_storage)) return
  RTDAELinearVelocityDataset = &
    this%realization%uniform_velocity_dataset%time_storage% &
      time_interpolation_method == INTERPOLATION_LINEAR

end function RTDAELinearVelocityDataset

! ************************************************************************** !

subroutine RTDAECaptureBoundaries(this,seg)
  !
  ! Flow state and chemistry of the boundary connections. Totals come from
  ! the model's own closure at the constraints' free-ion molalities, so a
  ! boundary and a cell of equal composition have equal totals. A
  ! connection of an inactive cell has no constraint chemistry
  ! (RTUpdateAuxVars skips it); its total stays zero.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Reactive_Transport_module, only : RTUpdateAuxVars

  implicit none

  type(rt_dae_forcing_type) :: this
  type(rt_dae_segment_type) :: seg

  type(patch_type), pointer :: patch
  PetscInt :: iconn, status

  patch => this%realization%patch
  allocate(seg%total_bc(this%naqcomp,max(this%conn%nconn_bc,1)))
  seg%total_bc = 0.d0
  if (this%conn%nconn_bc == 0) return

  ! free-ion molalities of the boundary constraints (possibly equilibrated
  ! at the boundary state)
  call RTUpdateAuxVars(this%realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)
  do iconn = 1, this%conn%nconn_bc
    call GlobalAuxVarCopy(patch%aux%Global%auxvars_bc(iconn), &
                          this%gaux_bc(iconn),this%option)
    if (patch%imat(this%conn%bc_ghost(iconn)) <= 0) cycle
    this%scratch%pri_molal = patch%aux%RT%auxvars_bc(iconn)%pri_molal
    call RTDAEAqueousClosure(this%scratch,this%gaux_bc(iconn), &
                             this%reaction,this%option,this%gam,status)
    if (status /= RT_DAE_CHEM_OK) then
      this%option%io_buffer = 'GIRT_DAE: the chemistry of a boundary &
        &constraint could not be evaluated.'
      call PrintErrMsgByRank(this%option)
    endif
    seg%total_bc(:,iconn) = this%scratch%total(:,1)
  enddo

end subroutine RTDAECaptureBoundaries

! ************************************************************************** !

subroutine RTDAECaptureSources(this,seg)
  !
  ! Chemistry of the source/sink connections: the totals of their current
  ! constraints.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Coupler_module
  use Transport_Constraint_RT_module, only : TranConstraintRTGetAuxVar

  implicit none

  type(rt_dae_forcing_type) :: this
  type(rt_dae_segment_type) :: seg

  type(coupler_type), pointer :: coupler
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar_ss
  PetscInt :: iconn, n

  allocate(seg%total_ss(this%naqcomp,max(this%conn%nconn_ss,1)))
  seg%total_ss = 0.d0
  n = 0
  coupler => this%realization%patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    rt_auxvar_ss => TranConstraintRTGetAuxVar(coupler%tran_condition% &
                                                cur_constraint_coupler)
    do iconn = 1, coupler%connection_set%num_connections
      n = n + 1
      if (associated(rt_auxvar_ss)) seg%total_ss(:,n) = rt_auxvar_ss%total(:,1)
    enddo
    coupler => coupler%next
  enddo

end subroutine RTDAECaptureSources

! ************************************************************************** !

subroutine RTDAECheckForcing(this)
  !
  ! Rejects forcing outside the model: an active cell whose liquid
  ! saturation is below the transport minimum.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module

  implicit none

  type(rt_dae_forcing_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscInt :: local_id, ghosted_id

  patch => this%realization%patch
  grid => patch%grid
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (this%gaux(ghosted_id)%sat(1) < rt_min_saturation) then
      this%option%io_buffer = 'MODE GIRT_DAE does not yet support cells &
        &that dry out (liquid saturation below the transport minimum).'
      call PrintErrMsgByRank(this%option)
    endif
  enddo

end subroutine RTDAECheckForcing

! ************************************************************************** !

function RTDAESegmentJumps(old,new,flow_coupled,option)
  !
  ! Whether the forcing jumps where segment new meets old: the boundary and
  ! source chemistry and, uncoupled, the face fluxes (the flow state is
  ! constant). Changes within round-off (a boundary's molarity is rebuilt
  ! from its molality at the current boundary density) are not jumps.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(rt_dae_segment_type) :: old
  type(rt_dae_segment_type) :: new
  PetscBool :: flow_coupled
  type(option_type) :: option
  PetscBool :: RTDAESegmentJumps

  PetscErrorCode :: ierr

  RTDAESegmentJumps = &
    Differ(new%total_bc,old%total_bc) .or. &
    Differ(new%total_ss,old%total_ss)
  if (.not.flow_coupled) RTDAESegmentJumps = RTDAESegmentJumps .or. &
    Differ(new%vel_int,old%vel_int_end) .or. &
    Differ(new%vel_bc,old%vel_bc_end)
  call MPI_Allreduce(MPI_IN_PLACE,RTDAESegmentJumps,ONE_INTEGER_MPI, &
                     MPI_C_BOOL,MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)

contains

  function Differ(a,b)
    PetscReal :: a(:,:), b(:,:)
    PetscBool :: Differ
    PetscReal, parameter :: rtol = 1.d-10
    Differ = any(abs(a-b) > rtol*max(abs(a),abs(b)))
  end function Differ

end function RTDAESegmentJumps

! ************************************************************************** !

subroutine RTDAEForcingStore(this,k)
  !
  ! Keeps the current segment as record k (see TS_Model_class).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(rt_dae_forcing_type) :: this
  PetscInt :: k

  type(rt_dae_segment_type), allocatable :: grown(:)

  if (.not.allocated(this%records)) allocate(this%records(16))
  if (k > size(this%records)) then
    allocate(grown(2*size(this%records)))
    grown(1:size(this%records)) = this%records
    call move_alloc(grown,this%records)
  endif
  this%records(k) = this%seg

end subroutine RTDAEForcingStore

! ************************************************************************** !

subroutine RTDAEForcingRestore(this,k)
  !
  ! Makes record k the current segment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(rt_dae_forcing_type) :: this
  PetscInt :: k

  this%seg = this%records(k)
  this%forcing_segment = -1
  this%forcing_time = -MAX_DOUBLE

end subroutine RTDAEForcingRestore

! ************************************************************************** !

subroutine RTDAEForcingInterpolate(this,time,maux,changed)
  !
  ! Uncoupled: the porosity (into maux) and face fluxes at time. A pure
  ! function of (segment, time), recomputed only when either changes;
  ! changed: they were.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(rt_dae_forcing_type) :: this
  PetscReal :: time
  type(material_auxvar_type) :: maux(:)
  PetscBool :: changed

  PetscInt :: ghosted_id
  PetscReal :: wv

  changed = PETSC_FALSE
  if (this%forcing_segment == this%seg%id) then
    if (.not.this%seg%transient_velocity) return
    if (time == this%forcing_time) return
  endif
  changed = PETSC_TRUE
  this%forcing_segment = this%seg%id
  this%forcing_time = time

  do ghosted_id = 1, this%ngmax
    maux(ghosted_id)%porosity = this%seg%por(ghosted_id)
  enddo
  wv = 0.d0
  if (this%seg%transient_velocity) then
    wv = (time-this%seg%t_start)/(this%seg%t_end-this%seg%t_start)
    wv = max(0.d0,min(1.d0,wv))
  endif
  this%vel_int = (1.d0-wv)*this%seg%vel_int + wv*this%seg%vel_int_end
  this%vel_bc = (1.d0-wv)*this%seg%vel_bc + wv*this%seg%vel_bc_end
  if (associated(this%ccv)) &
    this%ccv = (1.d0-wv)*this%seg%ccv + wv*this%seg%ccv_end

end subroutine RTDAEForcingInterpolate

! ************************************************************************** !

subroutine RTDAEForcingFromFlowState(this,fs,maux,changed)
  !
  ! Coupled: the flow state, porosity (into maux) and face fluxes of fs,
  ! recomputed only when it or the segment changes; changed: they were.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(rt_dae_forcing_type) :: this
  type(ts_flow_state_type) :: fs
  type(material_auxvar_type) :: maux(:)
  PetscBool :: changed

  PetscInt :: ghosted_id, iconn

  changed = PETSC_FALSE
  if (this%forcing_segment == this%seg%id .and. &
      this%forcing_flow_id == fs%id) return
  changed = PETSC_TRUE
  this%forcing_segment = this%seg%id
  this%forcing_flow_id = fs%id
  do ghosted_id = 1, this%ngmax
    this%gaux(ghosted_id)%sat(1) = fs%sat(ghosted_id)
    this%gaux(ghosted_id)%den(1) = fs%den(ghosted_id)
    this%gaux(ghosted_id)%den_kg(1) = fs%den_kg(ghosted_id)
    this%gaux(ghosted_id)%pres(1) = fs%pres(ghosted_id)
    maux(ghosted_id)%porosity = fs%por(ghosted_id)
  enddo
  do iconn = 1, this%conn%nconn_bc
    this%gaux_bc(iconn)%sat(1) = fs%sat_bc(iconn)
    this%gaux_bc(iconn)%den(1) = fs%den_bc(iconn)
    this%gaux_bc(iconn)%den_kg(1) = fs%den_kg_bc(iconn)
    this%gaux_bc(iconn)%pres(1) = fs%pres_bc(iconn)
  enddo
  this%vel_int = 0.d0
  this%vel_bc = 0.d0
  if (this%conn%nconn_int > 0) &
    this%vel_int(1,1:this%conn%nconn_int) = fs%v_int(1:this%conn%nconn_int)
  if (this%conn%nconn_bc > 0) &
    this%vel_bc(1,1:this%conn%nconn_bc) = fs%v_bc(1:this%conn%nconn_bc)

end subroutine RTDAEForcingFromFlowState

! ************************************************************************** !

subroutine RTDAEForcingDestroy(this)
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
  use Reactive_Transport_Aux_module, only : RTAuxVarStrip

  implicit none

  type(rt_dae_forcing_type) :: this

  PetscInt :: i

  if (associated(this%vel_int)) deallocate(this%vel_int)
  if (associated(this%vel_bc)) deallocate(this%vel_bc)
  if (associated(this%ccv)) deallocate(this%ccv)
  if (associated(this%gaux)) then
    do i = 1, size(this%gaux)
      call GlobalAuxVarStrip(this%gaux(i))
    enddo
    deallocate(this%gaux)
  endif
  if (associated(this%gaux_bc)) then
    do i = 1, this%conn%nconn_bc
      call GlobalAuxVarStrip(this%gaux_bc(i))
    enddo
    deallocate(this%gaux_bc)
  endif
  if (allocated(this%records)) deallocate(this%records)
  if (associated(this%scratch)) then
    call RTAuxVarStrip(this%scratch)
    deallocate(this%scratch)
  endif
  call RTDAEGammaDestroy(this%gam)

end subroutine RTDAEForcingDestroy

end module RT_DAE_Forcing_module
