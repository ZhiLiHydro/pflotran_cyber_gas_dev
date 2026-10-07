module Adjoint_RT_DAE_class
  !
  ! Discrete adjoint of MODE GIRT_DAE (RT_DAE_Model_class), alone or with
  ! RICHARDS_DAE flow (Flow_RT_DAE_Model_class), for the ADJOINT block's
  ! costs (functions of the final state) and parameters.
  !
  ! A reaction parameter's cell derivatives, at the model's converged
  ! chemistry, go into the model's rows:
  !   F_N = dN/dt + flux + R - Q  ->  dF_N/dp =  dR/dp
  !   F_K = dK/dt - j             ->  dF_K/dp = -dj/dp
  ! A material permeability (with flow) acts through the Darcy fluxes on
  ! the flow rows and, by advection and dispersion, on the transport rows.
  ! The gradient at the initial state is pulled back to the physical
  ! initial data and, with flow, the initial pressures.
  !
  ! Output (.out file; FULL_FIELD also writes binary Vecs):
  !   === ADJOINT SUMMARY ===
  !   Cost J_i = <name>
  !     J_i(T)              = <value>
  !     dJ_i/dp_j [<name>] = <value>
  !     max |dJ_i/du_phys(0)| = <value>     (primary species)
  !     max |dJ_i/dK(0)| = <value>          (kinetic inventories, if any)
  !     max |dJ_i/dP(0)| = <value>          (initial pressures, with flow)
  !   === END ADJOINT ===
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use PFLOTRAN_Constants_module
  use TS_Adjoint_class
  use TS_Model_class
  use RT_DAE_Model_class
  use RT_DAE_Flow_Coupling_module
  use Flow_RT_DAE_Model_class
  use Richards_DAE_Model_class
  use Adjoint_Cost_module, only : adjoint_cost_type
  use Adjoint_Param_module
  use Realization_Subsurface_class
  use Option_module

  implicit none

  private

  type, public, extends(ts_adjoint_type) :: adjoint_rt_dae_type
    class(realization_subsurface_type), pointer :: realization => null()
    ! -adjoint_debug_jacobianp: check the first dF/dp by central differences
    PetscBool :: debug_jacobianp = PETSC_FALSE
    PetscReal :: debug_eps = 1.d-6
    PetscBool :: debug_done = PETSC_FALSE
  contains
    procedure, public :: Setup => AdjRTDAESetup
    procedure, public :: CheckTerminalTime => AdjRTDAECheckTerminalTime
    procedure, public :: TerminalGradients => AdjRTDAETerminalGradients
    procedure, public :: ParameterJacobian => AdjRTDAEParameterJacobian
    procedure, public :: Report => AdjRTDAEReport
    procedure, public :: GetParameters => AdjRTDAEGetParameters
    procedure, public :: SetParameters => AdjRTDAESetParameters
    procedure, public :: Destroy => AdjRTDAEDestroy
  end type adjoint_rt_dae_type

  public :: AdjointRTDAECreate

contains

! ************************************************************************** !

function AdjointRTDAECreate(option)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  implicit none

  type(option_type), pointer :: option
  class(adjoint_rt_dae_type), pointer :: AdjointRTDAECreate

  allocate(AdjointRTDAECreate)
  call TSAdjInit(AdjointRTDAECreate,option)

end function AdjointRTDAECreate

! ************************************************************************** !

function AdjRTDAEModel(this,model)
  !
  ! The GIRT_DAE model: the model, or the transport part of a coupled one.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type), target :: model
  class(rt_dae_model_type), pointer :: AdjRTDAEModel

  nullify(AdjRTDAEModel)
  select type(model)
    class is(rt_dae_model_type)
      AdjRTDAEModel => model
    class is(flow_rt_dae_model_type)
      AdjRTDAEModel => model%rt
    class default
      call PrintErrMsg(this%option,'The GIRT_DAE adjoint was given another &
        &equation model.')
  end select

end function AdjRTDAEModel

! ************************************************************************** !

function AdjRTDAECoupled(model)
  !
  ! The coupled flow and transport model, or null.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(ts_model_type), target :: model
  class(flow_rt_dae_model_type), pointer :: AdjRTDAECoupled

  nullify(AdjRTDAECoupled)
  select type(model)
    class is(flow_rt_dae_model_type)
      AdjRTDAECoupled => model
  end select

end function AdjRTDAECoupled

! ************************************************************************** !

subroutine AdjRTDAESetup(this,model,final_time)
  !
  ! Resolves the costs and parameters; the model rejects unsupported physics.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Adjoint_Cost_module, only : adj_cost_list, AdjCostCount
  use Adjoint_Param_module, only : adj_param_list, AdjParamCount

  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type), target :: model
  PetscReal :: final_time

  class(rt_dae_model_type), pointer :: rt
  type(adjoint_cost_type), pointer :: cur_cost
  type(adjoint_param_type), pointer :: cur_param
  type(option_type), pointer :: option
  PetscBool :: flag, valid
  PetscErrorCode :: ierr

  option => this%option
  rt => AdjRTDAEModel(this,model)
  this%realization => rt%realization
  ! a BAIJ factor of block size 15 in natural ordering does a forward solve
  ! for a transposed one (PETSc 3.24.5); the condensed operator does not
  ! rely on it
  if (rt%assembled_jacobian .and. rt%ndof == 15) then
    call PrintErrMsg(option,'ADJOINT with -rt_dae_assembled_jacobian needs &
      &a block size (twice the components plus the kinetic inventories) &
      &other than 15: PETSc solves the transpose of such a factor wrongly.')
  endif

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-adjoint_debug_jacobianp',this%debug_jacobianp, &
                           ierr);CHKERRQ(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           '-adjoint_debug_jacobianp_eps',this%debug_eps, &
                           flag,ierr);CHKERRQ(ierr)
  if (.not.(this%debug_eps > 0.d0)) then
    call PrintErrMsg(option,'ADJOINT debug epsilon must be positive &
      &(-adjoint_debug_jacobianp_eps).')
  endif

  cur_cost => adj_cost_list
  do while (associated(cur_cost))
    call cur_cost%Setup(this%realization,option)
    call cur_cost%ValidateTime(final_time,valid)
    if (.not.valid) then
      call PrintErrMsg(option,'MISFIT_LEAST_SQUARES final observation time &
        &must match the simulation final time; multi-time objectives are not &
        &supported.')
    endif
    cur_cost => cur_cost%next
  enddo
  cur_param => adj_param_list
  do while (associated(cur_param))
    if (cur_param%kind == ADJ_PARAM_PERMEABILITY .and. &
        .not.associated(AdjRTDAECoupled(model))) then
      call PrintErrMsg(option,'PARAMETER PERMEABILITY requires flow in &
        &the equation model (MODE RICHARDS_DAE with GIRT_DAE); GIRT_DAE &
        &alone has no flow.')
    endif
    if (cur_param%kind == ADJ_PARAM_PERMEABILITY) &
      call AdjRTDAECheckScaledSources(this%realization%patch,option)
    call cur_param%Setup(this%realization,option)
    cur_param => cur_param%next
  enddo
  this%ncosts = AdjCostCount()
  this%nparams = AdjParamCount()

  ! the parameter Jacobian reads the chemistry of the model's last closure
  call RTDAEKeepActivities(rt)

end subroutine AdjRTDAESetup

! ************************************************************************** !

subroutine AdjRTDAECheckScaledSources(patch,option)
  !
  ! Stops at a source/sink whose rate is scaled by permeability
  ! (SCALED_MASS_RATE or SCALED_VOLUMETRIC_RATE, unless by VOLUME): dF/dk
  ! omits the derivative of the scaling.
  !
  ! Author: Hansol Suh
  ! Date: 10/13/26
  !
  use Patch_module
  use Coupler_module
  use Condition_module, only : SCALE_BY_VOLUME

  implicit none

  type(patch_type) :: patch
  type(option_type) :: option

  type(coupler_type), pointer :: coupler

  coupler => patch%source_sink_list%first
  do while (associated(coupler))
    if (associated(coupler%flow_condition)) then
      if (associated(coupler%flow_condition%rate)) then
        select case(coupler%flow_condition%rate%itype)
          case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
            if (coupler%flow_condition%rate%isubtype /= SCALE_BY_VOLUME) &
              call PrintErrMsg(option,'PARAMETER PERMEABILITY does not &
                &support source/sink ' // trim(coupler%name) // ': its &
                &rate is scaled by permeability (scale by VOLUME instead).')
        end select
      endif
    endif
    coupler => coupler%next
  enddo

end subroutine AdjRTDAECheckScaledSources

! ************************************************************************** !

subroutine AdjRTDAECheckTerminalTime(this,time,valid)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  use Adjoint_Cost_module, only : adj_cost_list

  implicit none

  class(adjoint_rt_dae_type) :: this
  PetscReal :: time
  PetscBool :: valid

  type(adjoint_cost_type), pointer :: cur_cost
  PetscBool :: cost_valid

  valid = PETSC_TRUE
  cur_cost => adj_cost_list
  do while (associated(cur_cost))
    call cur_cost%ValidateTime(time,cost_valid)
    valid = valid .and. cost_valid
    cur_cost => cur_cost%next
  enddo

end subroutine AdjRTDAECheckTerminalTime

! ************************************************************************** !

subroutine AdjRTDAETerminalGradients(this,model,time,Y)
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !

  use Adjoint_Cost_module, only : adj_cost_list

  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type), target :: model
  PetscReal :: time
  Vec :: Y

  type(adjoint_cost_type), pointer :: cur_cost
  class(rt_dae_model_type), pointer :: rt
  class(flow_rt_dae_model_type), pointer :: coupled
  Vec :: Y_f, Y_t, lam_f, lam_t
  PetscInt :: i
  PetscErrorCode :: ierr

  rt => AdjRTDAEModel(this,model)
  coupled => AdjRTDAECoupled(model)
  ! the costs are functions of the transport state
  if (associated(coupled)) call FlowRTDAEGetParts(coupled,Y,Y_f,Y_t)
  i = 0
  cur_cost => adj_cost_list
  do while (associated(cur_cost))
    i = i + 1
    call VecZeroEntries(this%lambda(i),ierr);CHKERRQ(ierr)
    if (associated(coupled)) then
      call FlowRTDAEGetParts(coupled,this%lambda(i),lam_f,lam_t)
      call cur_cost%FillTerminalLambda(Y_t,lam_t,rt%ndof,this%option)
      call FlowRTDAERestoreParts(coupled,this%lambda(i),lam_f,lam_t)
    else
      call cur_cost%FillTerminalLambda(Y,this%lambda(i),rt%ndof,this%option)
    endif
    this%cost_value(i) = cur_cost%J_value
    cur_cost => cur_cost%next
  enddo
  if (associated(coupled)) call FlowRTDAERestoreParts(coupled,Y,Y_f,Y_t)

end subroutine AdjRTDAETerminalGradients

! ************************************************************************** !

subroutine AdjRTDAEParameterJacobian(this,model,time,Y,Ydot,Jacp,status)
  !
  ! Jacp = dF/dp at (time,Y): reaction parameters' cell derivatives and,
  ! with flow, permeabilities' flux derivatives.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Adjoint_Param_module, only : adj_param_list
  use Patch_module
  use Grid_module

  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type), target :: model
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  Mat :: Jacp
  PetscInt :: status

  class(rt_dae_model_type), pointer :: rt
  class(flow_rt_dae_model_type), pointer :: coupled
  type(adjoint_param_type), pointer :: cur_param
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  Vec :: Y_f, Y_t
  PetscReal, allocatable :: dR(:), dj(:), values(:)
  PetscInt, allocatable :: rows(:)
  PetscInt :: local_id, ghosted_id, row_start, row_end, nc, nk, n, i
  PetscInt :: tran_start, col, cell_status, param_status
  PetscBool :: invalid(this%nparams)
  PetscErrorCode :: ierr

  rt => AdjRTDAEModel(this,model)
  coupled => AdjRTDAECoupled(model)
  patch => rt%realization%patch
  grid => patch%grid
  nc = rt%ncomp
  nk = rt%nkin
  n = nc + nk
  allocate(dR(rt%ndof),dj(rt%ndof),values(rt%ndof),rows(rt%ndof))
  invalid = PETSC_FALSE
  status = TS_MODEL_OK

  call MatZeroEntries(Jacp,ierr);CHKERRQ(ierr)
  call MatGetOwnershipRange(Jacp,row_start,row_end,ierr);CHKERRQ(ierr)
  ! the first transport row of this process
  tran_start = row_start
  if (associated(coupled)) then
    tran_start = row_start + coupled%nflow
    call FlowRTDAEGetParts(coupled,Y,Y_f,Y_t)
    call FlowRTDAESetFlowState(coupled,time,Y_f,status)
    call RTDAEPrepareReactions(rt,time,Y_t,cell_status)
  else
    call RTDAEPrepareReactions(rt,time,Y,cell_status)
  endif
  status = max(status,cell_status)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.rt%cell_ok(ghosted_id)) cycle
    call RTDAELoadCell(rt,local_id)
    ! rows N, then K, of the cell (0-based global)
    do i = 1, n
      rows(i) = tran_start + (local_id-1)*rt%ndof + nc + i - 1
    enddo
    col = 0
    cur_param => adj_param_list
    do while (associated(cur_param))
      ! a permeability is a flux parameter, below
      if (cur_param%kind /= ADJ_PARAM_PERMEABILITY) then
        dR(1:nc) = 0.d0
        dj(1:max(nk,1)) = 0.d0
        call AdjRTDAECellDerivative(rt,cur_param,local_id,dR(1:nc), &
                                    dj(1:nk),param_status)
        if (param_status /= 0) invalid(col+1) = PETSC_TRUE
        values(1:nc) = dR(1:nc)
        if (nk > 0) values(nc+1:n) = -dj(1:nk)
        call MatSetValues(Jacp,n,rows,ONE_INTEGER,[col],values, &
                          INSERT_VALUES,ierr);CHKERRQ(ierr)
      endif
      col = col + 1
      cur_param => cur_param%next
    enddo
  enddo

  if (associated(coupled)) then
    col = 0
    cur_param => adj_param_list
    do while (associated(cur_param))
      if (cur_param%kind == ADJ_PARAM_PERMEABILITY) &
        call AdjRTDAEPermeabilityColumn(coupled,cur_param,col,Jacp,row_start)
      col = col + 1
      cur_param => cur_param%next
    enddo
    call FlowRTDAERestoreParts(coupled,Y,Y_f,Y_t)
  endif
  call MatAssemblyBegin(Jacp,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Jacp,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! errors are collective: the cell loop is rank local
  call MPI_Allreduce(MPI_IN_PLACE,status,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,this%option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(MPI_IN_PLACE,invalid,this%nparams,MPI_C_BOOL, &
                     MPI_LOR,this%option%mycomm,ierr);CHKERRQ(ierr)
  col = 0
  cur_param => adj_param_list
  do while (associated(cur_param))
    col = col + 1
    if (invalid(col)) then
      this%option%io_buffer = 'ADJOINT: the derivative with respect to ' // &
        trim(cur_param%name) // ' is not defined at a state of the forward &
        &run (HALF_SATURATION_CONSTANT requires K + concentration - &
        &threshold > 0 in every active cell).'
      call PrintErrMsg(this%option)
    endif
    cur_param => cur_param%next
  enddo

  if (this%debug_jacobianp .and. .not.this%debug_done) then
    this%debug_done = PETSC_TRUE
    call AdjRTDAECheckParameterJacobian(this,model,time,Y,Ydot,Jacp)
  endif
  deallocate(dR,dj,values,rows)

end subroutine AdjRTDAEParameterJacobian

! ************************************************************************** !

subroutine AdjRTDAECellDerivative(rt,param,local_id,dR,dj,status)
  !
  ! Derivatives with respect to reaction parameter param of the reaction
  ! terms of active cell local_id, at the chemistry RTDAELoadCell left in
  ! rt%scratch:
  !   dR(ncomp)  component sinks [mol/s] (consumption positive)
  !   dj(nkin)   kinetic transfer rates into the kinetic inventories [mol/s]
  ! Both arrive zeroed. A rate linear in its rate constant has as
  ! derivative the rate at unit rate constant, so the model's rate laws are
  ! evaluated with unit_rate. status is nonzero where the derivative is
  ! undefined; the caller reports it collectively.
  !
  ! Author: Hansol Suh
  ! Date: 05/27/26
  !
  use Reaction_Mineral_Aux_module
  use Reaction_Mineral_module, only : ReactionMnrlTSTSimpleRate
  use Reaction_Immobile_module, only : ReactionImDecay

  implicit none

  class(rt_dae_model_type) :: rt
  type(adjoint_param_type) :: param
  PetscInt :: local_id
  PetscReal :: dR(:)
  PetscReal :: dj(:)
  PetscInt :: status

  type(mineral_type), pointer :: mineral
  PetscReal :: ln_act(rt%naqcomp), jac(size(dR),size(dR))
  PetscReal :: lnQK, QK, rate_constant, rate
  PetscInt :: ghosted_id, imnrl, i, icomp
  PetscBool :: may_dissolve

  status = 0
  ghosted_id = rt%realization%patch%grid%nL2G(local_id)
  select case(param%kind)
    case(ADJ_PARAM_RATE_CONSTANT)
      if (param%kinetics == ADJ_KINETICS_MINERAL) then
        ! as in RTDAECellReactions, at k = 1
        mineral => rt%reaction%mineral
        imnrl = param%imnrl
        ln_act = log(rt%scratch%pri_molal(1:rt%naqcomp)* &
                     rt%scratch%pri_act_coef(1:rt%naqcomp))
        may_dissolve = rt%mnrl_regime(imnrl,local_id) /= RT_DAE_MNRL_ABSENT
        call ReactionMnrlTSTSimpleRate(mineral,imnrl,ln_act, &
                                       rt%scratch%ln_act_h2o, &
                                       rt%mnrl_area(imnrl,local_id), &
                                       may_dissolve,lnQK,QK,rate_constant, &
                                       rate,unit_rate=PETSC_TRUE)
        dj(imnrl) = rt%maux(ghosted_id)%volume*rate
        do i = 1, mineral%kinmnrlspecid_in_residual(0,imnrl)
          icomp = mineral%kinmnrlspecid_in_residual(i,imnrl)
          dR(icomp) = dR(icomp) + &
            mineral%kinmnrlstoich_in_residual(i,imnrl)*dj(imnrl)
        enddo
      else
        call AdjRTDAEMicrobialDerivative(rt,param,ghosted_id,dR,status)
      endif
    case(ADJ_PARAM_DECAY_RATE)
      call ReactionImDecay(dR,jac,PETSC_FALSE,rt%scratch, &
                           rt%forcing%gaux(ghosted_id),rt%maux(ghosted_id), &
                           rt%reaction,rt%option, &
                           irxn_filter=param%reaction_index, &
                           unit_rate=PETSC_TRUE)
    case default
      call AdjRTDAEMicrobialDerivative(rt,param,ghosted_id,dR,status)
  end select

end subroutine AdjRTDAECellDerivative

! ************************************************************************** !

subroutine AdjRTDAEMicrobialDerivative(rt,param,ghosted_id,dR,status)
  !
  ! dR/dp of a microbial reaction parameter in a cell, at the chemistry in
  ! rt%scratch. ReactionMicrobRate evaluates the parameter's reaction with
  ! the differentiated factor (rate constant, inhibition term) set to one,
  ! then the result is scaled by that factor's derivative, so it stays
  ! defined where the factor is zero. A half-saturation constant instead
  ! scales the full rate by d ln(monod)/dK = -1/(K + C - C_th) (status = 1
  ! where that denominator is nonpositive); a yield, the biomass row's
  ! growth yield*rate, has the rate. The reaction data are read live, so
  ! AdjParamSetValue is seen.
  !
  ! Author: Hansol Suh
  ! Date: 06/23/26
  !
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use Reaction_Microbial_module, only : ReactionMicrobRate, &
                                        ReactionMicrobConcentrations
  use Reaction_Microbial_Aux_module
  use Reaction_Inhibition_Aux_module

  implicit none

  class(rt_dae_model_type) :: rt
  type(adjoint_param_type) :: param
  PetscInt :: ghosted_id
  PetscReal :: dR(:)
  PetscInt :: status

  type(microbial_type), pointer :: microbial
  PetscReal :: residual(size(dR)), jacobian(size(dR),size(dR))
  PetscReal :: conc(rt%reaction%naqcomp), dconc_dmolal(rt%reaction%naqcomp)
  PetscReal :: concentration, denom, scale, rate, a, threshold
  PetscReal :: da, dw
  PetscInt :: component, row, omit, irxn, term
  PetscBool :: unit_rate

  microbial => rt%reaction%microbial
  status = 0
  irxn = param%reaction_index
  term = param%term
  unit_rate = param%kind == ADJ_PARAM_RATE_CONSTANT
  omit = 0
  if (param%kind == ADJ_PARAM_INHIBITION_CONSTANT .or. &
      param%kind == ADJ_PARAM_SMOOTHSTEP_INTERVAL) omit = term
  scale = 1.d0
  if (param%kind == ADJ_PARAM_HALF_SATURATION .or. omit /= 0) then
    if (param%kind == ADJ_PARAM_HALF_SATURATION) then
      component = microbial%monod_specid(term)
    else
      component = microbial%inhibition_specid(term)
    endif
    ! the concentration the rate law reads
    call ReactionMicrobConcentrations(rt%scratch,rt%forcing%gaux(ghosted_id), &
                                      microbial,conc,dconc_dmolal)
    concentration = conc(component)
    if (param%kind == ADJ_PARAM_HALF_SATURATION) then
      denom = microbial%monod_K(term)+concentration-microbial%monod_Cth(term)
      if (.not.ieee_is_finite(denom) .or. denom <= 0.d0) then
        status = 1
        return
      endif
      scale = -1.d0/denom
    else
      threshold = microbial%inhibition_C(term)
      a = abs(threshold)
      if (microbial%inhibition_type(term) == INHIBITION_MONOD) then
        scale = -sign(1.d0,threshold)*(concentration/(a+concentration))/ &
                (a+concentration)
      else
        call ReactionInhibitionSmoothstepParameters(concentration,threshold, &
               microbial%inhibition_C2(term),da,dw)
        scale = da
        if (param%kind == ADJ_PARAM_SMOOTHSTEP_INTERVAL) scale = dw
      endif
    endif
  endif
  residual = 0.d0
  jacobian = 0.d0
  call ReactionMicrobRate(residual,jacobian,PETSC_FALSE,rt%scratch, &
         rt%forcing%gaux(ghosted_id),rt%maux(ghosted_id),rt%reaction, &
         rt%option, &
         irxn_filter=irxn,rate_value=rate,unit_rate=unit_rate, &
         omit_inhibition=omit)
  if (param%kind == ADJ_PARAM_YIELD) then
    ! the biomass row of the reaction: its growth is yield*rate
    residual = 0.d0
    row = microbial%biomassid(irxn)
    if (row < 0) row = rt%reaction%offset_immobile-row
    residual(row) = -rate
  endif
  dR = residual*scale

end subroutine AdjRTDAEMicrobialDerivative

! ************************************************************************** !

subroutine AdjRTDAEPermeabilityColumn(coupled,param,col,Jacp,row_start)
  !
  ! Column col of Jacp, dF/dk of a material permeability at the models'
  ! last prepared states: flow rows directly, transport rows through the
  ! Darcy fluxes.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  class(flow_rt_dae_model_type) :: coupled
  type(adjoint_param_type) :: param
  PetscInt :: col
  Mat :: Jacp
  PetscInt :: row_start

  Vec :: column, c_f, c_t
  PetscReal, allocatable :: dv_int(:), dv_bc(:)
  PetscReal, pointer :: c_p(:)
  PetscInt, allocatable :: rows(:)
  PetscInt :: nlocal, i
  PetscErrorCode :: ierr

  allocate(dv_int(max(coupled%flow%conn%nconn_int,1)))
  allocate(dv_bc(max(coupled%flow%conn%nconn_bc,1)))
  call MatCreateVecs(Jacp,PETSC_NULL_VEC,column,ierr);CHKERRQ(ierr)
  call VecZeroEntries(column,ierr);CHKERRQ(ierr)
  call FlowRTDAEGetParts(coupled,column,c_f,c_t)
  call RichardsDAEPermeabilityDerivative(coupled%flow,param%material_id, &
                                         param%value,c_f,dv_int,dv_bc)
  call RTDAEAddVelocityDerivative(coupled%rt,dv_int,dv_bc,c_t)
  call FlowRTDAERestoreParts(coupled,column,c_f,c_t)
  call VecGetLocalSize(column,nlocal,ierr);CHKERRQ(ierr)
  allocate(rows(nlocal))
  do i = 1, nlocal
    rows(i) = row_start + i - 1
  enddo
  call VecGetArrayRead(column,c_p,ierr);CHKERRQ(ierr)
  call MatSetValues(Jacp,nlocal,rows,ONE_INTEGER,[col],c_p,INSERT_VALUES, &
                    ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(column,c_p,ierr);CHKERRQ(ierr)
  call VecDestroy(column,ierr);CHKERRQ(ierr)
  deallocate(rows,dv_int,dv_bc)

end subroutine AdjRTDAEPermeabilityColumn

! ************************************************************************** !

subroutine AdjRTDAECheckParameterJacobian(this,model,time,Y,Ydot,Jacp)
  !
  ! -adjoint_debug_jacobianp: each column of dF/dp against
  ! (F(p+h) - F(p-h))/(2h), h = eps*|p| (eps if p = 0). p is restored and
  ! the model cache invalidated, so the run is unchanged. An error above
  ! 1e-5 of the column's largest entry stops the run.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Adjoint_Param_module, only : adj_param_list

  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type) :: model
  PetscReal :: time
  Vec :: Y
  Vec :: Ydot
  Mat :: Jacp

  type(adjoint_param_type), pointer :: cur_param
  type(option_type), pointer :: option
  Vec :: Fp, Fm, column
  PetscReal :: p0, h, abs_err, col_norm, rel_err
  PetscInt :: col, status_p, status_m
  PetscErrorCode :: ierr

  option => this%option
  call VecDuplicate(Y,Fp,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,Fm,ierr);CHKERRQ(ierr)
  call VecDuplicate(Y,column,ierr);CHKERRQ(ierr)
  col = 0
  cur_param => adj_param_list
  do while (associated(cur_param))
    p0 = cur_param%value
    h = this%debug_eps*abs(p0)
    if (p0 == 0.d0) h = this%debug_eps
    call cur_param%SetValue(this%realization,p0+h)
    call model%InvalidateCache()
    call model%Evaluate(time,Y,Ydot,Fp,status_p)
    call cur_param%SetValue(this%realization,p0-h)
    call model%InvalidateCache()
    call model%Evaluate(time,Y,Ydot,Fm,status_m)
    call cur_param%SetValue(this%realization,p0)
    call model%InvalidateCache()
    call TSModelReduceStatus(status_p,option%mycomm)
    call TSModelReduceStatus(status_m,option%mycomm)
    call VecAXPY(Fp,-1.d0,Fm,ierr);CHKERRQ(ierr)
    call VecScale(Fp,0.5d0/h,ierr);CHKERRQ(ierr)
    call MatGetColumnVector(Jacp,column,col,ierr);CHKERRQ(ierr)
    call VecNorm(column,NORM_INFINITY,col_norm,ierr);CHKERRQ(ierr)
    call VecAXPY(Fp,-1.d0,column,ierr);CHKERRQ(ierr)
    call VecNorm(Fp,NORM_INFINITY,abs_err,ierr);CHKERRQ(ierr)
    rel_err = abs_err/max(col_norm,tiny(1.d0))
    write(option%io_buffer,'(a,es10.3,a,es10.3,a,es10.3)') &
      ' ADJOINT debug dF/dp [' // trim(cur_param%name) // &
      ']: max |column| = ',col_norm,', max abs diff = ',abs_err, &
      ', max rel diff = ',rel_err
    call PrintMsg(option)
    if (status_p /= TS_MODEL_OK .or. status_m /= TS_MODEL_OK .or. &
        .not.(rel_err <= 1.d-5)) then
      option%io_buffer = 'ADJOINT debug check FAILED: dF/dp [' // &
        trim(cur_param%name) // '] differs from its central difference.'
      call PrintErrMsg(option)
    endif
    col = col + 1
    cur_param => cur_param%next
  enddo
  call VecDestroy(Fp,ierr);CHKERRQ(ierr)
  call VecDestroy(Fm,ierr);CHKERRQ(ierr)
  call VecDestroy(column,ierr);CHKERRQ(ierr)

end subroutine AdjRTDAECheckParameterJacobian

! ************************************************************************** !

subroutine AdjRTDAEReport(this,model,time0,Y0)
  !
  ! Prints each cost's summary. FULL_FIELD also writes PETSc binary Vecs:
  ! dJ/dc0, entry (natural id-1)*ncomp + i-1 for component i, and,
  ! with flow, dJ/dP0, entry natural id-1.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Adjoint_module, only : adjoint_output_full_field
  use Adjoint_Cost_module, only : adj_cost_list
  use Adjoint_Param_module, only : adj_param_list
  use Grid_module

  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type), target :: model
  PetscReal :: time0
  Vec :: Y0

  class(rt_dae_model_type), pointer :: rt
  class(flow_rt_dae_model_type), pointer :: coupled
  type(adjoint_cost_type), pointer :: cur_cost
  type(adjoint_param_type), pointer :: cur_param
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  PetscReal, allocatable :: dJ_dc(:,:), dJ_dK(:,:), dJ_dP(:)
  PetscReal, pointer :: mu_p(:)
  PetscReal :: mu_value, max_c, max_K
  PetscInt :: icost, iparam, col_start, col_end, nc, nk, local_id, i, status
  Vec :: Y0_f, Y0_t, lam_f, lam_t
  Vec :: field
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: name
  PetscErrorCode :: ierr

  option => this%option
  rt => AdjRTDAEModel(this,model)
  grid => rt%realization%patch%grid
  nc = rt%ncomp
  nk = rt%nkin
  allocate(dJ_dc(nc,grid%nlmax))
  allocate(dJ_dK(max(nk,1),grid%nlmax))
  allocate(dJ_dP(grid%nlmax))
  coupled => AdjRTDAECoupled(model)
  if (associated(coupled)) then
    call FlowRTDAEGetParts(coupled,Y0,Y0_f,Y0_t)
    call FlowRTDAESetFlowState(coupled,time0,Y0_f,status)
  endif

  call PrintMsg(option,'=== ADJOINT SUMMARY ===')
  icost = 0
  cur_cost => adj_cost_list
  do while (associated(cur_cost))
    icost = icost + 1
    write(option%io_buffer,'("Cost J_",i0," = ",a)') cur_cost%id, &
      trim(cur_cost%name)
    call PrintMsg(option)
    write(option%io_buffer,'("  J_",i0,"(T)              = ",es24.16)') &
      cur_cost%id,cur_cost%J_value
    call PrintMsg(option)
    write(name,'("J_",i0)') cur_cost%id
    call this%KeepResult(name,cur_cost%J_value)

    ! mu is distributed like the columns of Jacp
    call VecGetOwnershipRange(this%mu(icost),col_start,col_end, &
                              ierr);CHKERRQ(ierr)
    call VecGetArrayRead(this%mu(icost),mu_p,ierr);CHKERRQ(ierr)
    iparam = 0
    cur_param => adj_param_list
    do while (associated(cur_param))
      iparam = iparam + 1
      mu_value = 0.d0
      if (iparam > col_start .and. iparam <= col_end) &
        mu_value = mu_p(iparam-col_start)
      call MPI_Allreduce(MPI_IN_PLACE,mu_value,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm, &
                         ierr);CHKERRQ(ierr)
      write(option%io_buffer,'("  dJ_",i0,"/dp_",i0," [",a,"] = ",es24.16)') &
        cur_cost%id,iparam,trim(cur_param%name),mu_value
      call PrintMsg(option)
      write(name,'("dJ_",i0,"/dp_",i0)') cur_cost%id,iparam
      call this%KeepResult(name,mu_value)
      cur_param => cur_param%next
    enddo
    call VecRestoreArrayRead(this%mu(icost),mu_p,ierr);CHKERRQ(ierr)

    ! the initial state: the model's initialization map pulled back
    if (associated(coupled)) then
      call FlowRTDAEGetParts(coupled,this%lambda(icost),lam_f,lam_t)
      call RTDAEInitialSensitivity(rt,time0,Y0_t,lam_t,dJ_dc,dJ_dK)
      call AdjRTDAEPressureSensitivity(coupled,Y0_t,lam_f,lam_t,dJ_dc,dJ_dP)
      call FlowRTDAERestoreParts(coupled,this%lambda(icost),lam_f,lam_t)
    else
      call RTDAEInitialSensitivity(rt,time0,Y0,this%lambda(icost),dJ_dc, &
                                   dJ_dK)
    endif
    max_c = maxval(abs(dJ_dc))
    call MPI_Allreduce(MPI_IN_PLACE,max_c,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                       ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("  max |dJ_",i0,"/du_phys(0)| = ",es24.16)') &
      cur_cost%id,max_c
    call PrintMsg(option)
    write(name,'("max |dJ_",i0,"/du_phys(0)|")') cur_cost%id
    call this%KeepResult(name,max_c)
    if (nk > 0) then
      max_K = maxval(abs(dJ_dK(1:nk,:)))
      call MPI_Allreduce(MPI_IN_PLACE,max_K,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                         ierr);CHKERRQ(ierr)
      write(option%io_buffer,'("  max |dJ_",i0,"/dK(0)| = ",es24.16)') &
        cur_cost%id,max_K
      call PrintMsg(option)
      write(name,'("max |dJ_",i0,"/dK(0)|")') cur_cost%id
      call this%KeepResult(name,max_K)
    endif
    if (associated(coupled)) then
      max_c = maxval(abs(dJ_dP))
      call MPI_Allreduce(MPI_IN_PLACE,max_c,ONE_INTEGER_MPI, &
                         MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                         ierr);CHKERRQ(ierr)
      write(option%io_buffer,'("  max |dJ_",i0,"/dP(0)| = ",es24.16)') &
        cur_cost%id,max_c
      call PrintMsg(option)
      write(name,'("max |dJ_",i0,"/dP(0)|")') cur_cost%id
      call this%KeepResult(name,max_c)
    endif

    if (adjoint_output_full_field) then
      call VecCreateMPI(option%mycomm,grid%nlmax*nc,grid%nmax*nc,field, &
                        ierr);CHKERRQ(ierr)
      call VecSetBlockSize(field,nc,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax
        do i = 1, nc
          call VecSetValue(field,(grid%nG2A(grid%nL2G(local_id))-1)*nc+i-1, &
                           dJ_dc(i,local_id),INSERT_VALUES,ierr);CHKERRQ(ierr)
        enddo
      enddo
      call VecAssemblyBegin(field,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(field,ierr);CHKERRQ(ierr)
      write(filename,'(a,a,a,i0,a)') trim(option%global_prefix), &
        trim(option%group_prefix),'-adjoint-dJ',cur_cost%id,'-du0.bin'
      call PetscViewerBinaryOpen(option%mycomm,filename,FILE_MODE_WRITE, &
                                 viewer,ierr);CHKERRQ(ierr)
      call VecView(field,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
      call VecDestroy(field,ierr);CHKERRQ(ierr)
      call PrintMsg(option,'  Initial-condition sensitivity: ' // &
                    trim(filename))
      if (associated(coupled)) then
        call VecCreateMPI(option%mycomm,grid%nlmax,grid%nmax,field, &
                          ierr);CHKERRQ(ierr)
        do local_id = 1, grid%nlmax
          call VecSetValue(field,grid%nG2A(grid%nL2G(local_id))-1, &
                           dJ_dP(local_id),INSERT_VALUES,ierr);CHKERRQ(ierr)
        enddo
        call VecAssemblyBegin(field,ierr);CHKERRQ(ierr)
        call VecAssemblyEnd(field,ierr);CHKERRQ(ierr)
        write(filename,'(a,a,a,i0,a)') trim(option%global_prefix), &
          trim(option%group_prefix),'-adjoint-dJ',cur_cost%id,'-dP0.bin'
        call PetscViewerBinaryOpen(option%mycomm,filename,FILE_MODE_WRITE, &
                                   viewer,ierr);CHKERRQ(ierr)
        call VecView(field,viewer,ierr);CHKERRQ(ierr)
        call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        call VecDestroy(field,ierr);CHKERRQ(ierr)
        call PrintMsg(option,'  Initial-pressure sensitivity: ' // &
                      trim(filename))
      endif
    endif
    cur_cost => cur_cost%next
  enddo
  call PrintMsg(option,'=== END ADJOINT ===')
  if (associated(coupled)) call FlowRTDAERestoreParts(coupled,Y0,Y0_f,Y0_t)
  deallocate(dJ_dc,dJ_dK,dJ_dP)

end subroutine AdjRTDAEReport

! ************************************************************************** !

subroutine AdjRTDAEGetParameters(this,p)
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use Adjoint_Param_module, only : adj_param_list

  implicit none

  class(adjoint_rt_dae_type) :: this
  PetscReal :: p(:)

  type(adjoint_param_type), pointer :: cur_param
  PetscInt :: j

  j = 0
  cur_param => adj_param_list
  do while (associated(cur_param))
    j = j + 1
    p(j) = cur_param%value
    cur_param => cur_param%next
  enddo

end subroutine AdjRTDAEGetParameters

! ************************************************************************** !

subroutine AdjRTDAESetParameters(this,model,p)
  !
  ! Author: Hansol Suh
  ! Date: 10/07/26
  !
  use Adjoint_Param_module, only : adj_param_list

  implicit none

  class(adjoint_rt_dae_type) :: this
  class(ts_model_type) :: model
  PetscReal :: p(:)

  type(adjoint_param_type), pointer :: cur_param
  PetscInt :: j

  j = 0
  cur_param => adj_param_list
  do while (associated(cur_param))
    j = j + 1
    call cur_param%SetValue(this%realization,p(j))
    cur_param => cur_param%next
  enddo
  call model%InvalidateCache()

end subroutine AdjRTDAESetParameters

! ************************************************************************** !

subroutine AdjRTDAEPressureSensitivity(coupled,Y0_t,lam_f,lam_t,dJ_dc,dJ_dP)
  !
  ! dJ/dP0 through W0 = M(P0), N0 = A(c0,K0;P0) and, at the constraint's
  ! fixed molarity, c0 = molarity*1000/rho(P0) (a per-cell dataset
  ! constraint follows its own constraint types instead). With dJ/dc0 of
  ! RTDAEInitialSensitivity:
  !   dJ/dP0 = lambda_P + dM/dP lambda_W + (dA/dP)^T lambda_N
  !            - sum_aq dJ/dc0 c0 drho/dP / rho,
  ! at the last Prepare (the initial state). Zero in inactive cells.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Patch_module
  use Grid_module

  implicit none

  class(flow_rt_dae_model_type) :: coupled
  Vec :: Y0_t
  Vec :: lam_f
  Vec :: lam_t
  PetscReal :: dJ_dc(:,:)
  PetscReal :: dJ_dP(:)

  class(rt_dae_model_type), pointer :: rt
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  PetscReal, pointer :: lf_p(:), lt_p(:), y_p(:)
  PetscReal :: dA_dP(coupled%rt%ncomp), c0(coupled%rt%naqcomp)
  PetscInt :: local_id, ghosted_id, off_f, off_t, nc, naq
  PetscErrorCode :: ierr

  rt => coupled%rt
  patch => rt%realization%patch
  grid => patch%grid
  nc = rt%ncomp
  naq = rt%naqcomp
  dJ_dP = 0.d0
  call VecGetArrayRead(lam_f,lf_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(lam_t,lt_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(Y0_t,y_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.rt%cell_ok(ghosted_id)) cycle
    off_f = (local_id-1)*coupled%flow%ndof
    off_t = (local_id-1)*rt%ndof
    call RTDAEInventoryPressureDerivative(rt,ghosted_id,dA_dP)
    c0 = y_p(off_t+1:off_t+naq)
    if (rt%log_formulation) c0 = exp(c0)
    dJ_dP(local_id) = lf_p(off_f+1) + &
                      coupled%flow%dacc_dp(local_id)*lf_p(off_f+2) + &
                      dot_product(dA_dP,lt_p(off_t+nc+1:off_t+2*nc)) - &
                      dot_product(dJ_dc(1:naq,local_id),c0)* &
                      coupled%fs%dden_kg_dp(ghosted_id)/ &
                      coupled%fs%den_kg(ghosted_id)
  enddo
  call VecRestoreArrayRead(Y0_t,y_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(lam_t,lt_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(lam_f,lf_p,ierr);CHKERRQ(ierr)

end subroutine AdjRTDAEPressureSensitivity

! ************************************************************************** !

subroutine AdjRTDAEDestroy(this)
  !
  ! Also destroys the deck's cost and parameter lists.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  use Adjoint_module, only : AdjointDestroy

  implicit none

  class(adjoint_rt_dae_type) :: this

  nullify(this%realization)
  call AdjointDestroy()
  call TSAdjDestroy(this)

end subroutine AdjRTDAEDestroy

end module Adjoint_RT_DAE_class
