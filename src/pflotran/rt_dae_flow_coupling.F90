module RT_DAE_Flow_Coupling_module
  !
  ! Derivatives of the GIRT_DAE model (RT_DAE_Model_class) with respect to
  ! the flow state, for MODE RICHARDS_DAE with MODE GIRT_DAE as one model
  ! (Flow_RT_DAE_Model_class) and its adjoint: dF/dP of the transport rows
  ! (RTDAEAddFlowCoupling) and dF/dv for given changes of the Darcy fluxes
  ! (RTDAEAddVelocityDerivative).
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
#include "petsc/finclude/petscts.h"
  use petscts

  use TS_Model_class
  use PFLOTRAN_Constants_module
  use Option_module
  use Reactive_Transport_Aux_module
  use RT_DAE_Model_class
  use TS_Flow_State_module

  implicit none

  private

  public :: RTDAESetFlowCoupling, &
            RTDAEKeepActivities, &
            RTDAEAddFlowCoupling, &
            RTDAEInventoryPressureDerivative, &
            RTDAEAddVelocityDerivative

contains

! ************************************************************************** !

subroutine RTDAESetFlowCoupling(this,fs)
  !
  ! Couples the model to the flow state fs (see RTDAEForcingAt) and
  ! rejects transport physics whose flow derivatives RTDAEAddFlowCoupling
  ! lacks.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Reaction_Surface_Complexation_Aux_module

  implicit none

  class(rt_dae_model_type) :: this
  type(ts_flow_state_type), target :: fs

  type(reactive_transport_param_type), pointer :: rt_parameter
  type(surface_complexation_type), pointer :: srfcplx
  character(len=MAXSTRINGLENGTH) :: unsupported
  PetscBool :: rock_multirate
  PetscInt :: i

  rt_parameter => this%realization%patch%aux%RT%rt_parameter
  unsupported = ''
  if (rt_parameter%calculate_transverse_dispersion) &
    call TSModelListAdd(unsupported,'transverse dispersion')
  if (this%reaction%nradiodecay_rxn > 0) &
    call TSModelListAdd(unsupported,'RADIOACTIVE_DECAY_REACTION')
  ! with soil compressibility, sorption on rock sites follows the solid
  ! fraction 1 - phi(P): the equilibrium inventories carry it
  ! (RTDAEInventoryPressureDerivative), the multirate rates do not
  if (this%option%flow%transient_porosity) then
    srfcplx => this%reaction%surface_complexation
    rock_multirate = PETSC_FALSE
    do i = 1, srfcplx%nkinmrsrfcplxrxn
      if (srfcplx%srfcplxrxn_surf_type( &
            srfcplx%kinmrsrfcplxrxn_to_srfcplxrxn(i)) == ROCK_SURFACE) &
        rock_multirate = PETSC_TRUE
    enddo
    if (rock_multirate) call TSModelListAdd(unsupported,'multirate &
      &sorption on ROCK_DENSITY sites with soil compressibility')
    if (this%eq_srfcplx) then
      allocate(this%sorb_rock(this%naqcomp,this%ngmax))
      this%sorb_rock = 0.d0
    endif
  endif
  if (len_trim(unsupported) > 0) then
    this%option%io_buffer = 'MODE GIRT_DAE coupled to MODE RICHARDS_DAE &
      &does not support ' // trim(unsupported) // '.'
    call PrintErrMsg(this%option)
  endif
  this%fs => fs
  this%forcing%forcing_segment = -1
  this%prepared = PETSC_FALSE
  ! for these rates' pressure derivatives (RTDAEAddFlowCoupling)
  if (this%reaction%microbial%nrxn > 0 .or. &
      this%reaction%ngeneral_rxn > 0) call RTDAEKeepActivities(this)

end subroutine RTDAESetFlowCoupling

! ************************************************************************** !

subroutine RTDAEKeepActivities(this)
  !
  ! From now on keeps the activity coefficients of each cell's last
  ! closure, so RTDAELoadCell and RTDAELoadCellRates need not solve it
  ! again. The next Prepare fills them.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this

  if (allocated(this%pri_act_coef)) return
  allocate(this%pri_act_coef(this%naqcomp,this%ngmax))
  this%pri_act_coef = 1.d0
  this%prepared = PETSC_FALSE

end subroutine RTDAEKeepActivities

! ************************************************************************** !

subroutine RTDAEAddFlowCoupling(this,J,flow_ndof)
  !
  ! Adds dF/dP, the derivatives with respect to the coupled flow model's
  ! pressures at the last Prepare with derivatives, to J in local
  ! numbering: rows ndof per ghosted cell, columns the flow model's
  ! flow_ndof per ghosted cell (pressure first).
  !
  ! Pressure enters through fs (cell saturation, porosity and water
  ! density, Darcy fluxes, boundary state, source volume rates):
  !   F_c: A = phi*s*V*1000*T, T proportional to the water density at
  !        fixed free-ion molalities
  !   F_N: the transport coefficients (TFluxCoef, TDispersion,
  !        TDispersionBC, TSrcSinkCoef) times the totals, and rates per
  !        unit pore water (microbial rates without immobile biomass or in
  !        molarity, general rates)
  ! Boundary and source constraint totals are fixed in a segment.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module
  use Transport_module
  use Reaction_module, only : RGeneral
  use Reaction_Microbial_module, only : ReactionMicrobRate
  use Reaction_Microbial_Aux_module, only : MICROBIAL_MOLARITY

  implicit none

  class(rt_dae_model_type) :: this
  Mat :: J
  PetscInt :: flow_ndof

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(ts_flow_state_type), pointer :: fs
  PetscReal :: col(this%ndof), res(this%ncomp), jac(this%ncomp,this%ncomp)
  PetscReal :: dH_dsu(this%realization%patch%aux%RT%rt_parameter%ndiffcoef)
  PetscReal :: dH_dpu(size(dH_dsu)), dH_dsd(size(dH_dsu))
  PetscReal :: dH_dpd(size(dH_dsu)), dH_dq(size(dH_dsu))
  PetscReal :: dHu(this%naqcomp), dHd(this%naqcomp)
  PetscReal :: diff(this%naqcomp), adv(this%naqcomp)
  PetscReal :: dflux_u(this%naqcomp), dflux_d(this%naqcomp)
  PetscReal, pointer :: yl_p(:)
  PetscReal :: dlnw, dlnrho, a1000, q, dq
  PetscInt :: local_id, ghosted_id, iconn, up, dn, irxn, nc, naq
  PetscBool :: molarity
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  fs => this%fs
  nc = this%ncomp
  naq = this%naqcomp

  ! cells: inventory and reactions per unit pore water
  call VecGetArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    col = 0.d0
    dlnrho = fs%dden_kg_dp(ghosted_id)/fs%den_kg(ghosted_id)
    ! d ln(phi s) / dP
    dlnw = fs%dpor_dp(ghosted_id)/fs%por(ghosted_id) + &
           fs%dsat_dp(ghosted_id)/fs%sat(ghosted_id)
    call RTDAEInventoryPressureDerivative(this,ghosted_id,col(1:nc))
    col(1:nc) = -col(1:nc)
    if (allocated(this%pri_act_coef)) then
      call RTDAELoadCellRates(this,ghosted_id,yl_p)
      ! molarity: the concentrations are proportional to the density
      molarity = this%reaction%microbial%concentration_units == &
                 MICROBIAL_MOLARITY
      do irxn = 1, this%reaction%microbial%nrxn
        res = 0.d0
        jac = 0.d0
        call ReactionMicrobRate(res,jac,molarity,this%scratch, &
                                this%forcing%gaux(ghosted_id), &
                                this%maux(ghosted_id), &
                                this%reaction,this%option,irxn_filter=irxn)
        ! aqueous or no biomass: the rate is per unit pore water
        if (this%reaction%microbial%biomassid(irxn) >= 0) &
          col(nc+1:2*nc) = col(nc+1:2*nc) + res*dlnw
        if (molarity) &
          col(nc+1:2*nc) = col(nc+1:2*nc) + &
            matmul(jac(:,1:naq),this%scratch%pri_molal)*dlnrho
      enddo
      if (this%reaction%ngeneral_rxn > 0) then
        res = 0.d0
        call RGeneral(res,jac,PETSC_FALSE,this%scratch, &
                      this%forcing%gaux(ghosted_id),this%maux(ghosted_id), &
                      this%reaction,this%option)
        ! rates per kg of pore water
        col(nc+1:2*nc) = col(nc+1:2*nc) + res*(dlnw+dlnrho)
      endif
    endif
    call Add(ghosted_id,ghosted_id,1,2*nc,col(1:2*nc))
  enddo
  call VecRestoreArrayRead(this%Y_loc,yl_p,ierr);CHKERRQ(ierr)

  ! internal faces: flux from up to down, into rows N of up (+) and down (-)
  do iconn = 1, this%conn%nconn_int
    up = this%conn%int_up(iconn)
    dn = this%conn%int_dn(iconn)
    if (patch%imat(up) <= 0 .or. patch%imat(dn) <= 0) cycle
    if (grid%nG2L(up) <= 0 .and. grid%nG2L(dn) <= 0) cycle
    if (.not.(this%cell_ok(up) .and. this%cell_ok(dn))) cycle
    ! TFluxCoef: no flux out of or into a dry cell
    if (fs%sat(up) < rt_min_saturation .or. &
        fs%sat(dn) < rt_min_saturation) cycle
    q = this%forcing%vel_int(1,iconn)
    call TDispersionDerivative(this%forcing%gaux(up),this%maux(up), &
                               patch%material_property_array( &
                                 patch%imat(up))%ptr%dispersivity, &
                               this%forcing%gaux(dn),this%maux(dn), &
                               patch%material_property_array( &
                                 patch%imat(dn))%ptr%dispersivity, &
                               this%conn%int_dist(:,iconn),rt_parameter, &
                               this%option,this%forcing%vel_int(:,iconn), &
                               dH_dsu,dH_dpu,dH_dsd,dH_dpd,dH_dq)
    ! pressure derivatives of the harmonic coefficient, per component (one
    ! for all, or one each)
    if (size(dH_dsu) > 1) then
      dHu = dH_dsu*fs%dsat_dp(up) + dH_dpu*fs%dpor_dp(up) + &
            dH_dq*fs%dv_int(1,iconn)
      dHd = dH_dsd*fs%dsat_dp(dn) + dH_dpd*fs%dpor_dp(dn) + &
            dH_dq*fs%dv_int(2,iconn)
    else
      dHu = sum(dH_dsu*fs%dsat_dp(up) + dH_dpu*fs%dpor_dp(up) + &
                dH_dq*fs%dv_int(1,iconn))
      dHd = sum(dH_dsd*fs%dsat_dp(dn) + dH_dpd*fs%dpor_dp(dn) + &
                dH_dq*fs%dv_int(2,iconn))
    endif
    diff = this%tot(:,up) - this%tot(:,dn)
    if (q > 0.d0) then
      adv = this%tot(:,up)
    else
      adv = this%tot(:,dn)
    endif
    a1000 = this%conn%int_area(iconn)*1000.d0
    dflux_u = a1000*(dHu*diff + fs%dv_int(1,iconn)*adv) + &
              this%coef_up(:,iconn)*this%tot(:,up)* &
              fs%dden_kg_dp(up)/fs%den_kg(up)
    dflux_d = a1000*(dHd*diff + fs%dv_int(2,iconn)*adv) + &
              this%coef_dn(:,iconn)*this%tot(:,dn)* &
              fs%dden_kg_dp(dn)/fs%den_kg(dn)
    if (grid%nG2L(up) > 0) then
      call Add(up,up,nc+1,nc+naq,dflux_u)
      call Add(up,dn,nc+1,nc+naq,dflux_d)
    endif
    if (grid%nG2L(dn) > 0) then
      call Add(dn,up,nc+1,nc+naq,-dflux_u)
      call Add(dn,dn,nc+1,nc+naq,-dflux_d)
    endif
  enddo

  ! boundary faces: flux into the cell, out of rows N
  do iconn = 1, this%conn%nconn_bc
    ghosted_id = this%conn%bc_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    if (this%bc_itype(iconn) == MEMBRANE_BC) cycle
    if (fs%sat(ghosted_id) < rt_min_saturation) cycle
    q = this%forcing%vel_bc(1,iconn)
    call TDispersionBCDerivative(this%bc_itype(iconn), &
                                 this%forcing%gaux_bc(iconn), &
                                 this%forcing%gaux(ghosted_id), &
                                 this%maux(ghosted_id), &
                                 patch%material_property_array( &
                                   patch%imat(ghosted_id))%ptr%dispersivity, &
                                 this%conn%bc_dist(:,iconn),rt_parameter, &
                                 this%forcing%vel_bc(:,iconn), &
                                 dH_dsu,dH_dsd,dH_dpd, &
                                 dH_dq)
    ! (dH_dsu: the boundary's saturation)
    if (size(dH_dsu) > 1) then
      dHd = dH_dsu*fs%dsat_bc_dp(iconn) + dH_dsd*fs%dsat_dp(ghosted_id) + &
            dH_dpd*fs%dpor_dp(ghosted_id) + dH_dq*fs%dv_bc(iconn)
    else
      dHd = sum(dH_dsu*fs%dsat_bc_dp(iconn) + &
                dH_dsd*fs%dsat_dp(ghosted_id) + &
                dH_dpd*fs%dpor_dp(ghosted_id) + dH_dq*fs%dv_bc(iconn))
    endif
    diff = this%forcing%seg%total_bc(:,iconn) - this%tot(:,ghosted_id)
    if (q > 0.d0) then
      adv = this%forcing%seg%total_bc(:,iconn)
    else
      adv = this%tot(:,ghosted_id)
    endif
    a1000 = this%conn%bc_area(iconn)*1000.d0
    dflux_d = a1000*(dHd*diff + fs%dv_bc(iconn)*adv) + &
              this%coef_bc_dn(:,iconn)*this%tot(:,ghosted_id)* &
              fs%dden_kg_dp(ghosted_id)/fs%den_kg(ghosted_id)
    call Add(ghosted_id,ghosted_id,nc+1,nc+naq,-dflux_d)
  enddo

  ! sources: flux out of the cell (injection negative), into rows N
  do iconn = 1, this%conn%nconn_ss
    ghosted_id = this%conn%ss_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    if (fs%sat(ghosted_id) < rt_min_saturation) cycle
    q = fs%qvol_ss(iconn)
    dq = fs%dqvol_ss(iconn)
    dlnrho = fs%dden_kg_dp(ghosted_id)/fs%den_kg(ghosted_id)
    select case(this%ss_itype(iconn))
      case(EQUILIBRIUM_SS)
        dflux_d = this%coef_ss_in(iconn)*this%tot(:,ghosted_id)*dlnrho
      case(MASS_RATE_SS,MEMBRANE_BC)
        cycle
      case default
        if (q > 0.d0) then
          dflux_d = -1000.d0*dq*this%forcing%seg%total_ss(:,iconn)
        else
          dflux_d = -1000.d0*dq*this%tot(:,ghosted_id) + &
                    this%coef_ss_in(iconn)*this%tot(:,ghosted_id)*dlnrho
        endif
    end select
    call Add(ghosted_id,ghosted_id,nc+1,nc+naq,dflux_d)
  enddo

contains

  subroutine Add(ghosted_row,ghosted_col,first,last,values)
    ! rows first:last of cell ghosted_row, the pressure of ghosted_col
    PetscInt :: ghosted_row, ghosted_col, first, last
    PetscReal :: values(:)
    PetscInt :: rows(last-first+1), i
    do i = first, last
      rows(i-first+1) = (ghosted_row-1)*this%ndof + i - 1
    enddo
    call MatSetValuesLocal(J,last-first+1,rows,1, &
                           [(ghosted_col-1)*flow_ndof],values, &
                           ADD_VALUES,ierr);CHKERRQ(ierr)
  end subroutine Add

end subroutine RTDAEAddFlowCoupling

! ************************************************************************** !

subroutine RTDAEInventoryPressureDerivative(this,ghosted_id,dA_dP)
  !
  ! Coupled: dA/dP of an active cell's equilibrium inventory
  ! (RTDAEAccumulation) with respect to its pressure, at the last Prepare:
  ! the aqueous part phi*s*V*1000*T (T proportional to the water density)
  ! and, with soil compressibility, the sorbed part on rock sites (site
  ! density proportional to 1 - phi).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(rt_dae_model_type) :: this
  PetscInt :: ghosted_id
  PetscReal :: dA_dP(:)

  type(ts_flow_state_type), pointer :: fs

  fs => this%fs
  dA_dP = 0.d0
  dA_dP(1:this%naqcomp) = &
    RTDAEPoreLiters(fs%por(ghosted_id)*fs%sat(ghosted_id), &
                    this%maux(ghosted_id)%volume)*this%tot(:,ghosted_id)* &
    (fs%dpor_dp(ghosted_id)/fs%por(ghosted_id) + &
     fs%dsat_dp(ghosted_id)/fs%sat(ghosted_id) + &
     fs%dden_kg_dp(ghosted_id)/fs%den_kg(ghosted_id))
  if (allocated(this%sorb_rock)) &
    dA_dP(1:this%naqcomp) = dA_dP(1:this%naqcomp) - &
      this%maux(ghosted_id)%volume*this%sorb_rock(:,ghosted_id)* &
      fs%dpor_dp(ghosted_id)/(1.d0-fs%por(ghosted_id))

end subroutine RTDAEInventoryPressureDerivative

! ************************************************************************** !

subroutine RTDAEAddVelocityDerivative(this,dv_int,dv_bc,F)
  !
  ! Coupled: adds dF/dv * dv to F (this model's layout) for Darcy flux
  ! changes dv_int (internal, up to down) and dv_bc (boundary, into the
  ! cell) at fixed cell state and the last Prepare: advection and the
  ! velocity dependence of mechanical dispersion.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Patch_module
  use Grid_module
  use Transport_module

  implicit none

  class(rt_dae_model_type) :: this
  PetscReal :: dv_int(:)
  PetscReal :: dv_bc(:)
  Vec :: F

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(ts_flow_state_type), pointer :: fs
  PetscReal, pointer :: f_p(:)
  PetscReal :: dH_dsu(this%realization%patch%aux%RT%rt_parameter%ndiffcoef)
  PetscReal :: dH_dpu(size(dH_dsu)), dH_dsd(size(dH_dsu))
  PetscReal :: dH_dpd(size(dH_dsu)), dH_dq(size(dH_dsu))
  PetscReal :: dH(this%naqcomp), adv(this%naqcomp), dflux(this%naqcomp)
  PetscReal :: q
  PetscInt :: iconn, up, dn, lup, ldn, ghosted_id, off, nc, naq
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid
  rt_parameter => patch%aux%RT%rt_parameter
  fs => this%fs
  nc = this%ncomp
  naq = this%naqcomp
  call VecGetArray(F,f_p,ierr);CHKERRQ(ierr)

  do iconn = 1, this%conn%nconn_int
    if (dv_int(iconn) == 0.d0) cycle
    up = this%conn%int_up(iconn)
    dn = this%conn%int_dn(iconn)
    if (patch%imat(up) <= 0 .or. patch%imat(dn) <= 0) cycle
    lup = grid%nG2L(up)
    ldn = grid%nG2L(dn)
    if (lup <= 0 .and. ldn <= 0) cycle
    if (.not.(this%cell_ok(up) .and. this%cell_ok(dn))) cycle
    if (fs%sat(up) < rt_min_saturation .or. &
        fs%sat(dn) < rt_min_saturation) cycle
    q = this%forcing%vel_int(1,iconn)
    call TDispersionDerivative(this%forcing%gaux(up),this%maux(up), &
                               patch%material_property_array( &
                                 patch%imat(up))%ptr%dispersivity, &
                               this%forcing%gaux(dn),this%maux(dn), &
                               patch%material_property_array( &
                                 patch%imat(dn))%ptr%dispersivity, &
                               this%conn%int_dist(:,iconn),rt_parameter, &
                               this%option,this%forcing%vel_int(:,iconn), &
                               dH_dsu,dH_dpu,dH_dsd,dH_dpd,dH_dq)
    if (size(dH_dq) > 1) then
      dH = dH_dq
    else
      dH = dH_dq(1)
    endif
    if (q > 0.d0) then
      adv = this%tot(:,up)
    else
      adv = this%tot(:,dn)
    endif
    dflux = this%conn%int_area(iconn)*1000.d0*dv_int(iconn)* &
            (dH*(this%tot(:,up)-this%tot(:,dn)) + adv)
    if (lup > 0) then
      off = (lup-1)*this%ndof + nc
      f_p(off+1:off+naq) = f_p(off+1:off+naq) + dflux
    endif
    if (ldn > 0) then
      off = (ldn-1)*this%ndof + nc
      f_p(off+1:off+naq) = f_p(off+1:off+naq) - dflux
    endif
  enddo

  do iconn = 1, this%conn%nconn_bc
    if (dv_bc(iconn) == 0.d0) cycle
    ghosted_id = this%conn%bc_ghost(iconn)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (.not.this%cell_ok(ghosted_id)) cycle
    if (this%bc_itype(iconn) == MEMBRANE_BC) cycle
    if (fs%sat(ghosted_id) < rt_min_saturation) cycle
    q = this%forcing%vel_bc(1,iconn)
    call TDispersionBCDerivative(this%bc_itype(iconn), &
                                 this%forcing%gaux_bc(iconn), &
                                 this%forcing%gaux(ghosted_id), &
                                 this%maux(ghosted_id), &
                                 patch%material_property_array( &
                                   patch%imat(ghosted_id))%ptr%dispersivity, &
                                 this%conn%bc_dist(:,iconn),rt_parameter, &
                                 this%forcing%vel_bc(:,iconn), &
                                 dH_dsu,dH_dsd,dH_dpd, &
                                 dH_dq)
    if (size(dH_dq) > 1) then
      dH = dH_dq
    else
      dH = dH_dq(1)
    endif
    if (q > 0.d0) then
      adv = this%forcing%seg%total_bc(:,iconn)
    else
      adv = this%tot(:,ghosted_id)
    endif
    dflux = this%conn%bc_area(iconn)*1000.d0*dv_bc(iconn)* &
            (dH*(this%forcing%seg%total_bc(:,iconn)-this%tot(:,ghosted_id)) + &
             adv)
    off = (this%conn%bc_local(iconn)-1)*this%ndof + nc
    f_p(off+1:off+naq) = f_p(off+1:off+naq) - dflux
  enddo

  call VecRestoreArray(F,f_p,ierr);CHKERRQ(ierr)

end subroutine RTDAEAddVelocityDerivative

end module RT_DAE_Flow_Coupling_module
