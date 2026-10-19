module RT_DAE_Chem_module
  !
  ! Local chemistry closure of the GIRT_DAE equation model: from the
  ! free-ion molalities m of the aqueous primary species, the activity
  ! coefficients, secondary molalities, aqueous totals T(m) and dT/dm.
  ! With activity coefficients on, the ionic strength I is solved to
  ! convergence at every evaluation:
  !
  !   g(I) = I - h(I;m) = 0,
  !   h(I;m) = 1/2 [ sum_i z_i^2 m_i + sum_s z_s^2 m_s(m,I) ],
  !   ln m_s = ln K_s + sum_i nu_si (ln m_i + ln gamma_i(I)) - ln gamma_s(I)
  !
  ! with primary activities updated at every iterate (RActivityCoefficients'
  ! NEWTON algorithm holds them fixed). The sensitivity
  !
  !   dI/dm_k = (dh/dm_k)|_I / g'(I)
  !
  ! enters every derivative as a rank-one correction. The iteration starts
  ! from I at gamma = 1, so the result is a deterministic function of m.
  ! Active gases are in equilibrium with the solution (RTDAEGasClosure).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module

  implicit none

  private

  PetscInt, parameter, public :: RT_DAE_CHEM_OK = 0
  PetscInt, parameter, public :: RT_DAE_CHEM_FAILED = 1

  ! relative convergence tolerance on g(I) and iteration cap
  PetscReal, parameter :: GAMMA_RTOL = 1.d-13
  PetscInt, parameter :: GAMMA_MAX_IT = 100
  ! ionic strength below which the solution is treated as ideal
  PetscReal, parameter :: I_IDEAL = 1.d-30

  type, public :: rt_dae_gamma_type
    ! false when activity coefficients are off or the solution is ideal;
    ! the sensitivities below are then zero
    PetscBool :: active
    PetscReal :: ionic_strength
    PetscInt :: iterations
    PetscReal, allocatable :: dlngam_pri_dI(:)  ! d ln gamma_i / dI [kg/mol]
    PetscReal, allocatable :: dI_dm(:)          ! d I / d m_k      [-]
  end type rt_dae_gamma_type

  public :: RTDAEGammaInit, &
            RTDAEGammaDestroy, &
            RTDAEAqueousClosure, &
            RTDAEGasClosure, &
            RTDAEGammaChain

contains

! ************************************************************************** !

subroutine RTDAEGammaInit(gam,reaction)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  type(rt_dae_gamma_type) :: gam
  class(reaction_rt_type) :: reaction

  allocate(gam%dlngam_pri_dI(reaction%naqcomp))
  allocate(gam%dI_dm(reaction%naqcomp))
  gam%dlngam_pri_dI = 0.d0
  gam%dI_dm = 0.d0
  gam%active = PETSC_FALSE
  gam%ionic_strength = 0.d0
  gam%iterations = 0

end subroutine RTDAEGammaInit

! ************************************************************************** !

subroutine RTDAEGammaDestroy(gam)
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !

  implicit none

  type(rt_dae_gamma_type) :: gam

  if (allocated(gam%dlngam_pri_dI)) deallocate(gam%dlngam_pri_dI)
  if (allocated(gam%dI_dm)) deallocate(gam%dI_dm)

end subroutine RTDAEGammaDestroy

! ************************************************************************** !

subroutine RTDAEDebyeHuckel(Z,a0,reaction,I,lngam,dlngam_dI)
  !
  ! ln gamma and its derivative for one species by RActivityCoefficients'
  ! B-dot form, with gamma = 1 for neutral species.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Reaction_module, only : RDebyeHuckelLnGamma, RDebyeHuckelDLog10GammaDI

  implicit none

  PetscReal :: Z, a0
  class(reaction_rt_type) :: reaction
  PetscReal :: I
  PetscReal :: lngam, dlngam_dI

  PetscReal :: sqrt_I

  if (abs(Z) <= 1.d-10) then
    lngam = 0.d0
    dlngam_dI = 0.d0
    return
  endif
  sqrt_I = sqrt(max(I,I_IDEAL))
  lngam = RDebyeHuckelLnGamma(reaction,Z,a0,I,sqrt_I)
  dlngam_dI = RDebyeHuckelDLog10GammaDI(reaction,Z,a0,sqrt_I)*LOG_TO_LN

end subroutine RTDAEDebyeHuckel

! ************************************************************************** !

subroutine RTDAEIonicStrengthResidual(reaction,ln_m_ideal,strength,g,dg_dI)
  !
  ! g(I) = I - h(I) and g'(I) for fixed m. ln_m_ideal(0) holds
  ! 1/2 sum_i z_i^2 m_i, ln_m_ideal(s) ln K_s + sum_i nu_si ln m_i.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  class(reaction_rt_type) :: reaction
  PetscReal :: ln_m_ideal(0:)
  PetscReal :: strength
  PetscReal :: g, dg_dI

  PetscInt :: icplx, i, icomp
  PetscReal :: lngam, dlngam, lngam_s, dlngam_s
  PetscReal :: ln_ms, ms, sum_ln, sum_dln, h, dh_dI
  PetscReal :: lngam_pri(reaction%naqcomp), dlngam_pri(reaction%naqcomp)

  do icomp = 1, reaction%naqcomp
    call RTDAEDebyeHuckel(reaction%primary_spec_Z(icomp), &
                          reaction%primary_spec_a0(icomp),reaction,strength, &
                          lngam,dlngam)
    lngam_pri(icomp) = lngam
    dlngam_pri(icomp) = dlngam
  enddo

  h = ln_m_ideal(0)
  dh_dI = 0.d0
  do icplx = 1, reaction%neqcplx
    if (abs(reaction%eqcplx_Z(icplx)) <= 1.d-10) cycle
    call RTDAEDebyeHuckel(reaction%eqcplx_Z(icplx), &
                          reaction%eqcplx_a0(icplx),reaction,strength, &
                          lngam_s,dlngam_s)
    sum_ln = 0.d0
    sum_dln = 0.d0
    do i = 1, reaction%eqcplxspecid(0,icplx)
      icomp = reaction%eqcplxspecid(i,icplx)
      sum_ln = sum_ln + reaction%eqcplxstoich(i,icplx)*lngam_pri(icomp)
      sum_dln = sum_dln + reaction%eqcplxstoich(i,icplx)*dlngam_pri(icomp)
    enddo
    ln_ms = ln_m_ideal(icplx) + sum_ln - lngam_s
    ms = exp(ln_ms)
    h = h + 0.5d0*reaction%eqcplx_Z(icplx)**2*ms
    dh_dI = dh_dI + 0.5d0*reaction%eqcplx_Z(icplx)**2*ms*(sum_dln-dlngam_s)
  enddo

  g = strength - h
  dg_dI = 1.d0 - dh_dI

end subroutine RTDAEIonicStrengthResidual

! ************************************************************************** !

subroutine RTDAEAqueousClosure(rt_auxvar,global_auxvar,reaction,option, &
                               gam,status)
  !
  ! On entry rt_auxvar%pri_molal holds m and global_auxvar the liquid
  ! density. On exit rt_auxvar holds the converged activity coefficients,
  ! secondary molalities, totals [mol/L] and in aqueous%dtotal the complete
  ! dT/dm [kg/L], activity coefficient sensitivity included; gam holds the
  ! sensitivities for activity-based rate laws (RTDAEGammaChain).
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  use Option_module
  use Reaction_module, only : RTotalAqueous

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  type(rt_dae_gamma_type) :: gam
  PetscInt :: status

  PetscInt, parameter :: iphase = 1
  PetscInt :: icomp, jcomp, icplx, i, j, it, naqcomp, neqcplx
  PetscReal :: ln_m(reaction%naqcomp)
  PetscReal :: ln_m_ideal(0:reaction%neqcplx)
  PetscReal :: strength, g, dg_dI, dg_dI_step, strength_new
  PetscReal :: lngam, dlngam
  PetscReal :: G_s(reaction%neqcplx)
  PetscReal :: q(reaction%naqcomp)
  PetscReal :: dh_dm(reaction%naqcomp)
  PetscReal :: den_kg_per_L, xmass, ms, sum_dln
  PetscReal :: dlngam_sec
  PetscBool :: polish

  status = RT_DAE_CHEM_OK
  naqcomp = reaction%naqcomp
  neqcplx = reaction%neqcplx
  gam%active = PETSC_FALSE
  gam%ionic_strength = 0.d0
  gam%iterations = 0
  gam%dlngam_pri_dI = 0.d0
  gam%dI_dm = 0.d0

  ! zero molality (underflow) passes here; the finite checks below reject
  ! 0/0 in a complex's derivative
  do icomp = 1, naqcomp
    if (.not.(rt_auxvar%pri_molal(icomp) >= 0.d0) .or. &
        rt_auxvar%pri_molal(icomp) > huge(1.d0)) then
      status = RT_DAE_CHEM_FAILED
      return
    endif
  enddo

  rt_auxvar%pri_act_coef = 1.d0
  if (neqcplx > 0) rt_auxvar%sec_act_coef = 1.d0
  ! water activity 1 (ACTIVITY_H2O is not supported)
  rt_auxvar%ln_act_h2o = 0.d0

  if (reaction%act_coef_update_frequency /= ACT_COEF_FREQUENCY_OFF) then
    ln_m = log(rt_auxvar%pri_molal)
    ! primary part of h and the ideal log molality of each complex
    ln_m_ideal(0) = 0.d0
    do icomp = 1, naqcomp
      ln_m_ideal(0) = ln_m_ideal(0) + &
        0.5d0*reaction%primary_spec_Z(icomp)**2*rt_auxvar%pri_molal(icomp)
    enddo
    do icplx = 1, neqcplx
      ln_m_ideal(icplx) = -reaction%eqcplx_logK(icplx)*LOG_TO_LN
      do i = 1, reaction%eqcplxspecid(0,icplx)
        icomp = reaction%eqcplxspecid(i,icplx)
        ln_m_ideal(icplx) = ln_m_ideal(icplx) + &
          reaction%eqcplxstoich(i,icplx)*ln_m(icomp)
      enddo
    enddo

    ! deterministic start: strength at gamma = 1
    strength = ln_m_ideal(0)
    do icplx = 1, neqcplx
      strength = strength + &
        0.5d0*reaction%eqcplx_Z(icplx)**2*exp(ln_m_ideal(icplx))
    enddo

    if (strength > I_IDEAL) then
      it = 0
      polish = PETSC_FALSE
      do
        it = it + 1
        call RTDAEIonicStrengthResidual(reaction,ln_m_ideal,strength,g,dg_dI)
        if (.not.(abs(g) <= huge(1.d0))) then
          status = RT_DAE_CHEM_FAILED
          return
        endif
        ! once converged, one more Newton step reaches roundoff
        if (polish) exit
        if (abs(g) <= GAMMA_RTOL*strength) polish = PETSC_TRUE
        if (it > GAMMA_MAX_IT) then
          status = RT_DAE_CHEM_FAILED
          return
        endif
        ! g' > 0 near a well-posed root; the floor bounds the step where g'
        ! flattens
        dg_dI_step = max(dg_dI,0.1d0)
        strength_new = strength - g/dg_dI_step
        if (strength_new <= 0.d0) strength_new = 0.1d0*strength
        strength = strength_new
      enddo
      gam%iterations = it
      ! the sensitivity needs the true g'
      if (.not.(dg_dI > 1.d-8)) then
        status = RT_DAE_CHEM_FAILED
        return
      endif

      gam%active = PETSC_TRUE
      gam%ionic_strength = strength
      do icomp = 1, naqcomp
        call RTDAEDebyeHuckel(reaction%primary_spec_Z(icomp), &
                              reaction%primary_spec_a0(icomp),reaction, &
                              strength,lngam,dlngam)
        rt_auxvar%pri_act_coef(icomp) = exp(lngam)
        gam%dlngam_pri_dI(icomp) = dlngam
      enddo
      do icplx = 1, neqcplx
        call RTDAEDebyeHuckel(reaction%eqcplx_Z(icplx), &
                              reaction%eqcplx_a0(icplx),reaction,strength, &
                              lngam,dlngam_sec)
        rt_auxvar%sec_act_coef(icplx) = exp(lngam)
        sum_dln = 0.d0
        do i = 1, reaction%eqcplxspecid(0,icplx)
          icomp = reaction%eqcplxspecid(i,icplx)
          sum_dln = sum_dln + &
            reaction%eqcplxstoich(i,icplx)*gam%dlngam_pri_dI(icomp)
        enddo
        G_s(icplx) = sum_dln - dlngam_sec
      enddo
    endif
  endif

  ! totals and the activity-coefficient-fixed derivative
  call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
  if (.not.RTDAEAllFinite(rt_auxvar)) then
    status = RT_DAE_CHEM_FAILED
    return
  endif

  if (.not.gam%active) return

  ! rank-one correction of dT/dm for the ionic strength sensitivity
  xmass = 1.d0
  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  den_kg_per_L = global_auxvar%den_kg(iphase)*xmass*1.d-3
  do icomp = 1, naqcomp
    dh_dm(icomp) = 0.5d0*reaction%primary_spec_Z(icomp)**2
  enddo
  q = 0.d0
  do icplx = 1, neqcplx
    ms = rt_auxvar%sec_molal(icplx)
    do i = 1, reaction%eqcplxspecid(0,icplx)
      icomp = reaction%eqcplxspecid(i,icplx)
      q(icomp) = q(icomp) + reaction%eqcplxstoich(i,icplx)*ms*G_s(icplx)
      dh_dm(icomp) = dh_dm(icomp) + 0.5d0*reaction%eqcplx_Z(icplx)**2* &
        ms*reaction%eqcplxstoich(i,icplx)/rt_auxvar%pri_molal(icomp)
    enddo
  enddo
  q = q*den_kg_per_L
  ! g'(I) at the converged I
  call RTDAEIonicStrengthResidual(reaction,ln_m_ideal,gam%ionic_strength, &
                                  g,dg_dI)
  gam%dI_dm = dh_dm/dg_dI
  do jcomp = 1, naqcomp
    do j = 1, naqcomp
      rt_auxvar%aqueous%dtotal(j,jcomp,iphase) = &
        rt_auxvar%aqueous%dtotal(j,jcomp,iphase) + q(j)*gam%dI_dm(jcomp)
    enddo
  enddo
  if (.not.RTDAEAllFinite(rt_auxvar)) status = RT_DAE_CHEM_FAILED

end subroutine RTDAEAqueousClosure

! ************************************************************************** !

subroutine RTDAEGasClosure(rt_auxvar,global_auxvar,reaction,option,gam, &
                           status)
  !
  ! After RTDAEAqueousClosure: the active gases in equilibrium with the
  ! solution (ReactionGasTotalGas), their partial pressures [bar] in
  ! rt_auxvar%gas_pp, the gas-phase totals total(:,2) [mol/L gas] and in
  ! aqueous%dtotal(:,:,2) the complete dT/dm [kg water/L gas]. A partial
  ! pressure depends on the solution only through the primary activities,
  ! so the activity-coefficient sensitivity is RTDAEGammaChain's.
  !
  ! Author: Hansol Suh
  ! Date: 10/19/26
  !
  use Option_module
  use Reaction_Gas_module, only : ReactionGasTotalGas

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  type(rt_dae_gamma_type) :: gam
  PetscInt :: status

  PetscInt, parameter :: iphase = 2

  status = RT_DAE_CHEM_OK
  call ReactionGasTotalGas(rt_auxvar,global_auxvar,reaction,option)
  call RTDAEGammaChain(rt_auxvar%aqueous%dtotal(:,:,iphase), &
                       reaction%naqcomp,reaction%naqcomp, &
                       rt_auxvar%pri_molal,gam)
  ! zero molality in a gas's reaction gives 0/0 in its derivative
  if (.not.(all(abs(rt_auxvar%gas_pp) <= huge(1.d0)) .and. &
            all(abs(rt_auxvar%total(:,iphase)) <= huge(1.d0)) .and. &
            all(abs(rt_auxvar%aqueous%dtotal(:,:,iphase)) <= huge(1.d0)))) &
    status = RT_DAE_CHEM_FAILED

end subroutine RTDAEGasClosure

! ************************************************************************** !

function RTDAEAllFinite(rt_auxvar)
  !
  ! Whether the totals and their derivatives are finite.
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  PetscBool :: RTDAEAllFinite

  RTDAEAllFinite = all(abs(rt_auxvar%total(:,1)) <= huge(1.d0)) .and. &
                   all(abs(rt_auxvar%aqueous%dtotal(:,:,1)) <= huge(1.d0))

end function RTDAEAllFinite

! ************************************************************************** !

subroutine RTDAEGammaChain(J,nrow,naqcomp,pri_molal,gam)
  !
  ! Adds the activity-coefficient sensitivity to J(:,1:naqcomp), the
  ! fixed-gamma derivative of a quantity Q that depends on the chemistry
  ! only through a_i = gamma_i m_i. Then dQ/d ln gamma_i = m_i dQ/dm_i and
  ! the correction is rank one:
  !
  !   J(:,k) += [ sum_i J(:,i) m_i d ln gamma_i/dI ] dI/dm_k
  !
  ! Author: Hansol Suh
  ! Date: 09/22/26
  !
  implicit none

  PetscInt :: nrow, naqcomp
  PetscReal :: J(:,:)
  PetscReal :: pri_molal(naqcomp)
  type(rt_dae_gamma_type) :: gam

  PetscReal :: v(nrow)
  PetscInt :: i, k

  if (.not.gam%active) return
  v = 0.d0
  do i = 1, naqcomp
    v(:) = v(:) + J(1:nrow,i)*pri_molal(i)*gam%dlngam_pri_dI(i)
  enddo
  do k = 1, naqcomp
    J(1:nrow,k) = J(1:nrow,k) + v(:)*gam%dI_dm(k)
  enddo

end subroutine RTDAEGammaChain

end module RT_DAE_Chem_module
