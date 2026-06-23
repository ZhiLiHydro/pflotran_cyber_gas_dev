module THCFlow_Common_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Global_Aux_module
  use THCFlow_Aux_module

  implicit none

  private

  ! Liquid specific heat capacity c_p,l [J/(kg.K)] -- treated as constant over
  ! the 0-100 C range
  PetscReal, parameter :: thcflow_heat_capacity_liquid = 4184.d0

  ! Cutoff parameters
  PetscReal, parameter :: eps     = 1.d-8
  PetscReal, parameter :: floweps = 0.d0

  ! The branch below is implemented now so the energy BC Jacobian is complete.
  PetscInt, parameter :: THCFLOW_CONVECTIVE_BC = -2

  ! --- public interface ---------------------------------------------------
  public :: THCFlowAccumulation, &
            THCFlowFlux, &
            THCFlowBCFlux, &
            THCFlowSrcSink, &
            THCFlowAccumDerivative, &
            THCFlowFluxDerivative, &
            THCFlowBCFluxDerivative, &
            THCFlowSrcSinkDerivative

contains

! ************************************************************************** !

subroutine THCFlowAccumulation(thcflow_auxvar,global_auxvar,material_auxvar, &
                               option,Res,Jac,dResdparam,calculate_derivatives)
  !
  ! Computes the non-fixed portion of the accumulation term for the THCFLOW
  ! residual (flow, energy, solute) and its 3x3 accumulation Jacobian block.
  !
  ! Argument list mirrors ZFlowAccumulation. dResdparam is retained for
  ! interface compatibility but is a NO-OP here (zeroed, never filled) --
  ! THCFLOW does not yet support the adjoint-parameter sensitivities.
  !
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(thcflow_auxvar_type) :: thcflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jac(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparam(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscBool :: calculate_derivatives

  PetscReal :: porosity
  PetscReal :: saturation
  PetscReal :: por_sat
  PetscReal :: dpor_sat_dp
  PetscReal :: volume_over_dt
  PetscReal :: temperature
  PetscReal :: den_kg
  PetscReal :: den_kmol
  PetscReal :: c_p
  PetscReal :: rho_s_c_s
  PetscReal :: tempreal
  PetscReal, parameter :: L_per_m3 = 1.d3

  Res = 0.d0
  Jac = 0.d0
  dResdparam = 0.d0

  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use thcflow_auxvar%effective_porosity here as it enables numerical
  ! derivatives to be employed
  porosity    = thcflow_auxvar%effective_porosity
  saturation  = thcflow_auxvar%sat
  por_sat     = saturation * porosity
  ! d(por_sat)/dP = dS/dP * por + S * dpor/dP
  dpor_sat_dp = thcflow_auxvar%dsat_dp * porosity + &
                saturation * thcflow_auxvar%dpor_dp
  temperature = thcflow_auxvar%temp
  den_kg      = thcflow_auxvar%den_kg
  den_kmol    = thcflow_auxvar%den_kmol
  c_p         = thcflow_heat_capacity_liquid
  rho_s_c_s   = thcflow_auxvar%heat_cap_solid     ! rho_s * c_s [J/(m^3.K)]

  ! --- Flow (liquid mass) equation : units kmol/s -------------------------
  ! Res[kmol/sec] = sat[m^3 liq/m^3 void] * por[m^3 void/m^3 bulk] *
  !                 rho_kmol[kmol/m^3] * vol[m^3 bulk] / dt[sec]
  tempreal = volume_over_dt * den_kmol
  Res(thcflow_pressure_dof) = por_sat * tempreal

  ! --- Solute equation : units mol/s --------------------------------------
  ! Res[mole/sec] = c[mol/L] * 1000[L/m^3] * sat * por * vol / dt
  tempreal = L_per_m3 * volume_over_dt
  Res(thcflow_concentration_dof) = thcflow_auxvar%conc * tempreal * por_sat

  ! --- Energy equation : units J/s = W ------------------------------------
  ! Res[W] = [phi*S_l*rho_kg*c_p*T + (1-phi)*rho_s*c_s*T] * vol / dt
  !          \____ fluid thermal energy ___/  \__ solid thermal energy __/
  Res(thcflow_temperature_dof) = volume_over_dt * temperature * &
    (por_sat * den_kg * c_p + (1.d0 - porosity) * rho_s_c_s)

  if (calculate_derivatives) then

    ! ===================================================================== !
    ! Flow equation accumulation Jacobian
    !   A_flow = phi * S_l * rho_kmol * V/dt
    ! Note: rho_kmol derivatives are obtained from the stored mass-density
    ! derivatives via d(rho_kmol)/dX = d(rho_kg)/dX / FMWH2O.
    ! ===================================================================== !
    ! dA_flow/dP = V/dt * [rho_kmol * d(por_sat)/dP + por_sat * d(rho_kmol)/dP]
    Jac(thcflow_pressure_dof,thcflow_pressure_dof) = volume_over_dt * &
      (den_kmol * dpor_sat_dp + &
       por_sat * thcflow_auxvar%dden_dp / FMWH2O)
    ! dA_flow/dT = V/dt * por_sat * d(rho_kmol)/dT
    Jac(thcflow_pressure_dof,thcflow_temperature_dof) = volume_over_dt * &
      por_sat * thcflow_auxvar%dden_dT / FMWH2O
    ! dA_flow/dC = V/dt * por_sat * d(rho_kmol)/dC
    Jac(thcflow_pressure_dof,thcflow_concentration_dof) = volume_over_dt * &
      por_sat * thcflow_auxvar%dden_dC / FMWH2O

    ! ===================================================================== !
    ! Solute equation accumulation Jacobian
    !   A_sol = C * 1000 * phi * S_l * V/dt   (density-independent)
    ! ===================================================================== !
    ! dA_sol/dP = C * 1000 * V/dt * d(por_sat)/dP
    Jac(thcflow_concentration_dof,thcflow_pressure_dof) = &
      thcflow_auxvar%conc * L_per_m3 * volume_over_dt * dpor_sat_dp
    ! dA_sol/dT = 0    (no T dependence in solute accumulation)
    Jac(thcflow_concentration_dof,thcflow_temperature_dof) = 0.d0
    ! dA_sol/dC = 1000 * phi * S_l * V/dt
    Jac(thcflow_concentration_dof,thcflow_concentration_dof) = &
      L_per_m3 * volume_over_dt * por_sat

    ! ===================================================================== !
    ! Energy equation accumulation Jacobian
    !   A_energy = [phi*S_l*rho_kg*c_p*T + (1-phi)*rho_s*c_s*T] * V/dt
    ! ===================================================================== !
    ! dA_e/dP = V/dt * T * [ rho_kg*c_p*d(por_sat)/dP
    !                        + phi*S_l*c_p*d(rho_kg)/dP
    !                        - rho_s*c_s*dphi/dP ]
    Jac(thcflow_temperature_dof,thcflow_pressure_dof) = volume_over_dt * &
      temperature * (den_kg * c_p * dpor_sat_dp + &
                     por_sat * c_p * thcflow_auxvar%dden_dp - &
                     rho_s_c_s * thcflow_auxvar%dpor_dp)
    ! dA_e/dT = V/dt * [ phi*S_l*rho_kg*c_p + (1-phi)*rho_s*c_s
    !                    + phi*S_l*c_p*T*d(rho_kg)/dT ]
    !   first two terms = (rho c)_eff (effective volumetric heat capacity);
    !   third term = temperature-dependent-density correction.
    Jac(thcflow_temperature_dof,thcflow_temperature_dof) = volume_over_dt * &
      (por_sat * den_kg * c_p + (1.d0 - porosity) * rho_s_c_s + &
       por_sat * c_p * temperature * thcflow_auxvar%dden_dT)
    ! dA_e/dC = V/dt * phi*S_l*c_p*T*d(rho_kg)/dC
    Jac(thcflow_temperature_dof,thcflow_concentration_dof) = volume_over_dt * &
      por_sat * c_p * temperature * thcflow_auxvar%dden_dC

  endif

end subroutine THCFlowAccumulation

! ************************************************************************** !

subroutine THCFlowFlux(thcflow_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       thcflow_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       area, dist, &
                       thcflow_parameter, &
                       option,v_darcy,Res, &
                       Jup,Jdn, &
                       dResdparamup,dResdparamdn, &
                       calculate_derivatives, &
                       debug_connection)
  !
  ! Computes the internal flux terms (flow / energy / solute) for the THCFLOW
  ! residual and the full 3x3 (P,T,C) Jacobian blocks for a single internal
  ! connection.  Mirrors ZFlowFluxHarmonicPermOnly (zflow_common.F90), extended
  ! to variable density / viscosity and the mandatory 3-DOF system.
  !
  ! Harmonic averaging is used for BOTH the intrinsic permeability AND the
  ! effective thermal conductivity.  A single
  ! upwind direction -- the sign of the total potential difference -- governs
  ! the mobility, the transported density/enthalpy, and the advected
  ! concentration for all three equations.  The sign
  ! convention is identical to ZFLOW
  !
  ! dResdparamup / dResdparamdn are retained for interface compatibility but are
  ! a NO-OP here (zeroed, never filled) -- THCFLOW does not yet support the
  ! ZFLOW adjoint-parameter sensitivities.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  use Option_module
  use Material_Aux_module
  use Connection_module

  implicit none

  type(thcflow_auxvar_type) :: thcflow_auxvar_up, thcflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thcflow_parameter_type) :: thcflow_parameter
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jup(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: Jdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamup(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal, parameter :: L_per_m3 = 1.d3

  PetscReal :: dist_gravity      ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscReal :: g2                ! 0.5 * dist_gravity (face-average gravity wgt)

  ! permeability / mobility / volumetric flux
  PetscReal :: perm_up, perm_dn
  PetscReal :: perm_ave_over_dist
  PetscReal :: numerator, denominator
  PetscReal :: denom_D, denom_kappa
  PetscReal :: rho_avg, gravity_term, delta_pressure
  PetscReal :: kr, dkr_dpup, dkr_dpdn
  PetscReal :: vis_upw, mobility, kr_over_vis2
  PetscReal :: Gamma             ! geometric transmissibility = k_harm/d * area
  PetscReal :: q
  PetscReal :: wup, wdn          ! upwind weights (1/0) for transported scalars

  ! potential derivatives
  PetscReal :: dPi_dpup, dPi_dpdn, dPi_dTup, dPi_dTdn, dPi_dCup, dPi_dCdn
  ! mobility derivatives
  PetscReal :: dmob_dpup, dmob_dpdn, dmob_dTup, dmob_dTdn, dmob_dCup, dmob_dCdn
  ! volumetric flux derivatives
  PetscReal :: dq_dpup, dq_dpdn, dq_dTup, dq_dTdn, dq_dCup, dq_dCdn

  ! upwind-selected transported quantities
  PetscReal :: den_kmol_upw, den_kg_upw, temp_upw, conc_upw, h_upw
  ! molar-density derivative routing (kmol = mass / FMWH2O)
  PetscReal :: ddkmol_dp_up, ddkmol_dp_dn
  PetscReal :: ddkmol_dT_up, ddkmol_dT_dn
  PetscReal :: ddkmol_dC_up, ddkmol_dC_dn
  ! mass-density derivative routing (for enthalpy)
  PetscReal :: ddkg_dp_up, ddkg_dp_dn
  PetscReal :: ddkg_dT_up, ddkg_dT_dn
  PetscReal :: ddkg_dC_up, ddkg_dC_dn
  ! enthalpy derivatives
  PetscReal :: dh_dp_up, dh_dp_dn, dh_dT_up, dh_dT_dn, dh_dC_up, dh_dC_dn

  ! solute diffusion (hydrodynamic dispersion)
  PetscReal :: c_p
  PetscReal :: D_hyd_up, D_hyd_dn
  PetscReal :: Deff_over_dist, Gamma_D, delta_conc
  PetscReal :: dD_hyd_up_dp, dD_hyd_up_dT, dD_hyd_up_dC
  PetscReal :: dD_hyd_dn_dp, dD_hyd_dn_dT, dD_hyd_dn_dC
  PetscReal :: facD_up, facD_dn
  PetscReal :: dGamma_D_dpup, dGamma_D_dpdn
  PetscReal :: dGamma_D_dTup, dGamma_D_dTdn
  PetscReal :: dGamma_D_dCup, dGamma_D_dCdn
  PetscBool :: diffusion_on

  ! energy conduction
  PetscReal :: kappa_up, kappa_dn
  PetscReal :: kappa_over_dist, Gamma_kappa, delta_temp
  PetscReal :: dkappa_up_dp, dkappa_dn_dp
  PetscReal :: dGamma_kappa_dpup, dGamma_kappa_dpdn
  PetscBool :: conduction_on

  PetscReal :: tempreal

  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0
  dResdparamup = 0.d0
  dResdparamdn = 0.d0
  v_darcy = 0.d0

  q = 0.d0
  dq_dpup = 0.d0; dq_dpdn = 0.d0
  dq_dTup = 0.d0; dq_dTdn = 0.d0
  dq_dCup = 0.d0; dq_dCdn = 0.d0

  c_p = thcflow_heat_capacity_liquid

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  g2 = 0.5d0 * dist_gravity

  ! --------------------------------------------------------------------------
  ! Harmonic-mean intrinsic permeability
  !   k_harm / d_total = (k_up*k_dn) / (d_up*k_dn + d_dn*k_up)
  !   Gamma = (k_harm/d_total) * area
  ! --------------------------------------------------------------------------
  call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  numerator = perm_up * perm_dn
  denominator = dist_up*perm_dn + dist_dn*perm_up
  perm_ave_over_dist = numerator / denominator
  Gamma = perm_ave_over_dist * area

  ! --------------------------------------------------------------------------
  ! Total potential difference with face-averaged density gravity
  ! rho_avg = (rho_up + rho_dn)/2
  ! --------------------------------------------------------------------------
  rho_avg = 0.5d0 * (thcflow_auxvar_up%den_kg + thcflow_auxvar_dn%den_kg)
  gravity_term = rho_avg * dist_gravity
  delta_pressure = thcflow_auxvar_up%pres - &
                   thcflow_auxvar_dn%pres + &
                   gravity_term

  ! single upwind direction for all three equations:
  ! sign(q) = sign(delta_pressure) since Gamma, mobility >= 0.
  if (delta_pressure >= 0.d0) then
    wup = 1.d0
    wdn = 0.d0
  else
    wup = 0.d0
    wdn = 1.d0
  endif

  ! --------------------------------------------------------------------------
  ! Upwind mobility lambda = kr_up / vis_up  (variable viscosity)
  ! --------------------------------------------------------------------------
  kr = 0.d0
  dkr_dpup = 0.d0
  dkr_dpdn = 0.d0
  if (thcflow_auxvar_up%kr + thcflow_auxvar_dn%kr > floweps) then
    if (thcflow_tensorial_rel_perm) then
      if (delta_pressure >= 0.d0) then
        call THCFlowAuxTensorialRelPerm(thcflow_auxvar_up, &
                thcflow_parameter% &
                  tensorial_rel_perm_exponent(:,material_auxvar_up%id), &
                dist,kr,dkr_dpup,option)
      else
        call THCFlowAuxTensorialRelPerm(thcflow_auxvar_dn, &
                thcflow_parameter% &
                  tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                dist,kr,dkr_dpdn,option)
      endif
    else
      if (delta_pressure >= 0.d0) then
        kr = thcflow_auxvar_up%kr
        dkr_dpup = thcflow_auxvar_up%dkr_dp
      else
        kr = thcflow_auxvar_dn%kr
        dkr_dpdn = thcflow_auxvar_dn%dkr_dp
      endif
    endif
  endif

  ! upwind viscosity (mobility uses the SAME upwind cell as kr)
  vis_upw = wup*thcflow_auxvar_up%vis + wdn*thcflow_auxvar_dn%vis
  mobility = kr / vis_upw

  ! q[m^3/sec] = Gamma * mobility * delta_pressure
  q = Gamma * mobility * delta_pressure
  v_darcy = perm_ave_over_dist * mobility * delta_pressure

  ! --------------------------------------------------------------------------
  ! Upwind-selected transported quantities
  ! --------------------------------------------------------------------------
  den_kmol_upw = wup*thcflow_auxvar_up%den_kmol + wdn*thcflow_auxvar_dn%den_kmol
  den_kg_upw   = wup*thcflow_auxvar_up%den_kg   + wdn*thcflow_auxvar_dn%den_kg
  temp_upw     = wup*thcflow_auxvar_up%temp     + wdn*thcflow_auxvar_dn%temp
  conc_upw     = wup*thcflow_auxvar_up%conc     + wdn*thcflow_auxvar_dn%conc
  h_upw        = den_kg_upw * c_p * temp_upw

  ! --------------------------------------------------------------------------
  ! Solute hydrodynamic dispersion with harmonic-mean D_eff
  !   D_hyd = phi * S_l * tau * D_mol      (D_mol = D_mol(T,C), per-cell)
  ! --------------------------------------------------------------------------
  D_hyd_up = thcflow_auxvar_up%effective_porosity * thcflow_auxvar_up%sat * &
             material_auxvar_up%tortuosity * thcflow_auxvar_up%diff_mol
  D_hyd_dn = thcflow_auxvar_dn%effective_porosity * thcflow_auxvar_dn%sat * &
             material_auxvar_dn%tortuosity * thcflow_auxvar_dn%diff_mol
  numerator = D_hyd_up * D_hyd_dn
  denom_D = dist_up*D_hyd_dn + dist_dn*D_hyd_up
  diffusion_on = PETSC_TRUE
  if (denom_D == 0.d0) then
    ! turn off diffusion
    numerator = 0.d0
    denom_D = 1.d0
    diffusion_on = PETSC_FALSE
  endif
  Deff_over_dist = numerator / denom_D
  Gamma_D = area * Deff_over_dist
  delta_conc = thcflow_auxvar_up%conc - thcflow_auxvar_dn%conc

  ! --------------------------------------------------------------------------
  ! Energy conduction with harmonic-mean effective thermal conductivity
  ! kappa_face/d = (k_up*k_dn)/(d_up*k_dn+d_dn*k_up)
  ! --------------------------------------------------------------------------
  kappa_up = thcflow_auxvar_up%therm_cond_eff
  kappa_dn = thcflow_auxvar_dn%therm_cond_eff
  numerator = kappa_up * kappa_dn
  denom_kappa = dist_up*kappa_dn + dist_dn*kappa_up
  conduction_on = PETSC_TRUE
  if (denom_kappa == 0.d0) then
    ! turn off conduction
    numerator = 0.d0
    denom_kappa = 1.d0
    conduction_on = PETSC_FALSE
  endif
  kappa_over_dist = numerator / denom_kappa
  Gamma_kappa = area * kappa_over_dist
  delta_temp = thcflow_auxvar_up%temp - thcflow_auxvar_dn%temp

  ! --------------------------------------------------------------------------
  ! Residuals (same sign convention as ZFLOW)
  ! --------------------------------------------------------------------------
  ! Flow [kmol/s]
  Res(thcflow_pressure_dof) = den_kmol_upw * q
  ! Solute [mol/s] : advection + hydrodynamic dispersion
  Res(thcflow_concentration_dof) = (conc_upw * q + Gamma_D * delta_conc) * &
                                   L_per_m3
  ! Energy [W] : advected enthalpy + Fourier conduction
  Res(thcflow_temperature_dof) = h_upw * q + Gamma_kappa * delta_temp

  if (calculate_derivatives) then

    ! ===================================================================== !
    ! Building blocks
    ! ===================================================================== !
    ! Potential derivatives  dPi = d(P_up - P_dn + rho_avg*g*dz)
    dPi_dpup =  1.d0 + thcflow_auxvar_up%dden_dp * g2
    dPi_dpdn = -1.d0 + thcflow_auxvar_dn%dden_dp * g2
    dPi_dTup =         thcflow_auxvar_up%dden_dT * g2
    dPi_dTdn =         thcflow_auxvar_dn%dden_dT * g2
    dPi_dCup =         thcflow_auxvar_up%dden_dC * g2
    dPi_dCdn =         thcflow_auxvar_dn%dden_dC * g2

    ! Mobility derivatives  lambda = kr_upw / vis_upw  (upwind only)
    kr_over_vis2 = kr / (vis_upw*vis_upw)
    dmob_dpup = dkr_dpup / vis_upw
    dmob_dpdn = dkr_dpdn / vis_upw
    dmob_dTup = -kr_over_vis2 * (wup*thcflow_auxvar_up%dvis_dT)
    dmob_dTdn = -kr_over_vis2 * (wdn*thcflow_auxvar_dn%dvis_dT)
    dmob_dCup = -kr_over_vis2 * (wup*thcflow_auxvar_up%dvis_dC)
    dmob_dCdn = -kr_over_vis2 * (wdn*thcflow_auxvar_dn%dvis_dC)

    ! Volumetric flux derivatives  q = Gamma * lambda * dPi (product rule)
    dq_dpup = Gamma * (dmob_dpup*delta_pressure + mobility*dPi_dpup)
    dq_dpdn = Gamma * (dmob_dpdn*delta_pressure + mobility*dPi_dpdn)
    dq_dTup = Gamma * (dmob_dTup*delta_pressure + mobility*dPi_dTup)
    dq_dTdn = Gamma * (dmob_dTdn*delta_pressure + mobility*dPi_dTdn)
    dq_dCup = Gamma * (dmob_dCup*delta_pressure + mobility*dPi_dCup)
    dq_dCdn = Gamma * (dmob_dCdn*delta_pressure + mobility*dPi_dCdn)

    ! Upwind density derivative routing (kmol = mass/FMWH2O)
    ddkmol_dp_up = wup*thcflow_auxvar_up%dden_dp / FMWH2O
    ddkmol_dp_dn = wdn*thcflow_auxvar_dn%dden_dp / FMWH2O
    ddkmol_dT_up = wup*thcflow_auxvar_up%dden_dT / FMWH2O
    ddkmol_dT_dn = wdn*thcflow_auxvar_dn%dden_dT / FMWH2O
    ddkmol_dC_up = wup*thcflow_auxvar_up%dden_dC / FMWH2O
    ddkmol_dC_dn = wdn*thcflow_auxvar_dn%dden_dC / FMWH2O
    ddkg_dp_up   = wup*thcflow_auxvar_up%dden_dp
    ddkg_dp_dn   = wdn*thcflow_auxvar_dn%dden_dp
    ddkg_dT_up   = wup*thcflow_auxvar_up%dden_dT
    ddkg_dT_dn   = wdn*thcflow_auxvar_dn%dden_dT
    ddkg_dC_up   = wup*thcflow_auxvar_up%dden_dC
    ddkg_dC_dn   = wdn*thcflow_auxvar_dn%dden_dC

    ! Diffusive transmissibility derivatives (harmonic D_eff)
    dGamma_D_dpup = 0.d0; dGamma_D_dpdn = 0.d0
    dGamma_D_dTup = 0.d0; dGamma_D_dTdn = 0.d0
    dGamma_D_dCup = 0.d0; dGamma_D_dCdn = 0.d0
    if (diffusion_on) then
      ! cell-level D_eff derivatives
      dD_hyd_up_dp = (thcflow_auxvar_up%dsat_dp * &
                      thcflow_auxvar_up%effective_porosity + &
                      thcflow_auxvar_up%sat * thcflow_auxvar_up%dpor_dp) * &
                     material_auxvar_up%tortuosity * thcflow_auxvar_up%diff_mol
      dD_hyd_up_dT = thcflow_auxvar_up%effective_porosity * &
                     thcflow_auxvar_up%sat * material_auxvar_up%tortuosity * &
                     thcflow_auxvar_up%ddiff_dT
      dD_hyd_up_dC = thcflow_auxvar_up%effective_porosity * &
                     thcflow_auxvar_up%sat * material_auxvar_up%tortuosity * &
                     thcflow_auxvar_up%ddiff_dC
      dD_hyd_dn_dp = (thcflow_auxvar_dn%dsat_dp * &
                      thcflow_auxvar_dn%effective_porosity + &
                      thcflow_auxvar_dn%sat * thcflow_auxvar_dn%dpor_dp) * &
                     material_auxvar_dn%tortuosity * thcflow_auxvar_dn%diff_mol
      dD_hyd_dn_dT = thcflow_auxvar_dn%effective_porosity * &
                     thcflow_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thcflow_auxvar_dn%ddiff_dT
      dD_hyd_dn_dC = thcflow_auxvar_dn%effective_porosity * &
                     thcflow_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thcflow_auxvar_dn%ddiff_dC
      ! d(Deff_over_dist)/dD_hyd_up = D_hyd_dn^2 * dist_up / denom^2  (ZFLOW form)
      tempreal = denom_D * denom_D
      facD_up = area * D_hyd_dn * D_hyd_dn * dist_up / tempreal
      facD_dn = area * D_hyd_up * D_hyd_up * dist_dn / tempreal
      dGamma_D_dpup = facD_up * dD_hyd_up_dp
      dGamma_D_dpdn = facD_dn * dD_hyd_dn_dp
      dGamma_D_dTup = facD_up * dD_hyd_up_dT
      dGamma_D_dTdn = facD_dn * dD_hyd_dn_dT
      dGamma_D_dCup = facD_up * dD_hyd_up_dC
      dGamma_D_dCdn = facD_dn * dD_hyd_dn_dC
    endif

    ! Thermal transmissibility derivatives (harmonic kappa; kappa = kappa(S_l))
    dGamma_kappa_dpup = 0.d0
    dGamma_kappa_dpdn = 0.d0
    if (conduction_on) then
      dkappa_up_dp = thcflow_auxvar_up%dtherm_cond_dsat * &
                     thcflow_auxvar_up%dsat_dp
      dkappa_dn_dp = thcflow_auxvar_dn%dtherm_cond_dsat * &
                     thcflow_auxvar_dn%dsat_dp
      tempreal = denom_kappa * denom_kappa
      dGamma_kappa_dpup = area * kappa_dn * kappa_dn * dist_up / tempreal * &
                          dkappa_up_dp
      dGamma_kappa_dpdn = area * kappa_up * kappa_up * dist_dn / tempreal * &
                          dkappa_dn_dp
    endif

    ! Enthalpy derivatives  h = rho_kg_upw * c_p * T_upw
    dh_dp_up = c_p * temp_upw * ddkg_dp_up
    dh_dp_dn = c_p * temp_upw * ddkg_dp_dn
    dh_dT_up = c_p * (ddkg_dT_up*temp_upw + den_kg_upw*wup)
    dh_dT_dn = c_p * (ddkg_dT_dn*temp_upw + den_kg_upw*wdn)
    dh_dC_up = c_p * temp_upw * ddkg_dC_up
    dh_dC_dn = c_p * temp_upw * ddkg_dC_dn

    ! ===================================================================== !
    ! Flow equation flux Jacobian  F_flow = rho_kmol_upw * q
    ! ===================================================================== !
    Jup(thcflow_pressure_dof,thcflow_pressure_dof) = &
      ddkmol_dp_up*q + den_kmol_upw*dq_dpup
    Jup(thcflow_pressure_dof,thcflow_temperature_dof) = &
      ddkmol_dT_up*q + den_kmol_upw*dq_dTup
    Jup(thcflow_pressure_dof,thcflow_concentration_dof) = &
      ddkmol_dC_up*q + den_kmol_upw*dq_dCup
    Jdn(thcflow_pressure_dof,thcflow_pressure_dof) = &
      ddkmol_dp_dn*q + den_kmol_upw*dq_dpdn
    Jdn(thcflow_pressure_dof,thcflow_temperature_dof) = &
      ddkmol_dT_dn*q + den_kmol_upw*dq_dTdn
    Jdn(thcflow_pressure_dof,thcflow_concentration_dof) = &
      ddkmol_dC_dn*q + den_kmol_upw*dq_dCdn

    ! ===================================================================== !
    ! Solute equation flux Jacobian
    !   F_sol = (C_upw*q + Gamma_D*dC) * 1000
    !   d(C_upw)/dC_up = wup, d(C_upw)/dC_dn = wdn ; d(dC)/dC_up = +1, dn = -1
    ! ===================================================================== !
    Jup(thcflow_concentration_dof,thcflow_pressure_dof) = &
      (conc_upw*dq_dpup + dGamma_D_dpup*delta_conc) * L_per_m3
    Jup(thcflow_concentration_dof,thcflow_temperature_dof) = &
      (conc_upw*dq_dTup + dGamma_D_dTup*delta_conc) * L_per_m3
    Jup(thcflow_concentration_dof,thcflow_concentration_dof) = &
      (wup*q + conc_upw*dq_dCup + dGamma_D_dCup*delta_conc + Gamma_D) * L_per_m3
    Jdn(thcflow_concentration_dof,thcflow_pressure_dof) = &
      (conc_upw*dq_dpdn + dGamma_D_dpdn*delta_conc) * L_per_m3
    Jdn(thcflow_concentration_dof,thcflow_temperature_dof) = &
      (conc_upw*dq_dTdn + dGamma_D_dTdn*delta_conc) * L_per_m3
    Jdn(thcflow_concentration_dof,thcflow_concentration_dof) = &
      (wdn*q + conc_upw*dq_dCdn + dGamma_D_dCdn*delta_conc - Gamma_D) * L_per_m3

    ! ===================================================================== !
    ! Energy equation flux Jacobian
    !   F_energy = h_upw*q + Gamma_kappa*dT
    !   d(dT)/dT_up = +1, d(dT)/dT_dn = -1 ; Gamma_kappa indep of T,C
    ! ===================================================================== !
    Jup(thcflow_temperature_dof,thcflow_pressure_dof) = &
      dh_dp_up*q + h_upw*dq_dpup + dGamma_kappa_dpup*delta_temp
    Jup(thcflow_temperature_dof,thcflow_temperature_dof) = &
      dh_dT_up*q + h_upw*dq_dTup + Gamma_kappa
    Jup(thcflow_temperature_dof,thcflow_concentration_dof) = &
      dh_dC_up*q + h_upw*dq_dCup
    Jdn(thcflow_temperature_dof,thcflow_pressure_dof) = &
      dh_dp_dn*q + h_upw*dq_dpdn + dGamma_kappa_dpdn*delta_temp
    Jdn(thcflow_temperature_dof,thcflow_temperature_dof) = &
      dh_dT_dn*q + h_upw*dq_dTdn - Gamma_kappa
    Jdn(thcflow_temperature_dof,thcflow_concentration_dof) = &
      dh_dC_dn*q + h_upw*dq_dCdn

  endif

  if (debug_connection) then
    write(*,'(9x,"IFF(res_f,kr,dp): ",8es12.4)') Res(thcflow_pressure_dof), &
                                                 kr, delta_pressure
    write(*,'(9x,"IFE(res_e,Gk,dT): ",8es12.4)') Res(thcflow_temperature_dof), &
                                                 Gamma_kappa, delta_temp
    write(*,'(9x,"IFS(res_s,q,dC): ",8es12.4)') Res(thcflow_concentration_dof),&
                                                q, delta_conc
  endif

end subroutine THCFlowFlux

! ************************************************************************** !

subroutine THCFlowBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         thcflow_auxvar_up,global_auxvar_up, &
                         thcflow_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist, &
                         thcflow_parameter, &
                         option,v_darcy,Res,Jdn, &
                         dResdparamdn, &
                         calculate_derivatives, &
                         debug_connection)
  !
  ! Computes the boundary flux terms (flow / energy / solute) for the THCFLOW
  ! residual and the 3x3 boundary Jacobian block (interior cell only).
  ! Mirrors ZFlowBCFluxHarmonicPermOnly (zflow_common.F90), extended to the
  ! mandatory 3-DOF system with variable density / viscosity.
  !
  ! The "up" cell is the boundary ghost (its state is fixed by the BC), the
  ! "dn" cell is the interior cell.  Only Jdn (the interior cell) is filled.
  !
  ! Flow / solute BCs:
  !   DIRICHLET / SEEPAGE / CONDUCTANCE / HYDROSTATIC(+seepage/conductance) /
  !   PONDED_WATER  -- pressure-driven flux with variable-density hydrostatic
  !                    gravity term rho_avg = (rho_up + rho_dn)/2 ;
  !   NEUMANN       -- prescribed Darcy velocity.
  !
  ! Thermal BCs, keyed on the energy condition slot:
  !   DIRICHLET             -- prescribed T at ghost: advection + conduction ;
  !   NEUMANN               -- prescribed heat flux q_T [W/m^2] (no Jacobian) ;
  !   THCFLOW_CONVECTIVE_BC -- Robin / third-kind: h_conv*(T_ext - T_cell) .
  !
  ! Injection vs extraction for T and C is by the sign of the volumetric flux q
  ! inflow carries the boundary (ghost) state,
  ! outflow carries the interior cell state -- a single upwind direction shared
  ! by all three equations (consistent with THCFlowFlux and ZFLOW).
  !
  ! dResdparamdn is retained for interface compatibility but is a NO-OP here.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  use Option_module
  use Material_Aux_module
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(:)
  PetscInt :: auxvar_mapping(THCFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(thcflow_auxvar_type) :: thcflow_auxvar_up, thcflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thcflow_parameter_type) :: thcflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal, parameter :: L_per_m3 = 1.d3

  PetscInt :: bc_type, energy_bc_type
  PetscInt :: idof
  PetscReal :: perm_dn, perm_ave_over_dist
  PetscReal :: dist_gravity, g2
  PetscReal :: rho_avg, gravity_term, delta_pressure, ddelta_pressure_dpdn
  PetscReal :: boundary_pressure
  PetscReal :: kr, dkr_dpdn
  PetscReal :: vis_upw, mobility, kr_over_vis2
  PetscReal :: Gamma_bc, q
  PetscReal :: wup, wdn
  PetscBool :: derivative_toggle

  ! potential / mobility / flux derivatives (interior-cell only)
  PetscReal :: dPi_dpdn, dPi_dTdn, dPi_dCdn
  PetscReal :: dmob_dpdn, dmob_dTdn, dmob_dCdn
  PetscReal :: dq_dpdn, dq_dTdn, dq_dCdn

  ! upwind transported quantities
  PetscReal :: c_p
  PetscReal :: den_kmol_upw, den_kg_upw, temp_upw, conc_upw, h_upw
  PetscReal :: ddkmol_dp_dn, ddkmol_dT_dn, ddkmol_dC_dn
  PetscReal :: ddkg_dp_dn, ddkg_dT_dn, ddkg_dC_dn
  PetscReal :: dh_dp_dn, dh_dT_dn, dh_dC_dn

  ! solute diffusion
  PetscReal :: delta_conc, D_hyd_dn, Deff_over_dist, Gamma_D, dispersion_scale
  PetscReal :: dD_hyd_dn_dp, dD_hyd_dn_dT, dD_hyd_dn_dC
  PetscReal :: dGamma_D_dpdn, dGamma_D_dTdn, dGamma_D_dCdn

  ! energy conduction / prescribed-flux / convective
  PetscReal :: Gamma_kappa, delta_temp, dGamma_kappa_dpdn
  PetscReal :: heat_flux_bc, h_conv, temp_ext

  Res = 0.d0
  Jdn = 0.d0
  dResdparamdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  delta_pressure = 0.d0
  delta_conc = 0.d0
  ddelta_pressure_dpdn = 0.d0
  perm_ave_over_dist = 0.d0
  kr = 0.d0
  dkr_dpdn = 0.d0
  dq_dpdn = 0.d0; dq_dTdn = 0.d0; dq_dCdn = 0.d0
  wup = 0.d0; wdn = 0.d0
  Gamma_bc = 0.d0
  mobility = 0.d0
  vis_upw = thcflow_auxvar_dn%vis

  c_p = thcflow_heat_capacity_liquid

  ! dist(0) = scalar distance to boundary face; dist(1:3) = unit vector
  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
  g2 = 0.5d0 * dist_gravity

  ! ========================================================================= !
  ! Flow equation -- determine the volumetric flux q
  ! ========================================================================= !
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  bc_type = ibndtype(auxvar_mapping(THCFLOW_COND_WATER_INDEX))
  derivative_toggle = PETSC_TRUE
  select case(bc_type)
    case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
         HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
         PONDED_WATER_BC)
      if (thcflow_auxvar_up%kr + thcflow_auxvar_dn%kr > floweps) then

        if (bc_type == DIRICHLET_CONDUCTANCE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          ! conductance already includes 1/dist; viscosity applied via mobility
          idof = auxvar_mapping(THCFLOW_COND_WATER_AUX_INDEX)
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn / dist(0)
        endif

        ! variable-density hydrostatic gravity term
        rho_avg = 0.5d0 * (thcflow_auxvar_up%den_kg + thcflow_auxvar_dn%den_kg)
        gravity_term = rho_avg * dist_gravity
        boundary_pressure = thcflow_auxvar_up%pres
        delta_pressure = boundary_pressure - thcflow_auxvar_dn%pres + &
                         gravity_term
        ddelta_pressure_dpdn = -1.d0

        select case(bc_type)
          case(DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
               PONDED_WATER_BC)
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                thcflow_auxvar_up%pres - &
                  option%flow%reference_pressure < eps) then
              delta_pressure = 0.d0
              ddelta_pressure_dpdn = 0.d0
            endif
        end select

        ! upwind mobility direction (shared by all three equations)
        if (delta_pressure >= 0.d0) then
          wup = 1.d0
          wdn = 0.d0
        else
          wup = 0.d0
          wdn = 1.d0
        endif

        if (thcflow_tensorial_rel_perm) then
          if (delta_pressure >= 0.d0) then
            call THCFlowAuxTensorialRelPerm(thcflow_auxvar_up, &
                    thcflow_parameter% &
                      tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                    dist,kr,dkr_dpdn,option)
            ! upwind (boundary) pressure is fixed
            dkr_dpdn = 0.d0
          else
            call THCFlowAuxTensorialRelPerm(thcflow_auxvar_dn, &
                    thcflow_parameter% &
                      tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                    dist,kr,dkr_dpdn,option)
          endif
        else
          if (delta_pressure >= 0.d0) then
            kr = thcflow_auxvar_up%kr
            dkr_dpdn = 0.d0
          else
            kr = thcflow_auxvar_dn%kr
            dkr_dpdn = thcflow_auxvar_dn%dkr_dp
          endif
        endif

        ! upwind viscosity (boundary state fixed; interior varies)
        vis_upw = wup*thcflow_auxvar_up%vis + wdn*thcflow_auxvar_dn%vis
        mobility = kr / vis_upw
        ! v_darcy[m/sec] = perm/dist * kr/mu * dP
        v_darcy = perm_ave_over_dist * mobility * delta_pressure

        if (bc_type == PONDED_WATER_BC) then
          idof = auxvar_mapping(THCFLOW_COND_WATER_INDEX)
          if (v_darcy > auxvars(idof) / option%flow_dt) then
            v_darcy = auxvars(idof) / option%flow_dt
            ddelta_pressure_dpdn = 0.d0
            dkr_dpdn = 0.d0
            derivative_toggle = PETSC_FALSE
          endif
        endif
      endif

    case(NEUMANN_BC)
      idof = auxvar_mapping(THCFLOW_COND_WATER_INDEX)
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy = auxvars(idof)
      endif
      derivative_toggle = PETSC_FALSE

    case default
      option%io_buffer = &
        'Boundary condition type (' // trim(StringWrite(bc_type)) // &
        ') not recognized in THCFlowBCFlux phase loop.'
      call PrintErrMsg(option)
  end select

  ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
  q = v_darcy * area
  Gamma_bc = perm_ave_over_dist * area

  ! Upwind direction for transported scalars. For Dirichlet the wup/wdn set
  ! above (sign of delta_pressure) already matches sign(q); for Neumann (and
  ! the ponded cap) reset it from the prescribed velocity.
  if (.not. derivative_toggle .or. bc_type == NEUMANN_BC) then
    if (q >= 0.d0) then
      wup = 1.d0; wdn = 0.d0
    else
      wup = 0.d0; wdn = 1.d0
    endif
  endif

  ! Upwind transported quantities
  den_kmol_upw = wup*thcflow_auxvar_up%den_kmol + wdn*thcflow_auxvar_dn%den_kmol
  den_kg_upw   = wup*thcflow_auxvar_up%den_kg   + wdn*thcflow_auxvar_dn%den_kg
  temp_upw     = wup*thcflow_auxvar_up%temp     + wdn*thcflow_auxvar_dn%temp
  conc_upw     = wup*thcflow_auxvar_up%conc     + wdn*thcflow_auxvar_dn%conc
  h_upw        = den_kg_upw * c_p * temp_upw

  ! --- advective residual contributions (all three equations) -------------
  ! Flow [kmol/s]
  Res(thcflow_pressure_dof) = den_kmol_upw * q
  ! Solute advection [mol/s] (diffusion added below)
  Res(thcflow_concentration_dof) = conc_upw * q * L_per_m3
  ! Energy advection [W] (conduction / prescribed flux added below)
  Res(thcflow_temperature_dof) = h_upw * q

  ! ========================================================================= !
  ! Solute diffusion at the boundary (zero-gradient on outflow)
  ! ========================================================================= !
  delta_conc = thcflow_auxvar_up%conc - thcflow_auxvar_dn%conc
  if (q >= 0.d0) then
    dispersion_scale = 1.d0   ! inflow: full diffusive exchange with ghost
  else
    dispersion_scale = 0.d0   ! outflow: zero-gradient outlet BC
  endif
  D_hyd_dn = thcflow_auxvar_dn%effective_porosity * thcflow_auxvar_dn%sat * &
             material_auxvar_dn%tortuosity * thcflow_auxvar_dn%diff_mol
  Deff_over_dist = dispersion_scale * D_hyd_dn / dist(0)
  Gamma_D = area * Deff_over_dist
  Res(thcflow_concentration_dof) = Res(thcflow_concentration_dof) + &
                                   Gamma_D * delta_conc * L_per_m3

  ! ========================================================================= !
  ! Energy conduction / prescribed flux / convective BC
  ! ========================================================================= !
  energy_bc_type = ibndtype(auxvar_mapping(THCFLOW_COND_ENERGY_INDEX))
  Gamma_kappa = 0.d0
  dGamma_kappa_dpdn = 0.d0
  delta_temp = thcflow_auxvar_up%temp - thcflow_auxvar_dn%temp
  select case(energy_bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,DIRICHLET_SEEPAGE_BC, &
         DIRICHLET_CONDUCTANCE_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      ! prescribed T at ghost: Fourier conduction over the half-cell distance,
      ! using the interior effective conductivity.
      Gamma_kappa = thcflow_auxvar_dn%therm_cond_eff / dist(0) * area
      Res(thcflow_temperature_dof) = Res(thcflow_temperature_dof) + &
                                     Gamma_kappa * delta_temp
      dGamma_kappa_dpdn = thcflow_auxvar_dn%dtherm_cond_dsat * &
                          thcflow_auxvar_dn%dsat_dp / dist(0) * area
    case(NEUMANN_BC)
      ! prescribed heat flux q_T [W/m^2]; +q_T means heat enters the domain
      ! No Jacobian (value independent of solution).
      heat_flux_bc = auxvars(auxvar_mapping(THCFLOW_COND_ENERGY_INDEX))
      Res(thcflow_temperature_dof) = Res(thcflow_temperature_dof) + &
                                     heat_flux_bc * area
    case(THCFLOW_CONVECTIVE_BC)
      ! Robin / third-kind: heat flux INTO domain = h_conv*(T_ext - T_cell)
      ! h_conv stored in the water-aux slot, T_ext in the
      ! energy slot (placeholder plumbing -- see module header).
      temp_ext = auxvars(auxvar_mapping(THCFLOW_COND_ENERGY_INDEX))
      h_conv = auxvars(auxvar_mapping(THCFLOW_COND_WATER_AUX_INDEX))
      Res(thcflow_temperature_dof) = Res(thcflow_temperature_dof) + &
        h_conv * area * (temp_ext - thcflow_auxvar_dn%temp)
    case(ZERO_GRADIENT_BC,DIRICHLET_ZERO_GRADIENT_BC,NULL_CONDITION)
      ! no conductive boundary term -- advection only
    case default
      ! treat any other energy BC type as zero-gradient (advection only)
  end select

  ! ========================================================================= !
  ! Boundary Jacobian (interior cell only)
  ! ========================================================================= !
  if (calculate_derivatives) then

    ! --- potential / mobility / flux derivatives (Dirichlet family only) ---
    if (derivative_toggle) then
      dPi_dpdn = ddelta_pressure_dpdn + thcflow_auxvar_dn%dden_dp * g2
      dPi_dTdn = thcflow_auxvar_dn%dden_dT * g2
      dPi_dCdn = thcflow_auxvar_dn%dden_dC * g2

      kr_over_vis2 = kr / (vis_upw*vis_upw)
      dmob_dpdn = dkr_dpdn / vis_upw
      dmob_dTdn = -kr_over_vis2 * (wdn*thcflow_auxvar_dn%dvis_dT)
      dmob_dCdn = -kr_over_vis2 * (wdn*thcflow_auxvar_dn%dvis_dC)

      dq_dpdn = Gamma_bc * (dmob_dpdn*delta_pressure + mobility*dPi_dpdn)
      dq_dTdn = Gamma_bc * (dmob_dTdn*delta_pressure + mobility*dPi_dTdn)
      dq_dCdn = Gamma_bc * (dmob_dCdn*delta_pressure + mobility*dPi_dCdn)
    else
      dq_dpdn = 0.d0; dq_dTdn = 0.d0; dq_dCdn = 0.d0
    endif

    ! interior-cell upwind density derivatives (nonzero only on extraction)
    ddkmol_dp_dn = wdn*thcflow_auxvar_dn%dden_dp / FMWH2O
    ddkmol_dT_dn = wdn*thcflow_auxvar_dn%dden_dT / FMWH2O
    ddkmol_dC_dn = wdn*thcflow_auxvar_dn%dden_dC / FMWH2O
    ddkg_dp_dn   = wdn*thcflow_auxvar_dn%dden_dp
    ddkg_dT_dn   = wdn*thcflow_auxvar_dn%dden_dT
    ddkg_dC_dn   = wdn*thcflow_auxvar_dn%dden_dC

    ! enthalpy derivatives  h = rho_kg_upw * c_p * T_upw
    dh_dp_dn = c_p * temp_upw * ddkg_dp_dn
    dh_dT_dn = c_p * (ddkg_dT_dn*temp_upw + den_kg_upw*wdn)
    dh_dC_dn = c_p * temp_upw * ddkg_dC_dn

    ! diffusive transmissibility derivatives
    dGamma_D_dpdn = 0.d0; dGamma_D_dTdn = 0.d0; dGamma_D_dCdn = 0.d0
    if (dispersion_scale > 0.d0) then
      dD_hyd_dn_dp = (thcflow_auxvar_dn%dsat_dp * &
                      thcflow_auxvar_dn%effective_porosity + &
                      thcflow_auxvar_dn%sat * thcflow_auxvar_dn%dpor_dp) * &
                     material_auxvar_dn%tortuosity * thcflow_auxvar_dn%diff_mol
      dD_hyd_dn_dT = thcflow_auxvar_dn%effective_porosity * &
                     thcflow_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thcflow_auxvar_dn%ddiff_dT
      dD_hyd_dn_dC = thcflow_auxvar_dn%effective_porosity * &
                     thcflow_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thcflow_auxvar_dn%ddiff_dC
      dGamma_D_dpdn = area * dispersion_scale * dD_hyd_dn_dp / dist(0)
      dGamma_D_dTdn = area * dispersion_scale * dD_hyd_dn_dT / dist(0)
      dGamma_D_dCdn = area * dispersion_scale * dD_hyd_dn_dC / dist(0)
    endif

    ! ----------------------------------------------------------------------- !
    ! Flow equation  F_flow = rho_kmol_upw * q
    ! ----------------------------------------------------------------------- !
    Jdn(thcflow_pressure_dof,thcflow_pressure_dof) = &
      ddkmol_dp_dn*q + den_kmol_upw*dq_dpdn
    Jdn(thcflow_pressure_dof,thcflow_temperature_dof) = &
      ddkmol_dT_dn*q + den_kmol_upw*dq_dTdn
    Jdn(thcflow_pressure_dof,thcflow_concentration_dof) = &
      ddkmol_dC_dn*q + den_kmol_upw*dq_dCdn

    ! ----------------------------------------------------------------------- !
    ! Solute equation  F_sol = (C_upw*q + Gamma_D*dC)*1000
    !   d(C_upw)/dC_dn = wdn ; d(dC)/dC_dn = -1
    ! ----------------------------------------------------------------------- !
    Jdn(thcflow_concentration_dof,thcflow_pressure_dof) = &
      (conc_upw*dq_dpdn + dGamma_D_dpdn*delta_conc) * L_per_m3
    Jdn(thcflow_concentration_dof,thcflow_temperature_dof) = &
      (conc_upw*dq_dTdn + dGamma_D_dTdn*delta_conc) * L_per_m3
    Jdn(thcflow_concentration_dof,thcflow_concentration_dof) = &
      (wdn*q + conc_upw*dq_dCdn + dGamma_D_dCdn*delta_conc - Gamma_D) * L_per_m3

    ! ----------------------------------------------------------------------- !
    ! Energy equation  F_energy = h_upw*q + Gamma_kappa*dT  (+ Robin)
    !   d(dT)/dT_dn = -1 ; Gamma_kappa indep of T,C
    ! ----------------------------------------------------------------------- !
    Jdn(thcflow_temperature_dof,thcflow_pressure_dof) = &
      dh_dp_dn*q + h_upw*dq_dpdn + dGamma_kappa_dpdn*delta_temp
    Jdn(thcflow_temperature_dof,thcflow_temperature_dof) = &
      dh_dT_dn*q + h_upw*dq_dTdn - Gamma_kappa
    Jdn(thcflow_temperature_dof,thcflow_concentration_dof) = &
      dh_dC_dn*q + h_upw*dq_dCdn
    ! convective (Robin) thermal BC: d/dT_cell of h_conv*(T_ext - T_cell)
    if (energy_bc_type == THCFLOW_CONVECTIVE_BC) then
      h_conv = auxvars(auxvar_mapping(THCFLOW_COND_WATER_AUX_INDEX))
      Jdn(thcflow_temperature_dof,thcflow_temperature_dof) = &
        Jdn(thcflow_temperature_dof,thcflow_temperature_dof) - h_conv * area
    endif

  endif

  if (debug_connection) then
    write(*,'(9x,"BFF(res_f,kr,dp): ",8es12.4)') Res(thcflow_pressure_dof), &
                                                 kr, delta_pressure
    write(*,'(9x,"BFE(res_e,Gk,dT): ",8es12.4)') Res(thcflow_temperature_dof), &
                                                 Gamma_kappa, delta_temp
    write(*,'(9x,"BFS(res_s,q,dC): ",8es12.4)') Res(thcflow_concentration_dof),&
                                                q, delta_conc
  endif

end subroutine THCFlowBCFlux

! ************************************************************************** !

subroutine THCFlowSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                          flow_src_sink_type, &
                          thcflow_auxvar,global_auxvar,material_auxvar, &
                          ss_flow_vol_flux,Res,Jdn,dResdparamdn, &
                          calculate_derivatives)
  !
  ! Computes the source/sink terms (flow / energy / solute) for the THCFLOW
  ! residual and its 3x3 Jacobian.  Mirrors ZFlowSrcSink (zflow_common.F90),
  ! extended to the energy equation and to a solution-dependent density.
  !
  ! Injection vs extraction is by the sign of the volumetric rate q_vol
  !   q_vol >= 0 (injection):  solute at C_inj, energy at T_inj (both fixed) ;
  !   q_vol <  0 (extraction): solute at C_cell, energy at T_cell (interior) .
  !
  ! dResdparamdn is retained for interface compatibility but is a NO-OP here.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  type(thcflow_auxvar_type) :: thcflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscBool :: calculate_derivatives

  PetscReal, parameter :: L_per_m3 = 1.d3
  PetscReal :: tempreal
  PetscReal :: qsrc_m3, qsrc_L
  PetscReal :: c_p, den_kg, conc_used, temp_used
  PetscInt :: dof_index, index_

  Res = 0.d0
  Jdn = 0.d0
  dResdparamdn = 0.d0

  c_p = thcflow_heat_capacity_liquid
  den_kg = thcflow_auxvar%den_kg

  ! ========================================================================= !
  ! Flow (liquid mass) source/sink : units kmol/s
  !   Q_flow = q_vol * rho_kmol(P,T,C)   (density now solution-dependent)
  ! ========================================================================= !
  dof_index = flow_src_sink_mapping(THCFLOW_COND_WATER_INDEX)
  qsrc_m3 = 0.d0
  tempreal = flow_aux_real_var(dof_index)
  select case(flow_src_sink_type(dof_index))
    case(VOLUMETRIC_RATE_SS)
      qsrc_m3 = tempreal
    case(SCALED_VOLUMETRIC_RATE_SS)
      index_ = flow_src_sink_mapping(THCFLOW_COND_WATER_AUX_INDEX)
      qsrc_m3 = tempreal*flow_aux_real_var(index_)
    case default
      option%io_buffer = 'src_sink_type not supported in THCFlowSrcSink'
      call PrintErrMsg(option)
  end select
  ss_flow_vol_flux = qsrc_m3
  qsrc_L = qsrc_m3 * L_per_m3

  ! Res[kmol/sec]
  Res(thcflow_pressure_dof) = qsrc_m3 * thcflow_auxvar%den_kmol

  if (calculate_derivatives) then
    ! dQ_flow/dX = q_vol * d(rho_kmol)/dX  (rho_kmol = rho_kg / FMWH2O)
    Jdn(thcflow_pressure_dof,thcflow_pressure_dof) = &
      qsrc_m3 * thcflow_auxvar%dden_dp / FMWH2O
    Jdn(thcflow_pressure_dof,thcflow_temperature_dof) = &
      qsrc_m3 * thcflow_auxvar%dden_dT / FMWH2O
    Jdn(thcflow_pressure_dof,thcflow_concentration_dof) = &
      qsrc_m3 * thcflow_auxvar%dden_dC / FMWH2O
  endif

  ! ========================================================================= !
  ! Solute source/sink : units mol/s
  ! ========================================================================= !
  if (qsrc_m3 >= 0.d0) then
    ! injection at prescribed C_inj
    dof_index = flow_src_sink_mapping(THCFLOW_COND_SOLUTE_INDEX)
    conc_used = flow_aux_real_var(dof_index)
  else
    ! extraction carries interior cell concentration
    conc_used = thcflow_auxvar%conc
  endif
  Res(thcflow_concentration_dof) = qsrc_L * conc_used
  if (calculate_derivatives) then
    if (qsrc_m3 < 0.d0) then
      ! dQ_sol/dC = q_vol*1000   (C_cell depends on solution); 0 for injection
      Jdn(thcflow_concentration_dof,thcflow_concentration_dof) = qsrc_L
    endif
  endif

  ! ========================================================================= !
  ! Energy source/sink : units W
  !   Q_energy = q_vol * rho_kg(P,T,C) * c_p * T_used
  ! ========================================================================= !
  if (qsrc_m3 >= 0.d0) then
    ! injection at prescribed T_inj (fixed)
    dof_index = flow_src_sink_mapping(THCFLOW_COND_ENERGY_INDEX)
    temp_used = flow_aux_real_var(dof_index)
  else
    ! extraction carries interior cell temperature
    temp_used = thcflow_auxvar%temp
  endif
  Res(thcflow_temperature_dof) = qsrc_m3 * den_kg * c_p * temp_used
  if (calculate_derivatives) then
    ! density-coupling terms (present for both injection and extraction)
    Jdn(thcflow_temperature_dof,thcflow_pressure_dof) = &
      qsrc_m3 * c_p * temp_used * thcflow_auxvar%dden_dp
    Jdn(thcflow_temperature_dof,thcflow_concentration_dof) = &
      qsrc_m3 * c_p * temp_used * thcflow_auxvar%dden_dC
    if (qsrc_m3 >= 0.d0) then
      ! injection: T_inj fixed -> only density's T-dependence
      Jdn(thcflow_temperature_dof,thcflow_temperature_dof) = &
        qsrc_m3 * c_p * temp_used * thcflow_auxvar%dden_dT
    else
      ! extraction: T_used = T_cell -> product rule on rho_kg*T_cell
      Jdn(thcflow_temperature_dof,thcflow_temperature_dof) = &
        qsrc_m3 * c_p * (den_kg + temp_used * thcflow_auxvar%dden_dT)
    endif
  endif

end subroutine THCFlowSrcSink

! ************************************************************************** !

subroutine THCFlowAccumDerivative(thcflow_auxvar,global_auxvar, &
                                  material_auxvar, &
                                  option,Res,Jac,dResdparam)
  !
  ! Computes derivatives of the accumulation term for the Jacobian.
  !
  ! When thcflow_numerical_derivatives is .false. (default), the analytical 3x3
  ! accumulation block computed by THCFlowAccumulation is returned unchanged.
  ! When .true., the block is formed by first-order forward finite differences
  ! of the residual w.r.t. each perturbed primary variable (P,T,C):
  !   Jac(ieq,idof) = [Res(x+pert*e_idof) - Res(x)] / pert .
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(thcflow_auxvar_type) :: thcflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jac(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparam(THCFLOW_NDOF,THCFLOW_NDOF)

  PetscReal :: res_pert(THCFLOW_NDOF)
  PetscReal :: Jdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dJdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscInt :: idof, ieq

  call THCFlowAccumulation(thcflow_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar, &
                           option,Res,Jac,dResdparam, &
                           .not.thcflow_numerical_derivatives)

  if (thcflow_numerical_derivatives) then
    do idof = 1, option%nflowdof
      call THCFlowAccumulation(thcflow_auxvar(idof), &
                               global_auxvar, &
                               material_auxvar, &
                               option,res_pert,Jdum,dJdum, &
                               PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq))/thcflow_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine THCFlowAccumDerivative

! ************************************************************************** !

subroutine THCFlowFluxDerivative(thcflow_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 thcflow_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 area,dist,thcflow_parameter,option,v_darcy, &
                                 Res,Jup,Jdn,dResdparamup,dResdparamdn, &
                                 debug_connection)
  !
  ! Computes the derivatives of the internal flux terms for the Jacobian.
  !
  ! Analytical path (thcflow_numerical_derivatives = .false.): returns the
  ! analytical Jup/Jdn from THCFlowFlux.  Numerical path: perturbs upgradient
  ! DOFs (-> Jup) and then the downgradient DOFs (-> Jdn) one at a time & forms
  ! forward differences of the residual.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(thcflow_auxvar_type) :: thcflow_auxvar_up(0:), thcflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thcflow_parameter_type) :: thcflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jup(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: Jdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamup(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(THCFLOW_NDOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dJdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscInt :: idof, ieq

  Jup = 0.d0
  Jdn = 0.d0

  call THCFlowFlux(thcflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   thcflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   area,dist,thcflow_parameter,option,v_darcy, &
                   Res,Jup,Jdn,dResdparamup,dResdparamdn, &
                   .not.thcflow_numerical_derivatives, &
                   debug_connection)

  if (thcflow_numerical_derivatives) then
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call THCFlowFlux(thcflow_auxvar_up(idof),global_auxvar_up, &
                       material_auxvar_up, &
                       thcflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,thcflow_parameter,option,v_darcy_dum, &
                       res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jup(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thcflow_auxvar_up(idof)%pert
      enddo
    enddo
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call THCFlowFlux(thcflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                       material_auxvar_up, &
                       thcflow_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,thcflow_parameter,option,v_darcy_dum, &
                       res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thcflow_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine THCFlowFluxDerivative

! ************************************************************************** !

subroutine THCFlowBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   thcflow_auxvar_up, &
                                   global_auxvar_up, &
                                   thcflow_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   area,dist,thcflow_parameter,option,v_darcy, &
                                   Res,Jdn,dResdparamdn,debug_connection)
  !
  ! Computes the derivatives of the boundary flux terms for the Jacobian.
  ! Only the interior ("dn") cell DOFs are perturbed; the "up" cell is the fixed
  ! boundary ghost.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(:)
  PetscInt :: auxvar_mapping(THCFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(thcflow_auxvar_type) :: thcflow_auxvar_up, thcflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thcflow_parameter_type) :: thcflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparamdn(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(THCFLOW_NDOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dJdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscInt :: idof, ieq

  Jdn = 0.d0

  call THCFlowBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     thcflow_auxvar_up,global_auxvar_up, &
                     thcflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist,thcflow_parameter,option,v_darcy, &
                     Res,Jdn,dResdparamdn,.not.thcflow_numerical_derivatives, &
                     debug_connection)

  if (thcflow_numerical_derivatives) then
    ! downgradient (interior cell) derivatives
    do idof = 1, option%nflowdof
      call THCFlowBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         thcflow_auxvar_up,global_auxvar_up, &
                         thcflow_auxvar_dn(idof),global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist,thcflow_parameter,option,v_darcy_dum, &
                         res_pert,Jdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thcflow_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine THCFlowBCFluxDerivative

! ************************************************************************** !

subroutine THCFlowSrcSinkDerivative(option,flow_aux_real_var, &
                                    flow_src_sink_mapping, &
                                    flow_src_sink_type, &
                                    thcflow_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    ss_flow_vol_flux, &
                                    Res,Jac,dResdparam)
  !
  ! Computes the derivatives of the source/sink terms for the Jacobian.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  type(thcflow_auxvar_type) :: thcflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(THCFLOW_NDOF)
  PetscReal :: Jac(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dResdparam(THCFLOW_NDOF,THCFLOW_NDOF)

  PetscReal :: res_pert(THCFLOW_NDOF)
  PetscReal :: dummy_real
  PetscReal :: Jdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscReal :: dJdum(THCFLOW_NDOF,THCFLOW_NDOF)
  PetscInt :: idof, ieq

  Jac = 0.d0
  ! unperturbed thcflow_auxvar's value
  call THCFlowSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                      flow_src_sink_type, &
                      thcflow_auxvar(ZERO_INTEGER),global_auxvar, &
                      material_auxvar, &
                      ss_flow_vol_flux, &
                      Res,Jac,dResdparam,.not.thcflow_numerical_derivatives)

  if (thcflow_numerical_derivatives) then
    ! perturbed thcflow_auxvar's value
    do idof = 1, option%nflowdof
      call THCFlowSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                          flow_src_sink_type, &
                          thcflow_auxvar(idof),global_auxvar, &
                          material_auxvar, &
                          dummy_real, &
                          res_pert,Jdum,dJdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thcflow_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine THCFlowSrcSinkDerivative

! ************************************************************************** !

end module THCFlow_Common_module
