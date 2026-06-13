module Compositional_Flow_Common_module
#include "petsc/finclude/petscsys.h"

  use Compositional_Flow_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use petscsys

  implicit none

  private


#define SOLUTE_TRANSPORT

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 0.d0

  procedure(FluxDummy), pointer :: XXFlux => null()
  procedure(BCFluxDummy), pointer :: XXBCFlux => null()

  interface
    subroutine FluxDummy(compflow_auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         compflow_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area, dist, &
                         compflow_parameter, &
                         option,v_darcy,Res,Jup,Jdn, &
                         calculate_derivatives, &
                         debug_connection)
      use petscsys
      use Compositional_Flow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(compflow_auxvar_type) :: compflow_auxvar_up, compflow_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(compflow_parameter_type) :: compflow_parameter
      PetscReal :: Res(COMPFLOW_MAX_DOF)
      PetscReal :: Jup(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
      PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
      PetscBool :: calculate_derivatives
      PetscBool :: debug_connection
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                           compflow_auxvar_up,global_auxvar_up, &
                           compflow_auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           area,dist, &
                           compflow_parameter, &
                           option,v_darcy,Res,Jdn, &
                           calculate_derivatives, &
                           debug_connection)
      use petscsys
      use Compositional_Flow_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(option_type) :: option
      PetscInt :: ibndtype(1)
      PetscInt :: auxvar_mapping(COMPFLOW_MAX_INDEX)
      PetscReal :: auxvars(:) ! from aux_real_var array
      type(compflow_auxvar_type) :: compflow_auxvar_up, compflow_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_dn
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(compflow_parameter_type) :: compflow_parameter
      PetscReal :: v_darcy
      PetscReal :: Res(COMPFLOW_MAX_DOF)
      PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
      PetscBool :: calculate_derivatives
      PetscBool :: debug_connection
    end subroutine BCFluxDummy
  end interface

  public :: CompFlowAccumulation, &
            CompFlowFluxHarmonicPermOnly, &
            CompFlowBCFluxHarmonicPermOnly, &
            CompFlowAccumDerivative, &
            XXFluxDerivative, &
            XXBCFluxDerivative, &
            CompFlowSrcSinkDerivative

  public :: XXFlux, &
            XXBCFlux

contains

! ************************************************************************** !

subroutine CompFlowAccumulation(compflow_auxvar,global_auxvar,material_auxvar, &
                             option,Res,Jac,calculate_derivatives)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(compflow_auxvar_type) :: compflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jac(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: porosity
  PetscReal :: saturation
  PetscReal :: por_sat
  PetscReal :: volume_over_dt
  PetscReal :: tempreal
  PetscReal, parameter :: L_per_m3 = 1.d3

  Res = 0.d0
  Jac = 0.d0

    ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use compflow_auxvar%effective porosity here as it enables numerical
  ! derivatives to be employed
  porosity = compflow_auxvar%effective_porosity
  saturation = compflow_auxvar%sat
  por_sat = saturation * porosity

  if (compflow_liq_flow_eq > 0) then
    ! accumulation term units = kg/s
    ! Res[kg/sec] = den[kg/m^3] * sat[m^3 liquid/m^3 void] *
    !               por[m^3 void/m^3 bulk] * vol[m^3 bulk] / dt[sec]
    tempreal = volume_over_dt * compflow_auxvar%den
    Res(compflow_liq_flow_eq) = por_sat * tempreal

    if (calculate_derivatives) then
      Jac(compflow_liq_flow_eq,compflow_liq_flow_eq) = &
        volume_over_dt * (compflow_auxvar%den * &
          (compflow_auxvar%dsat_dp * porosity + &
           saturation * compflow_auxvar%dpor_dp) + &
          por_sat * compflow_auxvar%dden_dp)
    endif
  endif

#if defined(SOLUTE_TRANSPORT)
  if (compflow_sol_tran_eq > 0) then
    ! accumulation term units = mol/s
    ! Res[mole/sec] = c [mol/L] * 1000 [L/m^3 liquid]
    !                 sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
    !                 vol[m^3 bulk] / dt[sec]
    tempreal = L_per_m3 * volume_over_dt
    Res(compflow_sol_tran_eq) = compflow_auxvar%conc * tempreal * por_sat

    if (calculate_derivatives) then
      Jac(compflow_sol_tran_eq,compflow_sol_tran_eq) = tempreal * por_sat
      if (compflow_liq_flow_eq > 0) then
        Jac(compflow_sol_tran_eq,compflow_liq_flow_eq) = &
          compflow_auxvar%conc * tempreal * &
          (compflow_auxvar%dsat_dp * porosity + &
          saturation * compflow_auxvar%dpor_dp)
      endif
    endif
  endif
#else
  Res(compflow_sol_tran_eq) = 0.d0
  Jac(compflow_sol_tran_eq,compflow_sol_tran_eq) = 1.d0
#endif

end subroutine CompFlowAccumulation

! ************************************************************************** !

subroutine CompFlowFluxHarmonicPermOnly(compflow_auxvar_up,global_auxvar_up, &
                                     material_auxvar_up, &
                                     compflow_auxvar_dn,global_auxvar_dn, &
                                     material_auxvar_dn, &
                                     area, dist, &
                                     compflow_parameter, &
                                     option,v_darcy,Res, &
                                     Jup,Jdn, &
                                     calculate_derivatives, &
                                     debug_connection)
  !
  ! Computes the internal flux terms for the residual based on harmonic
  ! intrinsic permeability
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Option_module
  use Material_Aux_module
  use Connection_module

  implicit none

  type(compflow_auxvar_type) :: compflow_auxvar_up, compflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(compflow_parameter_type) :: compflow_parameter
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jup(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscReal :: perm_ave_over_dist
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: dkr_dpup, dkr_dpdn
  PetscReal :: dvis_dpup, dvis_dpdn
  PetscReal :: dvis_dtup, dvis_dtdn
  PetscReal :: vis_upwind
  PetscReal :: dkr_over_vis_dpup, dkr_over_vis_dpdn
  PetscReal :: tempreal
  PetscReal :: numerator, denominator
  PetscReal :: dperm_ave_dKup, dperm_ave_dKdn
  PetscReal :: delta_conc
  PetscReal :: dq_dpup, dq_dpdn
  PetscReal :: dq_dKup, dq_dKdn
  PetscReal :: density_ave
  PetscReal :: dden_ave_dpup, dden_ave_dpdn
  PetscReal :: dphi_dpup, dphi_dpdn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal, parameter :: L_per_m3 = 1.d3
  PetscReal :: conc_upwind, dconc_upwind_dup, dconc_upwind_ddn
  PetscReal :: D_hyd_up, D_hyd_dn
  PetscReal :: D_molecular
  PetscReal, parameter :: D_mech_up = 0.d0, D_mech_dn = 0.d0
  PetscReal :: Deff_over_dist
  PetscReal :: dD_mech_up_dpup, dD_mech_dn_dpdn
  PetscReal :: dDeff_over_dist_dpup, dDeff_over_dist_dpdn
  PetscReal :: dD_hyd_up_dpup, dD_hyd_dn_dpdn
  PetscReal :: dDeff_over_dist_dKup, dDeff_over_dist_dKdn

  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  dq_dpup = 0.d0
  dq_dpdn = 0.d0
  dq_dKup = 0.d0
  dq_dKdn = 0.d0
  delta_pressure = 0.d0
  delta_conc = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)

  if (compflow_liq_flow_eq > 0) then

    call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
    call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

    numerator = perm_up * perm_dn
    denominator = dist_up*perm_dn + dist_dn*perm_up
    perm_ave_over_dist = numerator / denominator

    if (compflow_auxvar_up%kr + compflow_auxvar_dn%kr > floweps) then

      density_ave = upweight*compflow_auxvar_up%den + &
                    (1.d0-upweight)*compflow_auxvar_dn%den
      dden_ave_dpup = upweight*compflow_auxvar_up%dden_dp
      dden_ave_dpdn = (1.d0-upweight)*compflow_auxvar_dn%dden_dp
      dgravity_dden_up = upweight*dist_gravity
      dgravity_dden_dn = (1.d0-upweight)*dist_gravity

      gravity_term = density_ave * dist_gravity
      delta_pressure = compflow_auxvar_up%pres - &
                       compflow_auxvar_dn%pres + &
                       gravity_term
      dphi_dpup = 1.d0 + dgravity_dden_up*compflow_auxvar_up%dden_dp
      dphi_dpdn = -1.d0 + dgravity_dden_dn*compflow_auxvar_dn%dden_dp

      dkr_dpup = 0.d0
      dkr_dpdn = 0.d0
      dvis_dpup = 0.d0
      dvis_dpdn = 0.d0
      dvis_dtup = 0.d0
      dvis_dtdn = 0.d0
      dkr_over_vis_dpup = 0.d0
      dkr_over_vis_dpdn = 0.d0
      if (delta_pressure >= 0.d0) then
        kr = compflow_auxvar_up%kr
        vis_upwind = compflow_auxvar_up%vis
        dkr_dpup = compflow_auxvar_up%dkr_dp
        dvis_dpup = compflow_auxvar_up%dvis_dp
        dvis_dtup = compflow_auxvar_up%dvis_dt
        dkr_over_vis_dpup = dkr_dpup/vis_upwind - &
          kr/(vis_upwind*vis_upwind)*dvis_dpup
      else
        kr = compflow_auxvar_dn%kr
        vis_upwind = compflow_auxvar_dn%vis
        dkr_dpdn = compflow_auxvar_dn%dkr_dp
        dvis_dpdn = compflow_auxvar_dn%dvis_dp
        dvis_dtdn = compflow_auxvar_dn%dvis_dt
        dkr_over_vis_dpdn = dkr_dpdn/vis_upwind - &
          kr/(vis_upwind*vis_upwind)*dvis_dpdn
      endif

      if (kr > floweps) then
        ! v_darcy[m/sec] = perm[m^2]/dist[m] * kr[-]/mu[Pa-s] * dP[Pa]
        v_darcy = perm_ave_over_dist * kr / vis_upwind * delta_pressure
        ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy * area
        ! Res[kg/sec]
        Res(compflow_liq_flow_eq) = Res(compflow_liq_flow_eq) + q * density_ave

        if (calculate_derivatives) then
          tempreal = perm_ave_over_dist * area
          if (delta_pressure >= 0.d0) then
            dq_dpup = tempreal * (dkr_over_vis_dpup*delta_pressure + &
              kr/vis_upwind*dphi_dpup)
            dq_dpdn = tempreal * kr/vis_upwind*dphi_dpdn
          else
            dq_dpup = tempreal * kr/vis_upwind*dphi_dpup
            dq_dpdn = tempreal * (dkr_over_vis_dpdn*delta_pressure + &
              kr/vis_upwind*dphi_dpdn)
          endif
          Jup(compflow_liq_flow_eq,compflow_liq_flow_eq) = &
            dq_dpup * density_ave + q * dden_ave_dpup
          Jdn(compflow_liq_flow_eq,compflow_liq_flow_eq) = &
            dq_dpdn * density_ave + q * dden_ave_dpdn
          if (compflow_heat_tran_eq > 0) then
            if (delta_pressure >= 0.d0) then
              tempreal = perm_ave_over_dist * area * &
                (-kr/(vis_upwind*vis_upwind)*dvis_dtup*delta_pressure + &
                 kr/vis_upwind*dgravity_dden_up*compflow_auxvar_up%dden_dt)
              Jup(compflow_liq_flow_eq,compflow_heat_tran_eq) = &
                tempreal * density_ave + q * upweight * &
                compflow_auxvar_up%dden_dt
            else
              tempreal = perm_ave_over_dist * area * &
                (-kr/(vis_upwind*vis_upwind)*dvis_dtdn*delta_pressure + &
                 kr/vis_upwind*dgravity_dden_dn*compflow_auxvar_dn%dden_dt)
              Jdn(compflow_liq_flow_eq,compflow_heat_tran_eq) = &
                tempreal * density_ave + q * (1.d0-upweight) * &
                compflow_auxvar_dn%dden_dt
            endif
          endif
        endif
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"IFF(res,kr,dp): ",8es12.4)') Res(compflow_liq_flow_eq), &
                                                 kr, delta_pressure
    endif
  endif

#if defined(SOLUTE_TRANSPORT)
  if (compflow_sol_tran_eq > 0) then
    if (q >= 0) then
      conc_upwind = compflow_auxvar_up%conc
      dconc_upwind_dup = 1.d0
      dconc_upwind_ddn = 0.d0
    else
      conc_upwind = compflow_auxvar_dn%conc
      dconc_upwind_dup = 0.d0
      dconc_upwind_ddn = 1.d0
    endif
    delta_conc = compflow_auxvar_up%conc - compflow_auxvar_dn%conc
    D_molecular = compflow_parameter%diffusion_coef
    D_hyd_up = D_mech_up + &
               compflow_auxvar_up%effective_porosity * compflow_auxvar_up%sat * &
               material_auxvar_up%tortuosity * D_molecular
    D_hyd_dn = D_mech_dn + &
               compflow_auxvar_dn%effective_porosity * compflow_auxvar_dn%sat * &
               material_auxvar_dn%tortuosity * D_molecular
    numerator = D_hyd_up * D_hyd_dn
    denominator = dist_up*D_hyd_dn + dist_dn*D_hyd_up
    if (denominator == 0.d0) then
      ! turn off diffusion
      numerator = 0.d0
      denominator = 1.d0
    endif
    Deff_over_dist = numerator / denominator
    ! Res[mol/sec]
    Res(compflow_sol_tran_eq) = Res(compflow_sol_tran_eq) + &
      (q * conc_upwind + & ! advection
       area * Deff_over_dist * delta_conc) * & ! hydrodynamic dispersion
      L_per_m3
    if (calculate_derivatives) then
      Jup(compflow_sol_tran_eq,compflow_sol_tran_eq) = &
        (q * dconc_upwind_dup + &
         area * Deff_over_dist * 1.d0) * &
        L_per_m3
      Jdn(compflow_sol_tran_eq,compflow_sol_tran_eq) = &
        (q * dconc_upwind_ddn + &
         area * Deff_over_dist * (-1.d0)) * &
        L_per_m3
      if (compflow_liq_flow_eq > 0) then
        if (D_mech_up > 0.d0) then
          ! implement derivatives here and in boundary flux
          option%io_buffer = 'Implement derivatives for D_mech_up'
          call PrintErrMsg(option)
        endif
        dD_mech_up_dpup = 0.d0
        dD_hyd_up_dpup = dD_mech_up_dpup + &
                         (compflow_auxvar_up%dsat_dp * &
                          compflow_auxvar_up%effective_porosity + &
                          compflow_auxvar_up%sat * compflow_auxvar_up%dpor_dp) * &
                         material_auxvar_up%tortuosity * D_molecular
        dD_mech_dn_dpdn = 0.d0
        dD_hyd_dn_dpdn = dD_mech_dn_dpdn + &
                         (compflow_auxvar_dn%dsat_dp * &
                          compflow_auxvar_dn%effective_porosity + &
                          compflow_auxvar_dn%sat * compflow_auxvar_dn%dpor_dp) * &
                         material_auxvar_dn%tortuosity * D_molecular
        tempreal = denominator * denominator
        dDeff_over_dist_dpup = dD_hyd_up_dpup * D_hyd_dn * D_hyd_dn * &
                               dist_up / tempreal
        dDeff_over_dist_dpdn = dD_hyd_dn_dpdn * D_hyd_up * D_hyd_up * &
                               dist_dn / tempreal
        Jup(compflow_sol_tran_eq,compflow_liq_flow_eq) = &
          (dq_dpup * conc_upwind + &
           area * dDeff_over_dist_dpup * delta_conc) * &
          L_per_m3
        Jdn(compflow_sol_tran_eq,compflow_liq_flow_eq) = &
          (dq_dpdn * conc_upwind + &
           area * dDeff_over_dist_dpdn * delta_conc) * &
          L_per_m3
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"IFS(res,q,dc): ",8es12.4)') Res(compflow_sol_tran_eq), q, &
                                                delta_conc
    endif
  endif
#endif

end subroutine CompFlowFluxHarmonicPermOnly

! ************************************************************************** !

subroutine CompFlowBCFluxHarmonicPermOnly(ibndtype,auxvar_mapping,auxvars, &
                                       compflow_auxvar_up,global_auxvar_up, &
                                       compflow_auxvar_dn,global_auxvar_dn, &
                                       material_auxvar_dn, &
                                       area,dist, &
                                       compflow_parameter, &
                                       option,v_darcy,Res,Jdn, &
                                       calculate_derivatives, &
                                       debug_connection)
  !
  ! Computes the boundary flux terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Option_module
  use Material_Aux_module
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1)
  PetscInt :: auxvar_mapping(COMPFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(compflow_auxvar_type) :: compflow_auxvar_up, compflow_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(compflow_parameter_type) :: compflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscInt :: bc_type
  PetscInt :: idof
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: dkr_dp, dvis_dp
  PetscReal :: vis_upwind
  PetscReal :: dkr_over_vis_dp
  PetscReal :: perm_dn
  PetscReal :: density_ave
  PetscReal :: dden_ave_dpdn
  PetscReal :: boundary_pressure
  PetscReal :: ddelta_pressure_dpdn
  PetscBool :: derivative_toggle
  PetscReal :: dperm_dK
  PetscReal :: delta_conc
  PetscReal :: dq_dpdn
  PetscReal :: dq_dKdn
  PetscReal, parameter :: L_per_m3 = 1.d3
  PetscReal :: conc_upwind, dconc_upwind_ddn
  PetscReal :: D_hyd_dn
  PetscReal :: D_molecular = 0.d0
  PetscReal, parameter :: D_mech_dn = 0.d0
  PetscReal :: Deff_over_dist
  PetscReal :: dispersion_scale
  PetscReal :: dD_mech_dn_dpdn
  PetscReal :: dDeff_over_dist_dpdn
  PetscReal :: dDeff_over_dist_dKdn
  PetscReal :: dD_hyd_dn_dpdn
  PetscReal :: tempreal

  Res = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  dq_dpdn = 0.d0
  dq_dKdn = 0.d0
  dperm_dK = 0.d0
  delta_pressure = 0.d0
  delta_conc = 0.d0

  if (compflow_liq_flow_eq > 0) then
    call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

    kr = 0.d0
    density_ave = compflow_auxvar_dn%den
    dden_ave_dpdn = compflow_auxvar_dn%dden_dp
    bc_type = ibndtype(auxvar_mapping(COMPFLOW_COND_WATER_INDEX))
    derivative_toggle = PETSC_TRUE
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
           HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
           PONDED_WATER_BC)
        if (compflow_auxvar_up%kr + compflow_auxvar_dn%kr > floweps) then
          ! dist(0) = scalar - magnitude of distance
          ! gravity = vector(3)
          ! dist(1:3) = vector(3) - unit vector
          dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

          boundary_pressure = compflow_auxvar_up%pres
          gravity_term = compflow_auxvar_dn%den * dist_gravity
          delta_pressure = boundary_pressure - &
                          compflow_auxvar_dn%pres + &
                          gravity_term
          ddelta_pressure_dpdn = -1.d0
          select case(bc_type)
            case(DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
                 HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
                 PONDED_WATER_BC)
                  ! flow in         ! boundary cell is <= pref
              if (delta_pressure > 0.d0 .and. &
                  compflow_auxvar_up%pres - &
                    option%flow%reference_pressure < eps) then
                delta_pressure = 0.d0
                ddelta_pressure_dpdn = 0.d0
              endif
          end select
          dkr_over_vis_dp = 0.d0
          if (delta_pressure >= 0.d0) then
            kr = compflow_auxvar_up%kr
            vis_upwind = compflow_auxvar_up%vis
            dkr_dp = 0.d0
            dvis_dp = compflow_auxvar_up%dvis_dp
            dkr_over_vis_dp = -kr/(vis_upwind*vis_upwind)*dvis_dp
            density_ave = compflow_auxvar_up%den
            dden_ave_dpdn = 0.d0
          else
            kr = compflow_auxvar_dn%kr
            vis_upwind = compflow_auxvar_dn%vis
            dkr_dp = compflow_auxvar_dn%dkr_dp
            dvis_dp = compflow_auxvar_dn%dvis_dp
            dkr_over_vis_dp = dkr_dp/vis_upwind - &
              kr/(vis_upwind*vis_upwind)*dvis_dp
            density_ave = compflow_auxvar_dn%den
            dden_ave_dpdn = compflow_auxvar_dn%dden_dp
          endif

          if (bc_type == DIRICHLET_CONDUCTANCE_BC .or. &
              bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
            idof = auxvar_mapping(COMPFLOW_COND_WATER_AUX_INDEX)
            perm_ave_over_dist = auxvars(idof) / vis_upwind
          else
            perm_ave_over_dist = perm_dn / (dist(0) * vis_upwind)
            dperm_dK = 1.d0 / (dist(0) * vis_upwind)
          endif

          ! v_darcy[m/sec] = perm[m^2]/dist[m] * kr[-]/mu[Pa-s] * dP[Pa]
          v_darcy = perm_ave_over_dist * kr / vis_upwind * delta_pressure

          if (bc_type == PONDED_WATER_BC) then
            idof = auxvar_mapping(COMPFLOW_COND_WATER_INDEX)
            tempreal = auxvars(auxvar_mapping(COMPFLOW_COND_WATER_INDEX)) / &
                       option%flow_dt
            if (v_darcy > tempreal) then
              v_darcy = tempreal
              ddelta_pressure_dpdn = 0.d0
              dkr_dp = 0.d0
              dkr_over_vis_dp = 0.d0
            endif
          endif
        endif

      case(NEUMANN_BC)
        idof = auxvar_mapping(COMPFLOW_COND_WATER_INDEX)
        if (dabs(auxvars(idof)) > floweps) then
          v_darcy = auxvars(idof)
        endif
        derivative_toggle = PETSC_FALSE
      case default
        option%io_buffer = &
          'Boundary condition type (' // trim(StringWrite(bc_type)) // &
          ') not recognized in CompFlowBCFlux phase loop.'
        call PrintErrMsg(option)
    end select
    if (dabs(v_darcy) > 0.d0 .or. kr > 0.d0) then
      ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy * area
      ! Res[kg/sec]
      Res(compflow_liq_flow_eq) = Res(compflow_liq_flow_eq) + q * density_ave
      if (calculate_derivatives .and. derivative_toggle) then
        ! derivative toggle takes care of the NEUMANN side
        dq_dpdn = area * perm_ave_over_dist * &
                  (dkr_over_vis_dp*delta_pressure + &
                   kr/vis_upwind*ddelta_pressure_dpdn)
        Jdn(compflow_liq_flow_eq,compflow_liq_flow_eq) = &
          dq_dpdn * density_ave + q * dden_ave_dpdn
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"BFF(res,kr,dp): ",8es12.4)') Res(compflow_liq_flow_eq), &
                                                 kr, delta_pressure
    endif
  endif

#if defined(SOLUTE_TRANSPORT)
  if (compflow_sol_tran_eq > 0) then
    if (q >= 0) then
      conc_upwind = compflow_auxvar_up%conc
      dconc_upwind_ddn = 0.d0
      dispersion_scale = 1.d0
    else
      conc_upwind = compflow_auxvar_dn%conc
      dconc_upwind_ddn = 1.d0
      dispersion_scale = 0.d0  ! creates zero-gradient outlet bc
    endif
    delta_conc = compflow_auxvar_up%conc - compflow_auxvar_dn%conc
    D_molecular = compflow_parameter%diffusion_coef
    D_hyd_dn = D_mech_dn + &
               compflow_auxvar_dn%effective_porosity * compflow_auxvar_dn%sat * &
               material_auxvar_dn%tortuosity * D_molecular
    Deff_over_dist = dispersion_scale * D_hyd_dn / dist(0)
    ! Res[mol/sec]
    Res(compflow_sol_tran_eq) = Res(compflow_sol_tran_eq) + &
      (q * conc_upwind + & ! advection
       area * Deff_over_dist * delta_conc) * & ! hydrodynamic dispersion
      L_per_m3
    if (calculate_derivatives) then
      Jdn(compflow_sol_tran_eq,compflow_sol_tran_eq) = &
        (q * dconc_upwind_ddn + &
         area * Deff_over_dist * (-1.d0)) * &
        L_per_m3
      if (compflow_liq_flow_eq > 0) then
        dD_mech_dn_dpdn = 0.d0
        dD_hyd_dn_dpdn = dD_mech_dn_dpdn + &
                         (compflow_auxvar_dn%dsat_dp * &
                          compflow_auxvar_dn%effective_porosity + &
                          compflow_auxvar_dn%sat * compflow_auxvar_dn%dpor_dp) * &
                         material_auxvar_dn%tortuosity * D_molecular
        dDeff_over_dist_dpdn = dD_hyd_dn_dpdn / dist(0)
        Jdn(compflow_sol_tran_eq,compflow_liq_flow_eq) = &
          (dq_dpdn * conc_upwind + &
           area * dDeff_over_dist_dpdn * (-1.d0)) * &
          L_per_m3
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"BFS(res,q,dc): ",8es12.4)') Res(compflow_sol_tran_eq), q, &
                                                delta_conc
    endif
  endif
#endif

end subroutine CompFlowBCFluxHarmonicPermOnly

! ************************************************************************** !

subroutine CompFlowSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                        flow_src_sink_type, &
                        compflow_auxvar,global_auxvar,material_auxvar, &
                        ss_flow_vol_flux,Res,Jdn, &
                        calculate_derivatives)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Option_module
  use Material_Aux_module
  use EOS_Water_module
  use EOS_Gas_module

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  type(compflow_auxvar_type) :: compflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: tempreal
  PetscReal :: qsrc_m3, qsrc_L
  PetscReal, parameter :: L_per_m3 = 1.d3
  PetscInt :: dof_index, index_

  Res = 0.d0
  Jdn = 0.d0

  if (compflow_liq_flow_eq > 0) then
    dof_index = flow_src_sink_mapping(COMPFLOW_COND_WATER_INDEX)
    qsrc_m3 = 0.d0
    tempreal = flow_aux_real_var(dof_index)
    select case(flow_src_sink_type(dof_index))
      case(VOLUMETRIC_RATE_SS)
        qsrc_m3 = tempreal
      case(SCALED_VOLUMETRIC_RATE_SS)
        index_ = flow_src_sink_mapping(COMPFLOW_COND_WATER_AUX_INDEX)
        qsrc_m3 = tempreal*flow_aux_real_var(index_)
      case default
        option%io_buffer = 'src_sink_type not supported in CompFlowSrcSink'
        call PrintErrMsg(option)
    end select
    ss_flow_vol_flux = qsrc_m3
    ! Res[kg/sec]
    Res(compflow_liq_flow_eq) = qsrc_m3 * compflow_auxvar%den

    if (calculate_derivatives) then
    endif
  endif

  if (compflow_sol_tran_eq > 0) then
    if (qsrc_m3 >= 0.d0) then
      ! injection
      dof_index = flow_src_sink_mapping(COMPFLOW_COND_SOLUTE_INDEX)
      tempreal = flow_aux_real_var(dof_index)
    else
      ! extraction
      tempreal = compflow_auxvar%conc
    endif
    qsrc_L = qsrc_m3 * L_per_m3
    Res(compflow_sol_tran_eq) = qsrc_L * tempreal
    if (calculate_derivatives) then
      if (qsrc_m3 < 0.d0) then
       Jdn(compflow_sol_tran_eq,compflow_sol_tran_eq) = qsrc_L
      endif
    endif
  endif

end subroutine CompFlowSrcSink

! ************************************************************************** !

subroutine CompFlowAccumDerivative(compflow_auxvar,global_auxvar, &
                                material_auxvar, &
                                material_auxvar_pert, &
                                option,Res,Jac)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(compflow_auxvar_type) :: compflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  type(option_type) :: option
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jac(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)

  PetscReal :: res_pert(COMPFLOW_MAX_DOF)
  PetscReal :: Jdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscReal :: dJdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  call CompFlowAccumulation(compflow_auxvar(ZERO_INTEGER), &
                         global_auxvar, &
                         material_auxvar, &
                         option,Res,Jac, &
                         .not.compflow_numerical_derivatives)

  if (compflow_numerical_derivatives) then
    do idof = 1, option%nflowdof
      call CompFlowAccumulation(compflow_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar, &
                             option,res_pert,Jdum, &
                             PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq))/compflow_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine CompFlowAccumDerivative

! ************************************************************************** !

subroutine XXFluxDerivative(compflow_auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            material_auxvar_pert_up, &
                            compflow_auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            material_auxvar_pert_dn, &
                            area, dist,compflow_parameter,option,v_darcy, &
                            Res,Jup,Jdn, &
                            debug_connection)
  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(compflow_auxvar_type) :: compflow_auxvar_up(0:), compflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_pert_up(:)
  type(material_auxvar_type) :: material_auxvar_pert_dn(:)
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(compflow_parameter_type) :: compflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jup(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(COMPFLOW_MAX_DOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscReal :: dJdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  Jup = 0.d0
  Jdn = 0.d0

  call XXFlux(compflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              compflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist,compflow_parameter,option,v_darcy, &
              Res,Jup,Jdn, &
              .not.compflow_numerical_derivatives, &
              debug_connection)

  if (compflow_numerical_derivatives) then
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(compflow_auxvar_up(idof),global_auxvar_up, &
                  material_auxvar_up, &
                  compflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,compflow_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,PETSC_FALSE,debug_connection)
      do ieq = 1, option%nflowdof
        Jup(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        compflow_auxvar_up(idof)%pert
      enddo
    enddo
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(compflow_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  compflow_auxvar_dn(idof),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,compflow_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,PETSC_FALSE,debug_connection)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        compflow_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine XXFluxDerivative

! ************************************************************************** !

subroutine XXBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                              compflow_auxvar_up, &
                              global_auxvar_up, &
                              compflow_auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              material_auxvar_pert_dn, &
                              area,dist,compflow_parameter,option,v_darcy, &
                              Res,Jdn,debug_connection)
  !
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:1)
  PetscInt :: auxvar_mapping(COMPFLOW_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(compflow_auxvar_type) :: compflow_auxvar_up, compflow_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_pert_dn(:)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(compflow_parameter_type) :: compflow_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jdn(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(COMPFLOW_MAX_DOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscReal :: dJdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  Jdn = 0.d0

  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                compflow_auxvar_up,global_auxvar_up, &
                compflow_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                material_auxvar_dn, &
                area,dist,compflow_parameter,option,v_darcy, &
                Res,Jdn,.not.compflow_numerical_derivatives, &
                debug_connection)

  if (compflow_numerical_derivatives) then
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                    compflow_auxvar_up,global_auxvar_up, &
                    compflow_auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    area,dist,compflow_parameter,option,v_darcy_dum, &
                    res_pert,Jdum,PETSC_FALSE,debug_connection)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        compflow_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine XXBCFluxDerivative

! ************************************************************************** !

subroutine CompFlowSrcSinkDerivative(option,flow_aux_real_var, &
                                  flow_src_sink_mapping, &
                                  flow_src_sink_type, &
                                  compflow_auxvar,global_auxvar, &
                                  material_auxvar, &
                                  material_auxvar_pert, &
                                  ss_flow_vol_flux, &
                                  Res,Jac)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  type(compflow_auxvar_type) :: compflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(COMPFLOW_MAX_DOF)
  PetscReal :: Jac(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)

  PetscReal :: res_pert(COMPFLOW_MAX_DOF)
  PetscReal :: dummy_real
  PetscReal :: Jdum(COMPFLOW_MAX_DOF,COMPFLOW_MAX_DOF)
  PetscInt :: idof, ieq

  Jac = 0.d0
  ! unperturbed compflow_auxvars value
  call CompFlowSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                    flow_src_sink_type, &
                    compflow_auxvar(ZERO_INTEGER),global_auxvar, &
                    material_auxvar, &
                    ss_flow_vol_flux, &
                    Res,Jac,.not.compflow_numerical_derivatives)

  if (compflow_numerical_derivatives) then
    ! perturbed compflow_auxvars values
    do idof = 1, option%nflowdof
      call CompFlowSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                        flow_src_sink_type, &
                        compflow_auxvar(idof),global_auxvar, &
                        material_auxvar, &
                        dummy_real, &
                        res_pert,Jdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
      Jac(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                      compflow_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine CompFlowSrcSinkDerivative

end module Compositional_Flow_Common_module
