module Immiscible_Common_module

  use Immiscible_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use petscsys

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#define SOLUTE_TRANSPORT

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 0.d0

  procedure(FluxDummy), pointer :: XXFlux => null()
  procedure(BCFluxDummy), pointer :: XXBCFlux => null()

  interface
    subroutine FluxDummy(immiscible_auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         immiscible_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area, dist, &
                         immiscible_parameter, &
                         option,v_darcy,Res,Jup,Jdn, &
                         dResdparamup,dResdparamdn, &
                         calculate_derivatives, &
                         debug_connection)
      use Immiscible_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(immiscible_parameter_type) :: immiscible_parameter
      PetscReal :: Res(IMMIS_MAX_DOF)
      PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscReal :: dResdparamup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscBool :: calculate_derivatives
      PetscBool :: debug_connection
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                           immiscible_auxvar_up,global_auxvar_up, &
                           immiscible_auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           area,dist, &
                           immiscible_parameter, &
                           option,v_darcy,Res,Jdn, &
                           dResdparamdn, &
                           calculate_derivatives, &
                           debug_connection)
      use Immiscible_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(option_type) :: option
      PetscInt :: ibndtype(1)
      PetscInt :: auxvar_mapping(IMMIS_MAX_INDEX)
      PetscReal :: auxvars(:) ! from aux_real_var array
      type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_dn
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(immiscible_parameter_type) :: immiscible_parameter
      PetscReal :: v_darcy
      PetscReal :: Res(IMMIS_MAX_DOF)
      PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscBool :: calculate_derivatives
      PetscBool :: debug_connection
    end subroutine BCFluxDummy
  end interface

  public :: ImmiscibleAccumulation, &
            ImmiscibleFluxHarmonicPermOnly, &
            ImmiscibleBCFluxHarmonicPermOnly, &
            ImmiscibleAccumDerivative, &
            XXFluxDerivative, &
            XXBCFluxDerivative, &
            ImmiscibleSrcSinkDerivative

  public :: XXFlux, &
            XXBCFlux

contains

! ************************************************************************** !

subroutine ImmiscibleAccumulation(immiscible_auxvar,global_auxvar,material_auxvar, &
                             option,Res,Jac,dResdparam,calculate_derivatives)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jac(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparam(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: porosity
  PetscReal :: saturation
  PetscReal :: por_sat
  PetscReal :: volume_over_dt
  PetscReal :: tempreal
  PetscReal, parameter :: L_per_m3 = 1.d3

  Res = 0.d0
  Jac = 0.d0
  dResdparam = 0.d0

    ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use immiscible_auxvar%effective porosity here as it enables numerical
  ! derivatives to be employed
  porosity = immiscible_auxvar%effective_porosity
  saturation = immiscible_auxvar%sat
  por_sat = saturation * porosity

  if (immis_liq_flow_eq > 0) then
    ! accumulation term units = m^3 liquid/s
    ! Res[m^3 liquid/sec] = sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
    !                       vol[m^3 bulk] / dt[sec]
    tempreal = volume_over_dt * immis_density_kmol
    Res(immis_liq_flow_eq) = por_sat * tempreal

    if (calculate_derivatives) then
      Jac(immis_liq_flow_eq,immis_liq_flow_eq) = &
        tempreal * (immiscible_auxvar%dsat_dp * porosity + &
                    saturation * immiscible_auxvar%dpor_dp)
      if (immis_calc_adjoint) then
        select case(immis_adjoint_parameter)
          case(IMMIS_ADJOINT_POROSITY)
            dResdparam(immis_liq_flow_eq,1) = &
              tempreal * immiscible_auxvar%sat
        end select
      endif
    endif
  endif

#if defined(SOLUTE_TRANSPORT)
  if (immis_sol_tran_eq > 0) then
    ! accumulation term units = mol/s
    ! Res[mole/sec] = c [mol/L] * 1000 [L/m^3 liquid]
    !                 sat[m^3 liquid/m^3 void] * por[m^3 void/m^3 bulk] *
    !                 vol[m^3 bulk] / dt[sec]
    tempreal = L_per_m3 * volume_over_dt
    Res(immis_sol_tran_eq) = immiscible_auxvar%conc * tempreal * por_sat

    if (calculate_derivatives) then
      Jac(immis_sol_tran_eq,immis_sol_tran_eq) = tempreal * por_sat
      if (immis_liq_flow_eq > 0) then
        Jac(immis_sol_tran_eq,immis_liq_flow_eq) = &
          immiscible_auxvar%conc * tempreal * &
          (immiscible_auxvar%dsat_dp * porosity + &
          saturation * immiscible_auxvar%dpor_dp)
      endif
      if (immis_calc_adjoint) then
        select case (immis_adjoint_parameter)
          case(IMMIS_ADJOINT_POROSITY)
            dResdparam(immis_sol_tran_eq,1) = &
              immiscible_auxvar%conc * tempreal * immiscible_auxvar%sat
        end select
      endif
    endif
  endif
#else
  Res(immis_sol_tran_eq) = 0.d0
  Jac(immis_sol_tran_eq,immis_sol_tran_eq) = 1.d0
#endif

end subroutine ImmiscibleAccumulation

! ************************************************************************** !

subroutine ImmiscibleFluxHarmonicPermOnly(immiscible_auxvar_up,global_auxvar_up, &
                                     material_auxvar_up, &
                                     immiscible_auxvar_dn,global_auxvar_dn, &
                                     material_auxvar_dn, &
                                     area, dist, &
                                     immiscible_parameter, &
                                     option,v_darcy,Res, &
                                     Jup,Jdn, &
                                     dResdparamup,dResdparamdn, &
                                     calculate_derivatives, &
                                     debug_connection)
  !
  ! Computes the internal flux terms for the residual based on harmonic
  ! intrinsic permeability
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module
  use Connection_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscReal :: perm_ave_over_dist_visc
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: dkr_dpup, dkr_dpdn
  PetscReal :: tempreal
  PetscReal :: numerator, denominator
  PetscReal :: dperm_ave_dKup, dperm_ave_dKdn
  PetscReal :: delta_conc
  PetscReal :: dq_dpup, dq_dpdn
  PetscReal :: dq_dKup, dq_dKdn
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
  dResdparamup = 0.d0
  dResdparamdn = 0.d0
  delta_pressure = 0.d0
  delta_conc = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)

  if (immis_liq_flow_eq > 0) then

    call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
    call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

    numerator = perm_up * perm_dn
    denominator = dist_up*perm_dn + dist_dn*perm_up
    perm_ave_over_dist_visc = numerator / (denominator * immis_viscosity)

    if (immiscible_auxvar_up%kr + immiscible_auxvar_dn%kr > floweps) then

      gravity_term = immis_density_kg * dist_gravity
      delta_pressure = immiscible_auxvar_up%pres - &
                       immiscible_auxvar_dn%pres + &
                       gravity_term
      dkr_dpup = 0.d0
      dkr_dpdn = 0.d0
      if (delta_pressure >= 0.d0) then
        kr = immiscible_auxvar_up%kr
        dkr_dpup = immiscible_auxvar_up%dkr_dp
      else
        kr = immiscible_auxvar_dn%kr
        dkr_dpdn = immiscible_auxvar_dn%dkr_dp
      endif

      if (kr > floweps) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy = perm_ave_over_dist_visc * kr * delta_pressure
        ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy * area
        ! Res[m^3 liquid/sec]
        Res(immis_liq_flow_eq) = Res(immis_liq_flow_eq) + &
                                                 q * immis_density_kmol

        if (calculate_derivatives) then
          tempreal = perm_ave_over_dist_visc * area
          dq_dpup = tempreal * (delta_pressure * dkr_dpup + kr)
          dq_dpdn = tempreal * (delta_pressure * dkr_dpdn - kr)
          Jup(immis_liq_flow_eq,immis_liq_flow_eq) = &
                                              dq_dpup * immis_density_kmol
          Jdn(immis_liq_flow_eq,immis_liq_flow_eq) = &
                                              dq_dpdn * immis_density_kmol
          if (immis_calc_adjoint) then
            select case(immis_adjoint_parameter)
              case(IMMIS_ADJOINT_PERMEABILITY)
                tempreal = denominator * denominator * immis_viscosity
                dperm_ave_dKup = perm_dn * perm_dn * dist_up / tempreal
                dperm_ave_dKdn = perm_up * perm_up * dist_dn / tempreal
                tempreal = kr * delta_pressure * area
                dq_dKup = tempreal * dperm_ave_dKup
                dq_dKdn = tempreal * dperm_ave_dKdn
                dResdparamup(immis_liq_flow_eq,1) = dq_dKup *immis_density_kmol
                dResdparamdn(immis_liq_flow_eq,1) = dq_dKdn *immis_density_kmol
              case(IMMIS_ADJOINT_POROSITY)
                ! not applicable
            end select
          endif
        endif
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"IFF(res,kr,dp): ",8es12.4)') Res(immis_liq_flow_eq), &
                                                 kr, delta_pressure
    endif
  endif

#if defined(SOLUTE_TRANSPORT)
  if (immis_sol_tran_eq > 0) then
    if (q >= 0) then
      conc_upwind = immiscible_auxvar_up%conc
      dconc_upwind_dup = 1.d0
      dconc_upwind_ddn = 0.d0
    else
      conc_upwind = immiscible_auxvar_dn%conc
      dconc_upwind_dup = 0.d0
      dconc_upwind_ddn = 1.d0
    endif
    delta_conc = immiscible_auxvar_up%conc - immiscible_auxvar_dn%conc
    D_molecular = immiscible_parameter%diffusion_coef
    D_hyd_up = D_mech_up + &
               immiscible_auxvar_up%effective_porosity * immiscible_auxvar_up%sat * &
               material_auxvar_up%tortuosity * D_molecular
    D_hyd_dn = D_mech_dn + &
               immiscible_auxvar_dn%effective_porosity * immiscible_auxvar_dn%sat * &
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
    Res(immis_sol_tran_eq) = Res(immis_sol_tran_eq) + &
      (q * conc_upwind + & ! advection
       area * Deff_over_dist * delta_conc) * & ! hydrodynamic dispersion
      L_per_m3
    if (calculate_derivatives) then
      Jup(immis_sol_tran_eq,immis_sol_tran_eq) = &
        (q * dconc_upwind_dup + &
         area * Deff_over_dist * 1.d0) * &
        L_per_m3
      Jdn(immis_sol_tran_eq,immis_sol_tran_eq) = &
        (q * dconc_upwind_ddn + &
         area * Deff_over_dist * (-1.d0)) * &
        L_per_m3
      if (immis_liq_flow_eq > 0) then
        if (D_mech_up > 0.d0) then
          ! implement derivatives here and in boundary flux
          option%io_buffer = 'Implement derivatives for D_mech_up'
          call PrintErrMsg(option)
        endif
        dD_mech_up_dpup = 0.d0
        dD_hyd_up_dpup = dD_mech_up_dpup + &
                         (immiscible_auxvar_up%dsat_dp * &
                          immiscible_auxvar_up%effective_porosity + &
                          immiscible_auxvar_up%sat * immiscible_auxvar_up%dpor_dp) * &
                         material_auxvar_up%tortuosity * D_molecular
        dD_mech_dn_dpdn = 0.d0
        dD_hyd_dn_dpdn = dD_mech_dn_dpdn + &
                         (immiscible_auxvar_dn%dsat_dp * &
                          immiscible_auxvar_dn%effective_porosity + &
                          immiscible_auxvar_dn%sat * immiscible_auxvar_dn%dpor_dp) * &
                         material_auxvar_dn%tortuosity * D_molecular
        tempreal = denominator * denominator
        dDeff_over_dist_dpup = dD_hyd_up_dpup * D_hyd_dn * D_hyd_dn * &
                               dist_up / tempreal
        dDeff_over_dist_dpdn = dD_hyd_dn_dpdn * D_hyd_up * D_hyd_up * &
                               dist_dn / tempreal
        Jup(immis_sol_tran_eq,immis_liq_flow_eq) = &
          (dq_dpup * conc_upwind + &
           area * dDeff_over_dist_dpup * delta_conc) * &
          L_per_m3
        Jdn(immis_sol_tran_eq,immis_liq_flow_eq) = &
          (dq_dpdn * conc_upwind + &
           area * dDeff_over_dist_dpdn * delta_conc) * &
          L_per_m3
      endif
      if (immis_calc_adjoint) then
        select case(immis_adjoint_parameter)
          case(IMMIS_ADJOINT_PERMEABILITY)
            dDeff_over_dist_dKup = 0.d0 ! for effect of dispersivity on Deff
            dDeff_over_dist_dKdn = 0.d0
            dResdparamup(immis_sol_tran_eq,1) = &
              (dq_dKup * conc_upwind + &
              area * dDeff_over_dist_dKup * delta_conc) * &
              L_per_m3
            dResdparamdn(immis_sol_tran_eq,1) = &
              (dq_dKdn * conc_upwind + &
              area * dDeff_over_dist_dKdn * delta_conc) * &
              L_per_m3
          case(IMMIS_ADJOINT_POROSITY)
            tempreal = D_mech_up + &
                       immiscible_auxvar_up%sat * material_auxvar_up%tortuosity * &
                       D_molecular
            numerator = tempreal * D_hyd_dn
            denominator = dist_up*D_hyd_dn + dist_dn*tempreal
            if (denominator == 0.d0) then
              ! turn off diffusion
              numerator = 0.d0
              denominator = 1.d0
            endif
            Deff_over_dist = numerator / denominator
            dResdparamup(immis_sol_tran_eq,1) = area * Deff_over_dist * &
                                                delta_conc * L_per_m3

            tempreal = D_mech_dn + &
                       immiscible_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                       D_molecular
            numerator = D_hyd_up * tempreal
            denominator = dist_up*tempreal + dist_dn*D_hyd_up
            if (denominator == 0.d0) then
              ! turn off diffusion
              numerator = 0.d0
              denominator = 1.d0
            endif
            Deff_over_dist = numerator / denominator
            dResdparamdn(immis_sol_tran_eq,1) = area * Deff_over_dist * &
                                                delta_conc * L_per_m3
        end select
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"IFS(res,q,dc): ",8es12.4)') Res(immis_sol_tran_eq), q, &
                                                delta_conc
    endif
  endif
#endif

end subroutine ImmiscibleFluxHarmonicPermOnly

! ************************************************************************** !

subroutine ImmiscibleBCFluxHarmonicPermOnly(ibndtype,auxvar_mapping,auxvars, &
                                       immiscible_auxvar_up,global_auxvar_up, &
                                       immiscible_auxvar_dn,global_auxvar_dn, &
                                       material_auxvar_dn, &
                                       area,dist, &
                                       immiscible_parameter, &
                                       option,v_darcy,Res,Jdn, &
                                       dResdparamdn, &
                                       calculate_derivatives, &
                                       debug_connection)
  !
  ! Computes the boundary flux terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1)
  PetscInt :: auxvar_mapping(IMMIS_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscInt :: bc_type
  PetscInt :: idof
  PetscReal :: perm_ave_over_dist_visc
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: dkr_dp
  PetscReal :: perm_dn
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
  dResdparamdn = 0.d0
  delta_pressure = 0.d0
  delta_conc = 0.d0

  if (immis_liq_flow_eq > 0) then
    call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

    kr = 0.d0
    bc_type = ibndtype(auxvar_mapping(IMMIS_COND_WATER_INDEX))
    derivative_toggle = PETSC_TRUE
    select case(bc_type)
      ! figure out the direction of flow
      case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
           HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        if (immiscible_auxvar_up%kr + immiscible_auxvar_dn%kr > floweps) then
          ! dist(0) = scalar - magnitude of distance
          ! gravity = vector(3)
          ! dist(1:3) = vector(3) - unit vector
          dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

          if (bc_type == DIRICHLET_CONDUCTANCE_BC .or. &
              bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
            idof = auxvar_mapping(IMMIS_COND_WATER_AUX_INDEX)
            perm_ave_over_dist_visc = auxvars(idof) / immis_viscosity
          else
            perm_ave_over_dist_visc = perm_dn / (dist(0) * immis_viscosity)
            dperm_dK = 1.d0 / (dist(0) * immis_viscosity)
          endif

          boundary_pressure = immiscible_auxvar_up%pres
          gravity_term = immis_density_kg * dist_gravity
          delta_pressure = boundary_pressure - &
                          immiscible_auxvar_dn%pres + &
                          gravity_term
          ddelta_pressure_dpdn = -1.d0
          select case(bc_type)
            case(DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
                 HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
                  ! flow in         ! boundary cell is <= pref
              if (delta_pressure > 0.d0 .and. &
                  immiscible_auxvar_up%pres - &
                    option%flow%reference_pressure < eps) then
                delta_pressure = 0.d0
                ddelta_pressure_dpdn = 0.d0
              endif
          end select
          if (delta_pressure >= 0.d0) then
            kr = immiscible_auxvar_up%kr
            dkr_dp = 0.d0
          else
            kr = immiscible_auxvar_dn%kr
            dkr_dp = immiscible_auxvar_dn%dkr_dp
          endif

          ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
          !                    dP[Pa]]
          v_darcy = perm_ave_over_dist_visc * kr * delta_pressure
        endif

      case(NEUMANN_BC)
        idof = auxvar_mapping(IMMIS_COND_WATER_INDEX)
        if (dabs(auxvars(idof)) > floweps) then
          v_darcy = auxvars(idof)
        endif
        derivative_toggle = PETSC_FALSE
      case default
        option%io_buffer = &
          'Boundary condition type (' // trim(StringWrite(bc_type)) // &
          ') not recognized in ImmiscibleBCFlux phase loop.'
        call PrintErrMsg(option)
    end select
    if (dabs(v_darcy) > 0.d0 .or. kr > 0.d0) then
      ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
      q = v_darcy * area
      ! Res[m^3 liquid/sec]
      Res(immis_liq_flow_eq) = Res(immis_liq_flow_eq) + q * immis_density_kmol
      if (calculate_derivatives .and. derivative_toggle) then
        ! derivative toggle takes care of the NEUMANN side
        dq_dpdn = area * perm_ave_over_dist_visc * &
                  (dkr_dp * delta_pressure + ddelta_pressure_dpdn * kr)
        Jdn(immis_liq_flow_eq,immis_liq_flow_eq) = dq_dpdn * immis_density_kmol
        if (immis_calc_adjoint) then
          select case(immis_adjoint_parameter)
            case(IMMIS_ADJOINT_PERMEABILITY)
              dq_dKdn = area * kr * delta_pressure * dperm_dK
              dResdparamdn(immis_liq_flow_eq,1) = dq_dKdn * immis_density_kmol
            case(IMMIS_ADJOINT_POROSITY)
              ! not applicable
          end select
        endif
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"BFF(res,kr,dp): ",8es12.4)') Res(immis_liq_flow_eq), &
                                                 kr, delta_pressure
    endif
  endif

#if defined(SOLUTE_TRANSPORT)
  if (immis_sol_tran_eq > 0) then
    if (q >= 0) then
      conc_upwind = immiscible_auxvar_up%conc
      dconc_upwind_ddn = 0.d0
      dispersion_scale = 1.d0
    else
      conc_upwind = immiscible_auxvar_dn%conc
      dconc_upwind_ddn = 1.d0
      dispersion_scale = 0.d0  ! creates zero-gradient outlet bc
    endif
    delta_conc = immiscible_auxvar_up%conc - immiscible_auxvar_dn%conc
    D_molecular = immiscible_parameter%diffusion_coef
    D_hyd_dn = D_mech_dn + &
               immiscible_auxvar_dn%effective_porosity * immiscible_auxvar_dn%sat * &
               material_auxvar_dn%tortuosity * D_molecular
    Deff_over_dist = dispersion_scale * D_hyd_dn / dist(0)
    ! Res[mol/sec]
    Res(immis_sol_tran_eq) = Res(immis_sol_tran_eq) + &
      (q * conc_upwind + & ! advection
       area * Deff_over_dist * delta_conc) * & ! hydrodynamic dispersion
      L_per_m3
    if (calculate_derivatives) then
      Jdn(immis_sol_tran_eq,immis_sol_tran_eq) = &
        (q * dconc_upwind_ddn + &
         area * Deff_over_dist * (-1.d0)) * &
        L_per_m3
      if (immis_liq_flow_eq > 0) then
        dD_mech_dn_dpdn = 0.d0
        dD_hyd_dn_dpdn = dD_mech_dn_dpdn + &
                         (immiscible_auxvar_dn%dsat_dp * &
                          immiscible_auxvar_dn%effective_porosity + &
                          immiscible_auxvar_dn%sat * immiscible_auxvar_dn%dpor_dp) * &
                         material_auxvar_dn%tortuosity * D_molecular
        dDeff_over_dist_dpdn = dD_hyd_dn_dpdn / dist(0)
        Jdn(immis_sol_tran_eq,immis_liq_flow_eq) = &
          (dq_dpdn * conc_upwind + &
           area * dDeff_over_dist_dpdn * (-1.d0)) * &
          L_per_m3
      endif
      if (immis_calc_adjoint) then
        select case (immis_adjoint_parameter)
          case(IMMIS_ADJOINT_PERMEABILITY)
            dDeff_over_dist_dKdn = 0.d0
            dResdparamdn(immis_sol_tran_eq,1) = &
              (dq_dKdn * conc_upwind + &
              area * dDeff_over_dist_dKdn * delta_conc) * &
              L_per_m3
          case(IMMIS_ADJOINT_POROSITY)
            tempreal = D_mech_dn + &
                       immiscible_auxvar_dn%sat * &
                       material_auxvar_dn%tortuosity * D_molecular
            Deff_over_dist = dispersion_scale * tempreal / dist(0)
            dResdparamdn(immis_sol_tran_eq,1) = &
              area * Deff_over_dist * delta_conc * L_per_m3
        end select
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"BFS(res,q,dc): ",8es12.4)') Res(immis_sol_tran_eq), q, &
                                                delta_conc
    endif
  endif
#endif

end subroutine ImmiscibleBCFluxHarmonicPermOnly

! ************************************************************************** !

subroutine ImmiscibleSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                        flow_src_sink_type, &
                        immiscible_auxvar,global_auxvar,material_auxvar, &
                        ss_flow_vol_flux,Res,Jdn,dResdparamdn, &
                        calculate_derivatives)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
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
  type(immiscible_auxvar_type) :: immiscible_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscReal :: tempreal
  PetscReal :: qsrc_m3, qsrc_L
  PetscReal, parameter :: L_per_m3 = 1.d3
  PetscInt :: dof_index, index_

  Res = 0.d0
  Jdn = 0.d0

  if (immis_liq_flow_eq > 0) then
    dof_index = flow_src_sink_mapping(IMMIS_COND_WATER_INDEX)
    qsrc_m3 = 0.d0
    tempreal = flow_aux_real_var(dof_index)
    select case(flow_src_sink_type(dof_index))
      case(VOLUMETRIC_RATE_SS)
        qsrc_m3 = tempreal
      case(SCALED_VOLUMETRIC_RATE_SS)
        index_ = flow_src_sink_mapping(IMMIS_COND_WATER_AUX_INDEX)
        qsrc_m3 = tempreal*flow_aux_real_var(index_)
      case default
        option%io_buffer = 'src_sink_type not supported in ImmiscibleSrcSink'
        call PrintErrMsg(option)
    end select
    ss_flow_vol_flux = qsrc_m3
    ! Res[m^3 liquid/sec]
    Res(immis_liq_flow_eq) = qsrc_m3 * immis_density_kmol

    if (calculate_derivatives) then
    endif
  endif

  if (immis_sol_tran_eq > 0) then
    if (qsrc_m3 >= 0.d0) then
      ! injection
      dof_index = flow_src_sink_mapping(IMMIS_COND_SOLUTE_INDEX)
      tempreal = flow_aux_real_var(dof_index)
    else
      ! extraction
      tempreal = immiscible_auxvar%conc
    endif
    qsrc_L = qsrc_m3 * L_per_m3
    Res(immis_sol_tran_eq) = qsrc_L * tempreal
    if (calculate_derivatives) then
      if (qsrc_m3 < 0.d0) then
       Jdn(immis_sol_tran_eq,immis_sol_tran_eq) = qsrc_L
      endif
    endif
  endif

end subroutine ImmiscibleSrcSink

! ************************************************************************** !

subroutine ImmiscibleAccumDerivative(immiscible_auxvar,global_auxvar, &
                                material_auxvar, &
                                material_auxvar_pert, &
                                option,Res,Jac,dResdparam)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  type(option_type) :: option
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jac(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparam(IMMIS_MAX_DOF,IMMIS_MAX_DOF)

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dJdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  call ImmiscibleAccumulation(immiscible_auxvar(ZERO_INTEGER), &
                         global_auxvar, &
                         material_auxvar, &
                         option,Res,Jac,dResdparam, &
                         .not.immis_numerical_derivatives)

  if (immis_numerical_derivatives) then
    do idof = 1, option%nflowdof
      call ImmiscibleAccumulation(immiscible_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar, &
                             option,res_pert,Jdum,dJdum, &
                             PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq))/immiscible_auxvar(idof)%pert
      enddo
    enddo
    if (immis_calc_adjoint) then
      call ImmiscibleAccumulation(immiscible_auxvar(option%nflowdof+1), &
                             global_auxvar, &
                             material_auxvar_pert(ONE_INTEGER), &
                             option,res_pert,Jdum,dJdum, &
                             PETSC_FALSE)
      idof = 1
      do ieq = 1, option%nflowdof
        dResdparam(ieq,idof) = (res_pert(ieq)-Res(ieq))/ &
                               immiscible_auxvar(ZERO_INTEGER)%mat_pert(idof)
      enddo
    endif
  endif

end subroutine ImmiscibleAccumDerivative

! ************************************************************************** !

subroutine XXFluxDerivative(immiscible_auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            material_auxvar_pert_up, &
                            immiscible_auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            material_auxvar_pert_dn, &
                            area, dist,immiscible_parameter,option,v_darcy, &
                            Res,Jup,Jdn,dResdparamup,dResdparamdn, &
                            debug_connection)
  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar_up(0:), immiscible_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_pert_up(:)
  type(material_auxvar_type) :: material_auxvar_pert_dn(:)
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dJdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  Jup = 0.d0
  Jdn = 0.d0

  call XXFlux(immiscible_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist,immiscible_parameter,option,v_darcy, &
              Res,Jup,Jdn,dResdparamup,dResdparamdn, &
              .not.immis_numerical_derivatives, &
              debug_connection)

  if (immis_numerical_derivatives) then
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(immiscible_auxvar_up(idof),global_auxvar_up, &
                  material_auxvar_up, &
                  immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,immiscible_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jup(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        immiscible_auxvar_up(idof)%pert
      enddo
    enddo
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(immiscible_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  immiscible_auxvar_dn(idof),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,immiscible_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        immiscible_auxvar_dn(idof)%pert
      enddo
    enddo
    if (immis_calc_adjoint) then
      idof = 1
      call XXFlux(immiscible_auxvar_up(option%nflowdof+1),global_auxvar_up, &
                  material_auxvar_pert_up(ONE_INTEGER), &
                  immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,immiscible_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        dResdparamup(ieq,idof) = &
          (res_pert(ieq)-Res(ieq))/ &
          immiscible_auxvar_up(ZERO_INTEGER)%mat_pert(idof)
      enddo
      call XXFlux(immiscible_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  immiscible_auxvar_dn(option%nflowdof+1),global_auxvar_dn, &
                  material_auxvar_pert_dn(ONE_INTEGER), &
                  area,dist,immiscible_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        dResdparamdn(ieq,idof) = &
          (res_pert(ieq)-Res(ieq))/ &
          immiscible_auxvar_dn(ZERO_INTEGER)%mat_pert(idof)
      enddo
    endif
  endif

end subroutine XXFluxDerivative

! ************************************************************************** !

subroutine XXBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                              immiscible_auxvar_up, &
                              global_auxvar_up, &
                              immiscible_auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              material_auxvar_pert_dn, &
                              area,dist,immiscible_parameter,option,v_darcy, &
                              Res,Jdn,dResdparamdn,debug_connection)
  !
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(1:1)
  PetscInt :: auxvar_mapping(IMMIS_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_pert_dn(:)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparamdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dJdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  Jdn = 0.d0

  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                immiscible_auxvar_up,global_auxvar_up, &
                immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                material_auxvar_dn, &
                area,dist,immiscible_parameter,option,v_darcy, &
                Res,Jdn,dResdparamdn,.not.immis_numerical_derivatives, &
                debug_connection)

  if (immis_numerical_derivatives) then
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                    immiscible_auxvar_up,global_auxvar_up, &
                    immiscible_auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    area,dist,immiscible_parameter,option,v_darcy_dum, &
                    res_pert,Jdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        immiscible_auxvar_dn(idof)%pert
      enddo
    enddo
    if (immis_calc_adjoint) then
      idof = 1
      call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                    immiscible_auxvar_up,global_auxvar_up, &
                    immiscible_auxvar_dn(option%nflowdof+1),global_auxvar_dn, &
                    material_auxvar_pert_dn(ONE_INTEGER), &
                    area,dist,immiscible_parameter,option,v_darcy_dum, &
                    res_pert,Jdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        dResdparamdn(ieq,idof) = &
          (res_pert(ieq)-Res(ieq))/ &
          immiscible_auxvar_dn(ZERO_INTEGER)%mat_pert(idof)
      enddo
    endif
  endif

end subroutine XXBCFluxDerivative

! ************************************************************************** !

subroutine ImmiscibleSrcSinkDerivative(option,flow_aux_real_var, &
                                  flow_src_sink_mapping, &
                                  flow_src_sink_type, &
                                  immiscible_auxvar,global_auxvar, &
                                  material_auxvar, &
                                  material_auxvar_pert, &
                                  ss_flow_vol_flux, &
                                  Res,Jac,dResdparam)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  type(immiscible_auxvar_type) :: immiscible_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jac(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: dResdparam(IMMIS_MAX_DOF,IMMIS_MAX_DOF)

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: dummy_real
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  Jac = 0.d0
  ! unperturbed immiscible_auxvars value
  call ImmiscibleSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                    flow_src_sink_type, &
                    immiscible_auxvar(ZERO_INTEGER),global_auxvar, &
                    material_auxvar, &
                    ss_flow_vol_flux, &
                    Res,Jac,dResdparam,.not.immis_numerical_derivatives)

  if (immis_numerical_derivatives) then
    ! perturbed immiscible_auxvars values
    do idof = 1, option%nflowdof
      call ImmiscibleSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                        flow_src_sink_type, &
                        immiscible_auxvar(idof),global_auxvar, &
                        material_auxvar, &
                        dummy_real, &
                        res_pert,Jdum,Jdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
      Jac(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                      immiscible_auxvar(idof)%pert
      enddo
    enddo
    if (immis_calc_adjoint) then
      idof = 1
      call ImmiscibleSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                        flow_src_sink_type, &
                        immiscible_auxvar(option%nflowdof+1),global_auxvar, &
                        material_auxvar_pert(ONE_INTEGER), &
                        dummy_real, &
                        res_pert,Jdum,Jdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        dResdparam(ieq,idof) = &
          (res_pert(ieq)-Res(ieq))/ &
          immiscible_auxvar(ZERO_INTEGER)%mat_pert(idof)
      enddo
    endif
  endif

end subroutine ImmiscibleSrcSinkDerivative

end module Immiscible_Common_module
