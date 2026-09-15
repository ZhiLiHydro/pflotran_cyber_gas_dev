module THC_EOS_Utils_module

! ----------------------------------------------------------------------------
! THC EOS utilities: liquid density, viscosity and their (P,T,C) derivatives
! via the EOS_Water extended interfaces, plus Somerton thermal conductivity.
!
! Author: Piyoosh Jaysaval
! Date:   06/09/26
! ----------------------------------------------------------------------------

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module, only : FMWH2O
  use EOS_Water_module, only : EOSWaterDensityExt, &
                               EOSWaterViscosityExt, &
                               EOSWaterDensity, &
                               EOSWaterViscosity, &
                               EOSWaterDensityIsBrine, &
                               EOSWaterSaturationPressure

  implicit none

  private

  ! solute molar mass [kg/mol]; NaCl default, set by SOLUTE_MOLAR_MASS
  PetscReal, public :: thc_molar_mass_solute = 0.05844d0
  ! reference density [kg/m^3] seeding the C -> salinity conversion
  PetscReal, parameter, public :: thc_density_reference = 1000.d0
  ! constituent thermal conductivities for the Somerton model [W/(m.K)]
  PetscReal, parameter, public :: thc_kappa_liquid = 0.6d0
  PetscReal, parameter, public :: thc_kappa_air    = 0.025d0

  ! finite-difference step in salinity mass fraction for d/dC
  PetscReal, parameter :: thc_salinity_perturbation = 1.d-6
  ! saturation floor keeping the Somerton derivative finite
  PetscReal, parameter, public :: thc_sat_floor = 1.d-6

  public :: THCConcToMassFraction, &
            THCDensityAndDerivs, &
            THCViscosityAndDerivs, &
            THCThermalConductivityEff

contains

! ****************************************************************************
subroutine THCConcToMassFraction(C, den_kg, M_s, s, ds_dC)
  !
  ! Converts concentration C [mol/L] to the salinity mass fraction s [kg/kg]
  ! required by the extended EOS routines, and returns ds/dC.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal, intent(in)  :: C       ! concentration [mol/L]
  PetscReal, intent(in)  :: den_kg  ! liquid density [kg/m^3]
  PetscReal, intent(in)  :: M_s     ! solute molar mass [kg/mol]
  PetscReal, intent(out) :: s       ! salinity mass fraction [kg/kg]
  PetscReal, intent(out) :: ds_dC   ! d(s)/d(C) [kg/kg / (mol/L)]

  ds_dC = M_s * 1.d3 / den_kg
  s     = C * ds_dC

end subroutine THCConcToMassFraction

! ****************************************************************************
subroutine THCDensityAndDerivs(T, P, C, den_kg, den_kmol, &
                                   dden_dp, dden_dT, dden_dC, ierr)
  !
  ! Liquid density rho_l(P,T,C) and its derivatives via EOSWaterDensityExt;
  ! d(rho)/dC by the chain rule through the salinity mass fraction.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal,      intent(in)    :: T         ! temperature [C]
  PetscReal,      intent(in)    :: P         ! pressure [Pa]
  PetscReal,      intent(in)    :: C         ! concentration [mol/L]
  PetscReal,      intent(out)   :: den_kg    ! rho_l [kg/m^3]
  PetscReal,      intent(out)   :: den_kmol  ! rho_l [kmol/m^3]
  PetscReal,      intent(out)   :: dden_dp   ! d(rho_kg)/dP [kg/(m^3.Pa)]
  PetscReal,      intent(out)   :: dden_dT   ! d(rho_kg)/dT [kg/(m^3.C)]
  PetscReal,      intent(out)   :: dden_dC   ! d(rho_kg)/dC [kg/(m^3.(mol/L))]
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: aux1(1)
  PetscReal :: dw, dwmol, dwp, dwt
  PetscReal :: dw_pert, dwmol_pert
  PetscReal :: s, ds_dC, dden_ds

  if (.not.EOSWaterDensityIsBrine()) then
    ! Pure-water correlation selected (IF97, IFC67, ...)
    call EOSWaterDensity(T, P, dw, dwmol, dwp, dwt, ierr)
    den_kg   = dw
    den_kmol = dwmol
    dden_dp  = dwp * FMWH2O
    dden_dT  = dwt * FMWH2O
    dden_dC  = 0.d0
    return
  endif

  ! C -> salinity mass fraction, one Picard pass on the EOS density
  call THCConcToMassFraction(C, thc_density_reference, &
                                 thc_molar_mass_solute, s, ds_dC)
  aux1(1) = s
  call EOSWaterDensityExt(T, P, aux1, dw, dwmol, ierr)
  call THCConcToMassFraction(C, dw, thc_molar_mass_solute, s, ds_dC)
  aux1(1) = s

  ! density with analytic P,T derivatives (molar basis)
  call EOSWaterDensityExt(T, P, aux1, dw, dwmol, dwp, dwt, ierr)
  den_kg   = dw
  den_kmol = dwmol

  ! molar -> mass density derivatives
  dden_dp = dwp * FMWH2O
  dden_dT = dwt * FMWH2O

  ! d(rho)/dC by chain rule, d(rho)/ds by finite difference
  aux1(1) = s + thc_salinity_perturbation
  call EOSWaterDensityExt(T, P, aux1, dw_pert, dwmol_pert, ierr)
  dden_ds = (dw_pert - dw) / thc_salinity_perturbation
  dden_dC = dden_ds * ds_dC

end subroutine THCDensityAndDerivs

! ****************************************************************************
subroutine THCViscosityAndDerivs(T, P, C, den_kg, vis, dvis_dT, dvis_dC, &
                                     ierr)
  !
  ! Liquid viscosity mu_l(T,C) and its derivatives via EOSWaterViscosityExt.
  ! Pressure dependence is neglected in the Jacobian, as in TH.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal,      intent(in)    :: T        ! temperature [C]
  PetscReal,      intent(in)    :: P        ! pressure [Pa]
  PetscReal,      intent(in)    :: C        ! concentration [mol/L]
  PetscReal,      intent(in)    :: den_kg   ! rho_l [kg/m^3]
  PetscReal,      intent(out)   :: vis      ! mu_l [Pa.s]
  PetscReal,      intent(out)   :: dvis_dT  ! d(mu_l)/dT [Pa.s/C]
  PetscReal,      intent(out)   :: dvis_dC  ! d(mu_l)/dC [Pa.s/(mol/L)]
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: aux1(1)
  PetscReal :: PS, dPS_dT
  PetscReal :: VW, dVW_dT, dVW_dP
  PetscReal :: VW_pert
  PetscReal :: s, ds_dC, dvis_ds

  ! saturation pressure is an interface input
  call EOSWaterSaturationPressure(T, PS, dPS_dT, ierr)

  if (.not.EOSWaterDensityIsBrine()) then
    ! Pure-water correlation
    call EOSWaterViscosity(T, P, PS, dPS_dT, VW, dVW_dT, dVW_dP, ierr)
    vis     = VW
    dvis_dT = dVW_dT
    dvis_dC = 0.d0
    return
  endif

  ! C -> salinity mass fraction at the current density
  call THCConcToMassFraction(C, den_kg, thc_molar_mass_solute, s, ds_dC)
  aux1(1) = s

  ! viscosity with analytic T derivative; dVW_dP discarded
  call EOSWaterViscosityExt(T, P, PS, dPS_dT, aux1, VW, dVW_dT, dVW_dP, &
                            ierr)
  vis     = VW
  dvis_dT = dVW_dT

  ! d(mu)/dC by chain rule, d(mu)/ds by finite difference
  aux1(1) = s + thc_salinity_perturbation
  call EOSWaterViscosityExt(T, P, PS, aux1, VW_pert, ierr)
  dvis_ds = (VW_pert - VW) / thc_salinity_perturbation
  dvis_dC = dvis_ds * ds_dC

end subroutine THCViscosityAndDerivs

! ****************************************************************************
subroutine THCThermalConductivityEff(S_l, phi, kappa_s, &
                                         kappa_eff, dkappa_dsat)
  !
  ! Effective thermal conductivity by the Somerton et al. (1974) mixing
  ! model and its saturation derivative.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal, intent(in)  :: S_l          ! liquid saturation [-]
  PetscReal, intent(in)  :: phi          ! porosity [-]
  PetscReal, intent(in)  :: kappa_s      ! solid thermal conductivity [W/(m.K)]
  PetscReal, intent(out) :: kappa_eff    ! effective conductivity [W/(m.K)]
  PetscReal, intent(out) :: dkappa_dsat  ! d(kappa_eff)/d(S_l) [W/(m.K)]

  PetscReal :: kappa_sat, kappa_dry, sqrt_sat, sat

  sat = max(S_l, thc_sat_floor)

  kappa_sat = kappa_s**(1.d0 - phi) * thc_kappa_liquid**phi
  kappa_dry = kappa_s**(1.d0 - phi) * thc_kappa_air**phi

  sqrt_sat    = sqrt(sat)
  kappa_eff   = kappa_dry + sqrt_sat * (kappa_sat - kappa_dry)
  dkappa_dsat = (kappa_sat - kappa_dry) / (2.d0 * sqrt_sat)

end subroutine THCThermalConductivityEff

end module THC_EOS_Utils_module
