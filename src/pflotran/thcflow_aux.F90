module THCFlow_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Material_Aux_module
  use EOS_Water_module
  use THCFlow_EOS_Utils_module

  implicit none

  private

  ! Solid / fluid thermal properties
  PetscReal, public :: thcflow_heat_capacity_fluid = 4184.d0  ! c_p,l [J/(kg.K)]
  PetscReal, public :: thcflow_heat_capacity_solid = 800.d0   ! c_s   [J/(kg.K)]
  PetscReal, public :: thcflow_density_solid       = 2650.d0  ! rho_s [kg/m^3]
  PetscReal, public :: thcflow_kappa_solid         = 3.d0     ! kappa_s [W/(m.K)]

  ! perturbation controls (mirror zflow_aux.F90)
  PetscReal, public :: thcflow_rel_pert = 1.d-8
  PetscReal, public :: thcflow_pres_min_pert = 1.d-2
  PetscReal, public :: thcflow_temp_min_pert = 1.d-6 ! not based on anything
  PetscReal, public :: thcflow_conc_min_pert = 1.d-6 ! not based on anything

  PetscReal, pointer, public :: thcflow_min_pert(:)

  PetscBool, public :: thcflow_calc_accum = PETSC_TRUE
  PetscBool, public :: thcflow_calc_flux = PETSC_TRUE
  PetscBool, public :: thcflow_calc_bcflux = PETSC_TRUE

  PetscBool, public :: thcflow_numerical_derivatives = PETSC_FALSE
  PetscBool, public :: thcflow_simult_function_evals = PETSC_TRUE
  PetscBool, public :: thcflow_tensorial_rel_perm = PETSC_FALSE
  PetscBool, public :: thcflow_acknowledge_no_compress = PETSC_FALSE

  ! debugging
  PetscInt, public :: thcflow_ni_count
  PetscInt, public :: thcflow_ts_cut_count
  PetscInt, public :: thcflow_ts_count
  PetscInt, public :: thcflow_debug_cell_id

  ! Fixed DOF indices -- all always active
  PetscInt, parameter, public :: thcflow_pressure_dof      = 1
  PetscInt, parameter, public :: thcflow_temperature_dof   = 2
  PetscInt, parameter, public :: thcflow_concentration_dof = 3
  PetscInt, parameter, public :: THCFLOW_NDOF = 3

  PetscInt, parameter, public :: THCFLOW_COND_WATER_INDEX = 1
  PetscInt, parameter, public :: THCFLOW_COND_ENERGY_INDEX = 2
  PetscInt, parameter, public :: THCFLOW_COND_SOLUTE_INDEX = 3
  PetscInt, parameter, public :: THCFLOW_COND_WATER_AUX_INDEX = 4
  PetscInt, parameter, public :: THCFLOW_MAX_INDEX = 4

  PetscInt, parameter, public :: THCFLOW_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: THCFLOW_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: THCFLOW_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: THCFLOW_UPDATE_FOR_BOUNDARY = 2

  PetscInt, parameter, public :: THCFLOW_LIQ_SAT_WRT_LIQ_PRES = 1
  PetscInt, parameter, public :: THCFLOW_LIQ_PRES_WRT_POROS = 2

  type, public :: thcflow_auxvar_type
    PetscReal :: pres ! liquid pressure P_l [Pa]
    PetscReal :: sat  ! liquid saturation S_l [-]
    PetscReal :: pc   ! capillary pressure [Pa]
    PetscReal :: kr   ! relative permeability [-]
    PetscReal :: effective_porosity   ! phi_eff [-]
    PetscReal :: dpor_dp              ! dphi/dP_l [1/Pa]
    PetscReal :: dsat_dp ! derivative of saturation wrt pressure [1/Pa]
    PetscReal :: dkr_dp  ! derivative of rel. perm. wrt pressure [1/Pa]
    PetscReal :: effective_saturation
    PetscReal :: deffsat_dp
    PetscReal :: conc ! concentration C [mol/L]
    PetscReal :: pert
    PetscReal :: mat_pert(1)

    PetscReal :: temp             ! T [C]
    PetscReal :: den_kg           ! rho_l [kg/m^3]
    PetscReal :: den_kmol         ! rho_l / M_w [kmol/m^3]
    PetscReal :: dden_dp          ! d(rho_l)/dP_l [kg/(m^3.Pa)]
    PetscReal :: dden_dT          ! d(rho_l)/dT [kg/(m^3.K)]
    PetscReal :: dden_dC          ! d(rho_l)/dC [kg/(m^3.(mol/L))]
    PetscReal :: vis              ! mu_l [Pa.s]
    PetscReal :: dvis_dT          ! d(mu_l)/dT [Pa.s/K]
    PetscReal :: dvis_dC          ! d(mu_l)/dC [Pa.s.L/mol]
    PetscReal :: diff_mol         ! D_mol [m^2/s]
    PetscReal :: ddiff_dT         ! d(D_mol)/dT [m^2/(s.K)]
    PetscReal :: ddiff_dC         ! d(D_mol)/dC [m^2.L/(s.mol)]
    PetscReal :: therm_cond_eff   ! kappa_eff [W/(m.K)]
    PetscReal :: dtherm_cond_dsat ! d(kappa_eff)/dS_l [W/(m.K)]
    PetscReal :: heat_cap_fluid   ! rho_l * c_p,l [J/(m^3.K)]
    PetscReal :: heat_cap_solid   ! rho_s * c_s   [J/(m^3.K)]
  end type thcflow_auxvar_type

  type, public :: thcflow_parameter_type
    PetscBool :: check_post_converged
    PetscReal, pointer :: tensorial_rel_perm_exponent(:,:)
    PetscReal :: diffusion_coef
  end type thcflow_parameter_type

  type, public :: thcflow_type
    PetscBool :: auxvars_up_to_date
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(thcflow_parameter_type), pointer :: thcflow_parameter
    type(thcflow_auxvar_type), pointer :: auxvars(:,:)
    type(thcflow_auxvar_type), pointer :: auxvars_bc(:)
    type(thcflow_auxvar_type), pointer :: auxvars_ss(:)
    type(material_auxvar_type), pointer :: material_auxvars_pert(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type thcflow_type

  interface THCFlowAuxVarDestroy
    module procedure THCFlowAuxVarSingleDestroy
    module procedure THCFlowAuxVarArray1Destroy
    module procedure THCFlowAuxVarArray2Destroy
  end interface THCFlowAuxVarDestroy

  interface THCFlowOutputAuxVars
    module procedure THCFlowOutputAuxVars1
  end interface THCFlowOutputAuxVars

  public :: THCFlowAuxCreate, &
            THCFlowAuxDestroy, &
            THCFlowAuxVarCompute, &
            THCFlowAuxVarInit, &
            THCFlowAuxVarCopy, &
            THCFlowAuxVarDestroy, &
            THCFlowAuxVarStrip, &
            THCFlowAuxVarPerturb, &
            THCFlowAuxMapConditionIndices, &
            THCFlowPrintAuxVars, &
            THCFlowOutputAuxVars, &
            THCFlowAuxTensorialRelPerm

contains

! ************************************************************************** !

function THCFlowAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(thcflow_type), pointer :: THCFlowAuxCreate

  type(thcflow_type), pointer :: aux

  nullify(thcflow_min_pert)

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%material_auxvars_pert)
  nullify(aux%matrix_zeroing)

  allocate(aux%thcflow_parameter)
  aux%thcflow_parameter%check_post_converged = PETSC_FALSE
  nullify(aux%thcflow_parameter%tensorial_rel_perm_exponent)
  aux%thcflow_parameter%diffusion_coef = 0.d0

  THCFlowAuxCreate => aux

end function THCFlowAuxCreate

! ************************************************************************** !

subroutine THCFlowAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  type(thcflow_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%kr = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%dpor_dp = 0.d0
  auxvar%dsat_dp = 0.d0
  auxvar%dkr_dp = 0.d0
  auxvar%effective_saturation = UNINITIALIZED_DOUBLE
  auxvar%deffsat_dp = UNINITIALIZED_DOUBLE
  auxvar%conc = 0.d0
  auxvar%pert = 0.d0
  auxvar%mat_pert = 0.d0

  auxvar%temp = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%den_kmol = 0.d0
  auxvar%dden_dp = 0.d0
  auxvar%dden_dT = 0.d0
  auxvar%dden_dC = 0.d0
  auxvar%vis = 0.d0
  auxvar%dvis_dT = 0.d0
  auxvar%dvis_dC = 0.d0
  auxvar%diff_mol = 0.d0
  auxvar%ddiff_dT = 0.d0
  auxvar%ddiff_dC = 0.d0
  auxvar%therm_cond_eff = 0.d0
  auxvar%dtherm_cond_dsat = 0.d0
  auxvar%heat_cap_fluid = 0.d0
  auxvar%heat_cap_solid = 0.d0

end subroutine THCFlowAuxVarInit

! ************************************************************************** !

subroutine THCFlowAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  type(thcflow_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%pc = auxvar%pc
  auxvar2%kr = auxvar%kr
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%dpor_dp = auxvar%dpor_dp
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%effective_saturation = auxvar%effective_saturation
  auxvar2%deffsat_dp = auxvar%deffsat_dp
  auxvar2%conc = auxvar%conc
  auxvar2%pert = auxvar%pert

  auxvar2%temp = auxvar%temp
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%den_kmol = auxvar%den_kmol
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dT = auxvar%dden_dT
  auxvar2%dden_dC = auxvar%dden_dC
  auxvar2%vis = auxvar%vis
  auxvar2%dvis_dT = auxvar%dvis_dT
  auxvar2%dvis_dC = auxvar%dvis_dC
  auxvar2%diff_mol = auxvar%diff_mol
  auxvar2%ddiff_dT = auxvar%ddiff_dT
  auxvar2%ddiff_dC = auxvar%ddiff_dC
  auxvar2%therm_cond_eff = auxvar%therm_cond_eff
  auxvar2%dtherm_cond_dsat = auxvar%dtherm_cond_dsat
  auxvar2%heat_cap_fluid = auxvar%heat_cap_fluid
  auxvar2%heat_cap_solid = auxvar%heat_cap_solid

end subroutine THCFlowAuxVarCopy

! ************************************************************************** !

subroutine THCFlowAuxVarCompute(x,thcflow_auxvar,global_auxvar, &
                                material_auxvar,characteristic_curves, &
                                natural_id,update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(:)
  type(thcflow_auxvar_type) :: thcflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscBool :: update_porosity
  PetscInt :: natural_id

  PetscBool :: saturated
  PetscReal :: dkr_dsat
  PetscReal :: deffsat_dsat
  PetscReal :: mobility, dmobility_dp, dmobility_dT, dmobility_dC
  PetscErrorCode :: ierr

  ! --------------------------------------------------------------------------
  ! Step 1: read primary variables (P,T,C) from solution vector
  ! --------------------------------------------------------------------------
  thcflow_auxvar%pres = x(thcflow_pressure_dof)
  thcflow_auxvar%temp = x(thcflow_temperature_dof)
  thcflow_auxvar%conc = x(thcflow_concentration_dof)

  ! --------------------------------------------------------------------------
  ! Step 2: effective porosity (soil compressibility)
  ! --------------------------------------------------------------------------
  ! %porosity should never be used. set to bogus value to catch misuse
  material_auxvar%porosity = -888.d0

  if (update_porosity .and. soil_compressibility_index > 0 .and. &
      thcflow_auxvar%pres > 0.d0) then
    call MaterialCompressSoil(material_auxvar,thcflow_auxvar%pres, &
                              thcflow_auxvar%effective_porosity, &
                              thcflow_auxvar%dpor_dp)
  else
    thcflow_auxvar%effective_porosity = material_auxvar%porosity_base
    thcflow_auxvar%dpor_dp = 0.d0
  endif
  material_auxvar%porosity = thcflow_auxvar%effective_porosity

  ! --------------------------------------------------------------------------
  ! Step 3: capillary pressure and saturation
  ! --------------------------------------------------------------------------
  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  thcflow_auxvar%pc = min(option%flow%reference_pressure - thcflow_auxvar%pres, &
                          characteristic_curves%saturation_function%pcmax)

  if (thcflow_auxvar%pc > 0.d0) then
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(thcflow_auxvar%pc, &
                                          thcflow_auxvar%sat, &
                                          thcflow_auxvar%dsat_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (thcflow_auxvar%dsat_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      ! ----------------------------------------------------------------------
      ! Step 4: relative permeability
      ! ----------------------------------------------------------------------
      dkr_dsat = 0.d0
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(thcflow_auxvar%sat, &
                                            thcflow_auxvar%kr, &
                                            dkr_dsat,option)
      thcflow_auxvar%dkr_dp = thcflow_auxvar%dsat_dp * dkr_dsat
      if (thcflow_tensorial_rel_perm) then
        call characteristic_curves%liq_rel_perm_function% &
          EffectiveSaturation(thcflow_auxvar%sat, &
                              thcflow_auxvar%effective_saturation, &
                              deffsat_dsat,option)
        thcflow_auxvar%deffsat_dp = deffsat_dsat * thcflow_auxvar%dsat_dp
      endif
    endif
  else
    saturated = PETSC_TRUE
  endif

  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when
  ! capillary pressure is very small
  if (saturated) then
    thcflow_auxvar%pc = 0.d0
    thcflow_auxvar%sat = 1.d0
    thcflow_auxvar%kr = 1.d0
    thcflow_auxvar%dsat_dp = 0.d0
    thcflow_auxvar%dkr_dp = 0.d0
    thcflow_auxvar%effective_saturation = 1.d0
    thcflow_auxvar%deffsat_dp = 0.d0
  endif

  if (option%iflag /= THCFLOW_UPDATE_FOR_DERIVATIVE) then
    global_auxvar%sat(1) = thcflow_auxvar%sat
    if (size(global_auxvar%sat) > 1) then
      global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
    endif
  endif

  ! --------------------------------------------------------------------------
  ! Step 5: liquid density rho_l(P,T,C) and its derivatives
  ! --------------------------------------------------------------------------
  ierr = 0
  call THCFlowDensityAndDerivs(thcflow_auxvar%temp,thcflow_auxvar%pres, &
                               thcflow_auxvar%conc,thcflow_auxvar%den_kg, &
                               thcflow_auxvar%den_kmol,thcflow_auxvar%dden_dp, &
                               thcflow_auxvar%dden_dT,thcflow_auxvar%dden_dC, &
                               ierr)
  if (ierr /= 0) then
    option%io_buffer = 'Error computing THCFLOW liquid density in &
      &THCFlowAuxVarCompute.'
    call PrintErrMsg(option)
  endif

  ! --------------------------------------------------------------------------
  ! Step 6: liquid viscosity mu_l(T,C) and its derivatives
  ! --------------------------------------------------------------------------
  call THCFlowViscosityAndDerivs(thcflow_auxvar%temp,thcflow_auxvar%pres, &
                                 thcflow_auxvar%conc,thcflow_auxvar%vis, &
                                 thcflow_auxvar%dvis_dT,thcflow_auxvar%dvis_dC, &
                                 ierr)
  if (ierr /= 0) then
    option%io_buffer = 'Error computing THCFLOW liquid viscosity in &
      &THCFlowAuxVarCompute.'
    call PrintErrMsg(option)
  endif

  ! --------------------------------------------------------------------------
  ! Step 7: effective porous-medium thermal conductivity (Somerton)
  ! --------------------------------------------------------------------------
  call THCFlowThermalConductivityEff(thcflow_auxvar%sat, &
                                     thcflow_auxvar%effective_porosity, &
                                     thcflow_kappa_solid, &
                                     thcflow_auxvar%therm_cond_eff, &
                                     thcflow_auxvar%dtherm_cond_dsat)

  ! --------------------------------------------------------------------------
  ! Step 8: derived quantities
  ! --------------------------------------------------------------------------
  ! volumetric heat capacities
  thcflow_auxvar%heat_cap_fluid = thcflow_auxvar%den_kg * &
                                  thcflow_heat_capacity_fluid
  thcflow_auxvar%heat_cap_solid = thcflow_density_solid * &
                                  thcflow_heat_capacity_solid

  ! derived mobility derivatives. These are recomputed
  ! by the flux kernel from the stored primitives (kr, dkr_dp, vis, dvis_dT,
  ! dvis_dC); evaluated here to document the chain rule.
  ! Assumes: no pressure dependence of viscosity mu_l, and
  !          no temperature or concentration dependence of kr. Then:
  !   mobility       = kr / mu_l
  !   d(mob)/dP_l    = dkr_dp / mu_l                      (through saturation)
  !   d(mob)/dT      = -kr / mu_l^2 * dmu_l/dT
  !   d(mob)/dC      = -kr / mu_l^2 * dmu_l/dC
  mobility     = thcflow_auxvar%kr / thcflow_auxvar%vis
  dmobility_dp = thcflow_auxvar%dkr_dp / thcflow_auxvar%vis
  dmobility_dT = -thcflow_auxvar%kr / thcflow_auxvar%vis**2 * &
                 thcflow_auxvar%dvis_dT
  dmobility_dC = -thcflow_auxvar%kr / thcflow_auxvar%vis**2 * &
                 thcflow_auxvar%dvis_dC

end subroutine THCFlowAuxVarCompute

! ************************************************************************** !

subroutine THCFlowAuxVarPerturb(x,thcflow_auxvar,global_auxvar, &
                                material_auxvar, &
                                material_auxvar_pert, &
                                characteristic_curves, &
                                natural_id, &
                                option)
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  PetscReal :: x(:)
  type(option_type) :: option
  PetscInt :: natural_id
  type(thcflow_auxvar_type) :: thcflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  class(characteristic_curves_type) :: characteristic_curves

  PetscInt :: idof
  PetscReal :: x_pert(THCFLOW_NDOF), pert

  ! THCFLOW_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = THCFLOW_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    pert = x(idof)*thcflow_rel_pert+thcflow_min_pert(idof)
    thcflow_auxvar(idof)%pert = pert
    x_pert(1:option%nflowdof) = x
    x_pert(idof) = x(idof) + pert
    call THCFlowAuxVarCompute(x_pert,thcflow_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves,natural_id, &
                              PETSC_TRUE,option)
  enddo

end subroutine THCFlowAuxVarPerturb

! ************************************************************************** !

subroutine THCFlowAuxTensorialRelPerm(auxvar,tensorial_rel_perm_exponent, &
                                      dist,rel_perm,drel_perm_dp,option)
  !
  ! Computes the tensorial (direction-dependent) relative permeability
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module
  use Utility_module

  implicit none

  type(thcflow_auxvar_type) :: auxvar
  PetscReal :: tensorial_rel_perm_exponent(3)
  PetscReal :: dist(-1:3)
  PetscReal :: rel_perm
  PetscReal :: drel_perm_dp
  type(option_type) :: option

  PetscReal :: exponent_
  PetscReal :: tensorial_scale

  exponent_ = UtilityTensorToScalar(dist,tensorial_rel_perm_exponent)

  ! remember that the default 0.5 was subtracted from the tensorial value
  ! in THCFlowSetup. If 0.5 is specified for the tensorial exponent in the
  ! input file, this value will be 0.
  tensorial_scale = auxvar%effective_saturation**exponent_
  rel_perm = auxvar%kr * tensorial_scale
  drel_perm_dp = auxvar%dkr_dp * tensorial_scale + &
                 exponent_ * rel_perm / &
                 auxvar%effective_saturation * &
                 auxvar%deffsat_dp

end subroutine THCFlowAuxTensorialRelPerm

! ************************************************************************** !

subroutine THCFlowPrintAuxVars(thcflow_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(thcflow_auxvar_type) :: thcflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', thcflow_auxvar%pres
  print *, '             temperature: ', thcflow_auxvar%temp
  print *, '           concentration: ', thcflow_auxvar%conc
  print *, '      capillary pressure: ', thcflow_auxvar%pc
  print *, '       liquid saturation: ', thcflow_auxvar%sat
  print *, '      liquid sat (deriv): ', thcflow_auxvar%dsat_dp
  print *, '         liquid rel perm: ', thcflow_auxvar%kr
  print *, ' liquid rel perm (deriv): ', thcflow_auxvar%dkr_dp
  print *, '      effective porosity: ', thcflow_auxvar%effective_porosity
  print *, '   eff. porosity (deriv): ', thcflow_auxvar%dpor_dp
  print *, '          liquid density: ', thcflow_auxvar%den_kg
  print *, '        liquid viscosity: ', thcflow_auxvar%vis
  print *, '  effective therm. cond.: ', thcflow_auxvar%therm_cond_eff
  print *, '--------------------------------------------------------'

end subroutine THCFlowPrintAuxVars

! ************************************************************************** !

subroutine THCFlowOutputAuxVars1(thcflow_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(thcflow_auxvar_type) :: thcflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=IUNIT_TEMP,file=string2,position='append')
  else
    open(unit=IUNIT_TEMP,file=string2)
  endif

  write(IUNIT_TEMP,*) '--------------------------------------------------------'
  write(IUNIT_TEMP,*) trim(string)
  write(IUNIT_TEMP,*) '                 cell id: ', natural_id
  write(IUNIT_TEMP,*) '         liquid pressure: ', thcflow_auxvar%pres
  write(IUNIT_TEMP,*) '             temperature: ', thcflow_auxvar%temp
  write(IUNIT_TEMP,*) '           concentration: ', thcflow_auxvar%conc
  write(IUNIT_TEMP,*) '      capillary pressure: ', thcflow_auxvar%pc
  write(IUNIT_TEMP,*) '       liquid saturation: ', thcflow_auxvar%sat
  write(IUNIT_TEMP,*) '      liquid sat (deriv): ', thcflow_auxvar%dsat_dp
  write(IUNIT_TEMP,*) '         liquid rel perm: ', thcflow_auxvar%kr
  write(IUNIT_TEMP,*) ' liquid rel perm (deriv): ', thcflow_auxvar%dkr_dp
  write(IUNIT_TEMP,*) '      effective porosity: ', thcflow_auxvar &
                                                      %effective_porosity
  write(IUNIT_TEMP,*) '   eff. porosity (deriv): ', thcflow_auxvar%dpor_dp
  write(IUNIT_TEMP,*) '          liquid density: ', thcflow_auxvar%den_kg
  write(IUNIT_TEMP,*) '        liquid viscosity: ', thcflow_auxvar%vis
  write(IUNIT_TEMP,*) '  effective therm. cond.: ', thcflow_auxvar &
                                                      %therm_cond_eff
  write(IUNIT_TEMP,*) '--------------------------------------------------------'

  close(IUNIT_TEMP)

end subroutine THCFlowOutputAuxVars1

! ************************************************************************** !

function THCFlowAuxMapConditionIndices(include_water_aux)
  !
  ! Maps indexing of conditions
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  PetscBool :: include_water_aux

  PetscInt, pointer :: THCFlowAuxMapConditionIndices(:)

  PetscInt, pointer :: mapping(:)

  PetscInt :: temp_int

  allocate(mapping(THCFLOW_MAX_INDEX))
  mapping = UNINITIALIZED_INTEGER

  ! THCFLOW always has all three equations active (fixed DOF), so every branch
  ! is taken unconditionally.
  temp_int = 0
  temp_int = temp_int + 1
  mapping(THCFLOW_COND_WATER_INDEX) = temp_int
  if (include_water_aux) then
    temp_int = temp_int + 1
    mapping(THCFLOW_COND_WATER_AUX_INDEX) = temp_int
  endif
  temp_int = temp_int + 1
  mapping(THCFLOW_COND_ENERGY_INDEX) = temp_int
  temp_int = temp_int + 1
  mapping(THCFLOW_COND_SOLUTE_INDEX) = temp_int

  THCFlowAuxMapConditionIndices => mapping

end function THCFlowAuxMapConditionIndices

! ************************************************************************** !

subroutine THCFlowAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(thcflow_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call THCFlowAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine THCFlowAuxVarSingleDestroy

! ************************************************************************** !

subroutine THCFlowAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(thcflow_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call THCFlowAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THCFlowAuxVarArray1Destroy

! ************************************************************************** !

subroutine THCFlowAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(thcflow_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call THCFlowAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THCFlowAuxVarArray2Destroy

! ************************************************************************** !

subroutine THCFlowMaterialAuxVarDestroy(auxvars)
  !
  ! Deallocates material auxiliary object for thcflow perturbation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(material_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call MaterialAuxVarStrip(auxvars(idof,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THCFlowMaterialAuxVarDestroy

! ************************************************************************** !

subroutine THCFlowAuxVarStrip(auxvar)
  !
  ! THCFlowAuxVarDestroy: Deallocates a thcflow auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(thcflow_auxvar_type) :: auxvar

end subroutine THCFlowAuxVarStrip

! ************************************************************************** !

subroutine THCFlowAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(thcflow_type), pointer :: aux

  call DeallocateArray(thcflow_min_pert)

  if (.not.associated(aux)) return

  call THCFlowAuxVarDestroy(aux%auxvars)
  call THCFlowAuxVarDestroy(aux%auxvars_bc)
  call THCFlowAuxVarDestroy(aux%auxvars_ss)
  call THCFlowMaterialAuxVarDestroy(aux%material_auxvars_pert)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%thcflow_parameter)) then
    call DeallocateArray(aux%thcflow_parameter%tensorial_rel_perm_exponent)
    deallocate(aux%thcflow_parameter)
  endif
  nullify(aux%thcflow_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine THCFlowAuxDestroy

end module THCFlow_Aux_module
