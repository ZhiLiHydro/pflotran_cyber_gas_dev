module Compositional_Flow_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Material_Aux_module

  implicit none

  private

  PetscReal, public :: compflow_rel_pert = 1.d-8
  PetscReal, public :: compflow_pres_min_pert = 1.d-2
  PetscReal, public :: compflow_temp_min_pert = 1.d-6 ! not based on anything
  PetscReal, public :: compflow_conc_min_pert = 1.d-6 ! not based on anything

  PetscReal, pointer, public :: compflow_min_pert(:)

  PetscBool, public :: compflow_calc_accum = PETSC_TRUE
  PetscBool, public :: compflow_calc_flux = PETSC_TRUE
  PetscBool, public :: compflow_calc_bcflux = PETSC_TRUE
  PetscBool, public :: compflow_calc_chem = PETSC_TRUE

  PetscBool, public :: compflow_numerical_derivatives = PETSC_FALSE
  PetscBool, public :: compflow_simult_function_evals = PETSC_TRUE
  PetscBool, public :: compflow_acknowledge_no_compress = PETSC_FALSE
  PetscBool, public :: compflow_shear_thinning = PETSC_FALSE
  PetscReal, public :: compflow_decay_rate_constant = 0.d0

  ! debugging
  PetscInt, public :: compflow_ni_count
  PetscInt, public :: compflow_ts_cut_count
  PetscInt, public :: compflow_ts_count
  PetscInt, public :: compflow_debug_cell_id

  ! process models
  PetscInt, public :: compflow_liq_flow_eq = UNINITIALIZED_INTEGER
  PetscInt, public :: compflow_heat_tran_eq = UNINITIALIZED_INTEGER
  PetscInt, public :: compflow_sol_tran_eq = UNINITIALIZED_INTEGER
  PetscInt, parameter, public :: COMPFLOW_MAX_DOF = 3

  PetscInt, parameter, public :: COMPFLOW_COND_WATER_INDEX = 1
  PetscInt, parameter, public :: COMPFLOW_COND_ENERGY_INDEX = 2
  PetscInt, parameter, public :: COMPFLOW_COND_SOLUTE_INDEX = 3
  PetscInt, parameter, public :: COMPFLOW_COND_WATER_AUX_INDEX = 4
  PetscInt, parameter, public :: COMPFLOW_MAX_INDEX = 4

  PetscInt, parameter, public :: COMPFLOW_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: COMPFLOW_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: COMPFLOW_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: COMPFLOW_UPDATE_FOR_BOUNDARY = 2

  PetscInt, parameter, public :: COMPFLOW_LIQ_SAT_WRT_LIQ_PRES = 1
  PetscInt, parameter, public :: COMPFLOW_LIQ_PRES_WRT_POROS = 2

  type, public :: compflow_auxvar_type
    PetscReal :: pres ! liquid pressure
    PetscReal :: sat  ! liquid saturation
    PetscReal :: pc   ! capillary pressure
    PetscReal :: kr   ! relative permeability
    PetscReal :: effective_porosity
    PetscReal :: dpor_dp
    PetscReal :: dsat_dp ! derivative of saturation wrt pressure
    PetscReal :: dkr_dp  ! derivative of rel. perm. wrt pressure
    PetscReal :: temp ! temperature
    PetscReal :: conc ! concentration
    PetscReal :: den
    PetscReal :: vis
    PetscReal :: dden_dp
    PetscReal :: dden_dt
    PetscReal :: dvis_dp
    PetscReal :: dvis_dt
    PetscReal :: dvis_dc
    PetscReal :: pert
    PetscReal :: mat_pert(1)
  end type compflow_auxvar_type

  type, public :: compflow_parameter_type
    PetscBool :: check_post_converged
    PetscReal :: diffusion_coef
  end type compflow_parameter_type

  type, public :: compflow_type
    PetscBool :: auxvars_up_to_date
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(compflow_parameter_type), pointer :: compflow_parameter
    type(compflow_auxvar_type), pointer :: auxvars(:,:)
    type(compflow_auxvar_type), pointer :: auxvars_bc(:)
    type(compflow_auxvar_type), pointer :: auxvars_ss(:)
    type(material_auxvar_type), pointer :: material_auxvars_pert(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type compflow_type

  interface CompFlowAuxVarDestroy
    module procedure CompFlowAuxVarSingleDestroy
    module procedure CompFlowAuxVarArray1Destroy
    module procedure CompFlowAuxVarArray2Destroy
  end interface CompFlowAuxVarDestroy

  interface CompFlowOutputAuxVars
    module procedure CompFlowOutputAuxVars1
  end interface CompFlowOutputAuxVars

  public :: CompFlowAuxCreate, &
            CompFlowAuxDestroy, &
            CompFlowAuxVarCompute, &
            CompFlowAuxVarInit, &
            CompFlowAuxVarCopy, &
            CompFlowAuxVarDestroy, &
            CompFlowAuxVarStrip, &
            CompFlowAuxVarPerturb, &
            CompFlowAuxMapConditionIndices, &
            CompFlowPrintAuxVars, &
            CompFlowOutputAuxVars

contains

! ************************************************************************** !

function CompFlowAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(compflow_type), pointer :: CompFlowAuxCreate

  type(compflow_type), pointer :: aux

  nullify(compflow_min_pert)

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

  allocate(aux%compflow_parameter)
  aux%compflow_parameter%check_post_converged = PETSC_FALSE
  aux%compflow_parameter%diffusion_coef = 0.d0

  CompFlowAuxCreate => aux

end function CompFlowAuxCreate

! ************************************************************************** !

subroutine CompFlowAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module

  implicit none

  type(compflow_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%kr = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%dpor_dp = 0.d0
  auxvar%dsat_dp = 0.d0
  auxvar%dkr_dp = 0.d0
  auxvar%temp = 0.d0
  auxvar%conc = 0.d0
  auxvar%den = UNINITIALIZED_DOUBLE
  auxvar%vis = UNINITIALIZED_DOUBLE
  auxvar%dden_dp = UNINITIALIZED_DOUBLE
  auxvar%dden_dt = UNINITIALIZED_DOUBLE
  auxvar%dvis_dp = UNINITIALIZED_DOUBLE
  auxvar%dvis_dt = UNINITIALIZED_DOUBLE
  auxvar%dvis_dc = UNINITIALIZED_DOUBLE

  auxvar%pert = 0.d0
  auxvar%mat_pert = 0.d0

end subroutine CompFlowAuxVarInit

! ************************************************************************** !

subroutine CompFlowAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module

  implicit none

  type(compflow_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%pc = auxvar%pc
  auxvar2%kr = auxvar%kr
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%dpor_dp = auxvar%dpor_dp
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%temp = auxvar%temp
  auxvar2%conc = auxvar%conc
  auxvar2%den = auxvar%den
  auxvar2%vis = auxvar%vis
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dt = auxvar%dden_dt
  auxvar2%dvis_dp = auxvar%dvis_dp
  auxvar2%dvis_dt = auxvar%dvis_dt
  auxvar2%dvis_dc = auxvar%dvis_dc
  auxvar2%pert = auxvar%pert

end subroutine CompFlowAuxVarCopy

! ************************************************************************** !

subroutine CompFlowAuxVarCompute(x,compflow_auxvar,global_auxvar, &
                              material_auxvar,characteristic_curves, &
                              natural_id,update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_module
  use EOS_Water_class

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(:)
  type(compflow_auxvar_type) :: compflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscBool :: update_porosity
  PetscInt :: natural_id

  PetscBool :: saturated
  PetscReal :: dkr_dsat
  PetscReal :: pw

  if (compflow_liq_flow_eq > 0) then
    compflow_auxvar%pres = x(compflow_liq_flow_eq)
  else
    compflow_auxvar%pres = option%flow%reference_pressure
  endif
  if (compflow_heat_tran_eq > 0) then
    compflow_auxvar%temp = x(compflow_heat_tran_eq)
  else
    compflow_auxvar%temp = option%flow%reference_temperature
  endif
  if (compflow_sol_tran_eq > 0) then
    compflow_auxvar%conc = x(compflow_sol_tran_eq)
  else
    compflow_auxvar%conc = 0.d0
  endif

  if (update_porosity .and. soil_compressibility_index > 0 .and. &
      compflow_auxvar%pres > 0.d0) then
    call MaterialCompressSoil(material_auxvar,compflow_auxvar%pres, &
                              compflow_auxvar%effective_porosity, &
                              compflow_auxvar%dpor_dp)
  else
    compflow_auxvar%effective_porosity = material_auxvar%porosity_base
    compflow_auxvar%dpor_dp = 0.d0
  endif
  material_auxvar%porosity = compflow_auxvar%effective_porosity

  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  compflow_auxvar%pc = min(option%flow%reference_pressure - compflow_auxvar%pres, &
                        characteristic_curves%saturation_function%pcmax)

  if (compflow_auxvar%pc > 0.d0) then
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(compflow_auxvar%pc, &
                                          compflow_auxvar%sat, &
                                          compflow_auxvar%dsat_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (compflow_auxvar%dsat_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      dkr_dsat = 0.d0
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(compflow_auxvar%sat, &
                                            compflow_auxvar%kr, &
                                            dkr_dsat,option)
      compflow_auxvar%dkr_dp = compflow_auxvar%dsat_dp * dkr_dsat
    endif
  else
    saturated = PETSC_TRUE
  endif

  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when
  ! capillary pressure is very small
  if (saturated) then
    compflow_auxvar%pc = 0.d0
    compflow_auxvar%sat = 1.d0
    compflow_auxvar%kr = 1.d0
    compflow_auxvar%dsat_dp = 0.d0
    compflow_auxvar%dkr_dp = 0.d0
  endif

  pw = option%flow%reference_pressure
  if (saturated) then
    pw = max(compflow_auxvar%pres, pw)
  endif

  compflow_auxvar%den = 1.d3
  compflow_auxvar%vis = 1.d-3
  compflow_auxvar%dden_dp = 0.d0
  compflow_auxvar%dden_dt = 0.d0
  compflow_auxvar%dvis_dp = 0.d0
  compflow_auxvar%dvis_dt = 0.d0
  compflow_auxvar%dvis_dc = 0.d0
  if (compflow_shear_thinning) then
    call CompFlowViscosity(compflow_auxvar,global_auxvar, &
                           material_auxvar,option)
  endif

  if (option%iflag /= COMPFLOW_UPDATE_FOR_DERIVATIVE) then
    global_auxvar%sat(1) = compflow_auxvar%sat
    if (size(global_auxvar%sat) > 1) then
      global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
    endif
  endif

end subroutine CompFlowAuxVarCompute

! ************************************************************************** !

subroutine CompFlowViscosity(auxvar,global_auxvar,material_auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Global_Aux_module
  use Option_module

  implicit none

  type(compflow_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal, parameter :: a = 1.16d-9
  PetscReal, parameter :: b = 2.769d0
  PetscReal, parameter :: alpha = 7.61d-9
  PetscReal, parameter :: beta = 2.737d0
  PetscReal, parameter :: n = 0.7d0
  PetscReal, parameter :: alpha_shear = 6.d0
  PetscReal, parameter :: mu_inf = 1.d-3
  PetscReal :: K
  PetscReal :: K_gamma
  PetscReal :: K_gamma_n
  PetscReal :: gamma
  PetscReal :: mu_0
  PetscReal :: one_over_one_plus_K_gamma_n
  PetscReal :: dmu_0_dc, dK_dc, df_dc

  if (auxvar%conc <= 0.d0) then
    auxvar%vis = mu_inf
    auxvar%dvis_dc = 0.d0
    return
  endif

  mu_0 = a * auxvar%conc ** b + mu_inf
  K = alpha * auxvar%conc ** beta
#if 0
  gamma = alpha_shear * global_auxvar%darcy_vel(1) / &  ! use Darcy velocity
          material_auxvar%soil_properties(mean_soil_grain_size_index)
#else
  gamma = alpha_shear * global_auxvar%darcy_vel(1) / &
          sqrt(material_auxvar%permeability(1)*material_auxvar%porosity)
#endif

  K_gamma = K * gamma
  K_gamma_n = K_gamma**n
  one_over_one_plus_K_gamma_n = 1.d0 / (1.d0 + K_gamma_n)
  auxvar%vis = mu_0 * one_over_one_plus_K_gamma_n

  if (compflow_numerical_derivatives) then
    auxvar%dvis_dc = UNINITIALIZED_DOUBLE
  else
    ! vis = mu_0 / (1 + (K*gamma)**n);  mu_0 = a*c^b + mu_inf
    ! dmu_0/dc = b*a*c^(b-1);  dK/dc = beta*alpha*c^(beta-1)
    ! gamma held fixed here
    dmu_0_dc = b * a * auxvar%conc**(b - 1.d0)
    if (K_gamma > 0.d0) then
      dK_dc = beta * alpha * auxvar%conc**(beta - 1.d0)
      df_dc = -n * one_over_one_plus_K_gamma_n**2 * &
              K_gamma**(n - 1.d0) * gamma * dK_dc
    else
      df_dc = 0.d0
    endif
    auxvar%dvis_dc = mu_0 * df_dc + &
                     one_over_one_plus_K_gamma_n * dmu_0_dc
  endif

end subroutine CompFlowViscosity

! ************************************************************************** !

subroutine CompFlowAuxVarPerturb(x,compflow_auxvar,global_auxvar, &
                              material_auxvar, &
                              characteristic_curves, &
                              natural_id,option)
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  PetscReal :: x(:)
  type(option_type) :: option
  PetscInt :: natural_id
  type(compflow_auxvar_type) :: compflow_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscInt :: idof
  PetscReal :: x_pert(COMPFLOW_MAX_DOF), pert

  ! COMPFLOW_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = COMPFLOW_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    pert = x(idof)*compflow_rel_pert+compflow_min_pert(idof)
    compflow_auxvar(idof)%pert = pert
    x_pert(1:option%nflowdof) = x
    x_pert(idof) = x(idof) + pert
    call CompFlowAuxVarCompute(x_pert,compflow_auxvar(idof),global_auxvar, &
                            material_auxvar, &
                            characteristic_curves, &
                            natural_id,PETSC_TRUE,option)
  enddo

end subroutine CompFlowAuxVarPerturb

! ************************************************************************** !

subroutine CompFlowPrintAuxVars(compflow_auxvar,global_auxvar,material_auxvar, &
                             natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(compflow_auxvar_type) :: compflow_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', compflow_auxvar%pres
  print *, '      capillary pressure: ', compflow_auxvar%pc
  print *, '       liquid saturation: ', compflow_auxvar%sat
  print *, '      liquid sat (deriv): ', compflow_auxvar%dsat_dp
  print *, '         liquid rel perm: ', compflow_auxvar%kr
  print *, ' liquid rel perm (deriv): ', compflow_auxvar%dkr_dp
  print *, '      effective porosity: ', compflow_auxvar%effective_porosity
  print *, '   eff. porosity (deriv): ', compflow_auxvar%dpor_dp
  print *, '--------------------------------------------------------'

end subroutine CompFlowPrintAuxVars

! ************************************************************************** !

subroutine CompFlowOutputAuxVars1(compflow_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(compflow_auxvar_type) :: compflow_auxvar
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
  write(IUNIT_TEMP,*) '             cell id: ', natural_id
  write(IUNIT_TEMP,*) '                 cell id: ', natural_id
  write(IUNIT_TEMP,*) '         liquid pressure: ', compflow_auxvar%pres
  write(IUNIT_TEMP,*) '      capillary pressure: ', compflow_auxvar%pc
  write(IUNIT_TEMP,*) '       liquid saturation: ', compflow_auxvar%sat
  write(IUNIT_TEMP,*) '      liquid sat (deriv): ', compflow_auxvar%dsat_dp
  write(IUNIT_TEMP,*) '         liquid rel perm: ', compflow_auxvar%kr
  write(IUNIT_TEMP,*) ' liquid rel perm (deriv): ', compflow_auxvar%dkr_dp
  write(IUNIT_TEMP,*) '      effective porosity: ', compflow_auxvar%effective_porosity
  write(IUNIT_TEMP,*) '   eff. porosity (deriv): ', compflow_auxvar%dpor_dp
  write(IUNIT_TEMP,*) '--------------------------------------------------------'

  close(IUNIT_TEMP)

end subroutine CompFlowOutputAuxVars1

! ************************************************************************** !

function CompFlowAuxMapConditionIndices(include_water_aux)
  !
  ! Maps indexing of conditions
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  use Option_module

  implicit none

  PetscBool :: include_water_aux

  PetscInt, pointer :: CompFlowAuxMapConditionIndices(:)

  PetscInt, pointer :: mapping(:)

  PetscInt :: temp_int

  allocate(mapping(COMPFLOW_MAX_INDEX))
  mapping = UNINITIALIZED_INTEGER

  temp_int = 0
  if (compflow_liq_flow_eq > 0) then
    temp_int = temp_int + 1
    mapping(COMPFLOW_COND_WATER_INDEX) = temp_int
    if (include_water_aux) then
      temp_int = temp_int + 1
      mapping(COMPFLOW_COND_WATER_AUX_INDEX) = temp_int
    endif
  endif
  if (compflow_heat_tran_eq > 0) then
    temp_int = temp_int + 1
    mapping(COMPFLOW_COND_ENERGY_INDEX) = temp_int
  endif
  if (compflow_sol_tran_eq > 0) then
    temp_int = temp_int + 1
    mapping(COMPFLOW_COND_SOLUTE_INDEX) = temp_int
  endif

  CompFlowAuxMapConditionIndices => mapping

end function CompFlowAuxMapConditionIndices

! ************************************************************************** !

subroutine CompFlowAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  implicit none

  type(compflow_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call CompFlowAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine CompFlowAuxVarSingleDestroy

! ************************************************************************** !

subroutine CompFlowAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  implicit none

  type(compflow_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call CompFlowAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine CompFlowAuxVarArray1Destroy

! ************************************************************************** !

subroutine CompFlowAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !

  implicit none

  type(compflow_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call CompFlowAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine CompFlowAuxVarArray2Destroy

! ************************************************************************** !

subroutine CompFlowMaterialAuxVarDestroy(auxvars)
  !
  ! Deallocates material auxiliary object for brine perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
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

end subroutine CompFlowMaterialAuxVarDestroy

! ************************************************************************** !

subroutine CompFlowAuxVarStrip(auxvar)
  !
  ! CompFlowAuxVarDestroy: Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(compflow_auxvar_type) :: auxvar

end subroutine CompFlowAuxVarStrip

! ************************************************************************** !

subroutine CompFlowAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/26
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(compflow_type), pointer :: aux

  call DeallocateArray(compflow_min_pert)

  if (.not.associated(aux)) return

  call CompFlowAuxVarDestroy(aux%auxvars)
  call CompFlowAuxVarDestroy(aux%auxvars_bc)
  call CompFlowAuxVarDestroy(aux%auxvars_ss)
  call CompFlowMaterialAuxVarDestroy(aux%material_auxvars_pert)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%compflow_parameter)) then
    deallocate(aux%compflow_parameter)
  endif
  nullify(aux%compflow_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine CompFlowAuxDestroy

end module Compositional_Flow_Aux_module
