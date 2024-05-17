module Immiscible_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Material_Aux_module

  implicit none

  private

  PetscReal, public :: immis_density_kg = 1000.d0
  PetscReal, public :: immis_density_kmol = UNINITIALIZED_DOUBLE
  PetscReal, public :: immis_viscosity = 1.d-3

  PetscReal, public :: immis_rel_pert = 1.d-8
  PetscReal, public :: immis_pres_min_pert = 1.d-2
  PetscReal, public :: immis_temp_min_pert = 1.d-6 ! not based on anything
  PetscReal, public :: immis_conc_min_pert = 1.d-6 ! not based on anything

  PetscReal, pointer, public :: immis_min_pert(:)

  PetscBool, public :: immis_calc_accum = PETSC_TRUE
  PetscBool, public :: immis_calc_flux = PETSC_TRUE
  PetscBool, public :: immis_calc_bcflux = PETSC_TRUE
  PetscBool, public :: immis_calc_chem = PETSC_TRUE

  PetscBool, public :: immis_numerical_derivatives = PETSC_FALSE
  PetscBool, public :: immis_simult_function_evals = PETSC_TRUE
  PetscBool, public :: immis_tensorial_rel_perm = PETSC_FALSE
  PetscBool, public :: immis_acknowledge_no_compress = PETSC_FALSE

  ! debugging
  PetscInt, public :: immis_ni_count
  PetscInt, public :: immis_ts_cut_count
  PetscInt, public :: immis_ts_count
  PetscInt, public :: immis_debug_cell_id

  ! process models
  PetscInt, public :: immis_liq_flow_eq = UNINITIALIZED_INTEGER
  PetscInt, public :: immis_heat_tran_eq = UNINITIALIZED_INTEGER
  PetscInt, public :: immis_sol_tran_eq = UNINITIALIZED_INTEGER
  PetscInt, parameter, public :: IMMIS_MAX_DOF = 3

  PetscInt, parameter, public :: IMMIS_COND_WATER_INDEX = 1
  PetscInt, parameter, public :: IMMIS_COND_ENERGY_INDEX = 2
  PetscInt, parameter, public :: IMMIS_COND_SOLUTE_INDEX = 3
  PetscInt, parameter, public :: IMMIS_COND_WATER_AUX_INDEX = 4
  PetscInt, parameter, public :: IMMIS_MAX_INDEX = 4

  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_BOUNDARY = 2

  PetscInt, parameter, public :: IMMIS_LIQ_SAT_WRT_LIQ_PRES = 1
  PetscInt, parameter, public :: IMMIS_LIQ_PRES_WRT_POROS = 2

  PetscInt, parameter, public :: IMMIS_ADJOINT_PERMEABILITY = 1
  PetscInt, parameter, public :: IMMIS_ADJOINT_POROSITY = 2
  PetscBool, public :: immis_calc_adjoint = PETSC_FALSE
  PetscInt, public :: immis_adjoint_parameter = IMMIS_ADJOINT_PERMEABILITY

  type, public :: immiscible_auxvar_type
    PetscReal :: pres ! liquid pressure
    PetscReal :: sat  ! liquid saturation
    PetscReal :: pc   ! capillary pressure
    PetscReal :: kr   ! relative permeability
    PetscReal :: effective_porosity
    PetscReal :: dpor_dp
    PetscReal :: dsat_dp ! derivative of saturation wrt pressure
    PetscReal :: dkr_dp  ! derivative of rel. perm. wrt pressure
    PetscReal :: effective_saturation
    PetscReal :: deffsat_dp
    PetscReal :: temp ! temperature
    PetscReal :: conc ! concentration
    PetscReal :: pert
    PetscReal :: mat_pert(1)
  end type immiscible_auxvar_type

  type, public :: immiscible_parameter_type
    PetscBool :: check_post_converged
    PetscReal, pointer :: tensorial_rel_perm_exponent(:,:)
    PetscReal :: diffusion_coef
  end type immiscible_parameter_type

  type, public :: immiscible_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(immiscible_parameter_type), pointer :: immiscible_parameter
    type(immiscible_auxvar_type), pointer :: auxvars(:,:)
    type(immiscible_auxvar_type), pointer :: auxvars_bc(:)
    type(immiscible_auxvar_type), pointer :: auxvars_ss(:)
    type(material_auxvar_type), pointer :: material_auxvars_pert(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type immiscible_type

  interface ImmiscibleAuxVarDestroy
    module procedure ImmiscibleAuxVarSingleDestroy
    module procedure ImmiscibleAuxVarArray1Destroy
    module procedure ImmiscibleAuxVarArray2Destroy
  end interface ImmiscibleAuxVarDestroy

  interface ImmiscibleOutputAuxVars
    module procedure ImmiscibleOutputAuxVars1
  end interface ImmiscibleOutputAuxVars

  public :: ImmiscibleAuxCreate, &
            ImmiscibleAuxDestroy, &
            ImmiscibleAuxVarCompute, &
            ImmiscibleAuxVarInit, &
            ImmiscibleAuxVarCopy, &
            ImmiscibleAuxVarDestroy, &
            ImmiscibleAuxVarStrip, &
            ImmiscibleAuxVarPerturb, &
            ImmiscibleAuxMapConditionIndices, &
            ImmisciblePrintAuxVars, &
            ImmiscibleOutputAuxVars, &
            ImmiscibleAuxTensorialRelPerm

contains

! ************************************************************************** !

function ImmiscibleAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(immiscible_type), pointer :: ImmiscibleAuxCreate

  type(immiscible_type), pointer :: aux

  nullify(immis_min_pert)

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%material_auxvars_pert)
  nullify(aux%matrix_zeroing)

  allocate(aux%immiscible_parameter)
  aux%immiscible_parameter%check_post_converged = PETSC_FALSE
  nullify(aux%immiscible_parameter%tensorial_rel_perm_exponent)
  aux%immiscible_parameter%diffusion_coef = 0.d0

  ImmiscibleAuxCreate => aux

end function ImmiscibleAuxCreate

! ************************************************************************** !

subroutine ImmiscibleAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar
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
  auxvar%temp = 0.d0
  auxvar%conc = 0.d0

  auxvar%pert = 0.d0
  auxvar%mat_pert = 0.d0

end subroutine ImmiscibleAuxVarInit

! ************************************************************************** !

subroutine ImmiscibleAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar, auxvar2
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
  auxvar2%temp = auxvar%temp
  auxvar2%conc = auxvar%conc
  auxvar2%pert = auxvar%pert

end subroutine ImmiscibleAuxVarCopy

! ************************************************************************** !

subroutine ImmiscibleAuxVarCompute(x,immis_auxvar,global_auxvar, &
                              material_auxvar,characteristic_curves, &
                              natural_id,update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(:)
  type(immiscible_auxvar_type) :: immis_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscBool :: update_porosity
  PetscInt :: natural_id

  PetscBool :: saturated
  PetscReal :: dkr_dsat
  PetscReal :: deffsat_dsat

  if (immis_liq_flow_eq > 0) then
    immis_auxvar%pres = x(immis_liq_flow_eq)
  else
    immis_auxvar%pres = option%flow%reference_pressure
  endif
  if (immis_heat_tran_eq > 0) then
    immis_auxvar%temp = x(immis_heat_tran_eq)
  else
    immis_auxvar%temp = option%flow%reference_temperature
  endif
  if (immis_sol_tran_eq > 0) then
    immis_auxvar%conc = x(immis_sol_tran_eq)
  else
    immis_auxvar%conc = 0.d0
  endif

  ! %porosity should never be used. set to bogus value to catch misuse
  material_auxvar%porosity = -888.d0

  if (update_porosity .and. soil_compressibility_index > 0 .and. &
      immis_auxvar%pres > 0.d0) then
    call MaterialCompressSoil(material_auxvar,immis_auxvar%pres, &
                              immis_auxvar%effective_porosity, &
                              immis_auxvar%dpor_dp)
  else
    immis_auxvar%effective_porosity = material_auxvar%porosity_base
    immis_auxvar%dpor_dp = 0.d0
  endif
  material_auxvar%porosity = immis_auxvar%effective_porosity

  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  immis_auxvar%pc = min(option%flow%reference_pressure - immis_auxvar%pres, &
                        characteristic_curves%saturation_function%pcmax)

  if (immis_auxvar%pc > 0.d0) then
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(immis_auxvar%pc, &
                                          immis_auxvar%sat, &
                                          immis_auxvar%dsat_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (immis_auxvar%dsat_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      dkr_dsat = 0.d0
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(immis_auxvar%sat, &
                                            immis_auxvar%kr, &
                                            dkr_dsat,option)
      immis_auxvar%dkr_dp = immis_auxvar%dsat_dp * dkr_dsat
      if (immis_tensorial_rel_perm) then
        call characteristic_curves%liq_rel_perm_function% &
          EffectiveSaturation(immis_auxvar%sat, &
                              immis_auxvar%effective_saturation, &
                              deffsat_dsat,option)
        immis_auxvar%deffsat_dp = deffsat_dsat * immis_auxvar%dsat_dp
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
    immis_auxvar%pc = 0.d0
    immis_auxvar%sat = 1.d0
    immis_auxvar%kr = 1.d0
    immis_auxvar%dsat_dp = 0.d0
    immis_auxvar%dkr_dp = 0.d0
    immis_auxvar%effective_saturation = 1.d0
    immis_auxvar%deffsat_dp = 0.d0
  endif

  if (option%iflag /= IMMIS_UPDATE_FOR_DERIVATIVE) then
    global_auxvar%sat(1) = immis_auxvar%sat
    if (size(global_auxvar%sat) > 1) then
      global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
    endif
  endif

end subroutine ImmiscibleAuxVarCompute

! ************************************************************************** !

subroutine ImmiscibleAuxVarPerturb(x,immis_auxvar,global_auxvar, &
                              material_auxvar, &
                              material_auxvar_pert, &
                              characteristic_curves, &
                              natural_id, &
                              option)
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  PetscReal :: x(:)
  type(option_type) :: option
  PetscInt :: natural_id
  type(immiscible_auxvar_type) :: immis_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  class(characteristic_curves_type) :: characteristic_curves

  PetscInt :: idof, i
  PetscReal :: x_pert(IMMIS_MAX_DOF), pert

  ! IMMIS_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = IMMIS_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    pert = x(idof)*immis_rel_pert+immis_min_pert(idof)
    immis_auxvar(idof)%pert = pert
    x_pert(1:option%nflowdof) = x
    x_pert(idof) = x(idof) + pert
    call ImmiscibleAuxVarCompute(x_pert,immis_auxvar(idof),global_auxvar, &
                            material_auxvar, &
                            characteristic_curves,natural_id, &
                            PETSC_TRUE,option)
  enddo

  if (immis_calc_adjoint) then
    idof = 1
    call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert(idof),option)
    if (immis_adjoint_parameter == IMMIS_ADJOINT_POROSITY) then
      pert = material_auxvar%porosity_base*immis_rel_pert
      material_auxvar_pert(idof)%porosity_base = &
        material_auxvar_pert(idof)%porosity_base + pert
    else if (immis_adjoint_parameter == IMMIS_ADJOINT_PERMEABILITY) then
      pert = material_auxvar%permeability(1)*immis_rel_pert
      do i = 1, size(material_auxvar%permeability)
        material_auxvar_pert(idof)%permeability(i) = &
          material_auxvar_pert(idof)%permeability(i) + pert
      enddo
    endif
    call ImmiscibleAuxVarCompute(x,immis_auxvar(option%nflowdof+1), &
                            global_auxvar, &
                            material_auxvar_pert(idof), &
                            characteristic_curves,natural_id, &
                            PETSC_TRUE,option)
    immis_auxvar(ZERO_INTEGER)%mat_pert(idof) = pert
  endif

end subroutine ImmiscibleAuxVarPerturb

! ************************************************************************** !

subroutine ImmiscibleAuxTensorialRelPerm(auxvar,tensorial_rel_perm_exponent, &
                                    dist,rel_perm,drel_perm_dp,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Utility_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar
  PetscReal :: tensorial_rel_perm_exponent(3)
  PetscReal :: dist(-1:3)
  PetscReal :: rel_perm
  PetscReal :: drel_perm_dp
  type(option_type) :: option

  PetscReal :: exponent_
  PetscReal :: tensorial_scale

  exponent_ = UtilityTensorToScalar(dist,tensorial_rel_perm_exponent)

  ! remember that the default 0.5 was subtracted from the tensorial value
  ! in ImmiscibleSetup. If 0.5 is specified for the tensorial exponent in the
  ! input file, this value will be 0.
  tensorial_scale = auxvar%effective_saturation**exponent_
  rel_perm = auxvar%kr * tensorial_scale
  drel_perm_dp = auxvar%dkr_dp * tensorial_scale + &
                 exponent_ * rel_perm / &
                 auxvar%effective_saturation * &
                 auxvar%deffsat_dp

end subroutine ImmiscibleAuxTensorialRelPerm

! ************************************************************************** !

subroutine ImmisciblePrintAuxVars(immis_auxvar,global_auxvar,material_auxvar, &
                             natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: immis_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', immis_auxvar%pres
  print *, '      capillary pressure: ', immis_auxvar%pc
  print *, '       liquid saturation: ', immis_auxvar%sat
  print *, '      liquid sat (deriv): ', immis_auxvar%dsat_dp
  print *, '         liquid rel perm: ', immis_auxvar%kr
  print *, ' liquid rel perm (deriv): ', immis_auxvar%dkr_dp
  print *, '      effective porosity: ', immis_auxvar%effective_porosity
  print *, '   eff. porosity (deriv): ', immis_auxvar%dpor_dp
  print *, '--------------------------------------------------------'

end subroutine ImmisciblePrintAuxVars

! ************************************************************************** !

subroutine ImmiscibleOutputAuxVars1(immis_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: immis_auxvar
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
  write(IUNIT_TEMP,*) '         liquid pressure: ', immis_auxvar%pres
  write(IUNIT_TEMP,*) '      capillary pressure: ', immis_auxvar%pc
  write(IUNIT_TEMP,*) '       liquid saturation: ', immis_auxvar%sat
  write(IUNIT_TEMP,*) '      liquid sat (deriv): ', immis_auxvar%dsat_dp
  write(IUNIT_TEMP,*) '         liquid rel perm: ', immis_auxvar%kr
  write(IUNIT_TEMP,*) ' liquid rel perm (deriv): ', immis_auxvar%dkr_dp
  write(IUNIT_TEMP,*) '      effective porosity: ', immis_auxvar%effective_porosity
  write(IUNIT_TEMP,*) '   eff. porosity (deriv): ', immis_auxvar%dpor_dp
  write(IUNIT_TEMP,*) '--------------------------------------------------------'

  close(IUNIT_TEMP)

end subroutine ImmiscibleOutputAuxVars1

! ************************************************************************** !

function ImmiscibleAuxMapConditionIndices(include_water_aux)
  !
  ! Maps indexing of conditions
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module

  implicit none

  PetscBool :: include_water_aux

  PetscInt, pointer :: ImmiscibleAuxMapConditionIndices(:)

  PetscInt, pointer :: mapping(:)

  PetscInt :: temp_int

  allocate(mapping(IMMIS_MAX_INDEX))
  mapping = UNINITIALIZED_INTEGER

  temp_int = 0
  if (immis_liq_flow_eq > 0) then
    temp_int = temp_int + 1
    mapping(IMMIS_COND_WATER_INDEX) = temp_int
    if (include_water_aux) then
      temp_int = temp_int + 1
      mapping(IMMIS_COND_WATER_AUX_INDEX) = temp_int
    endif
  endif
  if (immis_heat_tran_eq > 0) then
    temp_int = temp_int + 1
    mapping(IMMIS_COND_ENERGY_INDEX) = temp_int
  endif
  if (immis_sol_tran_eq > 0) then
    temp_int = temp_int + 1
    mapping(IMMIS_COND_SOLUTE_INDEX) = temp_int
  endif

  ImmiscibleAuxMapConditionIndices => mapping

end function ImmiscibleAuxMapConditionIndices

! ************************************************************************** !

subroutine ImmiscibleAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(immiscible_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call ImmiscibleAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine ImmiscibleAuxVarSingleDestroy

! ************************************************************************** !

subroutine ImmiscibleAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(immiscible_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call ImmiscibleAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ImmiscibleAuxVarArray1Destroy

! ************************************************************************** !

subroutine ImmiscibleAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(immiscible_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call ImmiscibleAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ImmiscibleAuxVarArray2Destroy

! ************************************************************************** !

subroutine ImmiscibleMaterialAuxVarDestroy(auxvars)
  !
  ! Deallocates material auxiliary object for immiscible perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
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

end subroutine ImmiscibleMaterialAuxVarDestroy

! ************************************************************************** !

subroutine ImmiscibleAuxVarStrip(auxvar)
  !
  ! ImmiscibleAuxVarDestroy: Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(immiscible_auxvar_type) :: auxvar

end subroutine ImmiscibleAuxVarStrip

! ************************************************************************** !

subroutine ImmiscibleAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(immiscible_type), pointer :: aux

  call DeallocateArray(immis_min_pert)

  if (.not.associated(aux)) return

  call ImmiscibleAuxVarDestroy(aux%auxvars)
  call ImmiscibleAuxVarDestroy(aux%auxvars_bc)
  call ImmiscibleAuxVarDestroy(aux%auxvars_ss)
  call ImmiscibleMaterialAuxVarDestroy(aux%material_auxvars_pert)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%immiscible_parameter)) then
    call DeallocateArray(aux%immiscible_parameter%tensorial_rel_perm_exponent)
    deallocate(aux%immiscible_parameter)
  endif
  nullify(aux%immiscible_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine ImmiscibleAuxDestroy

end module Immiscible_Aux_module
