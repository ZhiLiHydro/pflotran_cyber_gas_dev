module Immiscible_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Material_Aux_module

  implicit none

  private

  PetscReal, public :: immis_density_kg = 1000.d0
  PetscReal, public :: immis_viscosity = 1.d-3
  PetscReal, public :: immis_fmw(2) = [FMWH2O,FMWAIR]

  PetscReal, public :: immis_rel_pert = 1.d-8
  PetscReal, public :: immis_pres_min_pert = 1.d-2
  PetscReal, public :: immis_sat_min_pert = 1.d-10

  PetscReal, pointer, public :: immis_min_pert(:)

  PetscBool, public :: immis_calc_accum = PETSC_TRUE
  PetscBool, public :: immis_calc_flux = PETSC_TRUE
  PetscBool, public :: immis_calc_bcflux = PETSC_TRUE
  PetscBool, public :: immis_simultaneous_res_jac_calc = PETSC_TRUE

  PetscBool, public :: immis_numerical_derivatives = PETSC_FALSE

  ! debugging
  PetscInt, public :: immis_ni_count
  PetscInt, public :: immis_ts_cut_count
  PetscInt, public :: immis_ts_count
  PetscInt, public :: immis_debug_cell_id

  ! process models
  PetscInt, parameter, public :: IMMIS_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: IMMIS_GAS_SATURATION_DOF = 2
  PetscInt, parameter, public :: IMMIS_MAX_DOF = IMMIS_GAS_SATURATION_DOF

  PetscInt, parameter, public :: IMMIS_COND_LIQUID_INDEX = 1
  PetscInt, parameter, public :: IMMIS_COND_GAS_INDEX = 2
  PetscInt, parameter, public :: IMMIS_MAX_INDEX = 2

  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_BOUNDARY = 2

  PetscInt, parameter, public :: IMMIS_GAS_SAT_WRT_GAS_PRES = 1
  PetscInt, parameter, public :: IMMIS_GAS_PRES_WRT_POROS = 2

  type, public :: immiscible_auxvar_type
    PetscReal :: pres(3)
    PetscReal :: sat(2)
    PetscReal :: den_kg(2)
    PetscReal :: mu(2)
    PetscReal :: kr(2) ! relative permeability
    PetscReal :: pert
  end type immiscible_auxvar_type

  type, public :: immiscible_parameter_type
    PetscBool :: check_post_converged
  end type immiscible_parameter_type

  type, public :: immiscible_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(immiscible_parameter_type), pointer :: immiscible_parameter
    type(immiscible_auxvar_type), pointer :: auxvars(:,:)
    type(immiscible_auxvar_type), pointer :: auxvars_bc(:)
    type(immiscible_auxvar_type), pointer :: auxvars_ss(:)
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
            ImmisciblePrintAuxVars, &
            ImmiscibleOutputAuxVars

contains

! ************************************************************************** !

function ImmiscibleAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

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
  nullify(aux%matrix_zeroing)

  allocate(aux%immiscible_parameter)
  aux%immiscible_parameter%check_post_converged = PETSC_FALSE

  ImmiscibleAuxCreate => aux

end function ImmiscibleAuxCreate

! ************************************************************************** !

subroutine ImmiscibleAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mu = 0.d0
  auxvar%kr = 0.d0

  auxvar%pert = 0.d0

end subroutine ImmiscibleAuxVarInit

! ************************************************************************** !

subroutine ImmiscibleAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%mu = auxvar%mu
  auxvar2%kr = auxvar%kr

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

  PetscReal :: dkr_dsat

  immis_auxvar%pres = x(IMMIS_GAS_PRESSURE_DOF)
  immis_auxvar%sat = x(IMMIS_GAS_SATURATION_DOF)

#if 0
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
                                            immis_auxvar%kr(1), &
                                            dkr_dsat,option)
      immis_auxvar%dkr_dp = immis_auxvar%dsat_dp * dkr_dsat
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
  endif
#endif

end subroutine ImmiscibleAuxVarCompute

! ************************************************************************** !

subroutine ImmiscibleAuxVarPerturb(x,immis_auxvar,global_auxvar, &
                              material_auxvar, &
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
  class(characteristic_curves_type) :: characteristic_curves

  PetscInt :: idof
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

end subroutine ImmiscibleAuxVarPerturb

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

#if 0
  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', immis_auxvar%pres
  print *, '      capillary pressure: ', immis_auxvar%pc
  print *, '       liquid saturation: ', immis_auxvar%sat
  print *, '      liquid sat (deriv): ', immis_auxvar%dsat_dp
  print *, '         liquid rel perm: ', immis_auxvar%kr
  print *, ' liquid rel perm (deriv): ', immis_auxvar%dkr_dp
  print *, '--------------------------------------------------------'
#endif
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
#if 0
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
  write(IUNIT_TEMP,*) '--------------------------------------------------------'

  close(IUNIT_TEMP)
#endif
end subroutine ImmiscibleOutputAuxVars1

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

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%immiscible_parameter)) then
    deallocate(aux%immiscible_parameter)
  endif
  nullify(aux%immiscible_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine ImmiscibleAuxDestroy

end module Immiscible_Aux_module
