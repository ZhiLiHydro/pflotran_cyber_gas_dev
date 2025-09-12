module TH_Well_module

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module
  use TH_Aux_module
  use PM_Well_class

  implicit none

  private

  PetscReal :: liquid_density = 1.d3
  PetscReal :: thermal_diffusivity = 0.143d-6
  PetscReal :: volume = 1.d0

  public :: THWell

contains

! ************************************************************************** !

subroutine THWell(r,J,pm_well_base,calculate_derivatives)
  !
  ! Calculates contribution of well model to the TH residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !
  use Option_module
  use Petsc_Utility_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  Vec :: r
  Mat :: J
  class(pm_well_type), pointer :: pm_well_base
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(connection_set_type), pointer :: connection_set
  PetscReal :: Res(3)
  PetscReal :: Jac(3,3)
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

  do local_id = 1, size(pm_well%well_cells)
    ghosted_id = grid%nL2G(local_id)
    call THWellAccumulation(auxvars_well(ghosted_id),Res,Jac)
    if (calculate_derivatives) then
      call PUMSetValuesBlockedLocal(J,1,ghosted_id,1,ghosted_id, &
                                    Jac,ADD_VALUES,ierr)
    else
      i = local_id*option%nflowdof
      r_p(i) = r_p(i) - Res(th_well_dof)
    endif
  enddo

#if 0
  ! calculate heat transfer between well cells
  connection_set => this%connection_set
  do iconn = 1, connection_set%num_connections
    ! cells are ghosted in connection set
    ghosted_id_up = connection_set%id_up(iconn)
    ghosted_id_dn = connection_set%id_dn(iconn)
    dist = connection_set%dist(0,iconn)
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    tup = th_auxvar(ghosted_id_up)%pipe%temp
    tdn = th_auxvar(ghosted_id_dn)%pipe%temp

    ! convection
    if (velocity > 0) then
      enthalpy = th_auxvar(ghosted_id_up)%pipe%h
    else
      enthalpy = th_auxvar(ghosted_id_dn)%pipe%h
    endif

    q = velocity * area
    convective_energy_flux = q * density * enthalpy

    ! conduction
    diffusive_energy_flux = thermal_diffusivity * area * (tup-tdn) / dist

    totol_energy_flux = convective_energy_flux + diffusive_energy_flux
    residual(local_id_up) = residual(local_id_up) - total_energy_flux
    residual(local_id_up) = residual(local_id_up) + total_energy_flux

  enddo
#endif

  ! calculate heat transfer between well cells and subsurface cells

  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine THWell

! ************************************************************************** !

subroutine THWellAccumulation(auxvar_well,Res,Jac)
  !
  ! Calculates contribution of well model to the TH residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscReal :: Res(:)
  PetscReal :: Jac(:,:)

  PetscReal :: tempreal

  Res = 0.d0
  Jac = 0.d0
  tempreal = volume * liquid_density
  Res(1) = tempreal * auxvar_well%u
  Jac(3,3) = tempreal * auxvar_well%du_dT

end subroutine THWellAccumulation

end module TH_Well_module

