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
  PetscReal :: velocity = 1.1574074d-5

  public :: THWellSetup, &
            THWell

contains

! ************************************************************************** !

subroutine THWellSetup(pm_well_base)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: connection_set
  PetscReal :: dist_up, dist_dn
  PetscInt :: i, index, iconn
  PetscInt :: nx, ny, nz
  PetscInt :: offset_i, offset_j, offset_k
  PetscInt :: num_connections

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  nx = grid%structured_grid%nx
  ny = grid%structured_grid%ny
  nz = grid%structured_grid%nz
  num_connections = nx-1

  connection_set => &
    ConnectionCreate(num_connections,INTERNAL_FACE_CONNECTION_TYPE,NULL_GRID)
        allocate(connection_set%id_up(num_connections))
      allocate(connection_set%id_dn(num_connections))
      allocate(connection_set%dist(-1:3,num_connections))
      allocate(connection_set%area(num_connections))
      connection_set%id_up = 0
      connection_set%id_dn = 0
      connection_set%dist = 0.d0
      connection_set%area = 0.d0
  offset_i = 0
  offset_j = 1
  offset_k = 1
  do iconn = 1, num_connections
    index = nx*ny*offset_k + nx*offset_j + offset_i + iconn
    connection_set%id_up(iconn) = index
    connection_set%id_dn(iconn) = index + 1
    connection_set%dist(:,iconn) = 0.d0
    dist_up = 0.5d0*grid%structured_grid%dx(offset_i+1)
    dist_dn = 0.5d0*grid%structured_grid%dx(offset_i+2)
    connection_set%dist(0,iconn) = dist_up + dist_dn
    connection_set%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
    connection_set%dist(1,iconn) = 1.d0
    connection_set%area(iconn) = grid%structured_grid%dy(offset_j+1)* &
                                 grid%structured_grid%dz(offset_k+1)
  enddo
  pm_well%connection_set => connection_set
  allocate(pm_well%well_cells(nx))
  do i = 1, nx
    index = nx*ny*offset_k + nx*offset_j + offset_i + i
    pm_well%well_cells(i) = index
  enddo

end subroutine THWellSetup

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
  PetscInt :: i, iconn
  PetscReal :: dist, area
  PetscReal :: tup, tdn
  PetscReal :: enthalpy
  PetscReal :: q
  PetscReal :: convective_energy_flux, diffusive_energy_flux, total_energy_flux
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

  ! calculate heat transfer between well cells
  connection_set => pm_well%connection_set
  do iconn = 1, connection_set%num_connections
    ! cells are ghosted in connection set
    ghosted_id_up = connection_set%id_up(iconn)
    ghosted_id_dn = connection_set%id_dn(iconn)
    dist = connection_set%dist(0,iconn)
    area = connection_set%area(iconn)
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    tup = auxvars_well(ghosted_id_up)%temp
    tdn = auxvars_well(ghosted_id_dn)%temp

    ! convection
    if (velocity > 0) then
      enthalpy = auxvars_well(ghosted_id_up)%h
    else
      enthalpy = auxvars_well(ghosted_id_dn)%h
    endif

    q = velocity * area
    convective_energy_flux = q * liquid_density * enthalpy

    ! conduction
    diffusive_energy_flux = thermal_diffusivity * area * (tup-tdn) / dist

    total_energy_flux = convective_energy_flux + diffusive_energy_flux
    i = local_id_up*option%nflowdof
    r_p(i) = r_p(i) - total_energy_flux
    i = local_id_dn*option%nflowdof
    r_p(i) = r_p(i)  + total_energy_flux

  enddo

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

