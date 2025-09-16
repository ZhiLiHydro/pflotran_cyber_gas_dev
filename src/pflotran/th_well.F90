module TH_Well_module

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module
  use TH_Aux_module
  use PM_Well_class

  implicit none

  private

  PetscReal :: liquid_pressure = 1.d6
  PetscReal :: liquid_density_kmol = 1.d3 / FMWH2O
  PetscReal :: thermal_diffusivity = 0.143d-6
  PetscReal :: volume = 1.d0
  PetscReal :: velocity = 1.1574074d-5

  public :: THWellSetup, &
            THWellUpdateAuxVars, &
            THWellPerturb, &
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
  type(th_well_auxvar_type), pointer :: th_well_auxvars(:)
  PetscReal :: dist_up, dist_dn
  PetscInt :: i, index, iconn
  PetscInt :: nx, ny, nz
  PetscInt :: offset_i, offset_j, offset_k
  PetscInt :: num_connections
  PetscInt :: num_well_cells

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  nx = grid%structured_grid%nx
  ny = grid%structured_grid%ny
  nz = grid%structured_grid%nz
  num_well_cells = nx
  num_connections = num_well_cells-1

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

  allocate(pm_well%well_cells(num_well_cells))
  do i = 1, num_well_cells
    index = nx*ny*offset_k + nx*offset_j + offset_i + i
    pm_well%well_cells(i) = index
  enddo

  ! allocate auxvar data structures for all grid cells
  allocate(th_well_auxvars(num_well_cells))
  do i = 1, num_well_cells
    call THWellAuxVarInit(th_well_auxvars(i),option, &
                          th_numerical_derivatives)
  enddo
  patch%aux%TH%auxvars_well => th_well_auxvars
  patch%aux%TH%num_well_aux = num_well_cells
  allocate(th_well_auxvars(2))
  do i = 1, 2
    call THWellAuxVarInit(th_well_auxvars(i),option,PETSC_FALSE)
  enddo
  patch%aux%TH%auxvars_well_bc => th_well_auxvars
  patch%aux%TH%num_well_aux_bc = 2

end subroutine THWellSetup

! ************************************************************************** !

subroutine THWellUpdateAuxVars(xx,pm_well_base)
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

  PetscReal :: xx(:) ! ghosted solution vec
  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: i, iconn, ibc
  PetscInt :: istart, iend
  PetscInt :: num_well_cells
  PetscReal :: xxbc(3)

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  num_well_cells = size(pm_well%well_cells)
  do i = 1, num_well_cells
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    iend = ghosted_id * option%nflowdof
    istart = iend - option%nflowdof + 1
    call THWellAuxVarCompute(liquid_pressure,xx(istart:iend), &
                             auxvars_well(ghosted_id), &
                             natural_id,option)
  enddo

  xxbc = UNINITIALIZED_DOUBLE
  do ibc = 1, 2
    select case(ibc)
      case(1)
        xxbc(th_well_dof) = pm_well%inlet_temperature
      case(2)
        iconn = pm_well%connection_set%num_connections
        ! ghosted id of downwind cells on last connection
        ghosted_id = pm_well%connection_set%id_dn(iconn)
        ! to guarantee zero conduction at boundary
        xxbc(th_well_dof) = auxvars_well(ghosted_id)%temp
    end select
    call THWellAuxVarCompute(liquid_pressure,xxbc, &
                             auxvars_well_bc(ibc), &
                             UNINITIALIZED_INTEGER,option)
  enddo

end subroutine THWellUpdateAuxVars

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
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  type(connection_set_type), pointer :: connection_set
  PetscReal :: Res(3)
  PetscReal :: Jac(3,3)
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, iconn, ibc, offset
  PetscReal :: dist, area
  PetscReal :: tup, tdn
  PetscReal :: enthalpy
  PetscReal :: q
  PetscReal :: convective_energy_flux, diffusive_energy_flux, total_energy_flux
  PetscReal :: residual_scaling_factor
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

  do i = 1, size(pm_well%well_cells)
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    call THWellAccumulation(auxvars_well(ghosted_id),Res,Jac)
    if (calculate_derivatives) then
      call PUMSetValuesBlockedLocal(J,1,ghosted_id,1,ghosted_id, &
                                    Jac,ADD_VALUES,ierr)
    else
      offset = local_id*option%nflowdof
      r_p(offset) = r_p(offset) - Res(th_well_dof)
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
    convective_energy_flux = q * liquid_density_kmol * enthalpy

    ! conduction
    diffusive_energy_flux = thermal_diffusivity * area * (tup-tdn) / dist

    total_energy_flux = convective_energy_flux + diffusive_energy_flux

    if (calculate_derivatives) then
      call PUMSetValuesBlockedLocal(J,1,ghosted_id,1,ghosted_id, &
                                    Jac,ADD_VALUES,ierr)
      call PUMSetValuesBlockedLocal(J,1,ghosted_id,1,ghosted_id, &
                                    Jac,ADD_VALUES,ierr)
    else
      offset = local_id_up*option%nflowdof
      r_p(offset) = r_p(offset) - total_energy_flux
      offset = local_id_dn*option%nflowdof
      r_p(offset) = r_p(offset)  + total_energy_flux
    endif

  enddo

  do ibc = 1, 2
    select case(ibc)
      case(1)
        iconn = 1
        ! the upstream cell on first internal connection is downstream on bc
        ghosted_id_dn = connection_set%id_up(iconn)
        tup = auxvars_well_bc(ibc)%temp
        tdn = auxvars_well(ghosted_id_dn)%temp
        enthalpy = auxvars_well_bc(ibc)%temp
        residual_scaling_factor = 1.d0
      case(2)
        ! the downstram cell on last internal connection is downstream on bc
        iconn = connection_set%num_connections
        ghosted_id_dn = connection_set%id_dn(iconn)
        tup = auxvars_well(ghosted_id_dn)%temp
        tdn = tup ! to guarantee no conductive heat flux
!        tdn = auxvars_well_bc(ibc)%temp
        enthalpy = auxvars_well(ghosted_id_dn)%temp
        residual_scaling_factor = -1.d0
    end select
    local_id_dn = grid%nG2L(ghosted_id_dn)
    area = connection_set%area(iconn)
    ! boundary fluxes
    q = velocity * area
    convective_energy_flux = q * liquid_density_kmol * enthalpy
    ! conduction
    diffusive_energy_flux = thermal_diffusivity * area * (tup-tdn) / dist

    total_energy_flux = convective_energy_flux + diffusive_energy_flux
    offset = local_id_dn*option%nflowdof
    r_p(offset) = r_p(offset) - residual_scaling_factor * total_energy_flux
  enddo

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWell

! ************************************************************************** !

subroutine THWellDerivative(r,J,pm_well_base)
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

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  type(connection_set_type), pointer :: connection_set
  PetscReal :: Res(3)
  PetscReal :: Jac(3,3)
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, iconn, ibc, offset
  PetscReal :: dist, area
  PetscReal :: tup, tdn
  PetscReal :: enthalpy
  PetscReal :: q
  PetscReal :: convective_energy_flux, diffusive_energy_flux, total_energy_flux
  PetscReal :: residual_scaling_factor
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  if (th_numerical_derivatives) then
  else
  endif

end subroutine THWellDerivative

! ************************************************************************** !

subroutine THWellPerturb(pm_well_base)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: i
  PetscInt :: num_well_cells

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well

  num_well_cells = size(pm_well%well_cells)
  do i = 1, num_well_cells
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    call THWellAuxVarPerturb(liquid_pressure,auxvars_well(ghosted_id), &
                             natural_id,option)
  enddo

end subroutine THWellPerturb

! ************************************************************************** !

subroutine THWellAccumulation(auxvar_well,Res,Jac)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscReal :: Res(:)
  PetscReal :: Jac(:,:)

  PetscReal :: tempreal

  ! Res [MJ] = liquid_density_kmol [kmol liquid/m^3 liquid] *
  !            internal energy [MJ/kmol liquid] *
  !            volume [m^3 liquid]
  Res = 0.d0
  Jac = 0.d0
  tempreal = volume * liquid_density_kmol
  Res(1) = tempreal * auxvar_well%u
  Jac(3,3) = tempreal * auxvar_well%du_dT

end subroutine THWellAccumulation

! ************************************************************************** !

subroutine THWellFlux(auxvar_well_up,auxvar_welL_dn, &
                      dist,Res,Jup,Jdn,calculate_derivatives)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscReal :: Res(:)
  PetscReal :: Jup(:,:)
  PetscReal :: Jdn(:,:)
  PetscReal :: dist
  PetscBool :: calculate_derivatives

  PetscReal :: tup, tdn
  PetscReal :: dh_dTup, dh_dTdn
  PetscReal :: enthalpy, velocity, area, q_kmol
  PetscReal :: diffusive_energy_flux, convective_energy_flux
  PetscReal :: tempreal

  Res = 0.d0

  tup = auxvar_well_up%temp
  tdn = auxvar_well_dn%temp

  ! convection
  dh_dTup = 0.d0
  dh_dTdn = 0.d0
  if (velocity > 0) then
    enthalpy = auxvar_well_up%h
    dh_dTup = auxvar_well_up%dh_dT
  else
    enthalpy = auxvar_well_dn%h
    dh_dTdn = auxvar_well_dn%dh_dT
  endif

  q_kmol = velocity * area * liquid_density_kmol
  convective_energy_flux = q_kmol * enthalpy

  ! conduction
  ! [MJ/s] = [m^2 bulk/s] * [m^2 bulk] [C] / [m bulk]
  tempreal = thermal_diffusivity * area
  diffusive_energy_flux = tempreal * (tup-tdn)

  Res(th_well_dof) = convective_energy_flux + diffusive_energy_flux

  if (calculate_derivatives) then
    Jup = 0.d0
    Jup(3,3) = q_kmol * dh_dTup + tempreal
    Jdn = 0.d0
    Jdn(3,3) = q_kmol * dh_dTdn - tempreal
  endif

end subroutine THWellFlux

end module TH_Well_module

