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
            THWell, &
            THWellAccumulationTerms, &
            THWellMatrixZeroing

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
  use Matrix_Zeroing_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: connection_set
  type(th_well_auxvar_type), pointer :: th_well_auxvars(:)
  PetscReal :: dist_up, dist_dn
  PetscInt :: local_id, ghosted_id
  PetscInt :: i, iconn, idof
  PetscInt :: nx, ny, nz
  PetscInt :: offset_i, offset_j, offset_k
  PetscInt :: num_connections
  PetscInt :: num_well_cells
  PetscInt, allocatable :: inactive_cells(:)

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
  offset_j = int(ny/2)
  offset_k = int(nz/2)
  do iconn = 1, num_connections
    ghosted_id = nx*ny*offset_k + nx*offset_j + offset_i + iconn
    connection_set%id_up(iconn) = ghosted_id
    connection_set%id_dn(iconn) = ghosted_id + 1
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
    local_id = nx*ny*offset_k + nx*offset_j + offset_i + i
    pm_well%well_cells(i) = local_id
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

  ! generate matrix zeroing
  pm_well%matrix_zeroing => MatrixZeroingCreate()
  allocate(inactive_cells(grid%nlmax))
  inactive_cells = 1
  do i = 1, size(pm_well%well_cells)
    inactive_cells(pm_well%well_cells(i)) = 0
  enddo
  i = 0
  do local_id = 1, grid%nlmax
    if (inactive_cells(local_id) == 1) then
      i = i + 1
      inactive_cells(i) = i
    endif
  enddo
  call MatrixZeroingAllocateArray(pm_well%matrix_zeroing,i,option)
  do i = 1, pm_well%matrix_zeroing%n_zero_rows
    local_id = inactive_cells(i)
    ghosted_id = grid%nL2G(local_id)
    idof = (local_id-1)*option%nflowdof + th_well_dof
    pm_well%matrix_zeroing%zero_rows_local(i) = idof ! one-based
    idof = (ghosted_id-1)*option%nflowdof + th_well_dof
    pm_well%matrix_zeroing%zero_rows_local_ghosted(i) = idof-1 ! zero-based
  enddo
  deallocate(inactive_cells)

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

subroutine THWellAccumulationTerms(r,J,pm_well_base,calculate_derivatives)
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
  PetscReal :: Res(3)
  PetscReal :: Jblock(3,3)
  PetscInt :: ndof
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: i
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  ndof = option%nflowdof
  auxvars_well => patch%aux%TH%auxvars_well

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

  do i = 1, size(pm_well%well_cells)
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    call THWellAccumulation(auxvars_well(ghosted_id),Res,Jblock,ndof)
    Res = Res/option%flow_dt
    if (calculate_derivatives) then
      call PetUtilMatSVBL(J,ghosted_id,ghosted_id,Jblock,ndof)
    else
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
    endif
  enddo

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWellAccumulationTerms

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
  type(th_well_auxvar_type), pointer :: auxvar_up, auxvar_dn
  PetscReal :: Res(3)
  PetscReal :: Jup(3,3), Jdn(3,3)
  PetscInt :: ndof
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn, ibc
  PetscReal :: residual_scaling_factor
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  ndof = option%nflowdof
  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  call THWellAccumulationTerms(r,J,pm_well_base,calculate_derivatives)

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

  ! calculate heat transfer between well cells
  connection_set => pm_well%connection_set
  do iconn = 1, connection_set%num_connections
    ! cells are ghosted in connection set
    ghosted_id_up = connection_set%id_up(iconn)
    ghosted_id_dn = connection_set%id_dn(iconn)
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    call THWellFlux(auxvars_well(ghosted_id_up), &
                    auxvars_well(ghosted_id_dn), &
                    connection_set%dist(0,iconn), &
                    connection_set%area(iconn), &
                    Res,Jup,Jdn,ndof,calculate_derivatives)
    if (calculate_derivatives) then
      if (local_id_up > 0) then
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_up,Jup,ndof)
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_dn,Jdn,ndof)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_up,Jup,ndof)
      endif
    else
      if (local_id_up > 0) then
        call PetUtilVecSVBL(r_p,local_id_up,Res,ndof,PETSC_FALSE)
      endif
      if (local_id_dn > 0) then
        Res = -Res
        call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
      endif
    endif

  enddo

  do ibc = 1, 2
    select case(ibc)
      case(1)
        iconn = 1
        ! the upstream cell on first internal connection is downstream on bc
        ghosted_id_dn = connection_set%id_up(iconn)
        local_id_up = -1
        local_id_dn = grid%nG2L(ghosted_id_dn)
        auxvar_up => auxvars_well_bc(ibc)
        auxvar_dn => auxvars_well(ghosted_id_dn)
        residual_scaling_factor = 1.d0
      case(2)
        ! the downstream cell on last internal connection is upstream on bc
        iconn = connection_set%num_connections
        ghosted_id_dn = connection_set%id_dn(iconn)
        local_id_up = grid%nG2L(ghosted_id_dn)
        local_id_dn = -1
        auxvar_up => auxvars_well(ghosted_id_dn)
        ! to guarantee no conductive heat flux
        auxvar_dn => auxvars_well(ghosted_id_dn)
        residual_scaling_factor = -1.d0
    end select
    call THWellFlux(auxvars_well(ghosted_id_up), &
                    auxvars_well(ghosted_id_dn), &
                    connection_set%dist(0,iconn), &
                    connection_set%area(iconn), &
                    Res,Jup,Jdn,ndof,calculate_derivatives)
    if (calculate_derivatives) then
      if (local_id_up > 0) then
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_up,Jup,ndof)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
      endif
    else
      if (local_id_up > 0) then
        call PetUtilVecSVBL(r_p,local_id_up,Res,ndof,PETSC_FALSE)
      endif
      if (local_id_dn > 0) then
        Res = -Res
        call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
      endif
    endif
  enddo

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWell

! ************************************************************************** !

subroutine THWellMatrixZeroing(pm_well_base,r,J,calculate_derivatives)
  !
  ! Zeros dof for non-well grid cells
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/25
  !
  use Matrix_Zeroing_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base
  Vec :: r
  Mat :: J
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  if (calculate_derivatives) then
    call MatrixZeroingZeroMatEntries(pm_well%matrix_zeroing,J)
  else
    call MatrixZeroingZeroVecEntries(pm_well%matrix_zeroing,r)
  endif

end subroutine THWellMatrixZeroing

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

subroutine THWellAccumulation(auxvar_well,Res,Jac,ndof)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jac(ndof,ndof)

  PetscReal :: tempreal

  ! Res [MJ] = liquid_density_kmol [kmol liquid/m^3 liquid] *
  !            internal energy [MJ/kmol liquid] *
  !            volume [m^3 liquid]
  Res = 0.d0
  Jac = 0.d0
  tempreal = volume * liquid_density_kmol
  Res(ndof) = tempreal * auxvar_well%u
  Jac(ndof,ndof) = tempreal * auxvar_well%du_dT

end subroutine THWellAccumulation

! ************************************************************************** !

subroutine THWellAccumulationDerivative(auxvar_well,Res,Jac,ndof, &
                                        numerical_derivatives)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jac(ndof,ndof)
  PetscBool :: numerical_derivatives

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof)
  PetscInt :: idof, ieq

  call THWellAccumulation(auxvar_well,Res,Jac,ndof)
  if (numerical_derivatives) then
    do idof = 1, ndof
      call THWellAccumulation(auxvar_well%auxvar_pert(idof), &
                              Res_pert,Jdum,ndof)
      do ieq = 1, ndof
        Jac(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well%auxvar_pert(idof)%pert
      enddo
    enddo
  endif

end subroutine THWellAccumulationDerivative

! ************************************************************************** !

subroutine THWellFlux(auxvar_well_up,auxvar_well_dn, &
                      dist,area,Res,Jup,Jdn,ndof, &
                      calculate_derivatives)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jup(ndof,ndof)
  PetscReal :: Jdn(ndof,ndof)
  PetscReal :: dist
  PetscReal :: area
  PetscBool :: calculate_derivatives

  PetscReal :: tup, tdn
  PetscReal :: dh_dTup, dh_dTdn
  PetscReal :: enthalpy, q_kmol
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
    Jup(ndof,ndof) = q_kmol * dh_dTup + tempreal
    Jdn = 0.d0
    Jdn(ndof,ndof) = q_kmol * dh_dTdn - tempreal
  endif

end subroutine THWellFlux

! ************************************************************************** !

subroutine THWellFluxDerivative(auxvar_well_up,auxvar_well_dn, &
                                dist,area,Res,Jup,Jdn,ndof, &
                                numerical_derivatives)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jup(ndof,ndof)
  PetscReal :: Jdn(ndof,ndof)
  PetscReal :: dist
  PetscReal :: area
  PetscBool :: numerical_derivatives

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof), Jdum1(ndof,ndof)
  PetscInt :: idof, ieq

  call THWellFlux(auxvar_well_up,auxvar_well_dn, &
                  dist,area,Res,Jup,Jdn,ndof,.not.numerical_derivatives)
  if (numerical_derivatives) then
    Jup = 0.d0
    do idof = 1, ndof
      call THWellFlux(auxvar_well_up%auxvar_pert(idof),auxvar_well_dn, &
                      dist,area,Res_pert,Jdum,Jdum1,ndof,PETSC_FALSE)
      do ieq = 1, ndof
        Jup(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well_up%auxvar_pert(idof)%pert
      enddo
    enddo
    Jdn = 0.d0
    do idof = 1, ndof
      call THWellFlux(auxvar_well_up,auxvar_well_dn%auxvar_pert(idof), &
                      dist,area,Res_pert,Jdum,Jdum1,ndof,PETSC_FALSE)
      do ieq = 1, ndof
        Jdn(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well_dn%auxvar_pert(idof)%pert
      enddo
    enddo
  endif

end subroutine THWellFluxDerivative

end module TH_Well_module

