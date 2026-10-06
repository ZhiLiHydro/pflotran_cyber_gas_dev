module TS_Connections_module
  !
  ! The patch's internal, boundary and source/sink connections flattened
  ! into arrays, in list order: the geometry a TS model indexes its fluxes
  ! by. Built from the same lists, the tables of the models that share a
  ! TS (RICHARDS_DAE with GIRT_DAE) agree by construction.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: ts_connections_type
    PetscInt :: nconn_int = 0
    PetscInt :: nconn_bc = 0
    PetscInt :: nconn_ss = 0
    PetscInt, allocatable :: int_up(:)       ! (nconn_int) ghosted id
    PetscInt, allocatable :: int_dn(:)       ! (nconn_int) ghosted id
    PetscReal, allocatable :: int_area(:)    ! (nconn_int) [m^2]
    PetscReal, allocatable :: int_dist(:,:)  ! (-1:3,nconn_int)
    PetscInt, allocatable :: bc_local(:)     ! (nconn_bc) local id
    PetscInt, allocatable :: bc_ghost(:)     ! (nconn_bc) ghosted id
    PetscReal, allocatable :: bc_area(:)     ! (nconn_bc) [m^2]
    PetscReal, allocatable :: bc_dist(:,:)   ! (-1:3,nconn_bc)
    PetscInt, allocatable :: ss_local(:)     ! (nconn_ss) local id
    PetscInt, allocatable :: ss_ghost(:)     ! (nconn_ss) ghosted id
  end type ts_connections_type

  public :: TSConnectionsBuild

contains

! ************************************************************************** !

subroutine TSConnectionsBuild(conn,patch)
  !
  ! Flattens the patch's connection lists into conn.
  !
  ! Author: Hansol Suh
  ! Date: 10/05/26
  !
  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module

  implicit none

  type(ts_connections_type) :: conn
  type(patch_type), pointer :: patch

  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: coupler
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, n

  grid => patch%grid

  conn%nconn_int = 0
  cur_connection_set => grid%internal_connection_set_list%first
  do while (associated(cur_connection_set))
    conn%nconn_int = conn%nconn_int + cur_connection_set%num_connections
    cur_connection_set => cur_connection_set%next
  enddo
  conn%nconn_bc = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  conn%nconn_ss = CouplerGetNumConnectionsInList(patch%source_sink_list)

  allocate(conn%int_up(conn%nconn_int))
  allocate(conn%int_dn(conn%nconn_int))
  allocate(conn%int_area(conn%nconn_int))
  allocate(conn%int_dist(-1:3,conn%nconn_int))
  n = 0
  cur_connection_set => grid%internal_connection_set_list%first
  do while (associated(cur_connection_set))
    do iconn = 1, cur_connection_set%num_connections
      n = n + 1
      conn%int_up(n) = cur_connection_set%id_up(iconn)
      conn%int_dn(n) = cur_connection_set%id_dn(iconn)
      conn%int_area(n) = cur_connection_set%area(iconn)
      conn%int_dist(:,n) = cur_connection_set%dist(:,iconn)
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  allocate(conn%bc_local(conn%nconn_bc))
  allocate(conn%bc_ghost(conn%nconn_bc))
  allocate(conn%bc_area(conn%nconn_bc))
  allocate(conn%bc_dist(-1:3,conn%nconn_bc))
  n = 0
  coupler => patch%boundary_condition_list%first
  do while (associated(coupler))
    cur_connection_set => coupler%connection_set
    do iconn = 1, cur_connection_set%num_connections
      n = n + 1
      conn%bc_local(n) = cur_connection_set%id_dn(iconn)
      conn%bc_ghost(n) = grid%nL2G(conn%bc_local(n))
      conn%bc_area(n) = cur_connection_set%area(iconn)
      conn%bc_dist(:,n) = cur_connection_set%dist(:,iconn)
    enddo
    coupler => coupler%next
  enddo

  allocate(conn%ss_local(conn%nconn_ss))
  allocate(conn%ss_ghost(conn%nconn_ss))
  n = 0
  coupler => patch%source_sink_list%first
  do while (associated(coupler))
    cur_connection_set => coupler%connection_set
    do iconn = 1, cur_connection_set%num_connections
      n = n + 1
      conn%ss_local(n) = cur_connection_set%id_dn(iconn)
      conn%ss_ghost(n) = grid%nL2G(conn%ss_local(n))
    enddo
    coupler => coupler%next
  enddo

end subroutine TSConnectionsBuild

end module TS_Connections_module
