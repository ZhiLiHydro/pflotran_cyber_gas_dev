module Stokes_module

#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
  use petscsnes
  use petscsys
  use petscmat
  use PFLOTRAN_Constants_module
  use Global_Aux_module
  use Material_Aux_module
  use Utility_module, only : Equal
  use Stokes_Aux_module

  implicit none

  private

  public :: StokesSetup, &
            StokesInitializeTimestep, &
            StokesUpdateAuxVars, &
            StokesUpdateSolution, &
            StokesMaxChange, &
            StokesDestroy, &
            StokesResidual, &
            StokesJacobian, &
            StokesResidualInternalConn_structured, &
            StokesResidualBoundaryConn, &
            StokesResidualAccumulation, &
            StokesResidualSourceSink, &
            StokesJacobianInternalConn_structured, &
            StokesJacobianBoundaryConn, &
            StokesJacobianAccumulation, &
            StokesJacobianSourceSink

contains

! ************************************************************************** !

subroutine StokesSetup(realization,rho,mu)
  !
  ! Creates arrays for auxiliary variables for Stokes problem
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  type(patch_type),         pointer :: patch
  type(grid_type),          pointer :: grid
  type(material_auxvar_type), pointer :: mat_aux(:)
  type(stokes_type),    pointer :: stokes_aux
  type(option_type), pointer :: option

  PetscInt :: gid
  PetscBool :: dof_is_active(1)
  PetscInt :: iconn, sum_connection

  patch  => realization%patch
  grid   => patch%grid
  mat_aux => patch%aux%Material%auxvars
  option => realization%option

  if (.not.associated(patch%aux%Stokes)) then
    patch%aux%Stokes => StokesAuxCreate(realization%option)
  endif

  stokes_aux => patch%aux%Stokes
  if (.not.associated(stokes_aux%auxvars)) then
    allocate(stokes_aux%auxvars(grid%ngmax))
    do gid = 1, grid%ngmax
      call StokesAuxVarInit(stokes_aux%auxvars(gid),realization%option)
    enddo
    stokes_aux%num_aux = grid%ngmax
  endif

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(stokes_aux%auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call StokesAuxVarInit(stokes_aux%auxvars_bc(iconn),option)
    enddo
  endif
  stokes_aux%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(stokes_aux%auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call StokesAuxVarInit(stokes_aux%auxvars_ss(iconn),option)
    enddo
  endif
  stokes_aux%num_aux_ss = sum_connection

  dof_is_active = PETSC_TRUE
  call PatchCreateZeroArray(patch,dof_is_active, &
                            patch%aux%Stokes%matrix_zeroing,option)



end subroutine StokesSetup

! ************************************************************************** !

subroutine StokesInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  call StokesUpdateFixedAccum(realization)

end subroutine StokesInitializeTimestep

! ************************************************************************** !

subroutine StokesUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with the Stokes problem
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(patch_type),         pointer :: patch
  type(grid_type),          pointer :: grid
  type(option_type),        pointer :: option
  type(field_type),         pointer :: field
  type(stokes_type),        pointer :: stokes_aux
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: Ucell(STOKES_NUM_DOF)
  PetscInt :: gid, lid, istart
  PetscErrorCode :: ierr

  patch   => realization%patch
  grid    => patch%grid
  option  => realization%option
  field   => realization%field
  global_auxvars => patch%aux%Global%auxvars

  if (.not.associated(patch%aux%Stokes)) return

  stokes_aux => patch%aux%Stokes
  call VecGetArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  do gid = 1, grid%ngmax
    lid = grid%nG2L(gid)
    if (lid < 0) cycle
    if (patch%imat(gid) <= 0) cycle
    istart = (lid-1)*option%nflowdof + 1
    Ucell(1) = xx_loc_p(istart    )
    Ucell(2) = xx_loc_p(istart + 1)
    Ucell(3) = xx_loc_p(istart + 2)
    Ucell(4) = xx_loc_p(istart + 3)
    call StokesAuxVarCompute(Ucell, &
                             stokes_aux%auxvars(gid), &
                             global_auxvars(gid), &
                             grid%nG2A(gid), &
                             option)
  enddo

  call VecRestoreArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  stokes_aux%auxvars_up_to_date = PETSC_TRUE

end subroutine StokesUpdateAuxVars

! ************************************************************************** !

subroutine StokesUpdateSolution(realization)
  !
  ! Updates data in module after a successful time step
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

end subroutine StokesUpdateSolution

! ************************************************************************** !

subroutine StokesMaxChange(realization,du_max)
  !
  ! Computes maximum change in solution
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Base_class
  use Field_module

  implicit none

  class(realization_base_type) :: realization
  PetscReal :: du_max
  type(field_type), pointer :: field
  PetscErrorCode :: ierr

  field => realization%field
  du_max = 0.d0
  call VecNorm(field%flow_xx,NORM_INFINITY,du_max,ierr);CHKERRQ(ierr)

end subroutine StokesMaxChange

! ************************************************************************** !

subroutine StokesDestroy(realization)
  !
  ! Deallocates variables associated with Stokes
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

end subroutine StokesDestroy
! ************************************************************************** !

subroutine StokesResidual(snes,xx,r,realization,rho,mu,debug,ierr)
  !
  ! Computes Stokes residual for SNES solver
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec  :: xx, r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  type(debug_type) :: debug
  PetscErrorCode :: ierr
  PetscReal:: r_norm
  PetscMPIInt :: mpi_rank

  call StokesResidualPreliminaries(xx,r,realization,ierr)
  call StokesResidualInternalConn_structured(r,realization,mu,ierr)
  call StokesResidualBoundaryConn(r,realization,rho,mu,ierr)
  call StokesResidualAccumulation(r,realization,rho,ierr)
  call StokesResidualSourceSink(r,realization,rho,ierr)

  ! DIAGNOSTIC: Check residual vector
  call VecNorm(r,NORM_2,r_norm,ierr)
  write(*,'(a,e12.4)') '  Residual 2-norm=', r_norm
  write(*,'(a)') '  First 4 residual values:'
  call VecView(r, PETSC_VIEWER_STDOUT_WORLD, ierr);CHKERRQ(ierr)

  if (debug%vecview_residual) then
    call DebugVecView(debug,r,'Sresidual','',0,0,0,realization%option)
  endif
  if (debug%vecview_solution) then
    call DebugVecView(debug,xx,'Sxx','',0,0,0,realization%option)
  endif

end subroutine StokesResidual

! ************************************************************************** !

subroutine StokesResidualPreliminaries(xx,r,realization,ierr)
  !
  ! Prepares residual computation by zeroing entries and updating auxiliary variables
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Field_module
  use Discretization_module

  implicit none

  Vec :: xx, r
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  type(field_type),          pointer :: field
  type(discretization_type), pointer :: disc

  field => realization%field
  disc  => realization%discretization

  call VecZeroEntries(r,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(disc,xx,field%flow_xx_loc, &
                                   realization%option%nflowdof)
  call StokesUpdateAuxVars(realization)

end subroutine StokesResidualPreliminaries

! ************************************************************************** !

subroutine StokesResidualInternalConn_structured(r,realization,mu,ierr)
  !
  ! Internal tractions for Stokes residual on structured grids,
  ! consistent with the LaTeX FV discretization.
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Grid_module
  use Patch_module
  use Option_module
  use Connection_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscReal :: mu
  PetscErrorCode :: ierr
  type(grid_type),          pointer :: grid
  type(patch_type),         pointer :: patch
  type(option_type),        pointer :: option
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type),      pointer :: cset
  PetscReal, pointer :: r_p(:), u_loc(:)
  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: ibeg_up_loc, ibeg_dn_loc
  PetscReal :: dx, area
  PetscReal :: distvec(3)
  PetscReal :: n(3)
  PetscReal :: u_up(3), u_dn(3)
  PetscReal :: p_up, p_dn, p_face
  PetscReal :: du_dn(3)
  PetscReal :: Tvisc(3)
  PetscReal :: Tpx, Tpy, Tpz
  PetscReal :: u_n_up, u_n_dn
  PetscReal :: flux_div

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)

  conn_list => grid%internal_connection_set_list
  cset      => conn_list%first
  do
    if (.not.associated(cset)) exit
    do iconn = 1, cset%num_connections
      gid_up = cset%id_up(iconn)
      gid_dn = cset%id_dn(iconn)
      if (patch%imat(gid_up) <= 0 .or. patch%imat(gid_dn) <= 0) cycle
      lid_up = grid%nG2L(gid_up)
      lid_dn = grid%nG2L(gid_dn)
      if (lid_up <= 0 .and. lid_dn <= 0) cycle

      distvec(1:3) = cset%dist(1:3,iconn)
      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)

      n = 0.d0
      if (abs(distvec(1)) > 0.d0 .and. &
          abs(distvec(2)) == 0.d0 .and. &
          abs(distvec(3)) == 0.d0) then
        n(1) = distvec(1) / abs(distvec(1))
      else if (abs(distvec(2)) > 0.d0 .and. &
               abs(distvec(1)) == 0.d0 .and. &
               abs(distvec(3)) == 0.d0) then
        n(2) = distvec(2) / abs(distvec(2))
      else if (abs(distvec(3)) > 0.d0 .and. &
               abs(distvec(1)) == 0.d0 .and. &
               abs(distvec(2)) == 0.d0) then
        n(3) = distvec(3) / abs(distvec(3))
      else
        cycle
      endif

      ! local DOF starts use local ids
      if (lid_up > 0) then
        ibeg_up_loc = (lid_up-1)*option%nflowdof + 1
        u_up(1) = u_loc(ibeg_up_loc    )
        u_up(2) = u_loc(ibeg_up_loc + 1)
        u_up(3) = u_loc(ibeg_up_loc + 2)
        p_up    = u_loc(ibeg_up_loc + 3)
      endif

      if (lid_dn > 0) then
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        u_dn(1) = u_loc(ibeg_dn_loc    )
        u_dn(2) = u_loc(ibeg_dn_loc + 1)
        u_dn(3) = u_loc(ibeg_dn_loc + 2)
        p_dn    = u_loc(ibeg_dn_loc + 3)
      else
        ! if one side is not local, use ghost/global access via gid indices
        ibeg_dn_loc = (grid%nG2L(gid_dn)-1)*option%nflowdof + 1
        u_dn(1) = u_loc(ibeg_dn_loc    )
        u_dn(2) = u_loc(ibeg_dn_loc + 1)
        u_dn(3) = u_loc(ibeg_dn_loc + 2)
        p_dn    = u_loc(ibeg_dn_loc + 3)
      endif

      ! viscous traction: -mu (u_L - u_K)/dx * area
      du_dn(1) = (u_dn(1) - u_up(1)) / dx
      du_dn(2) = (u_dn(2) - u_up(2)) / dx
      du_dn(3) = (u_dn(3) - u_up(3)) / dx
      Tvisc(1) = -mu * du_dn(1) * area
      Tvisc(2) = -mu * du_dn(2) * area
      Tvisc(3) = -mu * du_dn(3) * area

      ! pressure traction: -p_F n |F|
      p_face = 0.5d0*(p_up + p_dn)
      Tpx = -p_face * n(1) * area
      Tpy = -p_face * n(2) * area
      Tpz = -p_face * n(3) * area

      ! continuity: u_n,F = 0.5 (u_K·n + u_L·n)
      u_n_up = u_up(1)*n(1) + u_up(2)*n(2) + u_up(3)*n(3)
      u_n_dn = u_dn(1)*n(1) + u_dn(2)*n(2) + u_dn(3)*n(3)
      flux_div = 0.5d0 * (u_n_up + u_n_dn) * area

      if (lid_up > 0) then
        ! momentum residual for cell K
        r_p(ibeg_up_loc    ) = r_p(ibeg_up_loc    ) + Tvisc(1) + Tpx
        r_p(ibeg_up_loc + 1) = r_p(ibeg_up_loc + 1) + Tvisc(2) + Tpy
        r_p(ibeg_up_loc + 2) = r_p(ibeg_up_loc + 2) + Tvisc(3) + Tpz
        ! continuity residual for cell K
        r_p(ibeg_up_loc + 3) = r_p(ibeg_up_loc + 3) + flux_div
      endif

      if (lid_dn > 0) then
        ! opposite sign in neighboring cell L
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        r_p(ibeg_dn_loc    ) = r_p(ibeg_dn_loc    ) - Tvisc(1) - Tpx
        r_p(ibeg_dn_loc + 1) = r_p(ibeg_dn_loc + 1) - Tvisc(2) - Tpy
        r_p(ibeg_dn_loc + 2) = r_p(ibeg_dn_loc + 2) - Tvisc(3) - Tpz
        r_p(ibeg_dn_loc + 3) = r_p(ibeg_dn_loc + 3) - flux_div
      endif
    enddo
    cset => cset%next
  enddo

  call VecRestoreArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesResidualInternalConn_structured
! ************************************************************************** !

subroutine StokesResidualBoundaryConn(r,realization,rho,mu,ierr)
  !
  ! Boundary conditions: pressure (Dirichlet), velocity (Dirichlet or natural)
  ! Natural BC: zero residual on unconstrained velocity at boundaries
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  PetscErrorCode :: ierr
  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type),pointer :: bc
  type(connection_set_type), pointer :: cset
  PetscReal, pointer :: r_p(:), u_loc(:)
  PetscInt :: iconn, local_id, gid, ibeg
  PetscReal :: u_bc(3), p_bc
  ! Track which cells are on boundary and which velocity components have Dirichlet BC
  PetscBool, allocatable :: on_boundary(:), vel_dirichlet(:,:)

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  ! Allocate tracking arrays
  allocate(on_boundary(grid%nlmax))
  allocate(vel_dirichlet(3, grid%nlmax))
  on_boundary = PETSC_FALSE
  vel_dirichlet = PETSC_FALSE

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)

  ! First pass: identify boundary cells and apply explicit Dirichlet BCs
  bc => patch%boundary_condition_list%first
  do
    if (.not.associated(bc)) exit

    if (.not.associated(bc%flow_condition)) then
      bc => bc%next
      cycle
    endif

    if (.not.associated(bc%flow_aux_real_var)) then
      bc => bc%next
      cycle
    endif

    cset => bc%connection_set
    if (.not.associated(cset)) then
      bc => bc%next
      cycle
    endif

    do iconn = 1, cset%num_connections
      local_id = cset%id_dn(iconn)
      gid      = grid%nL2G(local_id)
      if (patch%imat(gid) <= 0) cycle

      ! Mark this cell as being on a boundary
      on_boundary(local_id) = PETSC_TRUE

      ibeg = (local_id-1)*option%nflowdof + 1

      ! Pressure BC (Dirichlet)
      if (associated(bc%flow_condition%stokes%pressure)) then
        if (bc%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          if (size(bc%flow_aux_real_var,1) >= 4 .and. &
              size(bc%flow_aux_real_var,2) >= iconn) then
            p_bc = bc%flow_aux_real_var(4,iconn)
            r_p(ibeg+3) = u_loc(ibeg+3) - p_bc
            write(*,'(a,i3,a,e12.4)') '    Pressure BC: local_id=', local_id, ' p_bc=', p_bc
          endif
        endif
      endif

      ! Velocity BC (Dirichlet)
      if (associated(bc%flow_condition%stokes%velocity)) then
        if (bc%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          if (size(bc%flow_aux_real_var,1) >= 3 .and. &
              size(bc%flow_aux_real_var,2) >= iconn) then
            u_bc(1) = bc%flow_aux_real_var(1,iconn)
            u_bc(2) = bc%flow_aux_real_var(2,iconn)
            u_bc(3) = bc%flow_aux_real_var(3,iconn)
            r_p(ibeg    ) = u_loc(ibeg    ) - u_bc(1)
            r_p(ibeg + 1) = u_loc(ibeg + 1) - u_bc(2)
            r_p(ibeg + 2) = u_loc(ibeg + 2) - u_bc(3)
            ! Mark velocity as Dirichlet constrained
            vel_dirichlet(1, local_id) = PETSC_TRUE
            vel_dirichlet(2, local_id) = PETSC_TRUE
            vel_dirichlet(3, local_id) = PETSC_TRUE
          endif
        endif
      endif
    enddo
    bc => bc%next
  enddo

  ! Second pass: Apply natural BC (zero residual) ONLY on boundary cells
  ! with unconstrained velocity
  do local_id = 1, grid%nlmax
    gid = grid%nL2G(local_id)
    if (patch%imat(gid) <= 0) cycle

    ! Only apply natural BC if this cell is ON the boundary
    if (.not.on_boundary(local_id)) cycle

    ibeg = (local_id-1)*option%nflowdof + 1

    ! Zero velocity rows that don't have explicit Dirichlet BC (natural BC)
    if (.not.vel_dirichlet(1, local_id)) then
      r_p(ibeg    ) = 0.d0
    endif
    if (.not.vel_dirichlet(2, local_id)) then
      r_p(ibeg + 1) = 0.d0
    endif
    if (.not.vel_dirichlet(3, local_id)) then
      r_p(ibeg + 2) = 0.d0
    endif
  enddo

  call VecRestoreArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

  deallocate(on_boundary)
  deallocate(vel_dirichlet)

end subroutine StokesResidualBoundaryConn

! ************************************************************************** !

subroutine StokesResidualAccumulation(r,realization,rho,ierr)
  !
  ! Accumulation terms for Stokes residual
  ! r -= rho*vol*u^n, then r += rho*vol*u^{n+1}
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Grid_module
  use Patch_module
  use Option_module
  use Field_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  PetscErrorCode :: ierr
  type(grid_type),   pointer :: grid
  type(patch_type),  pointer :: patch
  type(option_type), pointer :: option
  type(field_type),  pointer :: field
  PetscReal, pointer :: r_p(:), xx_loc_p(:), accum_p(:)
  PetscInt :: lid, gid, ibeg
  PetscReal :: vol

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  field  => realization%field

  if (option%flow%steady_state) return

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

  ! Subtract previous accumulation: r -= accum_p (which contains rho*vol*u^n)
  r_p(1:grid%nlmax*option%nflowdof) = r_p(1:grid%nlmax*option%nflowdof) - &
                                      accum_p(1:grid%nlmax*option%nflowdof)

  ! Add current accumulation: r += rho*vol*u^{n+1}
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    ibeg = (lid-1)*option%nflowdof + 1
    vol  = patch%aux%Material%auxvars(gid)%volume

    ! Add rho*volume*u^{n+1} to residual
    r_p(ibeg    ) = r_p(ibeg    ) + rho * xx_loc_p(ibeg    ) * vol
    r_p(ibeg + 1) = r_p(ibeg + 1) + rho * xx_loc_p(ibeg + 1) * vol
    r_p(ibeg + 2) = r_p(ibeg + 2) + rho * xx_loc_p(ibeg + 2) * vol
  enddo

  call VecRestoreArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesResidualAccumulation

! ************************************************************************** !

subroutine StokesResidualSourceSink(r,realization,rho,ierr)
  !
  ! Source/sink terms for Stokes residual
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  PetscErrorCode :: ierr

end subroutine StokesResidualSourceSink
! ************************************************************************** !

subroutine StokesJacobian(snes,xx,A,B,realization,rho,mu,debug,ierr)
  !
  ! Computes Stokes Jacobian for SNES solver
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec  :: xx
  Mat  :: A, B
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  type(debug_type) :: debug
  PetscErrorCode :: ierr
  Mat :: J
  MatType :: mat_type_A, mat_type_B
  type(option_type), pointer :: option
  PetscReal :: mat_norm
  PetscInt :: mat_rows, mat_cols
  PetscMPIInt :: mpi_rank

  option => realization%option

  call MatGetType(A,mat_type_A,ierr);CHKERRQ(ierr)
  call MatGetType(B,mat_type_B,ierr);CHKERRQ(ierr)

  if (mat_type_A == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)
  call StokesJacobianInternalConn_structured(J,realization,mu,debug,ierr)
  call StokesJacobianAccumulation(J,realization,rho,debug,ierr)
  call StokesJacobianSourceSink(J,realization,rho,debug,ierr)

  ! Assembly happens once at the end
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! NOW apply BC constraints AFTER assembly
  call StokesJacobianBoundaryConn(J,realization,rho,mu,debug,ierr)

  ! Re-assemble after BC modifications
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr);CHKERRQ(ierr)

  if (mpi_rank == 0) then
    call MatGetSize(A, mat_rows, mat_cols, ierr)
    call MatNorm(A, NORM_FROBENIUS, mat_norm, ierr)
    write(*,'(a,i8,a,i8,a,e12.4)') &
      '  Matrix: ', mat_rows, 'x', mat_cols, ' Frobenius norm=', mat_norm
    write(*,'(a)') '  First 4 matrix diagonal blocks:'
    call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)
  endif

  if (A /= B .and. mat_type_A /= MATMFFD) then
    call MatConvert(J,mat_type_B,MAT_REUSE_MATRIX,B,ierr);CHKERRQ(ierr)
  endif

  if (debug%matview_Matrix) then
    call DebugMatView(debug,J,'Sjacobian','',0,0,0,option)
  endif

end subroutine StokesJacobian

! ************************************************************************** !

subroutine StokesJacobianInternalConn_structured(A,realization,mu,debug,ierr)
  !
  ! Internal connection Jacobian for structured grids
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Petsc_Utility_module
  use Debug_module

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: mu
  type(debug_type) :: debug
  PetscErrorCode :: ierr
  type(grid_type),          pointer :: grid
  type(patch_type),         pointer :: patch
  type(option_type),        pointer :: option
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type),      pointer :: cset
  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: istart_up, istart_dn
  PetscReal :: dx, area
  PetscReal :: distvec(3)
  PetscReal :: n(3)
  PetscReal :: J_KK(4,4), J_KL(4,4), J_LK(4,4), J_LL(4,4)
  PetscReal :: coef_mu, coef_p, coef_div

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  conn_list => grid%internal_connection_set_list
  cset      => conn_list%first
  do
    if (.not.associated(cset)) exit
    do iconn = 1, cset%num_connections
      gid_up = cset%id_up(iconn)
      gid_dn = cset%id_dn(iconn)
      if (patch%imat(gid_up) <= 0 .or. patch%imat(gid_dn) <= 0) cycle
      lid_up = grid%nG2L(gid_up)
      lid_dn = grid%nG2L(gid_dn)
      if (lid_up <= 0 .and. lid_dn <= 0) cycle

      distvec(1:3) = cset%dist(1:3,iconn)
      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)

      n = 0.d0
      if (abs(distvec(1)) > 0.d0 .and. &
          abs(distvec(2)) == 0.d0 .and. &
          abs(distvec(3)) == 0.d0) then
        n(1) = distvec(1) / abs(distvec(1))
      else if (abs(distvec(2)) > 0.d0 .and. &
               abs(distvec(1)) == 0.d0 .and. &
               abs(distvec(3)) == 0.d0) then
        n(2) = distvec(2) / abs(distvec(2))
      else if (abs(distvec(3)) > 0.d0 .and. &
               abs(distvec(1)) == 0.d0 .and. &
               abs(distvec(2)) == 0.d0) then
        n(3) = distvec(3) / abs(distvec(3))
      else
        cycle
      endif

      coef_mu  = mu * area / dx
      coef_p   = 0.5d0 * area
      coef_div = 0.5d0 * area

      J_KK = 0.d0
      J_KL = 0.d0
      J_LK = 0.d0
      J_LL = 0.d0

      ! Viscous contributions: R_K^u = -mu (u_L - u_K)/dx |F|
      J_KK(1,1) =  coef_mu
      J_KK(2,2) =  coef_mu
      J_KK(3,3) =  coef_mu
      J_KL(1,1) = -coef_mu
      J_KL(2,2) = -coef_mu
      J_KL(3,3) = -coef_mu

      ! R_L^u = -(R_K^u): derivative is the negative of above
      J_LL(1,1) =  coef_mu
      J_LL(2,2) =  coef_mu
      J_LL(3,3) =  coef_mu
      J_LK(1,1) = -coef_mu
      J_LK(2,2) = -coef_mu
      J_LK(3,3) = -coef_mu

      ! Pressure tractions: -p_F n |F| with p_F = 0.5(p_K + p_L)
      ! R_K^u components
      J_KK(1,4) = -coef_p * n(1)
      J_KK(2,4) = -coef_p * n(2)
      J_KK(3,4) = -coef_p * n(3)
      J_KL(1,4) = -coef_p * n(1)
      J_KL(2,4) = -coef_p * n(2)
      J_KL(3,4) = -coef_p * n(3)

      ! R_L^u = -R_K^u
      J_LK(1,4) =  coef_p * n(1)
      J_LK(2,4) =  coef_p * n(2)
      J_LK(3,4) =  coef_p * n(3)
      J_LL(1,4) =  coef_p * n(1)
      J_LL(2,4) =  coef_p * n(2)
      J_LL(3,4) =  coef_p * n(3)

      ! Continuity: R_K^p = 0.5(u_K·n + u_L·n)|F|
      J_KK(4,1) = coef_div * n(1)
      J_KK(4,2) = coef_div * n(2)
      J_KK(4,3) = coef_div * n(3)
      J_KL(4,1) = coef_div * n(1)
      J_KL(4,2) = coef_div * n(2)
      J_KL(4,3) = coef_div * n(3)

      ! R_L^p = -R_K^p
      J_LK(4,1) = -coef_div * n(1)
      J_LK(4,2) = -coef_div * n(2)
      J_LK(4,3) = -coef_div * n(3)
      J_LL(4,1) = -coef_div * n(1)
      J_LL(4,2) = -coef_div * n(2)
      J_LL(4,3) = -coef_div * n(3)

      if (lid_up > 0) then
        istart_up = lid_up - 1  ! Block index for BAIJ matrix
        if (lid_dn > 0) then
          istart_dn = lid_dn - 1  ! Block index for BAIJ matrix
          call PUMSetValuesBlockedLocal(A,1,istart_up,1,istart_up,J_KK,ADD_VALUES,ierr);CHKERRQ(ierr)
          call PUMSetValuesBlockedLocal(A,1,istart_up,1,istart_dn,J_KL,ADD_VALUES,ierr);CHKERRQ(ierr)
        else
          call PUMSetValuesBlockedLocal(A,1,istart_up,1,istart_up,J_KK,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      endif

      if (lid_dn > 0) then
        istart_dn = lid_dn - 1  ! Block index for BAIJ matrix
        if (lid_up > 0) then
          istart_up = lid_up - 1  ! Block index for BAIJ matrix
          call PUMSetValuesBlockedLocal(A,1,istart_dn,1,istart_dn,J_LL,ADD_VALUES,ierr);CHKERRQ(ierr)
          call PUMSetValuesBlockedLocal(A,1,istart_dn,1,istart_up,J_LK,ADD_VALUES,ierr);CHKERRQ(ierr)
        else
          call PUMSetValuesBlockedLocal(A,1,istart_dn,1,istart_dn,J_LL,ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      endif
    enddo
    cset => cset%next
  enddo

end subroutine StokesJacobianInternalConn_structured

! ************************************************************************** !
subroutine StokesJacobianBoundaryConn(A,realization,rho,mu,debug,ierr)
  !
  ! BC Jacobian for Stokes problem
  ! Collects all BC constraints and applies once per cell
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Petsc_Utility_module
  use Debug_module

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  type(debug_type) :: debug
  PetscErrorCode :: ierr
  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type),pointer :: bc
  type(connection_set_type), pointer :: cset
  PetscInt :: iconn, local_id, gid, iblock, lid
  PetscReal :: Jbc(4,4)
  ! Track which rows have been constrained for this cell
  PetscBool, allocatable :: bc_applied(:,:)
  PetscInt :: max_local_id

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  ! Allocate tracking array for which DOFs have BCs
  max_local_id = grid%nlmax
  allocate(bc_applied(4, max_local_id))
  bc_applied = PETSC_FALSE

  ! First pass: identify all BC constraints
  bc => patch%boundary_condition_list%first
  do
    if (.not.associated(bc)) exit
    if (.not.associated(bc%flow_condition)) then
      bc => bc%next
      cycle
    endif
    if (.not.associated(bc%flow_aux_real_var)) then
      bc => bc%next
      cycle
    endif
    cset => bc%connection_set
    if (.not.associated(cset)) then
      bc => bc%next
      cycle
    endif

    do iconn = 1, cset%num_connections
      local_id = cset%id_dn(iconn)
      gid      = grid%nL2G(local_id)
      if (patch%imat(gid) <= 0) cycle

      ! Mark which components have BCs
      if (associated(bc%flow_condition%stokes%velocity)) then
        if (bc%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          bc_applied(1, local_id) = PETSC_TRUE
          bc_applied(2, local_id) = PETSC_TRUE
          bc_applied(3, local_id) = PETSC_TRUE
        endif
      endif

      if (associated(bc%flow_condition%stokes%pressure)) then
        if (bc%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          bc_applied(4, local_id) = PETSC_TRUE
          write(*,'(a,i3,a)') '    Pressure BC applied to local_id=', local_id
        endif
      endif
    enddo
    bc => bc%next
  enddo

  ! Second pass: apply all accumulated BCs
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle

    iblock = lid - 1  ! Block index for BAIJ matrix

    Jbc = 0.d0

    ! ALWAYS apply velocity identity (Dirichlet or natural BC)
    Jbc(1,1) = 1.d0
    Jbc(2,2) = 1.d0
    Jbc(3,3) = 1.d0

    ! Apply pressure identity ONLY if Dirichlet BC
    if (bc_applied(4, lid)) then
      Jbc(4,4) = 1.d0
    endif

    ! Only call if at least one constraint applies
    if (bc_applied(1, lid) .or. bc_applied(4, lid)) then
      call PUMSetValuesBlockedLocal(A,1,iblock,1,iblock,Jbc,INSERT_VALUES,ierr);CHKERRQ(ierr)
    endif
  enddo

  deallocate(bc_applied)

end subroutine StokesJacobianBoundaryConn

! ************************************************************************** !

subroutine StokesJacobianAccumulation(A,realization,rho,debug,ierr)
  !
  ! Accumulation Jacobian for Stokes problem
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Petsc_Utility_module
  use Debug_module

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  type(debug_type) :: debug
  PetscErrorCode :: ierr
  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  PetscInt :: lid, gid, istart
  PetscReal :: Jacc(4,4), vol

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  if (option%flow%steady_state) return

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    istart = lid - 1  ! Block index for BAIJ matrix
    vol    = patch%aux%Material%auxvars(gid)%volume

    Jacc = 0.d0
    Jacc(1,1) = rho * vol
    Jacc(2,2) = rho * vol
    Jacc(3,3) = rho * vol

    call PUMSetValuesBlockedLocal(A,1,istart,1,istart,Jacc,ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

end subroutine StokesJacobianAccumulation

! ************************************************************************** !

subroutine StokesJacobianSourceSink(A,realization,rho,debug,ierr)
  !
  ! Source/sink Jacobian for Stokes problem
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Debug_module

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  type(debug_type) :: debug
  PetscErrorCode :: ierr

end subroutine StokesJacobianSourceSink

! ************************************************************************** !

subroutine StokesUpdateFixedAccum(realization)
  !
  ! Updates the fixed portion of the accumulation term
  ! Called at the beginning of each time step
  !
  ! Author: Satish Karra
  ! Date: 01/28/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  type(field_type),  pointer :: field
  PetscReal, pointer :: xx_p(:), accum_p(:)
  PetscInt :: lid, gid, istart
  PetscReal :: vol
  PetscErrorCode :: ierr

  option => realization%option
  field  => realization%field
  patch  => realization%patch
  grid   => patch%grid

  call VecGetArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

  ! Zero accumulation storage
  accum_p(:) = 0.d0

  ! Store u^n (velocity at previous time step) for accumulation
  ! This is the "fixed" part that doesn't change during Newton iterations
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    istart = (lid-1)*option%nflowdof + 1
    vol    = patch%aux%Material%auxvars(gid)%volume

    ! accum stores: rho * volume * u^n (previous velocity)
    ! This will be subtracted from residual to get: rho*vol*(u^{n+1} - u^n)/dt
    accum_p(istart    ) = xx_p(istart    ) * vol
    accum_p(istart + 1) = xx_p(istart + 1) * vol
    accum_p(istart + 2) = xx_p(istart + 2) * vol
  enddo

  call VecRestoreArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

end subroutine StokesUpdateFixedAccum

end module Stokes_module
