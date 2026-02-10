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
            StokesMapBCAuxVarsToGlobal, &
            StokesUpdateTransportVelocities, &
            StokesUpdateSolution, &
            StokesMaxChange, &
            StokesDestroy, &
            StokesResidual, &
            StokesJacobian, &
            StokesResidualInternalConn, &
            StokesResidualConvectionInternalConn, &
            StokesResidualAccumulation, &
            StokesIFunctionAccumulation, &
            StokesResidualSourceSink, &
            StokesResidualBoundaryConn, &
            StokesResidualBoundaryFinal, &
            StokesResidualPinPressure, &
            StokesJacobianInternalConn, &
            StokesJacobianBoundaryConn, &
            StokesJacobianBoundaryFacesConn, &
            StokesJacobianPinPressure, &
            StokesJacobianAccumulation, &
            StokesIJacobianAccumulation, &
            StokesJacobianSourceSink, &
            StokesApplyBoundaryConditions

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
  use Coupler_module
  use Connection_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  type(patch_type),         pointer :: patch
  type(grid_type),          pointer :: grid
  type(material_auxvar_type), pointer :: mat_aux(:)
  type(stokes_type),        pointer :: stokes_aux
  type(option_type),        pointer :: option

  PetscInt :: gid
  PetscBool :: dof_is_active(1)
  PetscInt :: iconn, sum_connection

  patch  => realization%patch
  grid   => patch%grid
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

subroutine StokesInitializeTimestep(realization,rho)
  !
  ! Update data in module prior to time step
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: rho

  call StokesApplyBoundaryConditions(realization)
  call StokesUpdateFixedAccum(realization,rho)

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
  call StokesMapBCAuxVarsToGlobal(realization)

end subroutine StokesUpdateAuxVars

! ************************************************************************** !

subroutine StokesMapBCAuxVarsToGlobal(realization)
  !
  ! Maps Stokes boundary condition auxvars to global equivalent.
  !
  ! Author: Satish Karra
  ! Date: 02/09/26
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),  pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)

  PetscInt :: sum_connection, iconn, iphase

  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return
  if (.not.associated(patch%aux%Global%auxvars_bc)) return

  global_auxvars_bc => patch%aux%Global%auxvars_bc

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%temp = option%flow%reference_temperature
      do iphase = 1, option%nphase
        global_auxvars_bc(sum_connection)%pres(iphase) = option%flow%reference_pressure
        global_auxvars_bc(sum_connection)%sat(iphase) = 0.d0
        global_auxvars_bc(sum_connection)%den_kg(iphase) = 0.d0
      enddo
      global_auxvars_bc(sum_connection)%sat(option%liquid_phase) = &
        option%flow%reference_saturation
      global_auxvars_bc(sum_connection)%den_kg(option%liquid_phase) = &
        option%flow%reference_density(option%liquid_phase)
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine StokesMapBCAuxVarsToGlobal

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

subroutine StokesResidual(snes,xx,r,realization,rho,mu,add_convection,debug,ierr)
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
  Vec :: xx, r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  PetscBool :: add_convection
  type(debug_type) :: debug
  PetscErrorCode :: ierr

  call StokesResidualPreliminaries(xx,r,realization,ierr)
  call StokesResidualInternalConn(r,realization,rho,mu, &
                                            add_convection,ierr)
  call StokesResidualAccumulation(r,realization,rho,ierr)
  call StokesResidualSourceSink(r,realization,rho,ierr)
  call StokesResidualPinPressure(r,xx,realization,0.d0,ierr)
  call StokesResidualBoundaryFinal(r,xx,realization,ierr)


  if (debug%vecview_residual) then
    call DebugVecView(debug,r,'Sresidual','',0,0,0,realization%option)
  endif
  if (debug%vecview_solution) then
    call DebugVecView(debug,xx,'Sxx','',0,0,0,realization%option)
  endif

end subroutine StokesResidual

! ************************************************************************** !

subroutine StokesApplyBoundaryConditions(realization)
  !
  ! Apply Dirichlet BCs to solution vector
  !
  ! Author: Satish Karra
  ! Date: 01/28/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Discretization_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id, iconn, offset, gid
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  option => realization%option

  if (.not.associated(patch%boundary_condition_list%first)) return

  call VecGetArray(realization%field%flow_xx, vec_ptr, ierr);CHKERRQ(ierr)

  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      offset = (local_id-1)*option%nflowdof

      ! Apply velocity BC (all 3 components)
      if (associated(boundary_condition%flow_condition%stokes%velocity)) then
        if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          vec_ptr(offset+1) = boundary_condition%flow_aux_real_var(1,iconn)
          vec_ptr(offset+2) = boundary_condition%flow_aux_real_var(2,iconn)
          vec_ptr(offset+3) = boundary_condition%flow_aux_real_var(3,iconn)
        endif
      endif

      ! Apply pressure BC
      if (associated(boundary_condition%flow_condition%stokes%pressure)) then
        if (boundary_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          vec_ptr(offset+4) = boundary_condition%flow_aux_real_var(4,iconn)
        endif
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArray(realization%field%flow_xx, vec_ptr, ierr);CHKERRQ(ierr)


  call DiscretizationGlobalToLocal(realization%discretization, &
                                   realization%field%flow_xx, &
                                   realization%field%flow_xx_loc, &
                                   realization%option%nflowdof)


end subroutine StokesApplyBoundaryConditions

! ************************************************************************** !

subroutine StokesUpdateTransportVelocities(realization)
  !
  ! Populate transport face velocities from Stokes velocity field.
  !
  ! Author: Satish Karra
  ! Date: 02/09/26
  !
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Stokes_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(stokes_type), pointer :: stokes_aux
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type), pointer :: cset

  PetscInt :: iconn, sum_connection
  PetscInt :: gid_up, gid_dn
  PetscReal :: u_face(3)
  PetscReal :: vdarcy

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  if (.not.associated(patch%internal_velocities)) return
  if (.not.associated(patch%boundary_velocities)) return
  if (.not.associated(patch%aux%Stokes)) return

  stokes_aux => patch%aux%Stokes

  patch%internal_velocities = 0.d0
  patch%boundary_velocities = 0.d0

  ! Internal flux terms
  conn_list => grid%internal_connection_set_list
  cset      => conn_list%first
  sum_connection = 0
  do
    if (.not.associated(cset)) exit
    do iconn = 1, cset%num_connections
      sum_connection = sum_connection + 1
      gid_up = cset%id_up(iconn)
      gid_dn = cset%id_dn(iconn)
      if (gid_up <= 0 .or. gid_dn <= 0) cycle
      if (patch%imat(gid_up) <= 0 .or. patch%imat(gid_dn) <= 0) cycle

      u_face(1) = 0.5d0*(stokes_aux%auxvars(gid_up)%u + &
                         stokes_aux%auxvars(gid_dn)%u)
      u_face(2) = 0.5d0*(stokes_aux%auxvars(gid_up)%v + &
                         stokes_aux%auxvars(gid_dn)%v)
      u_face(3) = 0.5d0*(stokes_aux%auxvars(gid_up)%w + &
                         stokes_aux%auxvars(gid_dn)%w)

      vdarcy = dot_product(u_face, cset%dist(1:3,iconn))
      patch%internal_velocities(1,sum_connection) = vdarcy
    enddo
    cset => cset%next
  enddo

  ! Boundary flux terms
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cset => boundary_condition%connection_set
    do iconn = 1, cset%num_connections
      sum_connection = sum_connection + 1

      if (associated(boundary_condition%flow_condition%stokes%velocity)) then
        if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          u_face(1) = boundary_condition%flow_aux_real_var(1,iconn)
          u_face(2) = boundary_condition%flow_aux_real_var(2,iconn)
          u_face(3) = boundary_condition%flow_aux_real_var(3,iconn)
        else
          gid_dn = cset%id_dn(iconn)
          if (gid_dn <= 0) cycle
          u_face(1) = stokes_aux%auxvars(gid_dn)%u
          u_face(2) = stokes_aux%auxvars(gid_dn)%v
          u_face(3) = stokes_aux%auxvars(gid_dn)%w
        endif
      else
        gid_dn = cset%id_dn(iconn)
        if (gid_dn <= 0) cycle
        u_face(1) = stokes_aux%auxvars(gid_dn)%u
        u_face(2) = stokes_aux%auxvars(gid_dn)%v
        u_face(3) = stokes_aux%auxvars(gid_dn)%w
      endif

      vdarcy = dot_product(u_face, cset%dist(1:3,iconn))
      patch%boundary_velocities(1,sum_connection) = vdarcy
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine StokesUpdateTransportVelocities

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

subroutine StokesResidualInternalConn(r,realization,rho,mu,add_convection,ierr)
  !
  ! Internal connection contributions for incompressible NS:
  !  - diffusion:     + mu A/d (u_K - u_L)
  !  - grad(p):       + (p_L - p_K)/d * S   with S = n*A
  !  - convection:    + rho*phi*u_upwind    (optional)
  !  - continuity:    + phi (cell K), -phi (cell L)
  !  - p-stab:        - tau (p_K - p_L) in continuity, tau = gamma A/d
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
  PetscReal :: mu, rho
  PetscBool :: add_convection
  PetscErrorCode :: ierr

  type(grid_type),          pointer :: grid
  type(patch_type),         pointer :: patch
  type(option_type),        pointer :: option
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type),      pointer :: cset

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: u_loc(:)

  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: ibeg_up_loc, ibeg_dn_loc

  PetscReal :: dx, area
  PetscReal :: n(3), S(3)
  PetscReal :: u_up(3), u_dn(3), u_face(3), u_upwind(3)
  PetscReal :: p_up, p_dn
  PetscReal :: a_mu, phi
  PetscReal :: Tvisc(3), Tpress(3), mconv(3)
  PetscReal :: flux_div
  PetscReal :: tau, gamma_pstab
  PetscReal :: pstab_K, pstab_L

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  gamma_pstab = 0.05  ! add this to input deck; start ~0.05

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

      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)
      n    = 0.d0
      n(1:3) = cset%dist(1:3,iconn)   ! outward normal from "up" cell (K)

      S(1) = n(1)*area
      S(2) = n(2)*area
      S(3) = n(3)*area

      ! Initialize
      u_up = 0.d0; u_dn = 0.d0
      p_up = 0.d0; p_dn = 0.d0

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
      endif

      ! Face velocity (central)
      u_face(1) = 0.5d0*(u_up(1) + u_dn(1))
      u_face(2) = 0.5d0*(u_up(2) + u_dn(2))
      u_face(3) = 0.5d0*(u_up(3) + u_dn(3))

      ! Volumetric flux phi = u_face · S (from up(K) to dn(L))
      phi = u_face(1)*S(1) + u_face(2)*S(2) + u_face(3)*S(3)

      ! Diffusion coefficient and stabilization
      a_mu = mu * area / dx
      tau  = gamma_pstab * area / dx

      ! Diffusion contribution: +a_mu (uK - uL)
      Tvisc(1) = a_mu*(u_up(1) - u_dn(1))
      Tvisc(2) = a_mu*(u_up(2) - u_dn(2))
      Tvisc(3) = a_mu*(u_up(3) - u_dn(3))

      ! Pressure gradient contribution: +((pL - pK)/dx)*S
      Tpress(1) = ((p_dn - p_up)/dx) * S(1)
      Tpress(2) = ((p_dn - p_up)/dx) * S(2)
      Tpress(3) = ((p_dn - p_up)/dx) * S(3)

      ! Convection: + rho*phi*u_upwind (optional)
      mconv = 0.d0
      if (add_convection) then
        if (phi >= 0.d0) then
          u_upwind = u_up
        else
          u_upwind = u_dn
        endif
        mconv(1) = rho * phi * u_upwind(1)
        mconv(2) = rho * phi * u_upwind(2)
        mconv(3) = rho * phi * u_upwind(3)
      endif

      ! Continuity flux (divergence): +phi for K, -phi for L
      flux_div = phi

      ! Pressure stabilization in continuity: -tau(pK - pL)
      pstab_K = -tau*(p_up - p_dn)
      pstab_L = -tau*(p_dn - p_up)

      ! Assemble into residual vector
      if (lid_up > 0) then
        ibeg_up_loc = (lid_up-1)*option%nflowdof + 1
        r_p(ibeg_up_loc    ) = r_p(ibeg_up_loc    ) + Tvisc(1) + Tpress(1) + mconv(1)
        r_p(ibeg_up_loc + 1) = r_p(ibeg_up_loc + 1) + Tvisc(2) + Tpress(2) + mconv(2)
        r_p(ibeg_up_loc + 2) = r_p(ibeg_up_loc + 2) + Tvisc(3) + Tpress(3) + mconv(3)
        r_p(ibeg_up_loc + 3) = r_p(ibeg_up_loc + 3) + flux_div + pstab_K
      endif

      if (lid_dn > 0) then
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        r_p(ibeg_dn_loc    ) = r_p(ibeg_dn_loc    ) - Tvisc(1) - Tpress(1) - mconv(1)
        r_p(ibeg_dn_loc + 1) = r_p(ibeg_dn_loc + 1) - Tvisc(2) - Tpress(2) - mconv(2)
        r_p(ibeg_dn_loc + 2) = r_p(ibeg_dn_loc + 2) - Tvisc(3) - Tpress(3) - mconv(3)
        r_p(ibeg_dn_loc + 3) = r_p(ibeg_dn_loc + 3) - flux_div + pstab_L
      endif

    enddo

    cset => cset%next
  enddo

  call VecRestoreArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesResidualInternalConn

! ************************************************************************** !

subroutine StokesResidualConvectionInternalConn(r,realization,rho,ierr)
  !
  ! Internal connection convection contributions for IMEX split
  !  - convection: + rho*phi*u_upwind
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
  PetscReal :: rho
  PetscErrorCode :: ierr

  type(grid_type),          pointer :: grid
  type(patch_type),         pointer :: patch
  type(option_type),        pointer :: option
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type),      pointer :: cset

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: u_loc(:)

  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: ibeg_up_loc, ibeg_dn_loc

  PetscReal :: area
  PetscReal :: n(3), S(3)
  PetscReal :: u_up(3), u_dn(3), u_face(3), u_upwind(3)
  PetscReal :: phi
  PetscReal :: mconv(3)

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

      area = cset%area(iconn)
      n    = 0.d0
      n(1:3) = cset%dist(1:3,iconn)

      S(1) = n(1)*area
      S(2) = n(2)*area
      S(3) = n(3)*area

      u_up = 0.d0; u_dn = 0.d0

      if (lid_up > 0) then
        ibeg_up_loc = (lid_up-1)*option%nflowdof + 1
        u_up(1) = u_loc(ibeg_up_loc    )
        u_up(2) = u_loc(ibeg_up_loc + 1)
        u_up(3) = u_loc(ibeg_up_loc + 2)
      endif

      if (lid_dn > 0) then
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        u_dn(1) = u_loc(ibeg_dn_loc    )
        u_dn(2) = u_loc(ibeg_dn_loc + 1)
        u_dn(3) = u_loc(ibeg_dn_loc + 2)
      endif

      u_face(1) = 0.5d0*(u_up(1) + u_dn(1))
      u_face(2) = 0.5d0*(u_up(2) + u_dn(2))
      u_face(3) = 0.5d0*(u_up(3) + u_dn(3))

      phi = u_face(1)*S(1) + u_face(2)*S(2) + u_face(3)*S(3)

      if (phi >= 0.d0) then
        u_upwind = u_up
      else
        u_upwind = u_dn
      endif

      mconv(1) = rho * phi * u_upwind(1)
      mconv(2) = rho * phi * u_upwind(2)
      mconv(3) = rho * phi * u_upwind(3)

      if (lid_up > 0) then
        ibeg_up_loc = (lid_up-1)*option%nflowdof + 1
        r_p(ibeg_up_loc    ) = r_p(ibeg_up_loc    ) + mconv(1)
        r_p(ibeg_up_loc + 1) = r_p(ibeg_up_loc + 1) + mconv(2)
        r_p(ibeg_up_loc + 2) = r_p(ibeg_up_loc + 2) + mconv(3)
      endif

      if (lid_dn > 0) then
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        r_p(ibeg_dn_loc    ) = r_p(ibeg_dn_loc    ) - mconv(1)
        r_p(ibeg_dn_loc + 1) = r_p(ibeg_dn_loc + 1) - mconv(2)
        r_p(ibeg_dn_loc + 2) = r_p(ibeg_dn_loc + 2) - mconv(3)
      endif

    enddo

    cset => cset%next
  enddo

  call VecRestoreArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesResidualConvectionInternalConn

! ************************************************************************** !

subroutine StokesResidualAccumulation(r,realization,rho,ierr)
  !
  ! Accumulation terms for Stokes residual
  ! r -= rho*vol*u^n, then r += rho*vol*u^{n+1}
  ! Pressure is NOT accumulated (stays determined by continuity)
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
  PetscReal, pointer :: r_p(:), xx_loc_p(:), accum_p(:), accum_p2(:)
  PetscInt :: lid, gid, ibeg
  PetscReal :: vol, vol_over_dt

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  field  => realization%field

  if (option%flow%steady_state) return

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)


  accum_p2 = 0.d0

  ! Subtract previous accumulation and add current accumulation
  ! Only for velocity components (1,2,3), NOT pressure (4)
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    ibeg = (lid-1)*option%nflowdof + 1
    vol  = patch%aux%Material%auxvars(gid)%volume
    vol_over_dt = vol / option%flow_dt

    ! Subtract old: r -= rho*vol*u^n (velocity only)
    r_p(ibeg    ) = r_p(ibeg    ) - accum_p(ibeg    )
    r_p(ibeg + 1) = r_p(ibeg + 1) - accum_p(ibeg + 1)
    r_p(ibeg + 2) = r_p(ibeg + 2) - accum_p(ibeg + 2)
    ! DO NOT touch pressure residual (ibeg + 3)

    ! Add new: r += rho*vol*u^{n+1} (velocity only)
    r_p(ibeg    ) = r_p(ibeg    ) + rho * xx_loc_p(ibeg    ) * vol_over_dt
    r_p(ibeg + 1) = r_p(ibeg + 1) + rho * xx_loc_p(ibeg + 1) * vol_over_dt
    r_p(ibeg + 2) = r_p(ibeg + 2) + rho * xx_loc_p(ibeg + 2) * vol_over_dt
    ! DO NOT touch pressure residual (ibeg + 3)

    accum_p2(ibeg    ) = rho * xx_loc_p(ibeg    ) * vol_over_dt
    accum_p2(ibeg + 1) = rho * xx_loc_p(ibeg + 1) * vol_over_dt
    accum_p2(ibeg + 2) = rho * xx_loc_p(ibeg + 2) * vol_over_dt
    accum_p2(ibeg + 3) = 0.d0
  enddo

  call VecRestoreArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesResidualAccumulation

! ************************************************************************** !

subroutine StokesIFunctionAccumulation(r,realization,rho,ierr)
  !
  ! Accumulation terms for TS IFunction
  ! r += rho*vol*udot (velocity only)
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
  PetscReal, pointer :: r_p(:), xxdot_p(:)
  PetscInt :: lid, gid, ibeg
  PetscReal :: vol

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  field  => realization%field

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_xxdot_loc,xxdot_p,ierr);CHKERRQ(ierr)

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    ibeg = (lid-1)*option%nflowdof + 1
    vol  = patch%aux%Material%auxvars(gid)%volume

    r_p(ibeg    ) = r_p(ibeg    ) + rho * xxdot_p(ibeg    ) * vol
    r_p(ibeg + 1) = r_p(ibeg + 1) + rho * xxdot_p(ibeg + 1) * vol
    r_p(ibeg + 2) = r_p(ibeg + 2) + rho * xxdot_p(ibeg + 2) * vol
  enddo

  call VecRestoreArrayRead(field%flow_xxdot_loc,xxdot_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesIFunctionAccumulation

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

subroutine StokesJacobian(snes,xx,A,B,realization,rho,mu,add_convection,debug,ierr)
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
  PetscBool :: add_convection


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
  call StokesJacobianInternalConn(J,realization,mu,rho, &
                                            add_convection, &
                                            ierr)
  call StokesJacobianAccumulation(J,realization,rho,ierr)
  call StokesJacobianSourceSink(J,realization,rho,ierr)

  ! Assembly happens once at the end
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  call StokesJacobianPinPressure(J,realization,ierr)

  ! NOW apply BC constraints AFTER assembly
  call StokesJacobianBoundaryConn(J,realization,rho,mu,ierr)

  call MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)


  if (A /= B .and. mat_type_A /= MATMFFD) then
    call MatConvert(J,mat_type_B,MAT_REUSE_MATRIX,B,ierr);CHKERRQ(ierr)
  endif

  if (debug%matview_Matrix) then
    call DebugMatView(debug,J,'Sjacobian','',0,0,0,option)
  endif

end subroutine StokesJacobian

! ************************************************************************** !

subroutine StokesJacobianInternalConn(A,realization,mu,rho,add_convection,ierr)
  !
  ! Internal connection Jacobian for structured grids
  !
  !  - diffusion:   +a_mu(uK-uL)
  !  - grad(p):     +((pL-pK)/d)*S
  !  - convection:  +rho*phi*u_upwind (optional, fully implicit, upwind)
  !  - continuity:  +phi (K), -phi (L), with phi = 0.5(uK+uL)·S
  !  - p-stab:      -tau(pK-pL) in continuity
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
  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: mu, rho
  PetscBool :: add_convection
  PetscErrorCode :: ierr

  type(grid_type),          pointer :: grid
  type(patch_type),         pointer :: patch
  type(option_type),        pointer :: option
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type),      pointer :: cset

  PetscReal, pointer :: u_loc(:)

  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: istart_up, istart_dn
  PetscReal :: dx, area
  PetscReal :: n(3), S(3)

  PetscReal :: J_KK(4,4), J_KL(4,4), J_LK(4,4), J_LL(4,4)
  PetscReal :: a_mu, tau, gamma_pstab

  ! convection Jacobian
  PetscReal :: uK(3), uL(3), u_face(3), phi
  PetscReal :: dm_dUK(3,3), dm_dUL(3,3)
  PetscReal :: I3(3,3)
  PetscInt :: i,j

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  gamma_pstab = 0.05

  I3 = 0.d0
  do i=1,3
    I3(i,i)=1.d0
  enddo

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

      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)
      n    = 0.d0
      n(1:3) = cset%dist(1:3,iconn)

      S(1)=n(1)*area
      S(2)=n(2)*area
      S(3)=n(3)*area

      a_mu = mu * area / dx
      tau  = gamma_pstab * area / dx

      ! Fetch uK/uL (for convection Jacobian)
      uK = 0.d0
      uL = 0.d0
      if (lid_up > 0) then
        uK(1) = u_loc((lid_up-1)*option%nflowdof + 1)
        uK(2) = u_loc((lid_up-1)*option%nflowdof + 2)
        uK(3) = u_loc((lid_up-1)*option%nflowdof + 3)
      endif
      if (lid_dn > 0) then
        uL(1) = u_loc((lid_dn-1)*option%nflowdof + 1)
        uL(2) = u_loc((lid_dn-1)*option%nflowdof + 2)
        uL(3) = u_loc((lid_dn-1)*option%nflowdof + 3)
      endif

      u_face(1)=0.5d0*(uK(1)+uL(1))
      u_face(2)=0.5d0*(uK(2)+uL(2))
      u_face(3)=0.5d0*(uK(3)+uL(3))
      phi = u_face(1)*S(1) + u_face(2)*S(2) + u_face(3)*S(3)

      ! zero blocks
      J_KK = 0.d0
      J_KL = 0.d0
      J_LK = 0.d0
      J_LL = 0.d0

      ! ---------------------------------------------------
      ! Diffusion in momentum: +a_mu*(uK-uL)
      ! ---------------------------------------------------
      do i=1,3
        J_KK(i,i) = J_KK(i,i) + a_mu
        J_KL(i,i) = J_KL(i,i) - a_mu
        J_LK(i,i) = J_LK(i,i) - a_mu
        J_LL(i,i) = J_LL(i,i) + a_mu
      enddo

      ! ---------------------------------------------------
      ! Pressure gradient in momentum: +((pL - pK)/dx) * S
      ! dR^u_K/dpK = +(1/dx)*S, dR^u_K/dpL = -(1/dx)*S
      ! ---------------------------------------------------
      J_KK(1,4) = J_KK(1,4) - S(1)/dx
      J_KK(2,4) = J_KK(2,4) - S(2)/dx
      J_KK(3,4) = J_KK(3,4) - S(3)/dx
      J_KL(1,4) = J_KL(1,4) + S(1)/dx
      J_KL(2,4) = J_KL(2,4) + S(2)/dx
      J_KL(3,4) = J_KL(3,4) + S(3)/dx

      ! L-row gets opposite face contribution
      J_LK(1,4) = J_LK(1,4) + S(1)/dx
      J_LK(2,4) = J_LK(2,4) + S(2)/dx
      J_LK(3,4) = J_LK(3,4) + S(3)/dx

      J_LL(1,4) = J_LL(1,4) - S(1)/dx
      J_LL(2,4) = J_LL(2,4) - S(2)/dx
      J_LL(3,4) = J_LL(3,4) - S(3)/dx
      ! ---------------------------------------------------
      ! Continuity (divergence): R^p_K += phi, R^p_L -= phi
      ! phi = 0.5(uK+uL)·S
      ! ---------------------------------------------------
      J_KK(4,1) = J_KK(4,1) + 0.5d0*S(1)
      J_KK(4,2) = J_KK(4,2) + 0.5d0*S(2)
      J_KK(4,3) = J_KK(4,3) + 0.5d0*S(3)

      J_KL(4,1) = J_KL(4,1) + 0.5d0*S(1)
      J_KL(4,2) = J_KL(4,2) + 0.5d0*S(2)
      J_KL(4,3) = J_KL(4,3) + 0.5d0*S(3)

      J_LK(4,1) = J_LK(4,1) - 0.5d0*S(1)
      J_LK(4,2) = J_LK(4,2) - 0.5d0*S(2)
      J_LK(4,3) = J_LK(4,3) - 0.5d0*S(3)

      J_LL(4,1) = J_LL(4,1) - 0.5d0*S(1)
      J_LL(4,2) = J_LL(4,2) - 0.5d0*S(2)
      J_LL(4,3) = J_LL(4,3) - 0.5d0*S(3)

      ! ---------------------------------------------------
      ! Pressure stabilization in continuity: -tau(pK-pL)
      ! R^p_K adds: -tau*pK + tau*pL
      ! R^p_L adds:  tau*pK - tau*pL
      ! ---------------------------------------------------
      J_KK(4,4) = J_KK(4,4) - tau
      J_KL(4,4) = J_KL(4,4) + tau
      J_LK(4,4) = J_LK(4,4) + tau
      J_LL(4,4) = J_LL(4,4) - tau

      ! ---------------------------------------------------
      ! Convection in momentum (optional): R^u_K += m, R^u_L -= m
      ! m = rho*phi*u_upwind
      ! fully implicit upwind Jacobian (piecewise in sign(phi))
      ! ---------------------------------------------------
      if (add_convection) then
        dm_dUK = 0.d0
        dm_dUL = 0.d0

        if (phi >= 0.d0) then
          ! m = rho*phi*uK
          do i=1,3
            do j=1,3
              dm_dUK(i,j) = rho * ( phi*I3(i,j) + uK(i)*0.5d0*S(j) )
              dm_dUL(i,j) = rho * (               uK(i)*0.5d0*S(j) )
            enddo
          enddo
        else
          ! m = rho*phi*uL
          do i=1,3
            do j=1,3
              dm_dUK(i,j) = rho * (               uL(i)*0.5d0*S(j) )
              dm_dUL(i,j) = rho * ( phi*I3(i,j) + uL(i)*0.5d0*S(j) )
            enddo
          enddo
        endif

        ! add to momentum rows: K gets +dm, L gets -dm
        do i=1,3
          do j=1,3
            J_KK(i,j) = J_KK(i,j) + dm_dUK(i,j)
            J_KL(i,j) = J_KL(i,j) + dm_dUL(i,j)
            J_LK(i,j) = J_LK(i,j) - dm_dUK(i,j)
            J_LL(i,j) = J_LL(i,j) - dm_dUL(i,j)
          enddo
        enddo
      endif

      ! ---------------------------------------------------
      ! Insert blocks into PETSc matrix (local indexing)
      ! ---------------------------------------------------
      if (lid_up > 0) then
        istart_up = (lid_up-1)*option%nflowdof
        istart_dn = (lid_dn-1)*option%nflowdof

        do i=1,4
          do j=1,4
            call PUMSetValuesLocal(A,1,istart_up+i-1,1,istart_up+j-1, J_KK(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
            call PUMSetValuesLocal(A,1,istart_up+i-1,1,istart_dn+j-1, J_KL(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
          enddo
        enddo
      endif

      if (lid_dn > 0) then
        istart_up = (lid_up-1)*option%nflowdof
        istart_dn = (lid_dn-1)*option%nflowdof
        do i=1,4
          do j=1,4
            call PUMSetValuesLocal(A,1,istart_dn+i-1,1,istart_up+j-1, J_LK(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
            call PUMSetValuesLocal(A,1,istart_dn+i-1,1,istart_dn+j-1, J_LL(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
          enddo
        enddo
      endif

    enddo

    cset => cset%next
  enddo

  call VecRestoreArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)

end subroutine StokesJacobianInternalConn

! ************************************************************************** !
subroutine StokesJacobianBoundaryConn(A,realization,rho,mu,ierr)
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

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  PetscErrorCode :: ierr
  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type),pointer :: bc
  type(connection_set_type), pointer :: cset
  PetscInt :: iconn, local_id, gid, iblock, lid, ndof
  PetscReal :: Jbc(4,4)
  ! Track which rows have been constrained for this cell
  PetscBool, allocatable :: bc_applied(:,:)
  PetscInt, allocatable :: bc_rows(:)
  PetscInt :: bc_row_count
  PetscInt :: max_local_id, irow

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  ndof = option%nflowdof

  ! Allocate tracking array for which DOFs have BCs
  max_local_id = grid%nlmax
  allocate(bc_applied(ndof, max_local_id))
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
        endif
      endif
    enddo
    bc => bc%next
  enddo

  ! Allocate array to store all BC rows
  allocate(bc_rows(grid%nlmax * ndof))  ! Max 4 rows per cell
  bc_row_count = 0

  ! Second pass: collect all constrained rows
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    if (.not.any(bc_applied(:, lid))) cycle

    if (bc_applied(1, lid)) then
      bc_row_count = bc_row_count + 1
      bc_rows(bc_row_count) = (grid%global_offset + lid - 1)*ndof + 0  ! u_x (global)
      bc_row_count = bc_row_count + 1
      bc_rows(bc_row_count) = (grid%global_offset + lid - 1)*ndof + 1  ! u_y (global)
      bc_row_count = bc_row_count + 1
      bc_rows(bc_row_count) = (grid%global_offset + lid - 1)*ndof + 2  ! u_z (global)
    endif

    if (bc_applied(4, lid)) then
      bc_row_count = bc_row_count + 1
      bc_rows(bc_row_count) = (grid%global_offset + lid - 1)*ndof + 3  ! p (global)
    endif
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatZeroRows(A,bc_row_count,bc_rows,1.d0, &
                   PETSC_NULL_VEC,PETSC_NULL_VEC,ierr);CHKERRQ(ierr)

  deallocate(bc_applied)
  deallocate(bc_rows)

end subroutine StokesJacobianBoundaryConn

! ************************************************************************** !

subroutine StokesJacobianAccumulation(A,realization,rho,ierr)
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

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  PetscErrorCode :: ierr
  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  PetscInt :: lid, gid, istart
  PetscReal :: Jacc(4,4), vol, vol_over_dt
  PetscInt :: i,j

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  if (option%flow%steady_state) return

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    istart = (lid - 1) * option%nflowdof + 1  !  index for BAIJ matrixBlock
    vol    = patch%aux%Material%auxvars(gid)%volume
    vol_over_dt = vol / option%flow_dt

    Jacc = 0.d0
    Jacc(1,1) = rho * vol_over_dt
    Jacc(2,2) = rho * vol_over_dt
    Jacc(3,3) = rho * vol_over_dt

    do i = 1,4
      do j = 1,4
        if (i /= j) cycle
        if (i == 4) cycle  ! No accumulation for pressure
          call PUMSetValuesLocal(A,1,istart+i-1,1,istart+j-1,Jacc(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    enddo

  enddo

end subroutine StokesJacobianAccumulation

! ************************************************************************** !

subroutine StokesIJacobianAccumulation(A,shift,realization,rho,ierr)
  !
  ! Accumulation Jacobian for TS IJacobian
  ! Adds shift*rho*vol on velocity DOFs
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Petsc_Utility_module

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: shift, rho
  PetscErrorCode :: ierr
  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(option_type), pointer :: option
  PetscInt :: lid, gid, istart
  PetscReal :: vol, val

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    istart = (lid - 1) * option%nflowdof
    vol    = patch%aux%Material%auxvars(gid)%volume
    val    = shift * rho * vol

    call PUMSetValuesLocal(A,1,istart+0,1,istart+0,val,ADD_VALUES,ierr);CHKERRQ(ierr)
    call PUMSetValuesLocal(A,1,istart+1,1,istart+1,val,ADD_VALUES,ierr);CHKERRQ(ierr)
    call PUMSetValuesLocal(A,1,istart+2,1,istart+2,val,ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

end subroutine StokesIJacobianAccumulation

! ************************************************************************** !

subroutine StokesJacobianSourceSink(A,realization,rho,ierr)
  !
  ! Source/sink Jacobian for Stokes problem
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !

  use Realization_Subsurface_class

  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  PetscErrorCode :: ierr

end subroutine StokesJacobianSourceSink

! ************************************************************************** !

subroutine StokesUpdateFixedAccum(realization,rho)
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
  PetscReal :: vol, vol_over_dt, rho
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
    vol_over_dt = vol / option%flow_dt

    ! accum stores: rho * volume * u^n (previous velocity)
    ! This will be subtracted from residual to get: rho*vol*(u^{n+1} - u^n)/dt
    accum_p(istart    ) = rho * xx_p(istart    ) * vol_over_dt
    accum_p(istart + 1) = rho * xx_p(istart + 1) * vol_over_dt
    accum_p(istart + 2) = rho * xx_p(istart + 2) * vol_over_dt
  enddo

  call VecRestoreArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

end subroutine StokesUpdateFixedAccum

! ************************************************************************** !
subroutine StokesResidualBoundaryConn(r,realization,rho,mu,ierr)
  !
  ! Boundary face contributions for flux-based BCs (velocity Dirichlet via viscous flux,
  ! pressure BC via traction).
  !
  ! Author: Satish Karra
  ! Date: 01/28/26
  !
  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho, mu
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cset

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: u_loc(:)

  PetscInt :: iconn, local_id, gid, istart, ndof
  PetscReal :: Ucell(4), Res(4)

  patch => realization%patch
  grid  => patch%grid
  option => realization%option
  ndof = option%nflowdof

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)

  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    cset => boundary_condition%connection_set

    do iconn = 1, cset%num_connections

      local_id = cset%id_dn(iconn)
      gid = grid%nL2G(local_id)

      if (patch%imat(gid) <= 0) cycle

      if (gid <= 0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      istart = (local_id-1)*ndof + 1

      Ucell(1) = u_loc(istart    )
      Ucell(2) = u_loc(istart + 1)
      Ucell(3) = u_loc(istart + 2)
      Ucell(4) = u_loc(istart + 3)

      call StokesBCFluxConn(boundary_condition,iconn, &
                            cset%area(iconn),cset%dist(:,iconn), &
                            rho,mu,option,Ucell,Res)

      r_p(istart    ) = r_p(istart    ) + Res(1)
      r_p(istart + 1) = r_p(istart + 1) + Res(2)
      r_p(istart + 2) = r_p(istart + 2) + Res(3)
      r_p(istart + 3) = r_p(istart + 3) + Res(4)

    enddo

    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

end subroutine StokesResidualBoundaryConn

! ************************************************************************** !

subroutine StokesBCFluxConn(boundary_condition,iconn,area,dist,rho,mu,option,Ucell,Res)
  !
  ! Compute boundary-face residual contribution for ONE boundary connection (ONE face).
  !
  ! Flux-based velocity Dirichlet:
  !   momentum:  +a_b (uK - u_b),  a_b = mu*A/d_b
  !   continuity: + (u_b · S)
  !
  ! Pressure traction (Option B):
  !   momentum:  - p_b * S
  !
  ! No pressure Dirichlet equation is imposed here.
  ! Author: Satish Karra
  ! Date: 01/24/26
  !
  !
  use Coupler_module
  use Option_module
  use PFLOTRAN_Constants_module
  implicit none

  type(coupler_type), intent(in) :: boundary_condition
  PetscInt, intent(in) :: iconn
  PetscReal, intent(in) :: area
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(in) :: rho, mu
  type(option_type), intent(in) :: option
  PetscReal, intent(in) :: Ucell(4)
  PetscReal, intent(out) :: Res(4)

  PetscReal :: n(3), S(3), d_b, a_b
  PetscReal :: uK(3), pK
  PetscReal :: u_b(3), p_b
  PetscReal :: phi_b
  PetscBool :: has_u_dir, has_p_trac

  Res = 0.d0
  u_b = 0.d0
  p_b = 0.d0

  n(1)=dist(1); n(2)=dist(2); n(3)=dist(3)
  S(1)=n(1)*area; S(2)=n(2)*area; S(3)=n(3)*area

  d_b = dist(0)

  if (d_b <= 0.d0) then
    print *, 'BAD boundary distance: dist(-1:0)=',dist(-1),dist(0)
    stop
  endif

  uK(1)=Ucell(1); uK(2)=Ucell(2); uK(3)=Ucell(3)
  pK = Ucell(4)

  has_u_dir  = PETSC_FALSE
  has_p_trac = PETSC_FALSE

  ! Velocity Dirichlet present?
  if (associated(boundary_condition%flow_condition%stokes%velocity)) then
    if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
      has_u_dir = PETSC_TRUE
      u_b(1) = boundary_condition%flow_aux_real_var(1,iconn)
      u_b(2) = boundary_condition%flow_aux_real_var(2,iconn)
      u_b(3) = boundary_condition%flow_aux_real_var(3,iconn)
    endif
  endif

  ! Pressure traction present? (treat LIQUID_PRESSURE DIRICHLET as traction value p_b)
  if (associated(boundary_condition%flow_condition%stokes%pressure)) then
    if (boundary_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
      if (.not.has_u_dir) then
        p_b = boundary_condition%flow_aux_real_var(4,iconn)
        has_p_trac = PETSC_TRUE
      endif
    endif
  endif

  ! ---- Flux-based velocity Dirichlet via viscous diffusion ----
  if (has_u_dir) then
    a_b = mu * area / d_b
    Res(1) = Res(1) + a_b*(uK(1) - u_b(1))
    Res(2) = Res(2) + a_b*(uK(2) - u_b(2))
    Res(3) = Res(3) + a_b*(uK(3) - u_b(3))

    ! continuity boundary flux uses prescribed u_b
    phi_b = u_b(1)*S(1) + u_b(2)*S(2) + u_b(3)*S(3)
    Res(4) = Res(4) - phi_b
  else
    ! No velocity Dirichlet: use "do-nothing/outflow" flux for continuity
    phi_b = uK(1)*S(1) + uK(2)*S(2) + uK(3)*S(3)
    Res(4) = Res(4) - phi_b
  endif

  ! ---- Pressure traction ----
  if (has_p_trac) then
    Res(1) = Res(1) - (pK - p_b)*S(1)/d_b
    Res(2) = Res(2) - (pK - p_b)*S(2)/d_b
    Res(3) = Res(3) - (pK - p_b)*S(3)/d_b
  endif

end subroutine StokesBCFluxConn

! ************************************************************************** !

subroutine StokesResidualBoundaryFinal(r,xx,realization,ierr)
  !
  ! Sets residuals for Dirichlet BC rows to (u - u_prescribed)
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !
  !
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Option_module
  use Connection_module
  use Field_module

  implicit none

  Vec :: r, xx
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id, iconn, offset
  PetscReal, pointer :: r_ptr(:), xx_ptr(:)
  PetscInt :: nflowdof

  patch => realization%patch
  option => realization%option
  nflowdof = option%nflowdof

  if (.not.associated(patch%boundary_condition_list%first)) return

  call VecGetArray(r, r_ptr, ierr);CHKERRQ(ierr)
  call VecGetArrayRead(xx, xx_ptr, ierr);CHKERRQ(ierr)

  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      offset = (local_id-1)*nflowdof

      ! Set velocity BC residual to (u - u_prescribed)
      if (associated(boundary_condition%flow_condition%stokes%velocity)) then
        if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          r_ptr(offset+1) = xx_ptr(offset+1) - boundary_condition%flow_aux_real_var(1,iconn)
          r_ptr(offset+2) = xx_ptr(offset+2) - boundary_condition%flow_aux_real_var(2,iconn)
          r_ptr(offset+3) = xx_ptr(offset+3) - boundary_condition%flow_aux_real_var(3,iconn)
        endif
      endif

      ! Set pressure BC residual to (p - p_prescribed)
      if (associated(boundary_condition%flow_condition%stokes%pressure)) then
        if (boundary_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          r_ptr(offset+4) = xx_ptr(offset+4) - boundary_condition%flow_aux_real_var(4,iconn)
        endif
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayRead(xx, xx_ptr, ierr);CHKERRQ(ierr)
  call VecRestoreArray(r, r_ptr, ierr);CHKERRQ(ierr)

end subroutine StokesResidualBoundaryFinal

! ************************************************************************** !
subroutine StokesJacobianBoundaryFacesConn(A,realization,mu,ierr)
  !
  ! Boundary face Jacobian contributions for flux-based BCs.
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !
  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Petsc_Utility_module
  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscReal :: mu
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cset

  PetscInt :: iconn, local_id, gid, istart, ndof
  PetscReal :: area, dist(-1:3), d_b, a_b
  PetscBool :: has_u_dir, has_p_trac
  PetscReal :: S(3)

  patch => realization%patch
  grid  => patch%grid
  option => realization%option
  ndof = option%nflowdof

  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    cset => boundary_condition%connection_set
    if (.not.associated(cset)) then
      boundary_condition => boundary_condition%next
      cycle
    endif

    do iconn = 1, cset%num_connections

      local_id = cset%id_dn(iconn)
      gid = grid%nL2G(local_id)


      if (patch%imat(gid) <= 0) cycle

      if (gid <= 0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      istart = (local_id-1)*ndof   ! 0-based for PUMSetValuesLocal

      area = cset%area(iconn)
      dist = cset%dist(:,iconn)

      d_b = dist(0)

      ! S = n*A (only needed if you use do-nothing flux Jacobian)
      S(1) = dist(1)*area
      S(2) = dist(2)*area
      S(3) = dist(3)*area

      has_u_dir = PETSC_FALSE
      has_p_trac = PETSC_FALSE
      if (associated(boundary_condition%flow_condition%stokes%velocity)) then
        if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          has_u_dir = PETSC_TRUE
        endif
      endif

      if (associated(boundary_condition%flow_condition%stokes%pressure)) then
        if (boundary_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          if (.not.has_u_dir) then
            has_p_trac = PETSC_TRUE
          endif
        endif
      endif

      ! If velocity Dirichlet: momentum has +a_b(uK-u_b) => diag += a_b
      if (has_u_dir) then
        a_b = mu * area / d_b
        call PUMSetValuesLocal(A,1,istart+0,1,istart+0,a_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,istart+1,1,istart+1,a_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,istart+2,1,istart+2,a_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        ! continuity uses prescribed flux u_b·S -> no Jacobian
      else
        ! do-nothing continuity flux: phi = -uK·S -> J(p,u) -= S^T

        call PUMSetValuesLocal(A,1,istart+3,1,istart+0,-S(1),ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,istart+3,1,istart+1,-S(2),ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,istart+3,1,istart+2,-S(3),ADD_VALUES,ierr);CHKERRQ(ierr)
      endif

      if (has_p_trac) then
        ! If pressure traction: momentum has -(pK - p_b)*S/d_b => diag -= S/d_b
        call PUMSetValuesLocal(A,1,istart+0,1,istart+3,-S(1)/d_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,istart+1,1,istart+3,-S(2)/d_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,istart+2,1,istart+3,-S(3)/d_b,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif

    enddo

    boundary_condition => boundary_condition%next
  enddo

end subroutine StokesJacobianBoundaryFacesConn

! ************************************************************************** !
subroutine StokesResidualPinPressure(r,xx,realization,p_ref,ierr)
  !
  ! Pins pressure at the first active local cell:
  !   r_p(row) = p(row) - p_ref
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  implicit none

  Vec :: r, xx
  class(realization_subsurface_type) :: realization
  PetscReal :: p_ref
  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(grid_type),  pointer :: grid
  type(option_type),pointer :: option

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: xx_p(:)

  PetscInt :: lid, gid, ndof, row
  PetscInt :: pin_gid_local, pin_gid_global

  patch => realization%patch
  grid  => patch%grid
  option => realization%option
  ndof = option%nflowdof

  call VecGetArray(r, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayRead(xx, xx_p, ierr);CHKERRQ(ierr)

  pin_gid_local = HUGE(1)
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (gid <= 0) cycle
    if (patch%imat(gid) <= 0) cycle
    pin_gid_local = grid%global_offset + lid
    exit
  enddo

  call MPI_Allreduce(pin_gid_local,pin_gid_global,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MIN,option%mycomm,ierr);CHKERRQ(ierr)

  row = -1
  if (pin_gid_global > 0 .and. pin_gid_global < HUGE(1)) then
    if (pin_gid_global > grid%global_offset .and. &
        pin_gid_global <= grid%global_offset + grid%nlmax) then
      row = (pin_gid_global - grid%global_offset - 1)*ndof + 4
    endif
  endif

  if (row > 0) then
    r_p(row) = xx_p(row) - p_ref
  endif

  call VecRestoreArrayRead(xx, xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArray(r, r_p, ierr);CHKERRQ(ierr)

end subroutine StokesResidualPinPressure

! ************************************************************************** !

subroutine StokesJacobianPinPressure(A,realization,ierr)
  !
  ! Pins pressure at the first active local cell by replacing that row with identity:
  ! A(row,:) = 0; A(row,row)=1
  !
  ! Author: Satish Karra
  ! Date: 01/24/26
  !
  !
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  implicit none

  Mat, intent(inout) :: A
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(grid_type),  pointer :: grid
  type(option_type),pointer :: option

  PetscInt :: lid, gid, ndof
  PetscInt :: row_global(1)
  PetscInt :: pin_gid_local, pin_gid_global
  PetscBool :: found

  patch => realization%patch
  grid  => patch%grid
  option => realization%option
  ndof = option%nflowdof

  pin_gid_local = HUGE(1)
  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (gid <= 0) cycle
    if (patch%imat(gid) <= 0) cycle
    pin_gid_local = grid%global_offset + lid
    exit
  enddo

  call MPI_Allreduce(pin_gid_local,pin_gid_global,ONE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MIN,option%mycomm,ierr);CHKERRQ(ierr)

  found = PETSC_FALSE
  row_global(1) = -1
  if (pin_gid_global > 0 .and. pin_gid_global < HUGE(1)) then
    if (pin_gid_global > grid%global_offset .and. &
        pin_gid_global <= grid%global_offset + grid%nlmax) then
      row_global(1) = (pin_gid_global - 1)*ndof + 3   ! global row for pressure dof
      found = PETSC_TRUE
    endif
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  if (found) then
    call MatZeroRows(A,1,row_global,1.d0,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  else
    call MatZeroRows(A,0,row_global,1.d0,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr);CHKERRQ(ierr)
  endif

end subroutine StokesJacobianPinPressure


end module Stokes_module
