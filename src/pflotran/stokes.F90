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
  use Matrix_Zeroing_module

  implicit none

  private

  public :: StokesSetup, &
            StokesCreateZeroArray, &
            StokesInitializeTimestep, &
            StokesUpdateAuxVars, &
            StokesMapBCAuxVarsToGlobal, &
            StokesUpdateTransportVelocities, &
            StokesUpdateSolution, &
            StokesMaxChange, &
            StokesDestroy, &
            StokesComputeViscousDiag, &
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
            StokesApplyBoundaryConditions, &
            StokesApplyPrescribedConditions

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
  use String_module
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
  PetscInt :: iconn, sum_connection
  type(coupler_type), pointer :: coupler

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

  call StokesCreateZeroArray(patch,patch%aux%Stokes%matrix_zeroing,option)

  ! ensure that prescribed_conditions are solely DIRICHLET type
  coupler => patch%prescribed_condition_list%first
  do
    if (.not.associated(coupler)) exit
    if (associated(coupler%flow_condition%stokes%velocity)) then
      select case(coupler%flow_condition%stokes%velocity%itype)
        case(DIRICHLET_BC)
        case default
          option%io_buffer = 'FLOW_CONDITION "' // &
            trim(coupler%flow_condition%name) // '" STOKES VELOCITY TYPE (' // &
            StringWrite(coupler%flow_condition%stokes%velocity%itype) // &
            ') not supported in STOKES PRESCRIBED_CONDITION "' // &
            coupler%name // '".'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(coupler%flow_condition%stokes%pressure)) then
      select case(coupler%flow_condition%stokes%pressure%itype)
        case(DIRICHLET_BC)
        case default
          option%io_buffer = 'FLOW_CONDITION "' // &
            trim(coupler%flow_condition%name) // '" STOKES PRESSURE TYPE (' // &
            StringWrite(coupler%flow_condition%stokes%pressure%itype) // &
            ') not supported in STOKES PRESCRIBED_CONDITION "' // &
            coupler%name // '".'
          call PrintErrMsg(option)
      end select
    endif
    coupler => coupler%next
  enddo

  ! Allocate cached viscous diagonal for Rhie-Chow stabilization
  if (.not.associated(stokes_aux%aP)) then
    allocate(stokes_aux%aP(grid%ngmax))
    stokes_aux%aP = 0.d0
  endif

end subroutine StokesSetup

! ************************************************************************** !

subroutine StokesCreateZeroArray(patch,matrix_zeroing,option)
  !
  ! Builds the matrix_zeroing object for Stokes with DOF-selective zeroing:
  !   - Inactive cells (imat <= 0): zero all ndof DOFs
  !   - Prescribed-condition cells: zero only the prescribed DOFs
  !       velocity DIRICHLET -> DOFs 1,2,3 (u,v,w)
  !       pressure DIRICHLET -> DOF  4     (p)
  !
  ! This is a Stokes-specific replacement for PatchCreateZeroArray, which
  ! incorrectly zeroes all DOFs for every prescribed-condition cell regardless
  ! of which DOFs are actually prescribed.
  !
  ! Author: Satish Karra
  ! Date: 03/18/26
  !

  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module

  implicit none

  type(patch_type), pointer :: patch
  type(matrix_zeroing_type), pointer :: matrix_zeroing
  type(option_type) :: option

  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: cur_coupler
  type(connection_set_type), pointer :: cur_connection_set
  ! per-local-cell bitmask: bit i set means DOF i is to be zeroed
  PetscBool, allocatable :: zero_dof(:,:)
  PetscInt :: local_id, ghosted_id, iconn, idof, ndof, ncount, n_zero_rows

  grid => patch%grid
  ndof = option%nflowdof

  allocate(zero_dof(ndof, grid%nlmax))
  zero_dof = PETSC_FALSE

  ! Inactive cells: zero all DOFs
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) then
      zero_dof(:, local_id) = PETSC_TRUE
    endif
  enddo

  ! Prescribed-condition cells: zero only the prescribed DOFs
  cur_coupler => patch%prescribed_condition_list%first
  do
    if (.not.associated(cur_coupler)) exit
    if (.not.associated(cur_coupler%flow_condition)) then
      cur_coupler => cur_coupler%next
      cycle
    endif
    cur_connection_set => cur_coupler%connection_set
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      if (local_id <= 0 .or. local_id > grid%nlmax) cycle
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle  ! already handled above

      ! Velocity DIRICHLET -> zero u, v, w (DOFs 1,2,3)
      if (associated(cur_coupler%flow_condition%stokes%velocity)) then
        if (cur_coupler%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          zero_dof(1, local_id) = PETSC_TRUE
          zero_dof(2, local_id) = PETSC_TRUE
          zero_dof(3, local_id) = PETSC_TRUE
        endif
      endif

      ! Pressure DIRICHLET -> zero p (DOF 4)
      if (associated(cur_coupler%flow_condition%stokes%pressure)) then
        if (cur_coupler%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          zero_dof(4, local_id) = PETSC_TRUE
        endif
      endif
    enddo
    cur_coupler => cur_coupler%next
  enddo

  ! Count total zero rows
  n_zero_rows = 0
  do local_id = 1, grid%nlmax
    do idof = 1, ndof
      if (zero_dof(idof, local_id)) n_zero_rows = n_zero_rows + 1
    enddo
  enddo

  call MatrixZeroingAllocateArray(matrix_zeroing, n_zero_rows, option)

  ncount = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    do idof = 1, ndof
      if (.not.zero_dof(idof, local_id)) cycle
      ncount = ncount + 1
      ! 1-based indexing
      matrix_zeroing%zero_rows_local(ncount) = (local_id-1)*ndof + idof
      ! 0-based indexing
      matrix_zeroing%zero_rows_local_ghosted(ncount) = &
        (ghosted_id-1)*ndof + idof - 1
    enddo
  enddo

  deallocate(zero_dof)

end subroutine StokesCreateZeroArray

! ************************************************************************** !

subroutine StokesComputeViscousDiag(realization,mu)
  !
  ! Computes the viscous diagonal aP(gid) = sum_faces (mu * A / dx)
  ! over internal connections.  Result is stored in
  ! patch%aux%Stokes%aP(:) (pre-allocated in StokesSetup).
  !
  ! Must be called once per residual/Jacobian evaluation BEFORE any
  ! routine that uses the Rhie-Chow pressure stabilization.
  !
  ! Author: Satish Karra
  ! Date: 03/17/26
  !

  use Realization_Subsurface_class
  use Grid_module
  use Patch_module
  use Connection_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal, intent(in) :: mu

  type(patch_type),  pointer :: patch
  type(grid_type),   pointer :: grid
  type(stokes_type), pointer :: stokes_aux
  type(connection_set_list_type), pointer :: conn_list
  type(connection_set_type),      pointer :: cset

  PetscInt  :: iconn, gid_up, gid_dn
  PetscReal :: area, dx, a_mu

  patch => realization%patch
  grid  => patch%grid
  if (.not.associated(patch%aux%Stokes)) return
  stokes_aux => patch%aux%Stokes

  stokes_aux%aP = 0.d0

  conn_list => grid%internal_connection_set_list
  cset      => conn_list%first
  do
    if (.not.associated(cset)) exit
    do iconn = 1, cset%num_connections
      gid_up = cset%id_up(iconn)
      gid_dn = cset%id_dn(iconn)
      if (patch%imat(gid_up) <= 0 .or. patch%imat(gid_dn) <= 0) cycle
      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)
      a_mu = mu * area / dx
      stokes_aux%aP(gid_up) = stokes_aux%aP(gid_up) + a_mu
      stokes_aux%aP(gid_dn) = stokes_aux%aP(gid_dn) + a_mu
    enddo
    cset => cset%next
  enddo

end subroutine StokesComputeViscousDiag

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
    istart = (gid-1)*option%nflowdof + 1
    Ucell(1) = xx_loc_p(istart    )
    Ucell(2) = xx_loc_p(istart + 1)
    Ucell(3) = xx_loc_p(istart + 2)
    Ucell(4) = xx_loc_p(istart + 3)
    call StokesAuxVarCompute(Ucell, &
                             stokes_aux%auxvars(gid), &
                             global_auxvars(gid), &
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
  ! Updates data in module after a successful time step.
  !
  ! Currently a no-op: transport velocities are updated via
  ! PMStokesUpdateAuxVars -> StokesUpdateTransportVelocities.
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
  ! Computes maximum change in solution between current (flow_xx) and
  ! previous (flow_yy) time step:  du_max = ||flow_yy - flow_xx||_inf
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
  ! flow_dxx = -1 * flow_xx + flow_yy  (i.e. flow_yy - flow_xx)
  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecNorm(field%flow_dxx,NORM_INFINITY,du_max,ierr);CHKERRQ(ierr)

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
  use Patch_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(patch_type), pointer :: patch

  patch => realization%patch
  if (associated(patch%aux%Stokes)) then
    call StokesAuxDestroy(patch%aux%Stokes)
  endif

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
  call StokesComputeViscousDiag(realization,mu)
  call StokesResidualInternalConn(r,realization,rho,mu, &
                                            add_convection,ierr)
  call StokesResidualBoundaryConn(r,realization,rho,mu,ierr)
  call StokesResidualAccumulation(r,realization,rho,ierr)
  call StokesResidualSourceSink(r,realization,rho,ierr)
  call StokesResidualPinPressure(r,xx,realization,0.d0,ierr)
  call StokesResidualBoundaryFinal(r,xx,realization,ierr)
  if (associated(realization%patch%aux%Stokes)) then
    call MatrixZeroingZeroVecEntries(realization%patch%aux%Stokes%matrix_zeroing,r)
  endif


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
  ! Backwards-compatible wrapper for prescribed conditions
  !
  use Realization_Subsurface_class
  implicit none

  class(realization_subsurface_type) :: realization

  call StokesApplyPrescribedConditions(realization)

end subroutine StokesApplyBoundaryConditions

! ************************************************************************** !

subroutine StokesApplyPrescribedConditions(realization)
  !
  ! Apply prescribed Dirichlet conditions to solution vector
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
  type(coupler_type), pointer :: prescribed_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id, iconn, offset, ghosted_id
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  option => realization%option

  if (.not.associated(patch%prescribed_condition_list%first)) return

  call VecGetArray(realization%field%flow_xx,vec_ptr,ierr);CHKERRQ(ierr)

  prescribed_condition => patch%prescribed_condition_list%first
  do
    if (.not.associated(prescribed_condition)) exit
    cur_connection_set => prescribed_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      if (local_id <= 0 .or. local_id > patch%grid%nlmax) cycle
      ghosted_id = patch%grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. ghosted_id > patch%grid%ngmax) cycle
      offset = (local_id-1)*option%nflowdof

      ! Apply prescribed velocity (all 3 components)
      if (associated(prescribed_condition%flow_condition%stokes%velocity)) then
        if (prescribed_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          vec_ptr(offset+1) = prescribed_condition%flow_aux_real_var(1,iconn)
          vec_ptr(offset+2) = prescribed_condition%flow_aux_real_var(2,iconn)
          vec_ptr(offset+3) = prescribed_condition%flow_aux_real_var(3,iconn)
        endif
      endif

      ! Apply prescribed pressure
      if (associated(prescribed_condition%flow_condition%stokes%pressure)) then
        if (prescribed_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          vec_ptr(offset+4) = prescribed_condition%flow_aux_real_var(4,iconn)
        endif
      endif

    enddo
    prescribed_condition => prescribed_condition%next
  enddo

  call VecRestoreArray(realization%field%flow_xx,vec_ptr,ierr);CHKERRQ(ierr)


  call DiscretizationGlobalToLocal(realization%discretization, &
                                   realization%field%flow_xx, &
                                   realization%field%flow_xx_loc, &
                                   realization%option%nflowdof)


end subroutine StokesApplyPrescribedConditions

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
  type(stokes_type),        pointer :: stokes_aux

  PetscReal, pointer :: r_p(:)

  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: ibeg_up_loc, ibeg_dn_loc

  PetscReal :: dx, area
  PetscReal :: n(3), S(3)
  PetscReal :: u_up(3), u_dn(3), u_face(3), u_upwind(3)
  PetscReal :: p_up, p_dn
  PetscReal :: a_mu, phi
  PetscReal :: Tvisc(3), Tpress(3), mconv(3)
  PetscReal :: flux_div
  PetscReal :: gamma_pstab
  PetscReal, pointer :: aP(:)
  PetscReal :: a_f, rc

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  gamma_pstab = 0.05  ! unused when Rhie-Chow is enabled

  if (.not.associated(patch%aux%Stokes)) return
  stokes_aux => patch%aux%Stokes

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

  ! Use pre-computed viscous diagonal aP for Rhie-Chow
  aP => stokes_aux%aP

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
      if ((lid_up <= 0 .or. lid_up > grid%nlmax) .and. &
          (lid_dn <= 0 .or. lid_dn > grid%nlmax)) cycle

      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)
      n    = 0.d0
      n(1:3) = cset%dist(1:3,iconn)   ! outward normal from "up" cell (K)

      S(1) = n(1)*area
      S(2) = n(2)*area
      S(3) = n(3)*area

      u_up(1) = stokes_aux%auxvars(gid_up)%u
      u_up(2) = stokes_aux%auxvars(gid_up)%v
      u_up(3) = stokes_aux%auxvars(gid_up)%w
      p_up    = stokes_aux%auxvars(gid_up)%pressure

      u_dn(1) = stokes_aux%auxvars(gid_dn)%u
      u_dn(2) = stokes_aux%auxvars(gid_dn)%v
      u_dn(3) = stokes_aux%auxvars(gid_dn)%w
      p_dn    = stokes_aux%auxvars(gid_dn)%pressure
      if (u_up(1) /= u_up(1) .or. u_up(2) /= u_up(2) .or. u_up(3) /= u_up(3) .or. &
          p_up /= p_up .or. u_dn(1) /= u_dn(1) .or. u_dn(2) /= u_dn(2) .or. &
          u_dn(3) /= u_dn(3) .or. p_dn /= p_dn) then
        write(option%io_buffer,'("StokesResidualInternalConn: NaN in auxvars ",2i8)') gid_up, gid_dn
        call PrintErrMsg(option)
      endif

      ! Face velocity (central)
      u_face(1) = 0.5d0*(u_up(1) + u_dn(1))
      u_face(2) = 0.5d0*(u_up(2) + u_dn(2))
      u_face(3) = 0.5d0*(u_up(3) + u_dn(3))

      ! Volumetric flux phi = u_face · S (from up(K) to dn(L))
      phi = u_face(1)*S(1) + u_face(2)*S(2) + u_face(3)*S(3)

      ! Diffusion coefficient
      a_mu = mu * area / dx

      ! Diffusion contribution: +a_mu (uK - uL)
      Tvisc(1) = a_mu*(u_up(1) - u_dn(1))
      Tvisc(2) = a_mu*(u_up(2) - u_dn(2))
      Tvisc(3) = a_mu*(u_up(3) - u_dn(3))

      ! Pressure contribution: +p_face * S, with p_face = 0.5(pK+pL)
      Tpress(1) = 0.5d0*(p_up + p_dn) * S(1)
      Tpress(2) = 0.5d0*(p_up + p_dn) * S(2)
      Tpress(3) = 0.5d0*(p_up + p_dn) * S(3)

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
      ! Rhie-Chow pressure correction: - (A^2/a_f) * (pL - pK)/dx
      a_f = 0.5d0*(aP(gid_up) + aP(gid_dn))
      if (a_f <= 0.d0) a_f = a_mu
      rc = (area*area/a_f) * (p_dn - p_up) / dx
      flux_div = phi - rc

      ! Assemble into residual vector
      if (lid_up > 0 .and. lid_up <= grid%nlmax) then
        ibeg_up_loc = (lid_up-1)*option%nflowdof + 1
        r_p(ibeg_up_loc    ) = r_p(ibeg_up_loc    ) + Tvisc(1) + Tpress(1) + mconv(1)
        r_p(ibeg_up_loc + 1) = r_p(ibeg_up_loc + 1) + Tvisc(2) + Tpress(2) + mconv(2)
        r_p(ibeg_up_loc + 2) = r_p(ibeg_up_loc + 2) + Tvisc(3) + Tpress(3) + mconv(3)
        r_p(ibeg_up_loc + 3) = r_p(ibeg_up_loc + 3) + flux_div
      endif

      if (lid_dn > 0 .and. lid_dn <= grid%nlmax) then
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        r_p(ibeg_dn_loc    ) = r_p(ibeg_dn_loc    ) - Tvisc(1) - Tpress(1) - mconv(1)
        r_p(ibeg_dn_loc + 1) = r_p(ibeg_dn_loc + 1) - Tvisc(2) - Tpress(2) - mconv(2)
        r_p(ibeg_dn_loc + 2) = r_p(ibeg_dn_loc + 2) - Tvisc(3) - Tpress(3) - mconv(3)
        r_p(ibeg_dn_loc + 3) = r_p(ibeg_dn_loc + 3) - flux_div
      endif

    enddo

    cset => cset%next
  enddo

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
  type(stokes_type),        pointer :: stokes_aux

  PetscReal, pointer :: r_p(:)

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

  if (.not.associated(patch%aux%Stokes)) return
  stokes_aux => patch%aux%Stokes

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

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
      if ((lid_up <= 0 .or. lid_up > grid%nlmax) .and. &
          (lid_dn <= 0 .or. lid_dn > grid%nlmax)) cycle

      area = cset%area(iconn)
      n    = 0.d0
      n(1:3) = cset%dist(1:3,iconn)

      S(1) = n(1)*area
      S(2) = n(2)*area
      S(3) = n(3)*area

      u_up(1) = stokes_aux%auxvars(gid_up)%u
      u_up(2) = stokes_aux%auxvars(gid_up)%v
      u_up(3) = stokes_aux%auxvars(gid_up)%w

      u_dn(1) = stokes_aux%auxvars(gid_dn)%u
      u_dn(2) = stokes_aux%auxvars(gid_dn)%v
      u_dn(3) = stokes_aux%auxvars(gid_dn)%w

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

      if (lid_up > 0 .and. lid_up <= grid%nlmax) then
        ibeg_up_loc = (lid_up-1)*option%nflowdof + 1
        r_p(ibeg_up_loc    ) = r_p(ibeg_up_loc    ) + mconv(1)
        r_p(ibeg_up_loc + 1) = r_p(ibeg_up_loc + 1) + mconv(2)
        r_p(ibeg_up_loc + 2) = r_p(ibeg_up_loc + 2) + mconv(3)
      endif

      if (lid_dn > 0 .and. lid_dn <= grid%nlmax) then
        ibeg_dn_loc = (lid_dn-1)*option%nflowdof + 1
        r_p(ibeg_dn_loc    ) = r_p(ibeg_dn_loc    ) - mconv(1)
        r_p(ibeg_dn_loc + 1) = r_p(ibeg_dn_loc + 1) - mconv(2)
        r_p(ibeg_dn_loc + 2) = r_p(ibeg_dn_loc + 2) - mconv(3)
      endif

    enddo

    cset => cset%next
  enddo

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
  type(stokes_auxvar_type), pointer :: stokes_auxvars(:)
  PetscReal, pointer :: r_p(:), xx_loc_p(:), accum_p(:), accum_p2(:)
  PetscInt :: lid, gid, ibeg_loc
  PetscReal :: vol, vol_over_dt

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  field  => realization%field
  stokes_auxvars => patch%aux%Stokes%auxvars

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
    ibeg_loc   = (lid-1)*option%nflowdof + 1
    vol  = patch%aux%Material%auxvars(gid)%volume
    vol_over_dt = vol / option%flow_dt

    ! Subtract old: r -= rho*vol*u^n (velocity only)
    r_p(ibeg_loc    ) = r_p(ibeg_loc    ) - accum_p(ibeg_loc    )
    r_p(ibeg_loc + 1) = r_p(ibeg_loc + 1) - accum_p(ibeg_loc + 1)
    r_p(ibeg_loc + 2) = r_p(ibeg_loc + 2) - accum_p(ibeg_loc + 2)
    ! DO NOT touch pressure residual (ibeg + 3)

    ! Add new: r += rho*vol*u^{n+1} (velocity only)
    r_p(ibeg_loc    ) = r_p(ibeg_loc    ) + rho * stokes_auxvars(gid)%u * vol_over_dt
    r_p(ibeg_loc + 1) = r_p(ibeg_loc + 1) + rho * stokes_auxvars(gid)%v * vol_over_dt
    r_p(ibeg_loc + 2) = r_p(ibeg_loc + 2) + rho * stokes_auxvars(gid)%w * vol_over_dt
    ! DO NOT touch pressure residual (ibeg + 3)

    accum_p2(ibeg_loc    ) = rho * stokes_auxvars(gid)%u * vol_over_dt
    accum_p2(ibeg_loc + 1) = rho * stokes_auxvars(gid)%v * vol_over_dt
    accum_p2(ibeg_loc + 2) = rho * stokes_auxvars(gid)%w * vol_over_dt
    accum_p2(ibeg_loc + 3) = 0.d0
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
  type(stokes_auxvar_type), pointer :: stokes_auxvars(:)

  PetscReal, pointer :: r_p(:)
  PetscInt :: lid, gid, istart, iend
  PetscReal :: vol


  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  field  => realization%field
  stokes_auxvars => patch%aux%Stokes%auxvars

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    iend = lid*option%nflowdof
    istart = iend-option%nflowdof+1
    vol  = patch%aux%Material%auxvars(gid)%volume

    r_p(istart    ) = r_p(istart    ) + rho * stokes_auxvars(gid)%du_dt(1) * vol
    r_p(istart + 1) = r_p(istart + 1) + rho * stokes_auxvars(gid)%du_dt(2) * vol
    r_p(istart + 2) = r_p(istart + 2) + rho * stokes_auxvars(gid)%du_dt(3) * vol
  enddo

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
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Stokes_Aux_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscReal :: rho
  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cset

  PetscReal, pointer :: r_p(:)
  PetscInt :: iconn, id_dn, local_id, ghosted_id
  PetscInt :: istart, iend
  PetscReal :: vol
  PetscReal :: force_u, force_v, force_w

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  if (.not.associated(patch%source_sink_list%first)) return

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

  source_sink => patch%source_sink_list%first
  do
    if (.not.associated(source_sink)) exit
    if (.not.associated(source_sink%flow_condition%rate)) then
      source_sink => source_sink%next
      cycle
    endif

    cset => source_sink%connection_set
    do iconn = 1, cset%num_connections
      id_dn = cset%id_dn(iconn)
      local_id = 0
      ghosted_id = 0

      if (id_dn > 0 .and. id_dn <= grid%nlmax) then
        local_id = id_dn
        ghosted_id = grid%nL2G(local_id)
      elseif (id_dn > 0 .and. id_dn <= grid%ngmax) then
        ghosted_id = id_dn
        local_id = grid%nG2L(ghosted_id)
      endif

      if (local_id <= 0 .or. local_id > grid%nlmax) cycle
      if (ghosted_id <= 0 .or. ghosted_id > grid%ngmax) cycle
      if (patch%imat(ghosted_id) <= 0) cycle

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      vol = patch%aux%Material%auxvars(ghosted_id)%volume

      force_u = source_sink%flow_aux_real_var(STOKES_U_DOF,iconn)
      force_v = source_sink%flow_aux_real_var(STOKES_V_DOF,iconn)
      force_w = source_sink%flow_aux_real_var(STOKES_W_DOF,iconn)

      r_p(istart    ) = r_p(istart    ) - rho*force_u*vol
      r_p(istart + 1) = r_p(istart + 1) - rho*force_v*vol
      r_p(istart + 2) = r_p(istart + 2) - rho*force_w*vol
    enddo

    source_sink => source_sink%next
  enddo

  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

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
  call StokesComputeViscousDiag(realization,mu)
  call StokesJacobianInternalConn(J,realization,mu,rho, &
                                            add_convection, &
                                            ierr)
  call StokesJacobianBoundaryFacesConn(J,realization,mu,ierr)
  call StokesJacobianAccumulation(J,realization,rho,ierr)
  call StokesJacobianSourceSink(J,realization,rho,ierr)

  ! Assembly happens once at the end
  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  if (associated(realization%patch%aux%Stokes)) then
    call MatrixZeroingZeroMatEntries(realization%patch%aux%Stokes%matrix_zeroing,J)
  endif

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
  use Matrix_Zeroing_module
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
  type(stokes_type),        pointer :: stokes_aux

  PetscInt :: iconn, gid_up, gid_dn, lid_up, lid_dn
  PetscInt :: rowstart_up, rowstart_dn, colstart_up, colstart_dn
  PetscReal :: dx, area
  PetscReal :: n(3), S(3)

  PetscReal :: J_KK(4,4), J_KL(4,4), J_LK(4,4), J_LL(4,4)
  PetscReal :: a_mu, gamma_pstab
  PetscReal, pointer :: aP(:)
  PetscReal :: a_f, rc

  ! convection Jacobian
  PetscReal :: uK(3), uL(3), u_face(3), phi
  PetscReal :: dm_dUK(3,3), dm_dUL(3,3)
  PetscReal :: I3(3,3)
  PetscInt :: i,j
  PetscInt :: gid, lid

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  gamma_pstab = 0.05  ! unused when Rhie-Chow is enabled

  I3 = 0.d0
  do i=1,3
    I3(i,i)=1.d0
  enddo

  if (.not.associated(patch%aux%Stokes)) return
  stokes_aux => patch%aux%Stokes

  ! Use pre-computed viscous diagonal aP for Rhie-Chow
  aP => stokes_aux%aP

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

      area = cset%area(iconn)
      dx   = cset%dist(0,iconn)
      n    = 0.d0
      n(1:3) = cset%dist(1:3,iconn)

      S(1)=n(1)*area
      S(2)=n(2)*area
      S(3)=n(3)*area

      a_mu = mu * area / dx

      ! Fetch uK/uL (for convection Jacobian)
      uK(1) = stokes_aux%auxvars(gid_up)%u
      uK(2) = stokes_aux%auxvars(gid_up)%v
      uK(3) = stokes_aux%auxvars(gid_up)%w

      uL(1) = stokes_aux%auxvars(gid_dn)%u
      uL(2) = stokes_aux%auxvars(gid_dn)%v
      uL(3) = stokes_aux%auxvars(gid_dn)%w
      if (uK(1) /= uK(1) .or. uK(2) /= uK(2) .or. uK(3) /= uK(3) .or. &
          uL(1) /= uL(1) .or. uL(2) /= uL(2) .or. uL(3) /= uL(3)) then
        write(option%io_buffer,'("StokesJacobianInternalConn: NaN in auxvars ",2i8)') gid_up, gid_dn
        call PrintErrMsg(option)
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
      ! Pressure contribution in momentum: +p_face * S, p_face=0.5(pK+pL)
      ! dR^u_K/dpK = +0.5*S, dR^u_K/dpL = +0.5*S
      ! dR^u_L/dpK = -0.5*S, dR^u_L/dpL = -0.5*S
      ! ---------------------------------------------------
      do i=1,3
        J_KK(i,4) = J_KK(i,4) + 0.5d0*S(i)
        J_KL(i,4) = J_KL(i,4) + 0.5d0*S(i)
        J_LK(i,4) = J_LK(i,4) - 0.5d0*S(i)
        J_LL(i,4) = J_LL(i,4) - 0.5d0*S(i)
      enddo
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
      ! Rhie-Chow pressure correction in continuity:
      ! flux = phi - rc*(pL - pK), rc = (A^2/a_f)/dx
      ! K row: +rc*pK - rc*pL ; L row: -rc*pK + rc*pL
      ! ---------------------------------------------------
      a_f = 0.5d0*(aP(gid_up) + aP(gid_dn))
      if (a_f <= 0.d0) a_f = a_mu
      rc = (area*area/a_f)/dx
      J_KK(4,4) = J_KK(4,4) + rc
      J_KL(4,4) = J_KL(4,4) - rc
      J_LK(4,4) = J_LK(4,4) - rc
      J_LL(4,4) = J_LL(4,4) + rc

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
        rowstart_up = (gid_up-1)*option%nflowdof
        colstart_up = (gid_up-1)*option%nflowdof
        colstart_dn = (gid_dn-1)*option%nflowdof

        do i=1,4
          do j=1,4
            call PUMSetValuesLocal(A,1,rowstart_up+i-1,1,colstart_up+j-1,J_KK(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
            call PUMSetValuesLocal(A,1,rowstart_up+i-1,1,colstart_dn+j-1,J_KL(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
          enddo
        enddo
      endif

      if (lid_dn > 0) then
        rowstart_dn = (gid_dn-1)*option%nflowdof
        colstart_up = (gid_up-1)*option%nflowdof
        colstart_dn = (gid_dn-1)*option%nflowdof

        do i=1,4
          do j=1,4
            call PUMSetValuesLocal(A,1,rowstart_dn+i-1,1,colstart_up+j-1,J_LK(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
            call PUMSetValuesLocal(A,1,rowstart_dn+i-1,1,colstart_dn+j-1,J_LL(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
          enddo
        enddo
      endif

    enddo

    cset => cset%next
  enddo

  end subroutine StokesJacobianInternalConn

! ************************************************************************** !
subroutine StokesJacobianBoundaryConn(A,realization,rho,mu,ierr)
  !
  ! Prescribed-condition Jacobian for Stokes problem
  ! Collects all prescribed constraints and applies once per cell
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
  PetscInt :: iconn, local_id, ghosted_id, gid, lid, ndof
  ! Track which rows have been constrained for this cell
  PetscBool, allocatable :: bc_applied(:,:)
  PetscInt, allocatable :: bc_rows(:)
  PetscInt :: bc_row_count
  PetscInt :: max_local_id

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option
  ndof = option%nflowdof

  if (.not.associated(patch%prescribed_condition_list%first)) then
    return
  endif

  ! Allocate tracking array for which DOFs have BCs
  max_local_id = grid%nlmax
  allocate(bc_applied(ndof, max_local_id))
  bc_applied = PETSC_FALSE

  ! First pass: identify all prescribed constraints
  bc => patch%prescribed_condition_list%first
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
      if (local_id <= 0 .or. local_id > grid%nlmax) cycle
      ghosted_id = grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. ghosted_id > grid%ngmax) cycle
      if (patch%imat(ghosted_id) <= 0) cycle

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
  PetscInt :: lid, gid, rowstart
  PetscReal :: Jacc(4,4), vol, vol_over_dt
  PetscInt :: i,j

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  if (option%flow%steady_state) return

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    rowstart = (gid - 1) * option%nflowdof
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
          call PUMSetValuesLocal(A,1,rowstart+i-1,1,rowstart+j-1,Jacc(i,j),ADD_VALUES,ierr);CHKERRQ(ierr)
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
  PetscInt :: lid, gid, rowstart, i, j
  PetscReal :: vol, val, val_p

  patch  => realization%patch
  grid   => patch%grid
  option => realization%option

  do lid = 1, grid%nlmax
    gid = grid%nL2G(lid)
    if (patch%imat(gid) <= 0) cycle
    rowstart = (gid - 1) * option%nflowdof
    vol    = patch%aux%Material%auxvars(gid)%volume
    val    = shift * rho * vol

    do i = 1,4
      do j = 1,4
        if (i /= j) cycle
        if (i == 4) cycle  ! No accumulation for pressure
          call PUMSetValuesLocal(A,1,rowstart+i-1,1,rowstart+j-1,val,ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
    enddo
  enddo

end subroutine StokesIJacobianAccumulation

! ************************************************************************** !

subroutine StokesJacobianSourceSink(A,realization,rho,ierr)
  !
  ! Source/sink Jacobian for Stokes problem.
  !
  ! Currently a no-op: only constant body forces (gravity) are supported
  ! as source terms, so the Jacobian contribution is zero.
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
  PetscReal, pointer :: aP(:)

  PetscInt :: iconn, local_id, ghosted_id, gid, istart, istart_ghost, ndof
  PetscReal :: Ucell(4), Res(4)

  patch => realization%patch
  grid  => patch%grid
  option => realization%option
  ndof = option%nflowdof

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(realization%field%flow_xx_loc,u_loc,ierr);CHKERRQ(ierr)

  ! Use pre-computed viscous diagonal aP for Rhie-Chow
  aP => patch%aux%Stokes%aP

  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    if (.not.associated(boundary_condition%flow_condition)) then
      boundary_condition => boundary_condition%next
      cycle
    endif
    cset => boundary_condition%connection_set
    if (.not.associated(cset)) then
      boundary_condition => boundary_condition%next
      cycle
    endif
    if (.not.associated(boundary_condition%flow_condition%stokes%velocity) .and. &
        .not.associated(boundary_condition%flow_condition%stokes%pressure)) then
      boundary_condition => boundary_condition%next
      cycle
    endif

    do iconn = 1, cset%num_connections

      local_id = cset%id_dn(iconn)
      if (local_id <= 0 .or. local_id > grid%nlmax) cycle
      ghosted_id = grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. ghosted_id > grid%ngmax) cycle
      if (patch%imat(ghosted_id) <= 0) cycle

      istart = (local_id-1)*ndof + 1
      istart_ghost = (ghosted_id-1)*ndof + 1

      Ucell(1) = u_loc(istart_ghost    )
      Ucell(2) = u_loc(istart_ghost + 1)
      Ucell(3) = u_loc(istart_ghost + 2)
      Ucell(4) = u_loc(istart_ghost + 3)

      call StokesBCFluxConn(boundary_condition,iconn, &
                            cset%area(iconn),cset%dist(:,iconn), &
                            rho,mu,option,Ucell,aP(ghosted_id),Res)

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

subroutine StokesBCFluxConn(boundary_condition,iconn,area,dist,rho,mu,option,Ucell,aP,Res)
  !
  ! Compute boundary-face residual contribution for ONE boundary connection (ONE face).
  !
  ! Flux-based velocity Dirichlet:
  !   momentum:  +a_b (uK - u_b),  a_b = mu*A/d_b
  !   continuity: + (u_b · S)
  !
  ! Pressure traction (Option B):
  !   momentum:  + p_b * S
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
  type(option_type), intent(inout) :: option
  PetscReal, intent(in) :: Ucell(4)
  PetscReal, intent(in) :: aP
  PetscReal, intent(out) :: Res(4)

  PetscReal :: n(3), S(3), d_b, a_b
  PetscReal :: uK(3), pK
  PetscReal :: u_b(3), p_b
  PetscReal :: phi_b
  PetscReal :: rc
  PetscBool :: has_u_dir, has_p_bc

  Res = 0.d0
  u_b = 0.d0
  p_b = 0.d0

  n(1)=dist(1); n(2)=dist(2); n(3)=dist(3)
  S(1)=n(1)*area; S(2)=n(2)*area; S(3)=n(3)*area

  d_b = dist(0)

  if (d_b <= 0.d0) then
    write(option%io_buffer,'("StokesBCFluxConn: bad boundary distance: dist(-1:0) =",2es12.4)') dist(-1),dist(0)
    call PrintErrMsg(option)
  endif

  uK(1)=Ucell(1); uK(2)=Ucell(2); uK(3)=Ucell(3)
  pK = Ucell(4)

  has_u_dir  = PETSC_FALSE
  has_p_bc   = PETSC_FALSE

  ! Velocity Dirichlet present?
  if (associated(boundary_condition%flow_condition%stokes%velocity)) then
    if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
      has_u_dir = PETSC_TRUE
      u_b(1) = boundary_condition%flow_aux_real_var(1,iconn)
      u_b(2) = boundary_condition%flow_aux_real_var(2,iconn)
      u_b(3) = boundary_condition%flow_aux_real_var(3,iconn)
    endif
  endif

  ! Pressure BC present? (treat LIQUID_PRESSURE DIRICHLET as traction value p_b)
  if (associated(boundary_condition%flow_condition%stokes%pressure)) then
    if (boundary_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
      p_b = boundary_condition%flow_aux_real_var(4,iconn)
      has_p_bc = PETSC_TRUE
    endif
  endif

  ! If no pressure BC, use interior pressure (no correction)
  if (.not.has_p_bc) p_b = pK

  ! ---- Flux-based velocity Dirichlet via viscous diffusion ----
  ! Boundary continuity uses the same Rhie-Chow correction as interior faces:
  !   phi_f = u_f·S - (A^2/aP) * (p_b - pK)/d_b
  ! If no pressure BC, set p_b = pK so the correction is zero.
  rc = 0.d0
  if (aP > 0.d0) then
    rc = (area*area/aP) * (p_b - pK) / d_b
  endif

  if (has_u_dir) then
    a_b = mu * area / d_b
    Res(1) = Res(1) + a_b*(uK(1) - u_b(1))
    Res(2) = Res(2) + a_b*(uK(2) - u_b(2))
    Res(3) = Res(3) + a_b*(uK(3) - u_b(3))

    ! continuity boundary flux uses prescribed u_b
    phi_b = u_b(1)*S(1) + u_b(2)*S(2) + u_b(3)*S(3)
    Res(4) = Res(4) + phi_b - rc
  else
    ! No velocity Dirichlet: use "do-nothing/outflow" flux for continuity
    phi_b = uK(1)*S(1) + uK(2)*S(2) + uK(3)*S(3)
    Res(4) = Res(4) + phi_b - rc
  endif

  ! ---- Pressure traction ----
  ! Sign: interior cell is id_dn. The outward normal n points from ghost (id_up)
  ! toward interior (id_dn). In the interior loop, id_dn gets -Tpress = -p_face*S.
  ! The same convention applies here: boundary interior cell gets -p_b*S.
  Res(1) = Res(1) - p_b*S(1)
  Res(2) = Res(2) - p_b*S(2)
  Res(3) = Res(3) - p_b*S(3)

end subroutine StokesBCFluxConn

! ************************************************************************** !

subroutine StokesResidualBoundaryFinal(r,xx,realization,ierr)
  !
  ! Sets residuals for prescribed rows to (u - u_prescribed)
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
  use Grid_module
  use Grid_Structured_module

  implicit none

  Vec :: r, xx
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(coupler_type), pointer :: prescribed_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id, iconn, offset, ghosted_id
  PetscReal, pointer :: r_ptr(:), xx_ptr(:)
  PetscInt :: nflowdof, istart, iend
  PetscReal :: p_bc, face_area, face_nx, face_ny, face_nz
  PetscReal :: cell_dx, cell_dy, cell_dz

  patch => realization%patch
  option => realization%option
  nflowdof = option%nflowdof

  if (.not.associated(patch%prescribed_condition_list%first)) return

  call VecGetArray(r,r_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(xx,xx_ptr,ierr);CHKERRQ(ierr)

  prescribed_condition => patch%prescribed_condition_list%first
  do
    if (.not.associated(prescribed_condition)) exit
    cur_connection_set => prescribed_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      local_id = cur_connection_set%id_dn(iconn)
      if (local_id <= 0 .or. local_id > patch%grid%nlmax) cycle
      ghosted_id = patch%grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. ghosted_id > patch%grid%ngmax) cycle
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*nflowdof
      istart = iend-nflowdof+1

      ! Set prescribed velocity residual to (u - u_prescribed)
      if (associated(prescribed_condition%flow_condition%stokes%velocity)) then
        if (prescribed_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          r_ptr(istart) = xx_ptr(istart) - prescribed_condition%flow_aux_real_var(1,iconn)
          r_ptr(istart+1) = xx_ptr(istart+1) - prescribed_condition%flow_aux_real_var(2,iconn)
          r_ptr(istart+2) = xx_ptr(istart+2) - prescribed_condition%flow_aux_real_var(3,iconn)
        endif
      endif

      ! Set prescribed pressure residual to (p - p_prescribed)
      ! Also add the missing boundary-face pressure traction -p_bc*S to the
      ! momentum residuals.  Interior cell is id_dn; the inward face normal S
      ! points from the ghost (id_up) into the domain, so the traction on the
      ! interior cell is -p_bc*S (same sign as id_dn in the interior loop).
      !
      ! PRESCRIBED_CONDITION connection sets use GENERIC_CONNECTION_TYPE, so
      ! dist/area are not populated by GridPopulateConnection.  Compute the
      ! face normal and area directly from the grid dimensions and coupler iface.
      !
      ! NOTE: If a velocity PRESCRIBED_CONDITION is also applied to this cell
      ! (corner cell), the velocity block above will have already set
      ! r(u) = u - u_bc, overwriting these momentum values.  Velocity couplers
      ! must therefore be listed AFTER pressure couplers in the input file.
      if (associated(prescribed_condition%flow_condition%stokes%pressure)) then
        if (prescribed_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          r_ptr(istart+3) = xx_ptr(istart+3) - prescribed_condition%flow_aux_real_var(4,iconn)
          p_bc = prescribed_condition%flow_aux_real_var(4,iconn)
          ! Compute face area and inward normal from structured grid dimensions
          cell_dx = patch%grid%structured_grid%dx(ghosted_id)
          cell_dy = patch%grid%structured_grid%dy(ghosted_id)
          cell_dz = patch%grid%structured_grid%dz(ghosted_id)
          face_nx = 0.d0; face_ny = 0.d0; face_nz = 0.d0
          select case(prescribed_condition%region%iface)
            case(WEST_FACE)
              face_nx  = 1.d0   ! inward normal: +x
              face_area = cell_dy * cell_dz
            case(EAST_FACE)
              face_nx  = -1.d0  ! inward normal: -x
              face_area = cell_dy * cell_dz
            case(SOUTH_FACE)
              face_ny  = 1.d0   ! inward normal: +y
              face_area = cell_dx * cell_dz
            case(NORTH_FACE)
              face_ny  = -1.d0  ! inward normal: -y
              face_area = cell_dx * cell_dz
            case(BOTTOM_FACE)
              face_nz  = 1.d0   ! inward normal: +z
              face_area = cell_dx * cell_dy
            case(TOP_FACE)
              face_nz  = -1.d0  ! inward normal: -z
              face_area = cell_dx * cell_dy
            case default
              face_area = 0.d0
          end select
          r_ptr(istart  ) = r_ptr(istart  ) - p_bc * face_nx * face_area
          r_ptr(istart+1) = r_ptr(istart+1) - p_bc * face_ny * face_area
          r_ptr(istart+2) = r_ptr(istart+2) - p_bc * face_nz * face_area
        endif
      endif

    enddo
    prescribed_condition => prescribed_condition%next
  enddo

  call VecRestoreArrayRead(xx,xx_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_ptr,ierr);CHKERRQ(ierr)

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

  PetscInt :: iconn, local_id, ghosted_id, gid, rowstart, ndof
  PetscReal :: area, dist(-1:3), d_b, a_b
  PetscReal, pointer :: aP(:)
  PetscReal :: rc
  PetscBool :: has_u_dir, has_p_bc
  PetscReal :: S(3)

  patch => realization%patch
  grid  => patch%grid
  option => realization%option
  ndof = option%nflowdof

  ! Use pre-computed viscous diagonal aP for Rhie-Chow
  aP => patch%aux%Stokes%aP

  boundary_condition => patch%boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    cset => boundary_condition%connection_set
    if (.not.associated(cset)) then
      boundary_condition => boundary_condition%next
      cycle
    endif
    if (.not.associated(boundary_condition%flow_condition)) then
      boundary_condition => boundary_condition%next
      cycle
    endif
    if (.not.associated(boundary_condition%flow_condition%stokes%velocity) .and. &
        .not.associated(boundary_condition%flow_condition%stokes%pressure)) then
      boundary_condition => boundary_condition%next
      cycle
    endif

    do iconn = 1, cset%num_connections

      local_id = cset%id_dn(iconn)
      if (local_id <= 0 .or. local_id > grid%nlmax) cycle
      ghosted_id = grid%nL2G(local_id)
      if (ghosted_id <= 0 .or. ghosted_id > grid%ngmax) cycle
      if (patch%imat(ghosted_id) <= 0) cycle

      rowstart = (ghosted_id-1)*ndof

      area = cset%area(iconn)
      dist = cset%dist(:,iconn)

      d_b = dist(0)

      ! S = n*A (only needed if you use do-nothing flux Jacobian)
      S(1) = dist(1)*area
      S(2) = dist(2)*area
      S(3) = dist(3)*area

      has_u_dir = PETSC_FALSE
      has_p_bc = PETSC_FALSE
      if (associated(boundary_condition%flow_condition%stokes%velocity)) then
        if (boundary_condition%flow_condition%stokes%velocity%itype == DIRICHLET_BC) then
          has_u_dir = PETSC_TRUE
        endif
      endif

      if (associated(boundary_condition%flow_condition%stokes%pressure)) then
        if (boundary_condition%flow_condition%stokes%pressure%itype == DIRICHLET_BC) then
          has_p_bc = PETSC_TRUE
        endif
      endif

      ! If velocity Dirichlet: momentum has +a_b(uK-u_b) => diag += a_b
      if (has_u_dir) then
        a_b = mu * area / d_b
        call PUMSetValuesLocal(A,1,rowstart+0,1,rowstart+0,a_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,rowstart+1,1,rowstart+1,a_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,rowstart+2,1,rowstart+2,a_b,ADD_VALUES,ierr);CHKERRQ(ierr)
        ! continuity uses prescribed flux u_b·S -> no Jacobian
      else
        ! do-nothing continuity flux: phi = uK·S -> J(p,u) += S^T
        call PUMSetValuesLocal(A,1,rowstart+3,1,rowstart+0,S(1),ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,rowstart+3,1,rowstart+1,S(2),ADD_VALUES,ierr);CHKERRQ(ierr)
        call PUMSetValuesLocal(A,1,rowstart+3,1,rowstart+2,S(3),ADD_VALUES,ierr);CHKERRQ(ierr)
      endif

      if (aP(ghosted_id) > 0.d0) then
        rc = (area*area/aP(ghosted_id)) / d_b
        call PUMSetValuesLocal(A,1,rowstart+3,1,rowstart+3,rc,ADD_VALUES,ierr);CHKERRQ(ierr)
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

  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(xx,xx_p,ierr);CHKERRQ(ierr)

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

  call VecRestoreArrayRead(xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

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
