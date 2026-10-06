module TS_Flow_State_module
  !
  ! Flow quantities a transport model reads, with their pressure
  ! derivatives, when flow and transport are one equation model
  ! (Flow_RT_DAE_Model_class). Connections follow the patch's connection
  ! lists, as in the models' connection tables.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  type, public :: ts_flow_state_type
    ! changes whenever the values change
    PetscInt :: id = 0
    ! ghosted cells: liquid saturation, density (in the units of
    ! global_auxvar_type%den), density [kg/m^3], pressure [Pa] and porosity
    PetscReal, allocatable :: sat(:)
    PetscReal, allocatable :: den(:)
    PetscReal, allocatable :: den_kg(:)
    PetscReal, allocatable :: pres(:)
    PetscReal, allocatable :: por(:)
    ! their derivatives with respect to the cell's pressure
    PetscReal, allocatable :: dsat_dp(:)
    PetscReal, allocatable :: dden_kg_dp(:)
    PetscReal, allocatable :: dpor_dp(:)
    ! internal connections: Darcy flux from up to down [m/s], and its
    ! derivatives with respect to the up and down pressures, (2,nconn_int)
    PetscReal, allocatable :: v_int(:)
    PetscReal, allocatable :: dv_int(:,:)
    ! boundary connections: Darcy flux into the cell [m/s], the boundary
    ! state, and their derivatives with respect to the cell's pressure
    ! (the state's are zero for a prescribed-pressure boundary)
    PetscReal, allocatable :: v_bc(:)
    PetscReal, allocatable :: dv_bc(:)
    PetscReal, allocatable :: sat_bc(:)
    PetscReal, allocatable :: den_bc(:)
    PetscReal, allocatable :: den_kg_bc(:)
    PetscReal, allocatable :: pres_bc(:)
    PetscReal, allocatable :: dsat_bc_dp(:)
    PetscReal, allocatable :: dden_kg_bc_dp(:)
    ! sources: volume rate into the cell [m^3/s], and its derivative with
    ! respect to the cell's pressure
    PetscReal, allocatable :: qvol_ss(:)
    PetscReal, allocatable :: dqvol_ss(:)
  end type ts_flow_state_type

  public :: TSFlowStateInit

contains

! ************************************************************************** !

subroutine TSFlowStateInit(fs,ngmax,nconn_int,nconn_bc,nconn_ss)
  !
  ! Allocates the arrays of a grid with ngmax ghosted cells and the given
  ! connection counts, zeroed.
  !
  ! Author: Hansol Suh
  ! Date: 09/23/26
  !
  implicit none

  type(ts_flow_state_type) :: fs
  PetscInt :: ngmax, nconn_int, nconn_bc, nconn_ss

  PetscInt :: nint, nbc, nss

  nint = max(nconn_int,1)
  nbc = max(nconn_bc,1)
  nss = max(nconn_ss,1)
  allocate(fs%sat(ngmax),fs%den(ngmax),fs%den_kg(ngmax),fs%pres(ngmax), &
           fs%por(ngmax))
  allocate(fs%dsat_dp(ngmax),fs%dden_kg_dp(ngmax),fs%dpor_dp(ngmax))
  allocate(fs%v_int(nint),fs%dv_int(2,nint))
  allocate(fs%v_bc(nbc),fs%dv_bc(nbc),fs%sat_bc(nbc),fs%den_bc(nbc), &
           fs%den_kg_bc(nbc),fs%pres_bc(nbc),fs%dsat_bc_dp(nbc), &
           fs%dden_kg_bc_dp(nbc))
  allocate(fs%qvol_ss(nss),fs%dqvol_ss(nss))
  fs%sat = 0.d0
  fs%den = 0.d0
  fs%den_kg = 0.d0
  fs%pres = 0.d0
  fs%por = 0.d0
  fs%dsat_dp = 0.d0
  fs%dden_kg_dp = 0.d0
  fs%dpor_dp = 0.d0
  fs%v_int = 0.d0
  fs%dv_int = 0.d0
  fs%v_bc = 0.d0
  fs%dv_bc = 0.d0
  fs%sat_bc = 0.d0
  fs%den_bc = 0.d0
  fs%den_kg_bc = 0.d0
  fs%pres_bc = 0.d0
  fs%dsat_bc_dp = 0.d0
  fs%dden_kg_bc_dp = 0.d0
  fs%qvol_ss = 0.d0
  fs%dqvol_ss = 0.d0

end subroutine TSFlowStateInit

end module TS_Flow_State_module
