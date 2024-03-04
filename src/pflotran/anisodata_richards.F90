module Anisotropy_Richards_Data_module

  ! Module for managing additional data needed in Richards' flux calculations to
  ! deal with full (off-diagonal) anisotropy terms
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  
  implicit none
  
  public :: PermInnerProd
  public :: PrintRichardsAnisoData
  
  ! -- 
  ! -- Data structures for calculating anisotropic gradient contributions in STRUCTURED GRIDS.
  ! --
  
  type, public :: aniso_richards_data_structured_type
    ! Interface anisotropy data for Richards eqn (structured grid case)
  
    ! Adjacent pressure data - working data, overwritten throughout simulation.  Only 8 are used.
    PetscReal :: pres_up_adjXm ! pressure in the upwind   cell adjacent in the X -ve direction
    PetscReal :: pres_up_adjXp ! pressure in the upwind   cell adjacent in the X +ve direction
    PetscReal :: pres_dn_adjXm ! pressure in the downwind cell adjacent in the X -ve direction
    PetscReal :: pres_dn_adjXp ! pressure in the downwind cell adjacent in the X +ve direction
    PetscReal :: pres_up_adjYm ! pressure in the upwind   cell adjacent in the Y -ve direction
    PetscReal :: pres_up_adjYp ! pressure in the upwind   cell adjacent in the Y +ve direction
    PetscReal :: pres_dn_adjYm ! pressure in the downwind cell adjacent in the Y -ve direction
    PetscReal :: pres_dn_adjYp ! pressure in the downwind cell adjacent in the Y +ve direction
    PetscReal :: pres_up_adjZm ! pressure in the upwind   cell adjacent in the Z -ve direction
    PetscReal :: pres_up_adjZp ! pressure in the upwind   cell adjacent in the Z +ve direction
    PetscReal :: pres_dn_adjZm ! pressure in the downwind cell adjacent in the Z -ve direction
    PetscReal :: pres_dn_adjZp ! pressure in the downwind cell adjacent in the Z +ve direction
  
    ! Jacobian terms from v_darcy anisotropic contributions - working data, overwritten throughout simulation
    PetscReal :: dv_darcy_aniso_dp_up_adjXm ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell adjacent in the X -ve direction)
    PetscReal :: dv_darcy_aniso_dp_dn_adjXm ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell adjacent in the X -ve direction)
    PetscReal :: dv_darcy_aniso_dp_up_adjXp ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell adjacent in the X +ve direction)
    PetscReal :: dv_darcy_aniso_dp_dn_adjXp ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell adjacent in the X +ve direction)
  
    PetscReal :: dv_darcy_aniso_dp_up_adjYm ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell adjacent in the Y -ve direction)
    PetscReal :: dv_darcy_aniso_dp_dn_adjYm ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell adjacent in the Y -ve direction)
    PetscReal :: dv_darcy_aniso_dp_up_adjYp ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell adjacent in the Y +ve direction)
    PetscReal :: dv_darcy_aniso_dp_dn_adjYp ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell adjacent in the Y +ve direction)
  
    PetscReal :: dv_darcy_aniso_dp_up_adjZm ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell adjacent in the Z -ve direction)
    PetscReal :: dv_darcy_aniso_dp_dn_adjZm ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell adjacent in the Z -ve direction)
    PetscReal :: dv_darcy_aniso_dp_up_adjZp ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell adjacent in the Z +ve direction)
    PetscReal :: dv_darcy_aniso_dp_dn_adjZp ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell adjacent in the Z +ve direction)
  
    PetscReal :: dv_darcy_aniso_dp_up       ! d(darcy aniso contribution from downwind cell)/d(pressure in the upwind   cell)
    PetscReal :: dv_darcy_aniso_dp_dn       ! d(darcy aniso contribution from downwind cell)/d(pressure in the downwind cell)
  
    PetscReal :: Jup_adjXm ! J contribution from upwind   cell adjacent in the X -ve direction
    PetscReal :: Jdn_adjXm ! J contribution from downwind cell adjacent in the X -ve direction
    PetscReal :: Jup_adjXp ! J contribution from upwind   cell adjacent in the X +ve direction
    PetscReal :: Jdn_adjXp ! J contribution from downwind cell adjacent in the X +ve direction
    PetscReal :: Jup_adjYm ! J contribution from upwind   cell adjacent in the Y -ve direction
    PetscReal :: Jdn_adjYm ! J contribution from downwind cell adjacent in the Y -ve direction
    PetscReal :: Jup_adjYp ! J contribution from upwind   cell adjacent in the Y +ve direction
    PetscReal :: Jdn_adjYp ! J contribution from downwind cell adjacent in the Y +ve direction
    PetscReal :: Jup_adjZm ! J contribution from upwind   cell adjacent in the Z -ve direction
    PetscReal :: Jdn_adjZm ! J contribution from downwind cell adjacent in the Z -ve direction
    PetscReal :: Jup_adjZp ! J contribution from upwind   cell adjacent in the Z +ve direction
    PetscReal :: Jdn_adjZp ! J contribution from downwind cell adjacent in the Z +ve direction
  end type aniso_richards_data_structured_type
  
  type, public :: aniso_richards_work_structured_type
    ! structure for storing work data for anisotropic derivative calculations - structured case
    PetscReal :: ddphi_x_dp_up_adjXm, ddphi_x_dp_dn_adjXm, ddphi_x_dp_up_adjXp, ddphi_x_dp_dn_adjXp
    PetscReal :: ddphi_y_dp_up_adjYm, ddphi_y_dp_dn_adjYm, ddphi_y_dp_up_adjYp, ddphi_y_dp_dn_adjYp
    PetscReal :: ddphi_z_dp_up_adjZm, ddphi_z_dp_dn_adjZm, ddphi_z_dp_up_adjZp, ddphi_z_dp_dn_adjZp
    PetscReal :: ddphi_x_dp_up, ddphi_x_dp_dn, ddphi_y_dp_up, ddphi_y_dp_dn, ddphi_z_dp_up, ddphi_z_dp_dn
  
    PetscReal :: dv_darcy_aniso_up_dp_up_adjXm
    PetscReal :: dv_darcy_aniso_up_dp_up_adjXp
    PetscReal :: dv_darcy_aniso_up_dp_up_adjYm
    PetscReal :: dv_darcy_aniso_up_dp_up_adjYp
    PetscReal :: dv_darcy_aniso_up_dp_up_adjZm
    PetscReal :: dv_darcy_aniso_up_dp_up_adjZp
    PetscReal :: dv_darcy_aniso_up_dp_up
  
    PetscReal :: dv_darcy_aniso_dn_dp_dn_adjXm
    PetscReal :: dv_darcy_aniso_dn_dp_dn_adjXp
    PetscReal :: dv_darcy_aniso_dn_dp_dn_adjYm
    PetscReal :: dv_darcy_aniso_dn_dp_dn_adjYp
    PetscReal :: dv_darcy_aniso_dn_dp_dn_adjZm
    PetscReal :: dv_darcy_aniso_dn_dp_dn_adjZp
    PetscReal :: dv_darcy_aniso_dn_dp_dn
  end type aniso_richards_work_structured_type
  
  ! -- 
  ! -- Data structures for calculating anisotropic gradient contributions in UNSTRUCTURED GRIDS.
  ! --
  
  type, public :: aniso_richards_unstructured_cell_type
    ! Permeability w.r.t. (u,v,w) coords in cells - working data, overwritten throughout simulation (in case timedep)
    PetscReal :: perm_uu
    PetscReal :: perm_vv
    PetscReal :: perm_ww
    PetscReal :: perm_uv
    PetscReal :: perm_uw
    PetscReal :: perm_vw
  
    ! Adjacent pressure data - working data, overwritten throughout simulation
    PetscReal, pointer :: pres_adj(:) ! pressure adjacent to upwind cell
  end type aniso_richards_unstructured_cell_type
  
  type, public :: aniso_richards_data_unstructured_type
    type(aniso_richards_unstructured_cell_type), pointer :: up => null()
    type(aniso_richards_unstructured_cell_type), pointer :: dn => null()
  
    ! derivative data
    PetscReal :: dv_darcy_aniso_dp_up
    PetscReal :: dv_darcy_aniso_dp_dn
    PetscReal, pointer :: dv_darcy_aniso_dp_up_adj(:)
    PetscReal, pointer :: dv_darcy_aniso_dp_dn_adj(:)
    PetscReal, pointer :: Jup_adj(:)
    PetscReal, pointer :: Jdn_adj(:)
  end type aniso_richards_data_unstructured_type
  
  type, public :: aniso_richards_work_unstructured_type
    ! structure for storing work data for anisotropic derivative calculations - unstructured case
    PetscReal :: ddphi_v_up_dp_up
    PetscReal :: ddphi_v_up_dp_dn
    PetscReal, pointer :: ddphi_v_up_dp_up_adj(:)
    PetscReal :: ddphi_w_up_dp_up
    PetscReal :: ddphi_w_up_dp_dn
    PetscReal, pointer :: ddphi_w_up_dp_up_adj(:)
    PetscReal :: ddphi_v_dn_dp_dn
    PetscReal :: ddphi_v_dn_dp_up
    PetscReal, pointer :: ddphi_v_dn_dp_dn_adj(:)
    PetscReal :: ddphi_w_dn_dp_dn
    PetscReal :: ddphi_w_dn_dp_up
    PetscReal, pointer :: ddphi_w_dn_dp_dn_adj(:)
  end type aniso_richards_work_unstructured_type

  ! --
  ! -- General data type containing structured and unstructured ptrs
  ! --

  type, public :: aniso_richards_data_type
    type( aniso_richards_data_structured_type), pointer :: struct => null()
    type( aniso_richards_data_unstructured_type), pointer :: unstruct => null()
  end type aniso_richards_data_type
  

  ! --
  ! -- Methods
  ! --
  
  
  CONTAINS
    function PermInnerProd( u, v, K_xx, K_yy, K_zz, K_xy, K_xz, K_yz )
      ! compute u^T K v, with K symmetric
  
      PetscReal :: u(3)
      PetscReal :: v(3)
      PetscReal :: K_xx,K_yy,K_zz,K_xy,K_xz,K_yz
      PetscReal :: PermInnerProd
  
      PermInnerProd =  u(1) * (K_xx*v(1) + K_xy*v(2) + K_xz*v(3)) &
                     + u(2) * (K_xy*v(1) + K_yy*v(2) + K_yz*v(3)) &
                     + u(3) * (K_xz*v(1) + K_yz*v(2) + K_zz*v(3))
    end function PermInnerProd
  
    ! for debugging:
    subroutine PrintRichardsAnisoData( aniso_rich, aniso_geom )
  
      use Anisotropy_Geom_Data_module

      type( aniso_richards_data_type ), pointer :: aniso_rich
      type( aniso_geom_data_type ), pointer :: aniso_geom
      integer ierr,rank,size,i
  
      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
  

      if(associated(aniso_rich%unstruct)) then
        ! transformed perm for face
        write(*, "('  perm up,dn')")
          write(*, "('    perm_uu = ',e10.3,',',e10.3)") &
            aniso_rich%unstruct%up%perm_uu, aniso_rich%unstruct%dn%perm_uu
          write(*, "('    perm_vv = ',e10.3,',',e10.3)") &
            aniso_rich%unstruct%up%perm_vv, aniso_rich%unstruct%dn%perm_vv
          write(*, "('    perm_ww = ',e10.3,',',e10.3)") &
            aniso_rich%unstruct%up%perm_ww, aniso_rich%unstruct%dn%perm_ww
          write(*, "('    perm_uv = ',e10.3,',',e10.3)") &
            aniso_rich%unstruct%up%perm_uv, aniso_rich%unstruct%dn%perm_uv
          write(*, "('    perm_uw = ',e10.3,',',e10.3)") &
            aniso_rich%unstruct%up%perm_uw, aniso_rich%unstruct%dn%perm_uw
          write(*, "('    perm_vw = ',e10.3,',',e10.3)") &
            aniso_rich%unstruct%up%perm_vw, aniso_rich%unstruct%dn%perm_vw
          ! pressures
          write(*, "('  Pressures - upstream')")
          do i=1,aniso_geom%unstruct%up%num_adj
            if (aniso_geom%unstruct%up%cell_id_adj(i) == 0) cycle
            write(*, "('    ',i2,'. pres_adj=',e10.3)") i,aniso_rich%unstruct%up%pres_adj(i)
          enddo
          write(*, "('  Pressures - downstream')")
          do i=1,aniso_geom%unstruct%dn%num_adj
            if (aniso_geom%unstruct%dn%cell_id_adj(i) == 0) cycle
            write(*, "('    ',i2,'. pres_adj=',e10.3)") i,aniso_rich%unstruct%dn%pres_adj(i)
          enddo
  
      endif
    end subroutine

end module Anisotropy_Richards_Data_module
  