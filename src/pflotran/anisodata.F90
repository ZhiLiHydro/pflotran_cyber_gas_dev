module Anisotropy_Data_module

! Module for managing additional data needed in flux calculations to
! deal with full (off-diagonal) anisotropy terms

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  public :: PermInnerProd
  public :: PrintAnisoData

  type, public :: aniso_data_structured_type

    !  For calculating flux across interface marked '---->'.
    !  Including data from cells 'up' and 'dn' only allows flux based on permeability components
    !  in the direction of the flux to be calculated.  Anisotropic permeablity components will
    !  be ignored.
    !  To include anisotropic components, the pressure gradient in the plane of interface must
    !  also be computed.  In a 2-D structured grid, this can be done by approximating the pressure
    !  gradient using the up/downstream adjacent cells marked A,B.
    !
    !  .......................
    !  :          :          :
    !  :   up_B   :   dn_B   :   <- y-adjacent up/dn cells (+ve dir)
    !  :          :          :
    !  -----------------------
    !  |          |          |
    !  |    up  --|--> dn    |
    !  |          |          |
    !  -----------------------
    !  :          :          :
    !  :   up_A   :   dn_A   :   <- y-adjacent up/dn cells (-ve dir)
    !  :          :          :
    !  .......................
    !
    !  Note. X_A, X_B should be ordered so that the coordinate X(=X,Y,orZ)is increasing the
    !  direction A->B, in order to make sense of direction of gradients.  In aniso_data the
    !  p(lus) and m(inus) suffix indicates this.
    !
    !  Note. residual equations do not need additional terms local_id_up_A, local_id_dn_A, ... etc.
    !  since flows are only between cells up/dn so these are the only source/sinks of mass.

    ! Orientation (of normal to interface)
    character (len=1) :: orientation ! 'X', 'Y' or 'Z'

    ! Grid data - fixed through simulation
    PetscInt :: id_up_adjXm ! id of upwind   cell adjacent in the X -ve direction (structured grids only)
    PetscInt :: id_up_adjXp ! id of upwind   cell adjacent in the X +ve direction (structured grids only)
    PetscInt :: id_dn_adjXm ! id of downwind cell adjacent in the X -ve direction (structured grids only)
    PetscInt :: id_dn_adjXp ! id of downwind cell adjacent in the X +ve direction (structured grids only)
    PetscInt :: id_up_adjYm ! id of upwind   cell adjacent in the Y -ve direction (structured grids only)
    PetscInt :: id_up_adjYp ! id of upwind   cell adjacent in the Y +ve direction (structured grids only)
    PetscInt :: id_dn_adjYm ! id of downwind cell adjacent in the Y -ve direction (structured grids only)
    PetscInt :: id_dn_adjYp ! id of downwind cell adjacent in the Y +ve direction (structured grids only)
    PetscInt :: id_up_adjZm ! id of upwind   cell adjacent in the Z -ve direction (structured grids only)
    PetscInt :: id_up_adjZp ! id of upwind   cell adjacent in the Z +ve direction (structured grids only)
    PetscInt :: id_dn_adjZm ! id of downwind cell adjacent in the Z -ve direction (structured grids only)
    PetscInt :: id_dn_adjZp ! id of downwind cell adjacent in the Z +ve direction (structured grids only)

    ! Gravity terms (usually = g.z where z is z coord of relevant cell centre
    PetscReal :: loc_gravity_up
    PetscReal :: loc_gravity_dn
    PetscReal :: loc_gravity_up_adjXm
    PetscReal :: loc_gravity_up_adjXp
    PetscReal :: loc_gravity_dn_adjXm
    PetscReal :: loc_gravity_dn_adjXp
    PetscReal :: loc_gravity_up_adjYm
    PetscReal :: loc_gravity_up_adjYp
    PetscReal :: loc_gravity_dn_adjYm
    PetscReal :: loc_gravity_dn_adjYp
    PetscReal :: loc_gravity_up_adjZm
    PetscReal :: loc_gravity_up_adjZp
    PetscReal :: loc_gravity_dn_adjZm
    PetscReal :: loc_gravity_dn_adjZp

    ! Distances from up/dn to adjacent cells (distances between cell centres)
    PetscReal :: dist_up_adjXm
    PetscReal :: dist_up_adjXp
    PetscReal :: dist_dn_adjXm
    PetscReal :: dist_dn_adjXp
    PetscReal :: dist_up_adjYm
    PetscReal :: dist_up_adjYp
    PetscReal :: dist_dn_adjYm
    PetscReal :: dist_dn_adjYp
    PetscReal :: dist_up_adjZm
    PetscReal :: dist_up_adjZp
    PetscReal :: dist_dn_adjZm
    PetscReal :: dist_dn_adjZp

    ! Adjacent pressure data - working data, overwritten throughout simulation
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

  end type aniso_data_structured_type


  type, public :: aniso_deriv_data_structured_type
    !structure for storing work data for anisotropic derivative calculations
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
  end type aniso_deriv_data_structured_type

  type, public :: aniso_data_unstructured_type
    ! not strictly necessary here - could be inferred from parent connection_set, but convenient to pass with single ptr
    PetscInt :: conn_face_id      ! id of connection face
    PetscInt :: cell_id_up        ! id of upstream cell
    PetscInt :: cell_id_dn        ! id of downstream cell
    ! adjacency info to up/dnstream cells
    PetscInt, pointer :: cell_id_up_adj(:) ! ids of cells adjacent to the upwind   cell that share a vertex with the connection face
    PetscInt, pointer :: cell_id_dn_adj(:) ! ids of cells adjacent to the downwind cell that share a vertex with the connection face
    PetscInt, pointer :: face_id_up_adj(:) ! ids of faces of the upwind   cell that share a vertex with the connection face
    PetscInt, pointer :: face_id_dn_adj(:) ! ids of faces of the downwind cell that share a vertex with the connection face
    ! (u,v,w) is a coord system relative to face.
    PetscReal, pointer :: uvec(:)      ! u-direction unit vector (normal to face, +ve dnstream) - in (x,y,z) coords
    PetscReal, pointer :: vvec(:)      ! v-direction unit vector (in plane of face) - in (x,y,z) coords
    PetscReal, pointer :: wvec(:)      ! w-direction unit vector (in plane of face) - in (x,y,z) coords

    ! Permeability w.r.t. (u,v,w) coords in up/dnstream cells - working data, overwritten throughout simulation (in case timedep)
    PetscReal :: perm_up_uu
    PetscReal :: perm_up_vv
    PetscReal :: perm_up_ww
    PetscReal :: perm_up_uv
    PetscReal :: perm_up_uw
    PetscReal :: perm_up_vw
    PetscReal :: perm_dn_uu
    PetscReal :: perm_dn_vv
    PetscReal :: perm_dn_ww
    PetscReal :: perm_dn_uv
    PetscReal :: perm_dn_uw
    PetscReal :: perm_dn_vw
  end type aniso_data_unstructured_type

  type, public :: aniso_data_type
    type( aniso_data_structured_type), pointer :: struct => null()
    type( aniso_data_unstructured_type), pointer :: unstruct => null()
  end type aniso_data_type



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
    subroutine PrintAnisoData( aniso_data )

      type( aniso_data_type ), pointer :: aniso_data
      integer ierr,rank,size,i

      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)

      if(associated(aniso_data%struct)) then
        write(*, "('aniso_data structured ids [',i3,'/',i3,'] = ', 12i4)") &
          rank,size, &
          aniso_data%struct%id_up_adjXm, &
          aniso_data%struct%id_up_adjXp, &
          aniso_data%struct%id_dn_adjXm, &
          aniso_data%struct%id_dn_adjXp, &
          aniso_data%struct%id_up_adjYm, &
          aniso_data%struct%id_up_adjYp, &
          aniso_data%struct%id_dn_adjYm, &
          aniso_data%struct%id_dn_adjYp, &
          aniso_data%struct%id_up_adjZm, &
          aniso_data%struct%id_up_adjZp, &
          aniso_data%struct%id_dn_adjZm, &
          aniso_data%struct%id_dn_adjZp
      endif

      if(associated(aniso_data%unstruct)) then
        write(*, "('aniso_data unstructured [',i3,'/',i3,'] face_id = ', i4)") &
          rank,size,aniso_data%unstruct%conn_face_id
        ! geometry
        write(*, "('  cell_id_up = ', i4)") aniso_data%unstruct%cell_id_up
        write(*, "('  cell_id_dn = ', i4)") aniso_data%unstruct%cell_id_dn
        write(*, "('  upstream')")
        do i=1,4 ! 4=MAX_VERT_PER_FACE
          if (aniso_data%unstruct%cell_id_up_adj(i) == 0) cycle
          write(*, "('    ',i2,'. adj_cell_id=',i4,' adj_face_id=',i4)") &
            i,aniso_data%unstruct%cell_id_up_adj(i), aniso_data%unstruct%face_id_up_adj(i)
        enddo
        write(*, "('  downstream')")
        do i=1,4 !4=MAX_VERT_PER_FACE
          if (aniso_data%unstruct%cell_id_up_adj(i) == 0) cycle
          write(*, "('    ',i2,'. adj_cell_id=',i4,' adj_face_id=',i4)") &
            i,aniso_data%unstruct%cell_id_dn_adj(i), aniso_data%unstruct%face_id_dn_adj(i)
        enddo
        ! transformed coords for face
        write(*, "('  local coords')")
        write(*, "('    u = (',e10.3,',',e10.3,',',e10.3,')')") &
          aniso_data%unstruct%uvec(1),aniso_data%unstruct%uvec(2),aniso_data%unstruct%uvec(3)
        write(*, "('    v = (',e10.3,',',e10.3,',',e10.3,')')") &
          aniso_data%unstruct%vvec(1),aniso_data%unstruct%vvec(2),aniso_data%unstruct%vvec(3)
        write(*, "('    w = (',e10.3,',',e10.3,',',e10.3,')')") &
          aniso_data%unstruct%wvec(1),aniso_data%unstruct%wvec(2),aniso_data%unstruct%wvec(3)
        ! transformed perm for face
        write(*, "('  perm up,dn')")
        write(*, "('    perm_uu = ',e10.3,',',e10.3)") &
          aniso_data%unstruct%perm_up_uu, aniso_data%unstruct%perm_dn_uu
        write(*, "('    perm_vv = ',e10.3,',',e10.3)") &
          aniso_data%unstruct%perm_up_vv, aniso_data%unstruct%perm_dn_vv
        write(*, "('    perm_ww = ',e10.3,',',e10.3)") &
          aniso_data%unstruct%perm_up_ww, aniso_data%unstruct%perm_dn_ww
        write(*, "('    perm_uv = ',e10.3,',',e10.3)") &
          aniso_data%unstruct%perm_up_uv, aniso_data%unstruct%perm_dn_uv
        write(*, "('    perm_uw = ',e10.3,',',e10.3)") &
          aniso_data%unstruct%perm_up_uw, aniso_data%unstruct%perm_dn_uw
        write(*, "('    perm_vw = ',e10.3,',',e10.3)") &
          aniso_data%unstruct%perm_up_vw, aniso_data%unstruct%perm_dn_vw

      endif
    end subroutine


end module Anisotropy_Data_module
