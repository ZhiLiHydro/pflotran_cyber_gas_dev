module Anisotropy_Geom_Data_module

! Module for managing additional data needed in flux calculations to
! deal with full (off-diagonal) anisotropy terms

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  ! List process model data types that use anisotropy geom data here
  use Anisotropy_Richards_Data_module

  implicit none

  public :: PrintAnisoData

  ! -- 
  ! -- Data structures for calculating anisotropic gradient contributions in STRUCTURED GRIDS.
  ! --

  type, public :: aniso_geom_data_structured_type

    !  For calculating the flux across interface marked '---->' in a structured grid setting.
    !  Including data from cells 'up' and 'dn' only allows flux based on, e.g., permeability components
    !  in the direction of the flux to be calculated.  Anisotropic permeablity components will
    !  be ignored.
    !  To include anisotropic components, the gradient in the plane of interface must
    !  also be computed.  In a 2-D structured grid, this can be done by approximating the
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
    !  Note. X_A, X_B should be ordered so that the coordinate X(=X,Y,or Z) is increasing the
    !  direction A->B, in order to make sense of direction of gradients.  In aniso_geom_data_structured_type
    !  the p(lus) and m(inus) suffix indicates this.

    ! Orientation (of normal to interface)
    character (len=1) :: orientation ! 'X', 'Y' or 'Z'

    ! Grid IDs - in 3-D, there are 4 adjacent cells in each of the up/dn directions that must be considered
    !            when calculating the head gradient, so 8 in total.  12 are stored here, but only 
    !            8 are used - depending on orientation of the interface.  E.g. for interface with
    !            normal in the x-direction, only the y,z terms will be used.
    PetscInt :: id_up_adjXm ! id of upwind   cell adjacent in the X -ve direction
    PetscInt :: id_up_adjXp ! id of upwind   cell adjacent in the X +ve direction
    PetscInt :: id_dn_adjXm ! id of downwind cell adjacent in the X -ve direction
    PetscInt :: id_dn_adjXp ! id of downwind cell adjacent in the X +ve direction
    PetscInt :: id_up_adjYm ! id of upwind   cell adjacent in the Y -ve direction
    PetscInt :: id_up_adjYp ! id of upwind   cell adjacent in the Y +ve direction
    PetscInt :: id_dn_adjYm ! id of downwind cell adjacent in the Y -ve direction
    PetscInt :: id_dn_adjYp ! id of downwind cell adjacent in the Y +ve direction
    PetscInt :: id_up_adjZm ! id of upwind   cell adjacent in the Z -ve direction
    PetscInt :: id_up_adjZp ! id of upwind   cell adjacent in the Z +ve direction
    PetscInt :: id_dn_adjZm ! id of downwind cell adjacent in the Z -ve direction
    PetscInt :: id_dn_adjZp ! id of downwind cell adjacent in the Z +ve direction

    ! Gravity terms (usually = g.z where z is z coord of relevant cell centre.  Only 8 are used.
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

    ! Distances from up/dn cell to adjacent cells (distances between cell centres).  Only 8 are used.
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
  end type aniso_geom_data_structured_type

  ! -- 
  ! -- Data structures for calculating anisotropic gradient contributions in UNSTRUCTURED GRIDS.
  ! --

  ! workarounds to allow array of pointers to be stored
  type, public :: aduc_ptr_container
    type(aniso_geom_unstructured_cell_type), pointer:: ptr => null()
  end type aduc_ptr_container

  ! represents adjacency data for one of the cells (up or dn) neighbouring an interface
  type, public :: aniso_geom_unstructured_cell_type
    PetscInt :: cell_id        ! id of cell (up or dn)
    PetscReal :: cell_pos(3)   ! (x,y,z coord of cell centre)
    ! adjacency info
    PetscInt :: num_adj ! number of adjacent cells (<=4)
    PetscInt, pointer :: cell_id_adj(:) ! ids of adjacent cells that share an edge with the connection face. -ve if ghost cell
    PetscInt, pointer :: face_id_adj(:) ! ids of faces of cell that share an edge with the connection face.  -ve if ghost-ghost face

    type(aduc_ptr_container), pointer :: adu_cell_adj(:) ! aniso data for cell_id_adj(:)

    ! g.c, where c is cell centre (this are differenced to get g.d, where d is distance between points)
    PetscReal :: loc_gravity

    ! (u,v,w) is coordinate system orthogonal to the common face.  u is normal to the face and (v,w) is in the plane of the face.    
    ! Interpolation machinery - 2x2 matrix inverse is required.
    PetscReal :: LSMatInv(2,2)
    PetscReal, pointer :: dists_u(:)    ! num_adj
    PetscReal, pointer :: dists_vw(:,:) ! 2 x num_adj
  end type aniso_geom_unstructured_cell_type

  type, public :: aniso_geom_data_unstructured_type
    ! not strictly necessary here - could be inferred from parent connection_set, but convenient to pass with single ptr
    PetscInt :: conn_face_id      ! id of connection face

    type(aniso_geom_unstructured_cell_type), pointer :: up => null()
    type(aniso_geom_unstructured_cell_type), pointer :: dn => null()

    ! (u,v,w) is a coord system relative to face.
    PetscReal, pointer :: uvec(:)      ! u-direction unit vector (normal to face, +ve dnstream) - in (x,y,z) coords
    PetscReal, pointer :: vvec(:)      ! v-direction unit vector (in plane of face) - in (x,y,z) coords
    PetscReal, pointer :: wvec(:)      ! w-direction unit vector (in plane of face) - in (x,y,z) coords
  end type aniso_geom_data_unstructured_type

  ! --
  ! -- General data type containing structured and unstructured ptrs
  ! --

  type, public :: aniso_geom_data_type
    type( aniso_geom_data_structured_type), pointer :: struct => null()
    type( aniso_geom_data_unstructured_type), pointer :: unstruct => null()

    ! List process model data types that use anisotropy geom data here
    type( aniso_richards_data_type), pointer :: aniso_rich => null()
  end type aniso_geom_data_type


  ! --
  ! -- Methods
  ! --
  
  CONTAINS
    ! for debugging:
    subroutine PrintAnisoData( aniso_geom )

      type( aniso_geom_data_type ), pointer :: aniso_geom
      integer ierr,rank,size,i

      call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)

      if(associated(aniso_geom%struct)) then
        write(*, "('aniso_geom structured ids [',i3,'/',i3,'] = ', 12i4)") &
          rank,size, &
          aniso_geom%struct%id_up_adjXm, &
          aniso_geom%struct%id_up_adjXp, &
          aniso_geom%struct%id_dn_adjXm, &
          aniso_geom%struct%id_dn_adjXp, &
          aniso_geom%struct%id_up_adjYm, &
          aniso_geom%struct%id_up_adjYp, &
          aniso_geom%struct%id_dn_adjYm, &
          aniso_geom%struct%id_dn_adjYp, &
          aniso_geom%struct%id_up_adjZm, &
          aniso_geom%struct%id_up_adjZp, &
          aniso_geom%struct%id_dn_adjZm, &
          aniso_geom%struct%id_dn_adjZp
      endif

      if(associated(aniso_geom%unstruct)) then
        write(*, "('aniso_geom unstructured [',i3,'/',i3,'] face_id = ', i4)") &
          rank,size,aniso_geom%unstruct%conn_face_id
        ! geometry
        write(*, "('  up cell_id = ', i4,'  (',e10.3,',',e10.3,',',e10.3,')')") aniso_geom%unstruct%up%cell_id, &
          aniso_geom%unstruct%up%cell_pos(1),aniso_geom%unstruct%up%cell_pos(2),aniso_geom%unstruct%up%cell_pos(3)
        write(*, "('  dn cell_id = ', i4,'  (',e10.3,',',e10.3,',',e10.3,')')") aniso_geom%unstruct%dn%cell_id, &
          aniso_geom%unstruct%dn%cell_pos(1),aniso_geom%unstruct%dn%cell_pos(2),aniso_geom%unstruct%dn%cell_pos(3)
        write(*, "('  upstream')")
        do i=1,aniso_geom%unstruct%up%num_adj
          if (aniso_geom%unstruct%up%cell_id_adj(i) == 0) cycle
          write(*, "('    ',i2,'. adj_cell_id=',i4,' adj_face_id=',i4,' (check=',i4,')')") &
            i,aniso_geom%unstruct%up%cell_id_adj(i), aniso_geom%unstruct%up%face_id_adj(i), &
            aniso_geom%unstruct%up%adu_cell_adj(i)%ptr%cell_id
        enddo
        write(*, "('  downstream')")
        do i=1,aniso_geom%unstruct%dn%num_adj
          if (aniso_geom%unstruct%dn%cell_id_adj(i) == 0) cycle
          write(*, "('    ',i2,'. adj_cell_id=',i4,' adj_face_id=',i4,' (check=',i4,',',i4,')')") &
            i,aniso_geom%unstruct%dn%cell_id_adj(i), aniso_geom%unstruct%dn%face_id_adj(i), &
            aniso_geom%unstruct%dn%adu_cell_adj(i)%ptr%cell_id
        enddo
        ! transformed coords for face
        write(*, "('  local coords')")
        write(*, "('    u = (',e10.3,',',e10.3,',',e10.3,')')") &
          aniso_geom%unstruct%uvec(1),aniso_geom%unstruct%uvec(2),aniso_geom%unstruct%uvec(3)
        write(*, "('    v = (',e10.3,',',e10.3,',',e10.3,')')") &
          aniso_geom%unstruct%vvec(1),aniso_geom%unstruct%vvec(2),aniso_geom%unstruct%vvec(3)
        write(*, "('    w = (',e10.3,',',e10.3,',',e10.3,')')") &
          aniso_geom%unstruct%wvec(1),aniso_geom%unstruct%wvec(2),aniso_geom%unstruct%wvec(3)
      endif
    end subroutine

end module Anisotropy_Geom_Data_module
