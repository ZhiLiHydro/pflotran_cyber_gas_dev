module Grid_Unstructured_Voronoi_module

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
  use petscsys
  use petscmat
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Explicit_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: UGridExplicitGenerateVoronoi

contains

! ************************************************************************** !

subroutine UGridExplicitGenerateVoronoi(unstructured_grid, filename, option)
  !
  ! Generate dual control-volume geometry from a primal mesh and route it
  ! into PFLOTRAN's existing explicit-grid reader path. This routine validates
  ! user options, dispatches format-specific primal mesh readers, computes
  ! Voronoi coefficients, writes a temporary *uge, and loads it as explicit-grid
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  type(unstructured_explicit_type), pointer :: explicit_grid
  character(len=MAXSTRINGLENGTH) :: uge_filename
  PetscInt :: status
  PetscInt :: ierr

  PetscInt :: num_pts
  PetscInt :: num_elem
  PetscInt, allocatable :: elem_conn(:,:)
  PetscReal, allocatable :: xyz(:,:)
  PetscReal, allocatable :: cell_vol(:)
  PetscInt, allocatable :: edge_i(:)
  PetscInt, allocatable :: edge_j(:)
  PetscReal, allocatable :: edge_area(:)
  PetscInt :: num_edges
  PetscInt :: io_rank

  explicit_grid => unstructured_grid%explicit_grid

  status = 0
  io_rank = option%comm%io_rank
  uge_filename = trim(filename) // '.pflotran_voronoi_tmp.uge'

  if (explicit_grid%explicit_grid_source /= EXPLICIT_GRID_SOURCE_PRIMAL_VORONOI) then
    option%io_buffer = 'UGridExplicitGenerateVoronoi called for non-primal source.'
    call PrintErrMsg(option)
  endif

  if (explicit_grid%control_volume_type /= CONTROL_VOLUME_TYPE_VORONOI) then
    option%io_buffer = 'Only CONTROL_VOLUME_TYPE VORONOI is currently supported.'
    call PrintErrMsg(option)
  endif

  if (OptionIsIORank(option)) then
    select case (explicit_grid%primal_mesh_format)
      case (PRIMAL_MESH_FORMAT_PFLOTRAN_UGI)
        call UGridVoronoiReadPrimalUGI(trim(filename),option,num_pts,num_elem,elem_conn,xyz)
      case (PRIMAL_MESH_FORMAT_AVS_UCD)
        call UGridVoronoiReadPrimalAVS(trim(filename),option,num_pts,num_elem,elem_conn,xyz)
      case (PRIMAL_MESH_FORMAT_TETGEN)
        call UGridVoronoiReadPrimalTetGen(trim(filename),option,num_pts,num_elem,elem_conn,xyz)
      case default
        status = 1
        option%io_buffer = 'PRIMAL_MESH_FORMAT must be PFLOTRAN_UGI, AVS_UCD, or TETGEN for PRIMAL_VORONOI.'
    end select

    if (status == 0) then
      ! fall through to global broadcast + parallel Voronoi assembly
    endif
  else
    num_pts = 0
    num_elem = 0
  endif

  call MPI_Bcast(status,ONE_INTEGER_MPI,MPI_INTEGER,io_rank, &
                 option%mycomm,ierr);CHKERRQ(ierr)
  if (status /= 0) then
    option%io_buffer = 'Failed generating Voronoi dual mesh internally.'
    call PrintErrMsg(option)
  endif

  ! Broadcast primal mesh metadata/arrays so Voronoi work can run on all ranks.
  call MPI_Bcast(num_pts,ONE_INTEGER_MPI,MPIU_INTEGER,io_rank, &
                 option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Bcast(num_elem,ONE_INTEGER_MPI,MPIU_INTEGER,io_rank, &
                 option%mycomm,ierr);CHKERRQ(ierr)

  if (.not.OptionIsIORank(option)) then
    allocate(elem_conn(4,num_elem))
    allocate(xyz(3,num_pts))
  endif

  call MPI_Bcast(elem_conn,4*num_elem,MPIU_INTEGER,io_rank, &
                 option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Bcast(xyz,3*num_pts,MPI_DOUBLE_PRECISION,io_rank, &
                 option%mycomm,ierr);CHKERRQ(ierr)

  ! Compute dual geometry in parallel, then write on IO rank for reuse by
  ! existing explicit-grid reader/decomposition paths.
  call UGridVoronoiComputeFromTets(num_pts,num_elem,elem_conn,xyz,cell_vol, &
                                   edge_i,edge_j,edge_area,num_edges,option)
  if (OptionIsIORank(option)) then
    call UGridVoronoiWriteUGE(trim(uge_filename),num_pts,xyz,cell_vol, &
                              num_edges,edge_i,edge_j,edge_area,option)
  endif

  if (allocated(elem_conn)) deallocate(elem_conn)
  if (allocated(xyz)) deallocate(xyz)
  if (allocated(cell_vol)) deallocate(cell_vol)
  if (allocated(edge_i)) deallocate(edge_i)
  if (allocated(edge_j)) deallocate(edge_j)
  if (allocated(edge_area)) deallocate(edge_area)

  call UGridExplicitRead(unstructured_grid,trim(uge_filename),option)

end subroutine UGridExplicitGenerateVoronoi

! ************************************************************************** !

subroutine UGridVoronoiReadPrimalUGI(ugi_filename, option, num_pts, num_elem, elem_conn, xyz)
  !
  ! Read a PFLOTRAN *ugi primal mesh into in-memory vertex coordinates and
  ! tetrahedral connectivity arrays used by the internal Voronoi generator.
  ! Only tets are supported.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module
  use Input_Aux_module

  implicit none

  character(len=*), intent(in) :: ugi_filename
  type(option_type) :: option
  PetscInt, intent(out) :: num_pts
  PetscInt, intent(out) :: num_elem
  PetscInt, allocatable, intent(out) :: elem_conn(:,:)
  PetscReal, allocatable, intent(out) :: xyz(:,:)

  type(input_type), pointer :: input
  PetscInt :: fid
  PetscInt :: i
  character(len=MAXWORDLENGTH) :: etype

  fid = 87
  input => InputCreate(fid,trim(ugi_filename),option)
  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,num_elem)
  call InputErrorMsg(input,option,'number of elements','UGI header')
  call InputReadInt(input,option,num_pts)
  call InputErrorMsg(input,option,'number of points','UGI header')

  allocate(elem_conn(4,num_elem))
  elem_conn = 0

  do i = 1, num_elem
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,etype)
    call InputErrorMsg(input,option,'element type','UGI element record')
    call InputReadInt(input,option,elem_conn(1,i))
    call InputErrorMsg(input,option,'vertex 1 id','UGI element record')
    call InputReadInt(input,option,elem_conn(2,i))
    call InputErrorMsg(input,option,'vertex 2 id','UGI element record')
    call InputReadInt(input,option,elem_conn(3,i))
    call InputErrorMsg(input,option,'vertex 3 id','UGI element record')
    call InputReadInt(input,option,elem_conn(4,i))
    call InputErrorMsg(input,option,'vertex 4 id','UGI element record')
    call UGridVoronoiToUpper(etype)
    if (trim(etype) /= 'T') then
      call InputDestroy(input)
      option%io_buffer = 'UGI primal mesh currently supports only T (tetrahedron) elements for PRIMAL_VORONOI.'
      call PrintErrMsg(option)
    endif
  enddo

  allocate(xyz(3,num_pts))
  do i = 1, num_pts
    call InputReadPflotranString(input,option)
    call InputReadDouble(input,option,xyz(1,i))
    call InputErrorMsg(input,option,'vertex x coordinate','UGI vertex record')
    call InputReadDouble(input,option,xyz(2,i))
    call InputErrorMsg(input,option,'vertex y coordinate','UGI vertex record')
    call InputReadDouble(input,option,xyz(3,i))
    call InputErrorMsg(input,option,'vertex z coordinate','UGI vertex record')
  enddo

  call InputDestroy(input)

end subroutine UGridVoronoiReadPrimalUGI

! ************************************************************************** !

subroutine UGridVoronoiReadPrimalAVS(avs_filename, option, num_pts, num_elem, elem_conn, xyz)
  !
  ! Read an AVS-UCD primal mesh into in-memory vertex coordinates and
  ! tetrahedral connectivity arrays used by the internal Voronoi generator.
  ! The routine enforces tet-only support for this initial implementation.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module
  use Input_Aux_module

  implicit none

  character(len=*), intent(in) :: avs_filename
  type(option_type) :: option
  PetscInt, intent(out) :: num_pts
  PetscInt, intent(out) :: num_elem
  PetscInt, allocatable, intent(out) :: elem_conn(:,:)
  PetscReal, allocatable, intent(out) :: xyz(:,:)

  type(input_type), pointer :: input
  PetscInt :: fid
  PetscInt :: i
  PetscInt :: n_node_atts, n_elem_atts, n_model_atts
  PetscInt :: id, mat
  character(len=MAXWORDLENGTH) :: etype

  fid = 86
  input => InputCreate(fid,trim(avs_filename),option)
  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,num_pts)
  call InputErrorMsg(input,option,'number of points','AVS header')
  call InputReadInt(input,option,num_elem)
  call InputErrorMsg(input,option,'number of elements','AVS header')
  call InputReadInt(input,option,n_node_atts)
  call InputErrorMsg(input,option,'number of node attributes','AVS header')
  call InputReadInt(input,option,n_elem_atts)
  call InputErrorMsg(input,option,'number of element attributes','AVS header')
  call InputReadInt(input,option,n_model_atts)
  call InputErrorMsg(input,option,'number of model attributes','AVS header')

  allocate(xyz(3,num_pts))
  do i = 1, num_pts
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,id)
    call InputErrorMsg(input,option,'node id','AVS node record')
    call InputReadDouble(input,option,xyz(1,i))
    call InputErrorMsg(input,option,'node x coordinate','AVS node record')
    call InputReadDouble(input,option,xyz(2,i))
    call InputErrorMsg(input,option,'node y coordinate','AVS node record')
    call InputReadDouble(input,option,xyz(3,i))
    call InputErrorMsg(input,option,'node z coordinate','AVS node record')
  enddo

  allocate(elem_conn(4,num_elem))
  elem_conn = 0
  do i = 1, num_elem
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,id)
    call InputErrorMsg(input,option,'element id','AVS element record')
    call InputReadInt(input,option,mat)
    call InputErrorMsg(input,option,'material id','AVS element record')
    call InputReadCard(input,option,etype)
    call InputErrorMsg(input,option,'element type','AVS element record')
    call InputReadInt(input,option,elem_conn(1,i))
    call InputErrorMsg(input,option,'vertex 1 id','AVS element record')
    call InputReadInt(input,option,elem_conn(2,i))
    call InputErrorMsg(input,option,'vertex 2 id','AVS element record')
    call InputReadInt(input,option,elem_conn(3,i))
    call InputErrorMsg(input,option,'vertex 3 id','AVS element record')
    call InputReadInt(input,option,elem_conn(4,i))
    call InputErrorMsg(input,option,'vertex 4 id','AVS element record')
    call UGridVoronoiToUpper(etype)
    if (trim(etype) /= 'TET') then
      call InputDestroy(input)
      option%io_buffer = 'AVS primal mesh currently supports only tet elements for PRIMAL_VORONOI.'
      call PrintErrMsg(option)
    endif
  enddo

  call InputDestroy(input)

end subroutine UGridVoronoiReadPrimalAVS

! ************************************************************************** !

subroutine UGridVoronoiReadPrimalTetGen(mesh_filename, option, num_pts, num_elem, elem_conn, xyz)
  !
  ! Read a TetGen primal mesh from .node/.ele files into in-memory vertex
  ! coordinates and tetrahedral connectivity arrays used by the internal
  ! Voronoi generator. The routine is for 4-node tet elements only.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module
  use Input_Aux_module

  implicit none

  character(len=*), intent(in) :: mesh_filename
  type(option_type) :: option
  PetscInt, intent(out) :: num_pts
  PetscInt, intent(out) :: num_elem
  PetscInt, allocatable, intent(out) :: elem_conn(:,:)
  PetscReal, allocatable, intent(out) :: xyz(:,:)

  character(len=MAXSTRINGLENGTH) :: node_filename, ele_filename
  type(input_type), pointer :: input
  PetscInt :: fid
  PetscInt :: i
  PetscInt :: num_dim, num_node_atts, node_has_markers
  PetscInt :: num_nodes_per_elem, num_elem_atts
  PetscInt :: id
  PetscInt :: id_min, id_max, id_shift
  PetscReal :: x, y, z
  PetscInt :: m1, m2, m3, m4

  call UGridVoronoiTetGenFilenames(trim(mesh_filename),node_filename,ele_filename)

  fid = 84
  input => InputCreate(fid,trim(node_filename),option)
  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,num_pts)
  call InputErrorMsg(input,option,'number of points','TetGen .node header')
  call InputReadInt(input,option,num_dim)
  call InputErrorMsg(input,option,'node dimension','TetGen .node header')
  call InputReadInt(input,option,num_node_atts)
  call InputErrorMsg(input,option,'number of node attributes','TetGen .node header')
  call InputReadInt(input,option,node_has_markers)
  call InputErrorMsg(input,option,'node marker flag','TetGen .node header')
  if (num_dim /= 3) then
    call InputDestroy(input)
    option%io_buffer = 'TetGen primal mesh must be 3D (node dimension = 3).'
    call PrintErrMsg(option)
  endif

  allocate(xyz(3,num_pts))
  id_shift = 0
  id_min = huge(id_min)
  id_max = -huge(id_max)
  do i = 1, num_pts
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,id)
    call InputErrorMsg(input,option,'node id','TetGen .node record')
    call InputReadDouble(input,option,x)
    call InputErrorMsg(input,option,'node x coordinate','TetGen .node record')
    call InputReadDouble(input,option,y)
    call InputErrorMsg(input,option,'node y coordinate','TetGen .node record')
    call InputReadDouble(input,option,z)
    call InputErrorMsg(input,option,'node z coordinate','TetGen .node record')
    id_min = min(id_min,id)
    id_max = max(id_max,id)
    if (id < 1 .or. id > num_pts) then
      ! id may be 0-based or sparse; remap below after reading all ids
    else
      xyz(1,id) = x
      xyz(2,id) = y
      xyz(3,id) = z
    endif
  enddo
  call InputDestroy(input)

  ! Handle 0-based TetGen ids by shifting to 1-based indexing.
  if (id_min == 0 .and. id_max == num_pts-1) then
    id_shift = 1
  elseif (id_min == 1 .and. id_max == num_pts) then
    id_shift = 0
  else
    option%io_buffer = 'TetGen node ids must be contiguous and either 0-based or 1-based.'
    call PrintErrMsg(option)
  endif

  if (id_shift /= 0) then
    input => InputCreate(fid,trim(node_filename),option)
    call InputReadPflotranString(input,option)
    xyz = 0.d0
    do i = 1, num_pts
      call InputReadPflotranString(input,option)
      call InputReadInt(input,option,id)
      call InputErrorMsg(input,option,'node id','TetGen .node record')
      call InputReadDouble(input,option,x)
      call InputErrorMsg(input,option,'node x coordinate','TetGen .node record')
      call InputReadDouble(input,option,y)
      call InputErrorMsg(input,option,'node y coordinate','TetGen .node record')
      call InputReadDouble(input,option,z)
      call InputErrorMsg(input,option,'node z coordinate','TetGen .node record')
      xyz(1,id+id_shift) = x
      xyz(2,id+id_shift) = y
      xyz(3,id+id_shift) = z
    enddo
    call InputDestroy(input)
  endif

  fid = 85
  input => InputCreate(fid,trim(ele_filename),option)
  call InputReadPflotranString(input,option)
  call InputReadInt(input,option,num_elem)
  call InputErrorMsg(input,option,'number of elements','TetGen .ele header')
  call InputReadInt(input,option,num_nodes_per_elem)
  call InputErrorMsg(input,option,'nodes per element','TetGen .ele header')
  call InputReadInt(input,option,num_elem_atts)
  call InputErrorMsg(input,option,'number of element attributes','TetGen .ele header')

  if (num_nodes_per_elem /= 4) then
    call InputDestroy(input)
    option%io_buffer = 'TetGen primal mesh must contain 4-node tetrahedra only for PRIMAL_VORONOI.'
    call PrintErrMsg(option)
  endif

  allocate(elem_conn(4,num_elem))
  do i = 1, num_elem
    call InputReadPflotranString(input,option)
    call InputReadInt(input,option,id)
    call InputErrorMsg(input,option,'element id','TetGen .ele record')
    call InputReadInt(input,option,m1)
    call InputErrorMsg(input,option,'vertex 1 id','TetGen .ele record')
    call InputReadInt(input,option,m2)
    call InputErrorMsg(input,option,'vertex 2 id','TetGen .ele record')
    call InputReadInt(input,option,m3)
    call InputErrorMsg(input,option,'vertex 3 id','TetGen .ele record')
    call InputReadInt(input,option,m4)
    call InputErrorMsg(input,option,'vertex 4 id','TetGen .ele record')
    elem_conn(1,i) = m1 + id_shift
    elem_conn(2,i) = m2 + id_shift
    elem_conn(3,i) = m3 + id_shift
    elem_conn(4,i) = m4 + id_shift
    if (minval(elem_conn(:,i)) < 1 .or. maxval(elem_conn(:,i)) > num_pts) then
      call InputDestroy(input)
      option%io_buffer = 'TetGen element connectivity references an out-of-range vertex id.'
      call PrintErrMsg(option)
    endif
  enddo
  call InputDestroy(input)

end subroutine UGridVoronoiReadPrimalTetGen

! ************************************************************************** !

subroutine UGridVoronoiTetGenFilenames(mesh_filename, node_filename, ele_filename)
  !
  ! Resolve TetGen .node/.ele filenames from an input mesh token. Supports
  ! base-name input or direct .node/.ele filename input.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  implicit none

  character(len=*), intent(in) :: mesh_filename
  character(len=*), intent(out) :: node_filename
  character(len=*), intent(out) :: ele_filename
  PetscInt :: n

  n = len_trim(mesh_filename)
  if (n >= 5 .and. mesh_filename(n-4:n) == '.node') then
    node_filename = trim(mesh_filename)
    ele_filename = trim(mesh_filename(1:n-5)) // '.ele'
  elseif (n >= 4 .and. mesh_filename(n-3:n) == '.ele') then
    ele_filename = trim(mesh_filename)
    node_filename = trim(mesh_filename(1:n-4)) // '.node'
  else
    node_filename = trim(mesh_filename) // '.node'
    ele_filename = trim(mesh_filename) // '.ele'
  endif

end subroutine UGridVoronoiTetGenFilenames

subroutine UGridVoronoiComputeFromTets(num_pts, num_elem, elem_conn, xyz, cell_vol, &
                                       edge_i, edge_j, edge_area, num_edges, option)
  !
  ! Compute vertex-centered Voronoi geometry from tetrahedral primal elements.
  ! For each tetrahedron, this routine computes sub-control-volume contributions
  ! and dual-face area contributions, then assembles global per-vertex volumes
  ! and per-edge connection areas. It also checks for invalid/non-positive
  ! resulting control-volume measures.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module

  implicit none

  PetscInt, intent(in) :: num_pts
  PetscInt, intent(in) :: num_elem
  PetscInt, intent(in) :: elem_conn(4,num_elem)
  PetscReal, intent(in) :: xyz(3,num_pts)
  PetscReal, allocatable, intent(out) :: cell_vol(:)
  PetscInt, allocatable, intent(out) :: edge_i(:)
  PetscInt, allocatable, intent(out) :: edge_j(:)
  PetscReal, allocatable, intent(out) :: edge_area(:)
  PetscInt, intent(out) :: num_edges
  type(option_type) :: option

  PetscInt :: i, j
  PetscErrorCode :: ierr
  PetscInt :: v(4)
  PetscReal :: v1(3), v2(3), v3(3), v4(3)
  PetscReal :: cc(3), c123(3), c124(3), c134(3), c234(3)
  PetscReal :: e12(3), e13(3), e14(3), e23(3), e24(3), e34(3)
  PetscReal :: area(6), vol(4), len(6)
  PetscInt :: iedge(6,2)
  PetscInt :: elem_start, elem_end
  PetscInt :: irank, nprocs
  PetscReal, allocatable :: cell_vol_local(:)
  Mat :: edge_area_mat, edge_conn_mat
  PetscInt :: rstart, rend, row, ncols
  PetscInt, pointer :: cols(:)
  PetscReal, pointer :: vals(:)
  PetscInt :: local_edges
  PetscInt, allocatable :: edge_i_local(:), edge_j_local(:)
  PetscReal, allocatable :: edge_area_local(:)
  PetscInt, allocatable :: recvcounts(:), displs(:)
  PetscReal, allocatable :: cell_vol_sum(:)
  PetscInt :: io_rank
  PetscInt :: total_edges
  PetscReal :: area_eps
  PetscInt :: row_idx(1), col_idx(1)
  PetscReal :: area_tmp(1)

  io_rank = option%comm%io_rank
  irank = option%comm%rank
  nprocs = option%comm%size
  area_eps = 1.d-30

  allocate(cell_vol_local(num_pts))
  cell_vol_local = 0.d0

  ! area matrix: accumulates total dual-face area per primal edge.
  call MatCreateAIJ(option%mycomm,PETSC_DECIDE,PETSC_DECIDE,num_pts,num_pts, &
                    12,PETSC_NULL_INTEGER_ARRAY,12,PETSC_NULL_INTEGER_ARRAY, &
                    edge_area_mat,ierr);CHKERRQ(ierr)
  call MatZeroEntries(edge_area_mat,ierr);CHKERRQ(ierr)
  ! connectivity matrix: stores edge existence, including zero-area edges.
  call MatCreateAIJ(option%mycomm,PETSC_DECIDE,PETSC_DECIDE,num_pts,num_pts, &
                    12,PETSC_NULL_INTEGER_ARRAY,12,PETSC_NULL_INTEGER_ARRAY, &
                    edge_conn_mat,ierr);CHKERRQ(ierr)
  call MatZeroEntries(edge_conn_mat,ierr);CHKERRQ(ierr)

  ! Static element partition by contiguous index range.
  elem_start = (irank*num_elem)/nprocs + 1
  elem_end = ((irank+1)*num_elem)/nprocs

  do i = elem_start, elem_end
    v = elem_conn(:,i)

    if (minval(v) < 1 .or. maxval(v) > num_pts) then
      option%io_buffer = 'Primal mesh connectivity references an out-of-range vertex id.'
      call PrintErrMsg(option)
    endif

    v1 = xyz(:,v(1))
    v2 = xyz(:,v(2))
    v3 = xyz(:,v(3))
    v4 = xyz(:,v(4))

    call UGridVoronoiTetGeom(v1,v2,v3,v4,cc,c123,c124,c134,c234,e12,e13,e14,e23,e24,e34,area,vol,len)

    cell_vol_local(v(1)) = cell_vol_local(v(1)) + abs(vol(1))
    cell_vol_local(v(2)) = cell_vol_local(v(2)) + abs(vol(2))
    cell_vol_local(v(3)) = cell_vol_local(v(3)) + abs(vol(3))
    cell_vol_local(v(4)) = cell_vol_local(v(4)) + abs(vol(4))

    iedge(1,1) = min(v(1),v(2)) - 1
    iedge(1,2) = max(v(1),v(2)) - 1
    iedge(2,1) = min(v(2),v(3)) - 1
    iedge(2,2) = max(v(2),v(3)) - 1
    iedge(3,1) = min(v(1),v(3)) - 1
    iedge(3,2) = max(v(1),v(3)) - 1
    iedge(4,1) = min(v(1),v(4)) - 1
    iedge(4,2) = max(v(1),v(4)) - 1
    iedge(5,1) = min(v(2),v(4)) - 1
    iedge(5,2) = max(v(2),v(4)) - 1
    iedge(6,1) = min(v(3),v(4)) - 1
    iedge(6,2) = max(v(3),v(4)) - 1
    do j = 1, 6
      call MatSetValue(edge_area_mat,iedge(j,1),iedge(j,2),abs(area(j)), &
                       ADD_VALUES,ierr);CHKERRQ(ierr)
      call MatSetValue(edge_conn_mat,iedge(j,1),iedge(j,2),1.d0, &
                       INSERT_VALUES,ierr);CHKERRQ(ierr)
    enddo
  enddo

  ! Sum per-vertex dual volumes across all ranks.
  allocate(cell_vol_sum(num_pts))
  call MPI_Allreduce(cell_vol_local,cell_vol_sum,num_pts,MPI_DOUBLE_PRECISION, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(edge_area_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(edge_area_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(edge_conn_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(edge_conn_mat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! Extract rank-owned rows from edge connectivity matrix, then query
  ! corresponding areas from edge_area_mat.
  call MatGetOwnershipRange(edge_conn_mat,rstart,rend,ierr);CHKERRQ(ierr)
  local_edges = 0
  do row = rstart, rend-1
    call MatGetRow(edge_conn_mat,row,ncols,cols,vals,ierr);CHKERRQ(ierr)
    do j = 1, ncols
      if (cols(j) > row .and. abs(vals(j)) > area_eps) local_edges = local_edges + 1
    enddo
    call MatRestoreRow(edge_conn_mat,row,ncols,cols,vals,ierr);CHKERRQ(ierr)
  enddo

  allocate(edge_i_local(max(local_edges,1)))
  allocate(edge_j_local(max(local_edges,1)))
  allocate(edge_area_local(max(local_edges,1)))
  edge_i_local = 0
  edge_j_local = 0
  edge_area_local = 0.d0

  local_edges = 0
  do row = rstart, rend-1
    call MatGetRow(edge_conn_mat,row,ncols,cols,vals,ierr);CHKERRQ(ierr)
    do j = 1, ncols
      if (cols(j) > row .and. abs(vals(j)) > area_eps) then
        local_edges = local_edges + 1
        edge_i_local(local_edges) = row + 1
        edge_j_local(local_edges) = cols(j) + 1
        row_idx(1) = row
        col_idx(1) = cols(j)
        call MatGetValues(edge_area_mat,1,row_idx,1,col_idx,area_tmp, &
                          ierr);CHKERRQ(ierr)
        edge_area_local(local_edges) = area_tmp(1)
      endif
    enddo
    call MatRestoreRow(edge_conn_mat,row,ncols,cols,vals,ierr);CHKERRQ(ierr)
  enddo

  if (OptionIsIORank(option)) then
    allocate(cell_vol(num_pts))
    cell_vol = cell_vol_sum
    if (minval(cell_vol) <= 0.d0) then
      option%io_buffer = 'Computed non-positive Voronoi cell volumes from primal mesh.'
      call PrintErrMsg(option)
    endif
    allocate(recvcounts(nprocs))
    allocate(displs(nprocs))
  else
    allocate(cell_vol(1))
    cell_vol(1) = 0.d0
    allocate(recvcounts(1))
    allocate(displs(1))
  endif

  ! Gather variable-length local edge lists to IO rank for UGE output.
  call MPI_Gather(local_edges,ONE_INTEGER_MPI,MPIU_INTEGER,recvcounts, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,io_rank,option%mycomm, &
                  ierr);CHKERRQ(ierr)

  if (OptionIsIORank(option)) then
    displs(1) = 0
    do i = 2, nprocs
      displs(i) = displs(i-1) + recvcounts(i-1)
    enddo
    total_edges = displs(nprocs) + recvcounts(nprocs)
    num_edges = total_edges
    allocate(edge_i(max(num_edges,1)))
    allocate(edge_j(max(num_edges,1)))
    allocate(edge_area(max(num_edges,1)))
  else
    num_edges = 0
    allocate(edge_i(1))
    allocate(edge_j(1))
    allocate(edge_area(1))
  endif

  call MPI_Gatherv(edge_i_local,local_edges,MPIU_INTEGER, &
                   edge_i,recvcounts,displs,MPIU_INTEGER,io_rank, &
                   option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Gatherv(edge_j_local,local_edges,MPIU_INTEGER, &
                   edge_j,recvcounts,displs,MPIU_INTEGER,io_rank, &
                   option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Gatherv(edge_area_local,local_edges,MPI_DOUBLE_PRECISION, &
                   edge_area,recvcounts,displs,MPI_DOUBLE_PRECISION,io_rank, &
                   option%mycomm,ierr);CHKERRQ(ierr)

  call MatDestroy(edge_area_mat,ierr);CHKERRQ(ierr)
  call MatDestroy(edge_conn_mat,ierr);CHKERRQ(ierr)
  deallocate(cell_vol_local)
  deallocate(cell_vol_sum)
  deallocate(edge_i_local)
  deallocate(edge_j_local)
  deallocate(edge_area_local)
  deallocate(recvcounts)
  deallocate(displs)

end subroutine UGridVoronoiComputeFromTets

! ************************************************************************** !

subroutine UGridVoronoiWriteUGE(uge_filename, num_pts, xyz, cell_vol, num_edges, edge_i, edge_j, edge_area, option)
  !
  ! Write internally generated Voronoi dual geometry to PFLOTRAN UGE format.
  ! Cell centroids are written at primal vertices with accumulated dual volumes,
  ! and connection entries are written per primal-edge pair with midpoint-based
  ! connection centroids and assembled dual-face areas.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module

  implicit none

  character(len=*), intent(in) :: uge_filename
  PetscInt, intent(in) :: num_pts
  PetscReal, intent(in) :: xyz(3,num_pts)
  PetscReal, intent(in) :: cell_vol(num_pts)
  PetscInt, intent(in) :: num_edges
  PetscInt, intent(in) :: edge_i(:)
  PetscInt, intent(in) :: edge_j(:)
  PetscReal, intent(in) :: edge_area(:)
  type(option_type) :: option

  integer :: fid
  integer :: ios
  PetscInt :: i
  PetscReal :: xc, yc, zc

  fid = 89
  open(unit=fid,file=trim(uge_filename),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing temporary internal UGE file: ' // trim(uge_filename)
    call PrintErrMsg(option)
  endif

  write(fid,'(a,1x,i0)') 'CELLS', num_pts
  do i = 1, num_pts
    write(fid,'(i0,1x,1pe20.12,1x,1pe20.12,1x,1pe20.12,1x,1pe20.12)') &
      i, xyz(1,i), xyz(2,i), xyz(3,i), cell_vol(i)
  enddo

  write(fid,'(a,1x,i0)') 'CONNECTIONS', num_edges
  do i = 1, num_edges
    xc = 0.5d0*(xyz(1,edge_i(i)) + xyz(1,edge_j(i)))
    yc = 0.5d0*(xyz(2,edge_i(i)) + xyz(2,edge_j(i)))
    zc = 0.5d0*(xyz(3,edge_i(i)) + xyz(3,edge_j(i)))
    write(fid,'(i0,1x,i0,1x,1pe20.12,1x,1pe20.12,1x,1pe20.12,1x,1pe20.12)') &
      edge_i(i), edge_j(i), xc, yc, zc, edge_area(i)
  enddo

  close(fid)

end subroutine UGridVoronoiWriteUGE

! ************************************************************************** !

subroutine UGridVoronoiAccumulateEdge(i_in, j_in, area_in, edge_i, edge_j, edge_area, num_edges)
  !
  ! Insert or update a primal-edge entry in the connection list.
  ! If an edge already exists, the routine accumulates area contribution;
  ! otherwise, it appends a new edge record.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  implicit none

  PetscInt, intent(in) :: i_in, j_in
  PetscReal, intent(in) :: area_in
  PetscInt, intent(inout) :: edge_i(:), edge_j(:)
  PetscReal, intent(inout) :: edge_area(:)
  PetscInt, intent(inout) :: num_edges

  PetscInt :: i, a, b

  a = min(i_in,j_in)
  b = max(i_in,j_in)

  do i = 1, num_edges
    if (edge_i(i) == a .and. edge_j(i) == b) then
      edge_area(i) = edge_area(i) + area_in
      return
    endif
  enddo

  num_edges = num_edges + 1
  edge_i(num_edges) = a
  edge_j(num_edges) = b
  edge_area(num_edges) = area_in

end subroutine UGridVoronoiAccumulateEdge

! ************************************************************************** !

subroutine UGridVoronoiTetGeom(v1,v2,v3,v4,cc,c123,c124,c134,c234,e12,e13,e14,e23,e24,e34,area,vol,len)
  !
  ! Compute tet-local geometric quantities required for Voronoi dual
  ! assembly: edge midpoints, tet circumcenter,
  ! dual-face area contributions for the six primal edges, and per-vertex
  ! sub-control-volume contributions within a tet.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  implicit none

  PetscReal, intent(in) :: v1(3), v2(3), v3(3), v4(3)
  PetscReal, intent(out) :: cc(3), c123(3), c124(3), c134(3), c234(3)
  PetscReal, intent(out) :: e12(3), e13(3), e14(3), e23(3), e24(3), e34(3)
  PetscReal, intent(out) :: area(6), vol(4), len(6)

  PetscReal :: a(6)

  e12 = 0.5d0*(v1 + v2)
  e13 = 0.5d0*(v1 + v3)
  e14 = 0.5d0*(v1 + v4)
  e23 = 0.5d0*(v2 + v3)
  e24 = 0.5d0*(v2 + v4)
  e34 = 0.5d0*(v3 + v4)

  call UGridVoronoiTriangleCircumcenter(v1,v2,v3,c123)
  call UGridVoronoiTriangleCircumcenter(v1,v2,v4,c124)
  call UGridVoronoiTriangleCircumcenter(v1,v3,v4,c134)
  call UGridVoronoiTriangleCircumcenter(v2,v3,v4,c234)

  call UGridVoronoiTetraCircumcenter(v1,v2,v3,v4,cc)

  area(1) = UGridVoronoiDualFaceArea(v2,v1,e12,c124,c123,cc)
  area(2) = UGridVoronoiDualFaceArea(v3,v2,e23,c234,c123,cc)
  area(3) = UGridVoronoiDualFaceArea(v3,v1,e13,c123,c134,cc)
  area(4) = UGridVoronoiDualFaceArea(v4,v1,e14,c134,c124,cc)
  area(5) = UGridVoronoiDualFaceArea(v4,v2,e24,c124,c234,cc)
  area(6) = UGridVoronoiDualFaceArea(v4,v3,e34,c234,c134,cc)

  a(1) = UGridVoronoiTetraVolume(v1,e12,c124,cc)
  a(2) = UGridVoronoiTetraVolume(v1,c124,e14,cc)
  a(3) = UGridVoronoiTetraVolume(v1,c134,e13,cc)
  a(4) = UGridVoronoiTetraVolume(v1,e14,c134,cc)
  a(5) = UGridVoronoiTetraVolume(v1,c123,e12,cc)
  a(6) = UGridVoronoiTetraVolume(v1,e13,c123,cc)
  vol(1) = sum(a)

  a(1) = UGridVoronoiTetraVolume(v2,c124,e12,cc)
  a(2) = UGridVoronoiTetraVolume(v2,c234,e24,cc)
  a(3) = UGridVoronoiTetraVolume(v2,e23,c234,cc)
  a(4) = UGridVoronoiTetraVolume(v2,e12,c123,cc)
  a(5) = UGridVoronoiTetraVolume(v2,c123,e23,cc)
  a(6) = UGridVoronoiTetraVolume(v2,e24,c124,cc)
  vol(2) = sum(a)

  a(1) = UGridVoronoiTetraVolume(v3,c134,e34,cc)
  a(2) = UGridVoronoiTetraVolume(v3,e34,c234,cc)
  a(3) = UGridVoronoiTetraVolume(v3,c234,e23,cc)
  a(4) = UGridVoronoiTetraVolume(v3,e13,c134,cc)
  a(5) = UGridVoronoiTetraVolume(v3,e23,c123,cc)
  a(6) = UGridVoronoiTetraVolume(v3,c123,e13,cc)
  vol(3) = sum(a)

  a(1) = UGridVoronoiTetraVolume(v4,c134,e14,cc)
  a(2) = UGridVoronoiTetraVolume(v4,e14,c124,cc)
  a(3) = UGridVoronoiTetraVolume(v4,c124,e24,cc)
  a(4) = UGridVoronoiTetraVolume(v4,e24,c234,cc)
  a(5) = UGridVoronoiTetraVolume(v4,c234,e34,cc)
  a(6) = UGridVoronoiTetraVolume(v4,e34,c134,cc)
  vol(4) = sum(a)

  len(1) = UGridVoronoiNorm(v1-v2)
  len(2) = UGridVoronoiNorm(v2-v3)
  len(3) = UGridVoronoiNorm(v1-v3)
  len(4) = UGridVoronoiNorm(v1-v4)
  len(5) = UGridVoronoiNorm(v2-v4)
  len(6) = UGridVoronoiNorm(v3-v4)

end subroutine UGridVoronoiTetGeom

! ************************************************************************** !

subroutine UGridVoronoiTriangleCircumcenter(v1,v2,v3,circ)
  !
  ! Compute the circumcenter of a triangle in 3D space. If the triangle is
  ! near-degenerate, fall back to the triangle centroid to keep coefficients
  ! finite and avoid numerical breakdown.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  implicit none

  PetscReal, intent(in) :: v1(3), v2(3), v3(3)
  PetscReal, intent(out) :: circ(3)
  PetscReal :: mm(3), denom, coef(3)

  call UGridVoronoiCross(mm, v1-v2, v2-v3)
  denom = 2.d0*UGridVoronoiNorm(mm)**2

  if (abs(denom) < 1.d-24) then
    circ = (v1+v2+v3)/3.d0
    return
  endif

  coef(1) = UGridVoronoiNorm(v2-v3)**2*dot_product(v1-v2,v1-v3)/denom
  coef(2) = UGridVoronoiNorm(v1-v3)**2*dot_product(v2-v1,v2-v3)/denom
  coef(3) = UGridVoronoiNorm(v1-v2)**2*dot_product(v3-v1,v3-v2)/denom
  circ = coef(1)*v1 + coef(2)*v2 + coef(3)*v3

end subroutine UGridVoronoiTriangleCircumcenter

! ************************************************************************** !

subroutine UGridVoronoiTetraCircumcenter(v1,v2,v3,v4,cent)
  !
  ! Compute the circumcenter of a tetrahedron using determinant-based formulas.
  ! If the tetrahedron is near-degenerate, fall back to the tet centroid.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  implicit none

  PetscReal, intent(in) :: v1(3), v2(3), v3(3), v4(3)
  PetscReal, intent(out) :: cent(3)

  PetscReal :: MatA(4,4), MatDx(4,4), MatDy(4,4), MatDz(4,4)
  PetscReal :: vsum(4)
  PetscReal :: a, dx, dy, dz

  MatA(1:3,1) = v1
  MatA(1:3,2) = v2
  MatA(1:3,3) = v3
  MatA(1:3,4) = v4
  MatA(4,:) = 1.d0

  vsum(1) = dot_product(v1,v1)
  vsum(2) = dot_product(v2,v2)
  vsum(3) = dot_product(v3,v3)
  vsum(4) = dot_product(v4,v4)

  MatDx(1,:) = vsum
  MatDx(2,:) = MatA(2,:)
  MatDx(3,:) = MatA(3,:)
  MatDx(4,:) = 1.d0

  MatDy(1,:) = vsum
  MatDy(2,:) = MatA(1,:)
  MatDy(3,:) = MatA(3,:)
  MatDy(4,:) = 1.d0

  MatDz(1,:) = vsum
  MatDz(2,:) = MatA(1,:)
  MatDz(3,:) = MatA(2,:)
  MatDz(4,:) = 1.d0

  a  = UGridVoronoiDet4x4(MatA(:,1),  MatA(:,2),  MatA(:,3),  MatA(:,4))
  dx = UGridVoronoiDet4x4(MatDx(:,1), MatDx(:,2), MatDx(:,3), MatDx(:,4))
  dy = -UGridVoronoiDet4x4(MatDy(:,1), MatDy(:,2), MatDy(:,3), MatDy(:,4))
  dz = UGridVoronoiDet4x4(MatDz(:,1), MatDz(:,2), MatDz(:,3), MatDz(:,4))

  if (abs(a) < 1.d-24) then
    cent = (v1+v2+v3+v4)/4.d0
  else
    cent(1) = dx/(2.d0*a)
    cent(2) = dy/(2.d0*a)
    cent(3) = dz/(2.d0*a)
  endif

end subroutine UGridVoronoiTetraCircumcenter

! ************************************************************************** !

PetscReal function UGridVoronoiDualFaceArea(vj,vi,eij,c1,c2,cc)

  implicit none

  PetscReal, intent(in) :: vj(3), vi(3), eij(3), c1(3), c2(3), cc(3)
  PetscReal :: a1(3), a2(3), nn(3)

  nn = (vj-vi)/UGridVoronoiNorm(vj-vi)
  call UGridVoronoiCross(a1, c1-eij, cc-eij)
  call UGridVoronoiCross(a2, cc-eij, c2-eij)
  UGridVoronoiDualFaceArea = dot_product(0.5d0*(a1+a2), nn)

end function UGridVoronoiDualFaceArea

! ************************************************************************** !

PetscReal function UGridVoronoiTetraVolume(a,b,c,d)

  implicit none

  PetscReal, intent(in) :: a(3), b(3), c(3), d(3)
  PetscReal :: v21(3), v31(3), v41(3), temp(3)

  v21 = b-a
  v31 = c-a
  v41 = d-a
  call UGridVoronoiCross(temp, v31, v41)
  UGridVoronoiTetraVolume = dot_product(v21,temp)/6.d0

end function UGridVoronoiTetraVolume

! ************************************************************************** !

subroutine UGridVoronoiCross(v12,v1,v2)
  !
  ! Compute 3D vector cross product used by circumcenter, area, and volume
  ! helper calculations in the Voronoi geometry kernel.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  implicit none

  PetscReal, intent(out) :: v12(3)
  PetscReal, intent(in) :: v1(3), v2(3)

  v12(1) = v1(2)*v2(3) - v2(2)*v1(3)
  v12(2) = v2(1)*v1(3) - v1(1)*v2(3)
  v12(3) = v1(1)*v2(2) - v2(1)*v1(2)

end subroutine UGridVoronoiCross

! ************************************************************************** !

PetscReal function UGridVoronoiNorm(v)

  implicit none

  PetscReal, intent(in) :: v(3)
  UGridVoronoiNorm = sqrt(sum(v**2))

end function UGridVoronoiNorm

! ************************************************************************** !

PetscReal function UGridVoronoiDet3x3(v1,v2,v3)

  implicit none

  PetscReal, intent(in) :: v1(3), v2(3), v3(3)
  PetscReal :: c1, c2, c3

  c1 = v1(1)*(v2(2)*v3(3)-v3(2)*v2(3))
  c2 = v1(2)*(v2(1)*v3(3)-v3(1)*v2(3))
  c3 = v1(3)*(v2(1)*v3(2)-v3(1)*v2(2))
  UGridVoronoiDet3x3 = c1 - c2 + c3

end function UGridVoronoiDet3x3

! ************************************************************************** !

PetscReal function UGridVoronoiDet4x4(v1,v2,v3,v4)

  implicit none

  PetscReal, intent(in) :: v1(4), v2(4), v3(4), v4(4)
  PetscReal :: c1, c2, c3, c4

  c1 = v1(1)*UGridVoronoiDet3x3(v2(2:4), v3(2:4), v4(2:4))
  c2 = v1(2)*UGridVoronoiDet3x3(v2((/1,3,4/)), v3((/1,3,4/)), v4((/1,3,4/)))
  c3 = v1(3)*UGridVoronoiDet3x3(v2((/1,2,4/)), v3((/1,2,4/)), v4((/1,2,4/)))
  c4 = v1(4)*UGridVoronoiDet3x3(v2(1:3), v3(1:3), v4(1:3))
  UGridVoronoiDet4x4 = c1 - c2 + c3 - c4

end function UGridVoronoiDet4x4

! ************************************************************************** !

subroutine UGridVoronoiToUpper(s)

  implicit none

  character(len=*), intent(inout) :: s
  integer :: i, ch

  do i = 1, len_trim(s)
    ch = iachar(s(i:i))
    if (ch >= iachar('a') .and. ch <= iachar('z')) s(i:i) = achar(ch-32)
  enddo

end subroutine UGridVoronoiToUpper

end module Grid_Unstructured_Voronoi_module
