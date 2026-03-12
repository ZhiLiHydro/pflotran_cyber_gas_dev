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
  PetscInt, allocatable :: elem_conn_local(:,:)
  PetscInt :: elem_start, elem_end, nlocal
  PetscInt :: num_edges
  PetscInt :: io_rank
  integer :: fid, ios

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
  call UGridVoronoiDistributeElemConn(num_elem, elem_conn, &
                elem_conn_local, elem_start, elem_end, option)

  nlocal = max(0, elem_end-elem_start+1)

  ! Write debug visualization files only when VORONOI_DEBUG_VTK is set.
  if (explicit_grid%voronoi_debug_vtk) then
    ! Write VTU pieces in parallel (all ranks)
    call UGridVoronoiWriteVTUSubTetsPiece(trim(filename)//'.voronoi_subtets', option%comm%rank, &
                                         num_pts, elem_conn_local, nlocal, xyz, option)
    call MPI_Barrier(option%mycomm, ierr)

    ! Rank0 writes PVTU
    if (OptionIsIORank(option)) then
      call UGridVoronoiWritePVTUSubTets(trim(filename)//'.voronoi_subtets', option)
    endif

    call UGridVoronoiWriteVTUDualFacePatchesPiece( trim(filename)//'.dual_patches', &
        option%comm%rank, num_pts, elem_start, elem_conn_local, nlocal, xyz, option )

    call MPI_Barrier(option%mycomm, ierr)
    if (OptionIsIORank(option)) then
      call UGridVoronoiWritePVTUDualFacePatches( trim(filename)//'.dual_patches', option )
    endif

    ! Write stitched Voronoi-face polygons per primal edge as VTU/PVTU.
    ! Only IO rank performs this since it requires global connectivity.
    if (OptionIsIORank(option)) then
      call UGridVoronoiWriteVTUDualEdgePolygons(trim(filename)//'.voronoi_polygons', &
          num_pts, num_elem, elem_conn, xyz, option)
    endif
    call MPI_Barrier(option%mycomm, ierr)

    call UGridVoronoiWriteVTUTetMarkersPiece( trim(filename)//'.tet_markers', &
        option%comm%rank, num_pts, elem_start, elem_conn_local, nlocal, xyz, option )

    call MPI_Barrier(option%mycomm, ierr)
    if (OptionIsIORank(option)) then
      call UGridVoronoiWritePVTUTetMarkers( trim(filename)//'.tet_markers', option )
    endif
  endif

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

  ! Remove the temporary UGE file now that it has been loaded.
  if (OptionIsIORank(option)) then
    open(newunit=fid,file=trim(uge_filename),status='old',iostat=ios)
    if (ios == 0) close(fid,status='delete')
  endif

end subroutine UGridExplicitGenerateVoronoi


subroutine UGridVoronoiDistributeElemConn(num_elem, elem_conn_global, elem_conn_local, &
                                         elem_start, elem_end, option)
  use Option_module
  implicit none
  PetscInt, intent(in) :: num_elem
  PetscInt, allocatable, intent(inout) :: elem_conn_global(:,:) ! only allocated on IO
  PetscInt, allocatable, intent(out) :: elem_conn_local(:,:)    ! (4,nlocal)
  PetscInt, intent(out) :: elem_start, elem_end
  type(option_type) :: option

  PetscInt :: irank, nprocs, io_rank, ierr
  PetscInt :: nlocal, i
  PetscInt, allocatable :: sendcounts(:), displs(:)
  PetscInt :: e0, e1

  irank  = option%comm%rank
  nprocs = option%comm%size
  io_rank = option%comm%io_rank

  ! same contiguous partition you already use
  e0 = (irank*num_elem)/nprocs + 1
  e1 = ((irank+1)*num_elem)/nprocs
  elem_start = e0
  elem_end   = e1
  nlocal = max(0, e1 - e0 + 1)

  allocate(elem_conn_local(4, max(nlocal,1)))
  elem_conn_local = 0

  allocate(sendcounts(nprocs), displs(nprocs))
  if (OptionIsIORank(option)) then
    do i=0,nprocs-1
      e0 = (i*num_elem)/nprocs + 1
      e1 = ((i+1)*num_elem)/nprocs
      sendcounts(i+1) = 4*max(0, e1-e0+1)
    enddo
    displs(1)=0
    do i=2,nprocs
      displs(i)=displs(i-1)+sendcounts(i-1)
    enddo
  else
    sendcounts = 0
    displs = 0
  endif

  call MPI_Scatterv(elem_conn_global, sendcounts, displs, MPIU_INTEGER, &
                    elem_conn_local, 4*nlocal, MPIU_INTEGER, io_rank, option%mycomm, ierr)

  deallocate(sendcounts, displs)

end subroutine UGridVoronoiDistributeElemConn

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
  PetscReal :: area(6), vol(4)
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
  call MatSetOption(edge_area_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr);CHKERRQ(ierr)
  call MatZeroEntries(edge_area_mat,ierr);CHKERRQ(ierr)
  ! connectivity matrix: stores edge existence, including zero-area edges.
  call MatCreateAIJ(option%mycomm,PETSC_DECIDE,PETSC_DECIDE,num_pts,num_pts, &
                    12,PETSC_NULL_INTEGER_ARRAY,12,PETSC_NULL_INTEGER_ARRAY, &
                    edge_conn_mat,ierr);CHKERRQ(ierr)
  call MatSetOption(edge_conn_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr);CHKERRQ(ierr)
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

    call UGridVoronoiTetGeom(v1,v2,v3,v4,cc,c123,c124,c134,c234,e12,e13,e14,e23,e24,e34,area,vol)

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
      call MatSetValue(edge_area_mat,iedge(j,1),iedge(j,2),area(j), &
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
        edge_area_local(local_edges) = abs(area_tmp(1))
        if (area_tmp(1) < 0.d0) then
          write(option%io_buffer, '("WARNING: Negative net dual-face area ",ES12.4, &
            &" for primal edge (",I0,",",I0,"). Mesh may not be Delaunay.")') &
            area_tmp(1), row+1, cols(j)+1
          call PrintMsg(option)
        endif
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

  open(newunit=fid,file=trim(uge_filename),status='replace',action='write',iostat=ios)
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

subroutine UGridVoronoiTetGeom(v1,v2,v3,v4,cc,c123,c124,c134,c234,e12,e13,e14,e23,e24,e34,area,vol)
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
  PetscReal, intent(out) :: area(6), vol(4)

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

subroutine UGridVoronoiWriteVTUSubTetsPiece(base, irank, num_pts, elem_conn_local, nlocal, xyz, option)
  use Option_module
  implicit none
  character(len=*), intent(in) :: base
  PetscInt, intent(in) :: irank, num_pts, nlocal
  PetscInt, intent(in) :: elem_conn_local(4, max(nlocal,1))
  PetscReal, intent(in) :: xyz(3, num_pts)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: fn
  integer :: fid, ios
  PetscInt :: npt, ncells
  PetscInt :: it, basep, c
  PetscInt :: v(4)

  PetscReal :: cc(3), c123(3), c124(3), c134(3), c234(3)
  PetscReal :: e12(3), e13(3), e14(3), e23(3), e24(3), e34(3)
  PetscReal :: area(6), vol(4)

  if (nlocal <= 0) return

  write(fn,'(A,"_",I4.4,".vtu")') trim(base), irank
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing VTU piece: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  npt    = 15*nlocal
  ncells = 24*nlocal

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <UnstructuredGrid>'
  write(fid,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', npt, '" NumberOfCells="', ncells, '">'

  ! Points
  write(fid,'(A)') '      <Points>'
  write(fid,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do it = 1, nlocal
    v = elem_conn_local(:,it)
    call UGridVoronoiTetGeom(xyz(:,v(1)), xyz(:,v(2)), xyz(:,v(3)), xyz(:,v(4)), &
                             cc, c123, c124, c134, c234, &
                             e12, e13, e14, e23, e24, e34, area, vol)

    call WriteXYZ(fid, xyz(:,v(1)))
    call WriteXYZ(fid, xyz(:,v(2)))
    call WriteXYZ(fid, xyz(:,v(3)))
    call WriteXYZ(fid, xyz(:,v(4)))
    call WriteXYZ(fid, cc)
    call WriteXYZ(fid, c123)
    call WriteXYZ(fid, c124)
    call WriteXYZ(fid, c134)
    call WriteXYZ(fid, c234)
    call WriteXYZ(fid, e12)
    call WriteXYZ(fid, e13)
    call WriteXYZ(fid, e14)
    call WriteXYZ(fid, e23)
    call WriteXYZ(fid, e24)
    call WriteXYZ(fid, e34)
  enddo
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </Points>'

  ! Cells connectivity/offsets/types
  write(fid,'(A)') '      <Cells>'

  ! connectivity: 4 ints per cell
  write(fid,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do it = 1, nlocal
    basep = (it-1)*15
    call Emit24SubTetsConnectivity(fid, basep)
  enddo
  write(fid,'(A)') '        </DataArray>'

  ! offsets: cumulative number of points per cell (4,8,12,...)
  write(fid,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  c = 0
  do it = 1, ncells
    c = c + 4
    write(fid,'(I0,1x)', advance='no') c
    if (mod(it,16)==0) write(fid,*)
  enddo
  if (mod(ncells,16)/=0) write(fid,*)
  write(fid,'(A)') '        </DataArray>'

  ! types: 10 for tetra
  write(fid,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do it = 1, ncells
    write(fid,'(I0,1x)', advance='no') 10
    if (mod(it,32)==0) write(fid,*)
  enddo
  if (mod(ncells,32)/=0) write(fid,*)
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '      </Cells>'

  ! CellData: vertex_id
  write(fid,'(A)') '      <CellData Scalars="vertex_id">'
  write(fid,'(A)') '        <DataArray type="Int32" Name="vertex_id" NumberOfComponents="1" format="ascii">'
  do it = 1, nlocal
    v = elem_conn_local(:,it)
    call Write6(fid, v(1))
    call Write6(fid, v(2))
    call Write6(fid, v(3))
    call Write6(fid, v(4))
  enddo
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </CellData>'

  write(fid,'(A)') '    </Piece>'
  write(fid,'(A)') '  </UnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

contains
  subroutine WriteXYZ(fid, p)
    integer, intent(in) :: fid
    PetscReal, intent(in) :: p(3)
    write(fid,'(3ES24.16,1x)') p(1), p(2), p(3)
  end subroutine

  subroutine Write6(fid, val)
    integer, intent(in) :: fid
    PetscInt, intent(in) :: val
    integer :: k
    do k=1,6
      write(fid,'(I0,1x)', advance='no') val
    enddo
    write(fid,*)
  end subroutine

  subroutine EmitTet4(fid, a,b,c,d)
    integer, intent(in) :: fid
    integer, intent(in) :: a,b,c,d
    write(fid,'(4(I0,1x))') a,b,c,d
  end subroutine

  subroutine Emit24SubTetsConnectivity(fid, basep)
    integer, intent(in) :: fid
    integer, intent(in) :: basep
    ! local ids in this tet-block:
    ! 0..3 v1..v4, 4 cc, 5..8 c123,c124,c134,c234, 9..14 e12,e13,e14,e23,e24,e34

    ! v1 group
    call EmitTet4(fid, basep+0, basep+9,  basep+6,  basep+4)
    call EmitTet4(fid, basep+0, basep+6,  basep+11, basep+4)
    call EmitTet4(fid, basep+0, basep+7,  basep+10, basep+4)
    call EmitTet4(fid, basep+0, basep+11, basep+7,  basep+4)
    call EmitTet4(fid, basep+0, basep+5,  basep+9,  basep+4)
    call EmitTet4(fid, basep+0, basep+10, basep+5,  basep+4)

    ! v2 group
    call EmitTet4(fid, basep+1, basep+6,  basep+9,  basep+4)
    call EmitTet4(fid, basep+1, basep+8,  basep+13, basep+4)
    call EmitTet4(fid, basep+1, basep+12, basep+8,  basep+4)
    call EmitTet4(fid, basep+1, basep+9,  basep+5,  basep+4)
    call EmitTet4(fid, basep+1, basep+5,  basep+12, basep+4)
    call EmitTet4(fid, basep+1, basep+13, basep+6,  basep+4)

    ! v3 group
    call EmitTet4(fid, basep+2, basep+7,  basep+14, basep+4)
    call EmitTet4(fid, basep+2, basep+14, basep+8,  basep+4)
    call EmitTet4(fid, basep+2, basep+8,  basep+12, basep+4)
    call EmitTet4(fid, basep+2, basep+10, basep+7,  basep+4)
    call EmitTet4(fid, basep+2, basep+12, basep+5,  basep+4)
    call EmitTet4(fid, basep+2, basep+5,  basep+10, basep+4)

    ! v4 group
    call EmitTet4(fid, basep+3, basep+7,  basep+11, basep+4)
    call EmitTet4(fid, basep+3, basep+11, basep+6,  basep+4)
    call EmitTet4(fid, basep+3, basep+6,  basep+13, basep+4)
    call EmitTet4(fid, basep+3, basep+13, basep+8,  basep+4)
    call EmitTet4(fid, basep+3, basep+8,  basep+14, basep+4)
    call EmitTet4(fid, basep+3, basep+14, basep+7,  basep+4)
  end subroutine

end subroutine UGridVoronoiWriteVTUSubTetsPiece

subroutine UGridVoronoiWritePVTUSubTets(base, option)
  use Option_module
  implicit none
  character(len=*), intent(in) :: base
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: fn, piece
  integer :: fid, ios
  PetscInt :: r, nprocs

  if (.not.OptionIsIORank(option)) return
  nprocs = option%comm%size

  fn = trim(base)//'.pvtu'
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing PVTU: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <PUnstructuredGrid GhostLevel="0">'
  write(fid,'(A)') '    <PPoints>'
  write(fid,'(A)') '      <PDataArray type="Float64" NumberOfComponents="3"/>'
  write(fid,'(A)') '    </PPoints>'
  write(fid,'(A)') '    <PCellData Scalars="vertex_id">'
  write(fid,'(A)') '      <PDataArray type="Int32" Name="vertex_id" NumberOfComponents="1"/>'
  write(fid,'(A)') '    </PCellData>'

  do r = 0, nprocs-1
    write(piece,'(A,"_",I4.4,".vtu")') trim(base), r
    write(fid,'(A,A,A)') '    <Piece Source="', trim(piece), '"/>'
  enddo

  write(fid,'(A)') '  </PUnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

  option%io_buffer = ' --> wrote PVTU: ' // trim(fn)
  call PrintMsg(option)
end subroutine UGridVoronoiWritePVTUSubTets

subroutine UGridVoronoiWriteVTUDualFacePatchesPiece(base, irank, num_pts, elem_start, &
                                                     elem_conn_local, nlocal, xyz, option)
  use Option_module
  implicit none
  character(len=*), intent(in) :: base
  PetscInt, intent(in) :: irank, num_pts, elem_start, nlocal
  PetscInt, intent(in) :: elem_conn_local(4, nlocal)
  PetscReal, intent(in) :: xyz(3, num_pts)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: fn
  integer :: fid, ios
  PetscInt :: npt, ncells
  PetscInt :: it, basep, tet_id
  PetscInt :: v(4)

  PetscReal :: cc(3), c123(3), c124(3), c134(3), c234(3)
  PetscReal :: e12(3), e13(3), e14(3), e23(3), e24(3), e34(3)
  PetscReal :: area(6), vol(4)

  if (nlocal <= 0) return

  write(fn,'(A,"_",I4.4,".vtu")') trim(base), irank
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing VTU piece: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  ! 11 points per tet (duplicated), 12 triangles per tet
  npt    = 11*nlocal
  ncells = 12*nlocal

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <UnstructuredGrid>'
  write(fid,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', npt, '" NumberOfCells="', ncells, '">'

  ! ---------------- Points ----------------
  write(fid,'(A)') '      <Points>'
  write(fid,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do it = 1, nlocal
    v = elem_conn_local(:,it)

    call UGridVoronoiTetGeom(xyz(:,v(1)), xyz(:,v(2)), xyz(:,v(3)), xyz(:,v(4)), &
                             cc, c123, c124, c134, c234, &
                             e12, e13, e14, e23, e24, e34, area, vol)

    ! point order in this 11-point block:
    ! 0 cc
    ! 1 c123, 2 c124, 3 c134, 4 c234
    ! 5 e12, 6 e13, 7 e14, 8 e23, 9 e24, 10 e34
    call WriteXYZ(fid, cc)
    call WriteXYZ(fid, c123)
    call WriteXYZ(fid, c124)
    call WriteXYZ(fid, c134)
    call WriteXYZ(fid, c234)
    call WriteXYZ(fid, e12)
    call WriteXYZ(fid, e13)
    call WriteXYZ(fid, e14)
    call WriteXYZ(fid, e23)
    call WriteXYZ(fid, e24)
    call WriteXYZ(fid, e34)
  enddo
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </Points>'

  ! ---------------- Cells ----------------
  write(fid,'(A)') '      <Cells>'

  ! connectivity: 3 ints per triangle
  write(fid,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do it = 1, nlocal
    basep = (it-1)*11
    call Emit12DualTrisConnectivity(fid, basep)
  enddo
  write(fid,'(A)') '        </DataArray>'

  ! offsets: cumulative 3,6,9,...,3*ncells
  write(fid,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  call WriteOffsets(fid, ncells, 3)
  write(fid,'(A)') '        </DataArray>'

  ! types: 5 for VTK_TRIANGLE
  write(fid,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  call WriteTypes(fid, ncells, 5)
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '      </Cells>'

  ! ---------------- CellData ----------------
  write(fid,'(A)') '      <CellData Scalars="edge_a">'
  write(fid,'(A)') '        <DataArray type="Int32" Name="parent_tet_id" NumberOfComponents="1" format="ascii">'
  do it = 1, nlocal
    tet_id = elem_start + it - 1
    call WriteRepeat(fid, tet_id, 12)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Int32" Name="edge_a" NumberOfComponents="1" format="ascii">'
  do it = 1, nlocal
    v = elem_conn_local(:,it)
    ! 12 triangles correspond to edges in this order:
    ! (1,2)x2, (2,3)x2, (1,3)x2, (1,4)x2, (2,4)x2, (3,4)x2
    call WriteRepeat(fid, v(1), 2)  ! edge 1-2
    call WriteRepeat(fid, v(2), 2)  ! edge 2-3
    call WriteRepeat(fid, v(1), 2)  ! edge 1-3
    call WriteRepeat(fid, v(1), 2)  ! edge 1-4
    call WriteRepeat(fid, v(2), 2)  ! edge 2-4
    call WriteRepeat(fid, v(3), 2)  ! edge 3-4
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Int32" Name="edge_b" NumberOfComponents="1" format="ascii">'
  do it = 1, nlocal
    v = elem_conn_local(:,it)
    call WriteRepeat(fid, v(2), 2)  ! edge 1-2
    call WriteRepeat(fid, v(3), 2)  ! edge 2-3
    call WriteRepeat(fid, v(3), 2)  ! edge 1-3
    call WriteRepeat(fid, v(4), 2)  ! edge 1-4
    call WriteRepeat(fid, v(4), 2)  ! edge 2-4
    call WriteRepeat(fid, v(4), 2)  ! edge 3-4
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '      </CellData>'

  write(fid,'(A)') '    </Piece>'
  write(fid,'(A)') '  </UnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

contains
  subroutine WriteXYZ(fid, p)
    integer, intent(in) :: fid
    PetscReal, intent(in) :: p(3)
    write(fid,'(3ES24.16,1x)') p(1), p(2), p(3)
  end subroutine

  subroutine WriteOffsets(fid, nc, stride)
    integer, intent(in) :: fid
    PetscInt, intent(in) :: nc, stride
    PetscInt :: i, off
    off = 0
    do i = 1, nc
      off = off + stride
      write(fid,'(I0,1x)', advance='no') off
      if (mod(i,16)==0) write(fid,*)
    enddo
    if (mod(nc,16)/=0) write(fid,*)
  end subroutine

  subroutine WriteTypes(fid, nc, t)
    integer, intent(in) :: fid
    PetscInt, intent(in) :: nc
    integer, intent(in) :: t
    PetscInt :: i
    do i = 1, nc
      write(fid,'(I0,1x)', advance='no') t
      if (mod(i,32)==0) write(fid,*)
    enddo
    if (mod(nc,32)/=0) write(fid,*)
  end subroutine

  subroutine WriteRepeat(fid, val, nrep)
    integer, intent(in) :: fid
    PetscInt, intent(in) :: val, nrep
    PetscInt :: i
    do i = 1, nrep
      write(fid,'(I0,1x)', advance='no') val
    enddo
    write(fid,*)
  end subroutine

  subroutine EmitTri(fid, a,b,c)
    integer, intent(in) :: fid
    integer, intent(in) :: a,b,c
    write(fid,'(3(I0,1x))') a,b,c
  end subroutine

  subroutine Emit12DualTrisConnectivity(fid, b)
    integer, intent(in) :: fid
    integer, intent(in) :: b
    ! indices in 11-point block:
    ! 0 cc
    ! 1 c123, 2 c124, 3 c134, 4 c234
    ! 5 e12, 6 e13, 7 e14, 8 e23, 9 e24, 10 e34

    ! edge 1-2: (e12,c124,cc) and (e12,cc,c123)
    call EmitTri(fid, b+5,  b+2,  b+0)
    call EmitTri(fid, b+5,  b+0,  b+1)

    ! edge 2-3: (e23,c234,cc) and (e23,cc,c123)
    call EmitTri(fid, b+8,  b+4,  b+0)
    call EmitTri(fid, b+8,  b+0,  b+1)

    ! edge 1-3: (e13,c123,cc) and (e13,cc,c134)
    call EmitTri(fid, b+6,  b+1,  b+0)
    call EmitTri(fid, b+6,  b+0,  b+3)

    ! edge 1-4: (e14,c134,cc) and (e14,cc,c124)
    call EmitTri(fid, b+7,  b+3,  b+0)
    call EmitTri(fid, b+7,  b+0,  b+2)

    ! edge 2-4: (e24,c124,cc) and (e24,cc,c234)
    call EmitTri(fid, b+9,  b+2,  b+0)
    call EmitTri(fid, b+9,  b+0,  b+4)

    ! edge 3-4: (e34,c234,cc) and (e34,cc,c134)
    call EmitTri(fid, b+10, b+4,  b+0)
    call EmitTri(fid, b+10, b+0,  b+3)
  end subroutine

end subroutine UGridVoronoiWriteVTUDualFacePatchesPiece

subroutine UGridVoronoiWritePVTUDualFacePatches(base, option)
  use Option_module
  implicit none
  character(len=*), intent(in) :: base
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: fn, piece
  integer :: fid, ios
  PetscInt :: r, nprocs

  if (.not.OptionIsIORank(option)) return
  nprocs = option%comm%size

  fn = trim(base)//'.pvtu'
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing PVTU: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <PUnstructuredGrid GhostLevel="0">'
  write(fid,'(A)') '    <PPoints>'
  write(fid,'(A)') '      <PDataArray type="Float64" NumberOfComponents="3"/>'
  write(fid,'(A)') '    </PPoints>'
  write(fid,'(A)') '    <PCellData Scalars="edge_a">'
  write(fid,'(A)') '      <PDataArray type="Int32" Name="parent_tet_id" NumberOfComponents="1"/>'
  write(fid,'(A)') '      <PDataArray type="Int32" Name="edge_a" NumberOfComponents="1"/>'
  write(fid,'(A)') '      <PDataArray type="Int32" Name="edge_b" NumberOfComponents="1"/>'
  write(fid,'(A)') '    </PCellData>'

  do r = 0, nprocs-1
    write(piece,'(A,"_",I4.4,".vtu")') trim(base), r
    write(fid,'(A,A,A)') '    <Piece Source="', trim(piece), '"/>'
  enddo

  write(fid,'(A)') '  </PUnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

  option%io_buffer = ' --> wrote PVTU: ' // trim(fn)
  call PrintMsg(option)
end subroutine UGridVoronoiWritePVTUDualFacePatches

subroutine UGridVoronoiWriteVTUTetMarkersPiece(base, irank, num_pts, elem_start, &
                                                elem_conn_local, nlocal, xyz, option)
  use Option_module
  implicit none
  character(len=*), intent(in) :: base
  PetscInt, intent(in) :: irank, num_pts, elem_start, nlocal
  PetscInt, intent(in) :: elem_conn_local(4, nlocal)
  PetscReal, intent(in) :: xyz(3, num_pts)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: fn
  integer :: fid, ios
  PetscInt :: it, tet_id
  PetscInt :: npt
  PetscInt :: v(4)

  PetscReal :: cc(3), c123(3), c124(3), c134(3), c234(3)
  PetscReal :: e12(3), e13(3), e14(3), e23(3), e24(3), e34(3)
  PetscReal :: area(6), vol(4)

  if (nlocal <= 0) return

  ! 1 cc + 4 face cc + 6 edge midpoints = 11 points per primal tet
  npt = 11*nlocal

  write(fn,'(A,"_",I4.4,".vtu")') trim(base), irank
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing VTU markers piece: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <UnstructuredGrid>'
  write(fid,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', npt, '" NumberOfCells="0">'

  ! Points
  write(fid,'(A)') '      <Points>'
  write(fid,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do it = 1, nlocal
    v = elem_conn_local(:,it)
    call UGridVoronoiTetGeom(xyz(:,v(1)), xyz(:,v(2)), xyz(:,v(3)), xyz(:,v(4)), &
                             cc, c123, c124, c134, c234, &
                             e12, e13, e14, e23, e24, e34, area, vol)

    ! order: cc, c123,c124,c134,c234, e12,e13,e14,e23,e24,e34
    call WriteXYZ(fid, cc)
    call WriteXYZ(fid, c123)
    call WriteXYZ(fid, c124)
    call WriteXYZ(fid, c134)
    call WriteXYZ(fid, c234)
    call WriteXYZ(fid, e12)
    call WriteXYZ(fid, e13)
    call WriteXYZ(fid, e14)
    call WriteXYZ(fid, e23)
    call WriteXYZ(fid, e24)
    call WriteXYZ(fid, e34)
  enddo
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </Points>'

  ! Empty Cells section required by some readers
  write(fid,'(A)') '      <Cells>'
  write(fid,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii"></DataArray>'
  write(fid,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii"></DataArray>'
  write(fid,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii"></DataArray>'
  write(fid,'(A)') '      </Cells>'

  ! PointData: point_type and parent_tet_id
  write(fid,'(A)') '      <PointData Scalars="point_type">'

  write(fid,'(A)') '        <DataArray type="Int32" Name="point_type" NumberOfComponents="1" format="ascii">'
  do it = 1, nlocal
    ! cc
    call WriteRepeat(fid, 1, 1)
    ! 4 face circumcenters
    call WriteRepeat(fid, 2, 4)
    ! 6 edge midpoints
    call WriteRepeat(fid, 3, 6)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Int32" Name="parent_tet_id" NumberOfComponents="1" format="ascii">'
  do it = 1, nlocal
    tet_id = elem_start + it - 1
    call WriteRepeat(fid, tet_id, 11)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '      </PointData>'

  write(fid,'(A)') '    </Piece>'
  write(fid,'(A)') '  </UnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

contains
  subroutine WriteXYZ(fid, p)
    integer, intent(in) :: fid
    PetscReal, intent(in) :: p(3)
    write(fid,'(3ES24.16,1x)') p(1), p(2), p(3)
  end subroutine

  subroutine WriteRepeat(fid, val, nrep)
    integer, intent(in) :: fid
    PetscInt, intent(in) :: val, nrep
    PetscInt :: i
    do i = 1, nrep
      write(fid,'(I0,1x)', advance='no') val
    enddo
    write(fid,*)
  end subroutine

end subroutine UGridVoronoiWriteVTUTetMarkersPiece

subroutine UGridVoronoiWritePVTUTetMarkers(base, option)
  use Option_module
  implicit none
  character(len=*), intent(in) :: base
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: fn, piece
  integer :: fid, ios
  PetscInt :: r, nprocs

  if (.not.OptionIsIORank(option)) return
  nprocs = option%comm%size

  fn = trim(base)//'.pvtu'
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing PVTU markers: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <PUnstructuredGrid GhostLevel="0">'
  write(fid,'(A)') '    <PPoints>'
  write(fid,'(A)') '      <PDataArray type="Float64" NumberOfComponents="3"/>'
  write(fid,'(A)') '    </PPoints>'
  write(fid,'(A)') '    <PPointData Scalars="point_type">'
  write(fid,'(A)') '      <PDataArray type="Int32" Name="point_type" NumberOfComponents="1"/>'
  write(fid,'(A)') '      <PDataArray type="Int32" Name="parent_tet_id" NumberOfComponents="1"/>'
  write(fid,'(A)') '    </PPointData>'

  do r = 0, nprocs-1
    write(piece,'(A,"_",I4.4,".vtu")') trim(base), r
    write(fid,'(A,A,A)') '    <Piece Source="', trim(piece), '"/>'
  enddo

  write(fid,'(A)') '  </PUnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

  option%io_buffer = ' --> wrote PVTU markers: ' // trim(fn)
  call PrintMsg(option)
end subroutine UGridVoronoiWritePVTUTetMarkers

! ************************************************************************** !

subroutine UGridVoronoiWriteVTUDualEdgePolygons(base, num_pts, num_elem, elem_conn, xyz, option)
  !
  ! Build stitched Voronoi-face polygons (one polygon per primal edge) from
  ! the full primal mesh and write a single VTU file.  Runs on IO rank only.
  !
  ! Author: Satish Karra
  ! Date: 02/12/2026
  !

  use Option_module
  implicit none

  character(len=*), intent(in) :: base
  PetscInt, intent(in) :: num_pts, num_elem
  PetscInt, intent(in) :: elem_conn(4,num_elem)
  PetscReal, intent(in) :: xyz(3, num_pts)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: fn
  integer :: fid, ios
  PetscInt :: max_occ
  PetscInt :: it, eid, c, v
  PetscInt :: nocc, nedge
  PetscInt :: nverts, point_cursor
  PetscInt :: num_points, num_cells, offset
  PetscReal :: xmin(3), xmax(3), diag, margin, area_max
  PetscInt, allocatable :: occ_edge(:), occ_tet(:)
  PetscReal, allocatable :: occ_ct(:,:), occ_e(:,:), occ_f1(:,:), occ_f2(:,:)
  PetscInt, allocatable :: edge_i(:), edge_j(:)
  PetscInt, allocatable :: nverts_edge(:)
  PetscReal, allocatable :: area_edge(:), centroid_edge(:,:)
  PetscReal, allocatable :: poly(:,:)
  PetscReal :: cc(3), c123(3), c124(3), c134(3), c234(3)
  PetscReal :: e12(3), e13(3), e14(3), e23(3), e24(3), e34(3)
  PetscReal :: area6(6), vol(4)
  PetscInt :: tet_id
  PetscInt :: va, vb, vc, vd

  max_occ = max(6*num_elem,1)
  allocate(occ_edge(max_occ))
  allocate(occ_tet(max_occ))
  allocate(occ_ct(3,max_occ))
  allocate(occ_e(3,max_occ))
  allocate(occ_f1(3,max_occ))
  allocate(occ_f2(3,max_occ))
  allocate(edge_i(max_occ))
  allocate(edge_j(max_occ))
  allocate(nverts_edge(max_occ))
  allocate(area_edge(max_occ))
  allocate(centroid_edge(3,max_occ))
  allocate(poly(3,max_occ))

  occ_edge = 0
  occ_tet = 0
  occ_ct = 0.d0
  occ_e = 0.d0
  occ_f1 = 0.d0
  occ_f2 = 0.d0
  edge_i = 0
  edge_j = 0
  nverts_edge = 0
  area_edge = 0.d0
  centroid_edge = 0.d0
  poly = 0.d0
  nocc = 0
  nedge = 0
  xmin(1) = minval(xyz(1,:)); xmax(1) = maxval(xyz(1,:))
  xmin(2) = minval(xyz(2,:)); xmax(2) = maxval(xyz(2,:))
  xmin(3) = minval(xyz(3,:)); xmax(3) = maxval(xyz(3,:))
  diag = sqrt( (xmax(1)-xmin(1))**2 + (xmax(2)-xmin(2))**2 + (xmax(3)-xmin(3))**2 )
  ! Visualization guardrails: reject faces that escape far outside the
  ! primal domain or have unrealistically large area for this mesh scale.
  margin = 1.d0*max(diag,1.d-12)
  area_max = 10.d0*max(diag*diag,1.d-18)

  ! Collect one circumcenter sample per (tet,edge).
  do it = 1, num_elem
    va = elem_conn(1,it)
    vb = elem_conn(2,it)
    vc = elem_conn(3,it)
    vd = elem_conn(4,it)
    tet_id = it

    call UGridVoronoiTetGeom(xyz(:,va), xyz(:,vb), xyz(:,vc), xyz(:,vd), &
                             cc, c123, c124, c134, c234, &
                             e12, e13, e14, e23, e24, e34, area6, vol)

    call AddOccurrence(min(va,vb), max(va,vb), cc, e12, c123, c124, tet_id)
    call AddOccurrence(min(vb,vc), max(vb,vc), cc, e23, c123, c234, tet_id)
    call AddOccurrence(min(va,vc), max(va,vc), cc, e13, c123, c134, tet_id)
    call AddOccurrence(min(va,vd), max(va,vd), cc, e14, c134, c124, tet_id)
    call AddOccurrence(min(vb,vd), max(vb,vd), cc, e24, c124, c234, tet_id)
    call AddOccurrence(min(vc,vd), max(vc,vd), cc, e34, c234, c134, tet_id)
  enddo

  ! First pass: count all stitched polygons and points.
  num_cells = 0
  num_points = 0
  do eid = 1, nedge
    call BuildEdgePolygon(eid, nverts, poly, area_edge(eid), centroid_edge(:,eid))
    if (.not.PolygonIsReasonable(nverts,poly,area_edge(eid))) nverts = 0
    if (nverts >= 3) then
      nverts_edge(eid) = nverts
      num_cells = num_cells + 1
      num_points = num_points + nverts
    else
      nverts_edge(eid) = 0
    endif
  enddo

  fn = trim(base)//'.vtu'
  open(newunit=fid,file=trim(fn),status='replace',action='write',iostat=ios)
  if (ios /= 0) then
    option%io_buffer = 'Failed writing Voronoi polygons VTU: ' // trim(fn)
    call PrintErrMsg(option)
  endif

  write(fid,'(A)') '<?xml version="1.0"?>'
  write(fid,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(fid,'(A)') '  <UnstructuredGrid>'
  write(fid,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', num_points, '" NumberOfCells="', num_cells, '">'

  ! Points
  write(fid,'(A)') '      <Points>'
  write(fid,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    call BuildEdgePolygon(eid, nverts, poly, area_edge(eid), centroid_edge(:,eid))
    do v = 1, nverts
      write(fid,'(3ES24.16,1x)') poly(1,v), poly(2,v), poly(3,v)
    enddo
  enddo
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </Points>'

  ! Cells
  write(fid,'(A)') '      <Cells>'
  write(fid,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  point_cursor = 0
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    do v = 1, nverts_edge(eid)
      write(fid,'(I0,1x)', advance='no') point_cursor + (v-1)
    enddo
    write(fid,*)
    point_cursor = point_cursor + nverts_edge(eid)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  c = 0
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    offset = offset + nverts_edge(eid)
    c = c + 1
    write(fid,'(I0,1x)', advance='no') offset
    if (mod(c,16) == 0) write(fid,*)
  enddo
  if (mod(c,16) /= 0) write(fid,*)
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  c = 0
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    c = c + 1
    write(fid,'(I0,1x)', advance='no') 7
    if (mod(c,32) == 0) write(fid,*)
  enddo
  if (mod(c,32) /= 0) write(fid,*)
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </Cells>'

  ! CellData
  write(fid,'(A)') '      <CellData Scalars="edge_a">'
  write(fid,'(A)') '        <DataArray type="Int32" Name="edge_a" NumberOfComponents="1" format="ascii">'
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    write(fid,'(I0,1x)') edge_i(eid)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Int32" Name="edge_b" NumberOfComponents="1" format="ascii">'
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    write(fid,'(I0,1x)') edge_j(eid)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Int32" Name="nverts" NumberOfComponents="1" format="ascii">'
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    write(fid,'(I0,1x)') nverts_edge(eid)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Float64" Name="area" NumberOfComponents="1" format="ascii">'
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    write(fid,'(ES24.16,1x)') area_edge(eid)
  enddo
  write(fid,'(A)') '        </DataArray>'

  write(fid,'(A)') '        <DataArray type="Float64" Name="centroid" NumberOfComponents="3" format="ascii">'
  do eid = 1, nedge

    if (nverts_edge(eid) <= 0) cycle
    write(fid,'(3ES24.16,1x)') centroid_edge(1,eid), centroid_edge(2,eid), centroid_edge(3,eid)
  enddo
  write(fid,'(A)') '        </DataArray>'
  write(fid,'(A)') '      </CellData>'

  write(fid,'(A)') '    </Piece>'
  write(fid,'(A)') '  </UnstructuredGrid>'
  write(fid,'(A)') '</VTKFile>'
  close(fid)

  option%io_buffer = ' --> wrote Voronoi polygons VTU: ' // trim(fn)
  call PrintMsg(option)

  deallocate(occ_edge,occ_tet,occ_ct,occ_e,occ_f1,occ_f2,edge_i,edge_j,nverts_edge,area_edge,centroid_edge,poly)

contains

  subroutine AddOccurrence(a,b,cc_in,e_in,f1_in,f2_in,tet_in)
    PetscInt, intent(in) :: a,b,tet_in
    PetscReal, intent(in) :: cc_in(3), e_in(3), f1_in(3), f2_in(3)
    PetscInt :: k, eid_local

    eid_local = 0
    do k = 1, nedge
      if (edge_i(k) == a .and. edge_j(k) == b) then
        eid_local = k
        exit
      endif
    enddo
    if (eid_local == 0) then
      nedge = nedge + 1
      if (nedge > max_occ) then
        option%io_buffer = 'Internal error: edge list overflow while building Voronoi polygons.'
        call PrintErrMsg(option)
      endif
      edge_i(nedge) = a
      edge_j(nedge) = b
      eid_local = nedge
    endif

    nocc = nocc + 1
    if (nocc > max_occ) then
      option%io_buffer = 'Internal error: occurrence list overflow while building Voronoi polygons.'
      call PrintErrMsg(option)
    endif
    occ_edge(nocc) = eid_local
    occ_tet(nocc) = tet_in
    occ_ct(:,nocc) = cc_in
    occ_e(:,nocc) = e_in
    occ_f1(:,nocc) = f1_in
    occ_f2(:,nocc) = f2_in
  end subroutine AddOccurrence

  subroutine BuildEdgePolygon(eid_in, nout, poly_out, area_out, centroid_out)
    PetscInt, intent(in) :: eid_in
    PetscInt, intent(out) :: nout
    PetscReal, intent(out) :: poly_out(3,max_occ)
    PetscReal, intent(out) :: area_out
    PetscReal, intent(out) :: centroid_out(3)

    PetscReal, parameter :: tol2 = 1.d-20
    PetscReal :: uniq(3,max_occ), p(3), d2
    PetscReal :: xi(3), xj(3), u(3), e1(3), e2(3), ref(3), ctr(3), q(3), cr(3)
    PetscReal :: angle(max_occ), tmpa
    PetscInt :: idx(max_occ), tmpi
    PetscInt :: i1, i2, j1, m, nuniq
    PetscReal :: un

    nout = 0
    area_out = 0.d0
    centroid_out = 0.d0
    poly_out = 0.d0
    uniq = 0.d0
    nuniq = 0

    do i1 = 1, nocc
      if (occ_edge(i1) /= eid_in) cycle
      p = occ_ct(:,i1)
      do j1 = 1, nuniq
        d2 = (uniq(1,j1)-p(1))*(uniq(1,j1)-p(1)) + &
             (uniq(2,j1)-p(2))*(uniq(2,j1)-p(2)) + &
             (uniq(3,j1)-p(3))*(uniq(3,j1)-p(3))
        if (d2 <= tol2) exit
      enddo
      if (j1 > nuniq) then
        nuniq = nuniq + 1
        uniq(:,nuniq) = p
      endif
    enddo

    if (nuniq < 3) then
      ! Boundary/low-connectivity fallback: use one local quad patch
      ! [edge-midpoint, face-cc-1, tet-cc, face-cc-2] so VTU is non-empty.
      do i1 = 1, nocc
        if (occ_edge(i1) == eid_in) then
          poly_out(:,1) = occ_e(:,i1)
          poly_out(:,2) = occ_f1(:,i1)
          poly_out(:,3) = occ_ct(:,i1)
          poly_out(:,4) = occ_f2(:,i1)
          nout = 4
          centroid_out = 0.25d0*(poly_out(:,1)+poly_out(:,2)+poly_out(:,3)+poly_out(:,4))
          call UGridVoronoiCross(cr, poly_out(:,2)-poly_out(:,1), poly_out(:,3)-poly_out(:,1))
          area_out = 0.5d0*UGridVoronoiNorm(cr)
          call UGridVoronoiCross(cr, poly_out(:,3)-poly_out(:,1), poly_out(:,4)-poly_out(:,1))
          area_out = area_out + 0.5d0*UGridVoronoiNorm(cr)
          return
        endif
      enddo
      nout = 0
      return
    endif

    xi = xyz(:,edge_i(eid_in))
    xj = xyz(:,edge_j(eid_in))
    u = xj - xi
    un = UGridVoronoiNorm(u)
    if (un <= 1.d-30) then
      return
    endif
    u = u/un

    ref = (/1.d0,0.d0,0.d0/)
    if (abs(dot_product(u,ref)) > 0.9d0) ref = (/0.d0,1.d0,0.d0/)
    call UGridVoronoiCross(e1,u,ref)
    un = UGridVoronoiNorm(e1)
    if (un <= 1.d-30) then
      ref = (/0.d0,0.d0,1.d0/)
      call UGridVoronoiCross(e1,u,ref)
      un = UGridVoronoiNorm(e1)
      if (un <= 1.d-30) return
    endif
    e1 = e1/un
    call UGridVoronoiCross(e2,u,e1)
    un = UGridVoronoiNorm(e2)
    if (un <= 1.d-30) return
    e2 = e2/un

    ctr = 0.d0
    do m = 1, nuniq
      ctr = ctr + uniq(:,m)
    enddo
    ctr = ctr/dble(nuniq)

    do m = 1, nuniq
      q = uniq(:,m) - ctr
      angle(m) = atan2(dot_product(q,e2), dot_product(q,e1))
      idx(m) = m
    enddo

    ! insertion sort by angle
    do i1 = 2, nuniq
      tmpa = angle(i1)
      tmpi = idx(i1)
      j1 = i1 - 1
      do while (j1 >= 1 .and. angle(j1) > tmpa)
        angle(j1+1) = angle(j1)
        idx(j1+1) = idx(j1)
        j1 = j1 - 1
      enddo
      angle(j1+1) = tmpa
      idx(j1+1) = tmpi
    enddo

    do m = 1, nuniq
      poly_out(:,m) = uniq(:,idx(m))
    enddo
    nout = nuniq

    centroid_out = 0.d0
    do m = 1, nout
      centroid_out = centroid_out + poly_out(:,m)
    enddo
    centroid_out = centroid_out/dble(nout)

    area_out = 0.d0
    do m = 1, nout
      i2 = m + 1
      if (i2 > nout) i2 = 1
      call UGridVoronoiCross(cr, poly_out(:,m)-centroid_out, poly_out(:,i2)-centroid_out)
      area_out = area_out + 0.5d0*abs(dot_product(cr,u))
    enddo
  end subroutine BuildEdgePolygon

  logical function PolygonIsReasonable(nv, polyv, area_in)
    PetscInt, intent(in) :: nv
    PetscReal, intent(in) :: polyv(3,max_occ)
    PetscReal, intent(in) :: area_in
    PetscInt :: k
    PolygonIsReasonable = .false.
    if (nv < 3) return
    if (area_in <= 1.d-30 .or. area_in > area_max) return
    if (centroid_outside_bounds(polv=polyv, nv=nv)) return
    do k = 1, nv
      if (polyv(1,k) < xmin(1)-margin .or. polyv(1,k) > xmax(1)+margin) return
      if (polyv(2,k) < xmin(2)-margin .or. polyv(2,k) > xmax(2)+margin) return
      if (polyv(3,k) < xmin(3)-margin .or. polyv(3,k) > xmax(3)+margin) return
    enddo
    PolygonIsReasonable = .true.
  end function PolygonIsReasonable

  logical function centroid_outside_bounds(polv, nv)
    PetscReal, intent(in) :: polv(3,max_occ)
    PetscInt, intent(in) :: nv
    PetscReal :: c(3)
    PetscInt :: k
    centroid_outside_bounds = .true.
    c = 0.d0
    do k = 1, nv
      c = c + polv(:,k)
    enddo
    c = c/dble(nv)
    if (c(1) < xmin(1)-margin .or. c(1) > xmax(1)+margin) return
    if (c(2) < xmin(2)-margin .or. c(2) > xmax(2)+margin) return
    if (c(3) < xmin(3)-margin .or. c(3) > xmax(3)+margin) return
    centroid_outside_bounds = .false.
  end function centroid_outside_bounds

end subroutine UGridVoronoiWriteVTUDualEdgePolygons

! ************************************************************************** !

end module Grid_Unstructured_Voronoi_module
