module Grid_Unstructured_Poly_HDF5_module
  !
  ! Distributed input of polyhedra meshes stored in the HDF5 mixed
  ! XDMF-code schema:
  !
  !   /Domain/Vertices                                 float64 [Nv, 3]
  !   /Domain/CellsPolyhedra/ConnectivityWithXDMFCodes int64   [L]
  !   /Domain/CellsPolyhedra/Offsets                   int64   [Nc + 1]
  !
  ! with attributes SchemaVersion = 1, CellTypeConvention = "XDMF",
  ! NodeIndexBase = 1, and OffsetIndexBase = 1 on /Domain/CellsPolyhedra.
  ! Cell i (natural id i) is ConnectivityWithXDMFCodes(Offsets(i):
  ! Offsets(i+1)-1); see Grid_Unstructured_XDMF_module for the records.
  !
  ! No rank reads or assembles the global mesh. Each rank reads a
  ! contiguous block of offsets, the matching connectivity slice, and a
  ! contiguous block of vertex rows with collective hyperslab reads. Faces
  ! are matched through a distributed directory keyed by the sorted node
  ! ids of each face, the cell graph is partitioned in parallel, cell
  ! records are migrated point to point, the natural-to-PETSc map uses a
  ! memory-scalable AO, and only the ghost cells and vertices a rank needs
  ! are fetched. Geometry is computed from the vertices after
  ! redistribution. Replicated data are limited to scalars and
  ! O(number of ranks) counts/offsets.
  !
  ! Point-to-point exchanges use MPI_Alltoallv with 32-bit counts; an
  ! exchange whose per-rank volume exceeds 2^31-1 entries is reported as an
  ! error (bounded batching is not implemented yet).
  !
#include "petsc/finclude/petscmat.h"
  use petscmat
  use hdf5
  use, intrinsic :: iso_c_binding

  use Geometry_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module
  use Grid_Unstructured_XDMF_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: PETSC_INT_ZERO = 0
  PetscInt64, parameter :: MAX_PETSC_INT = huge(PETSC_INT_ZERO)
  PetscMPIInt, parameter :: MPI_INT_ZERO = 0
  PetscInt64, parameter :: MAX_MPI_COUNT = huge(MPI_INT_ZERO)
  ! number of errors printed by each rank before messages are suppressed
  PetscInt, parameter :: MAX_ERRORS_PRINTED = 10
  ! tolerances for the geometric checks (relative)
  PetscReal, parameter :: ZERO_AREA_TOL = 1.d-12
  PetscReal, parameter :: ZERO_VOLUME_TOL = 1.d-12
  PetscReal, parameter :: SHARED_FACE_TOL = 1.d-10

  public :: UGridPolyhedraReadHDF5, &
            UGridPolyhedraDecomposeDistributed, &
            UGridPolyhedraWriteXDMFMesh, &
            UGridPolyhedraXDMFStreamLength

contains

! ************************************************************************** !

subroutine UGridPolyhedraReadHDF5(ugrid,filename,option)
  !
  ! Reads this rank's contiguous I/O block of a polyhedra mesh in the HDF5
  ! mixed XDMF-code schema, validates it, and stores it in
  ! ugrid%polyhedra_grid for UGridPolyhedraDecomposeDistributed().
  !
  ! For rank r of P and q = Nc/P, s = mod(Nc,P), the block holds
  ! nlocal = q + (1 if r < s) cells starting at natural id
  ! 1 + r*q + min(r,s). The rank reads Offsets(first:last+1) and then the
  ! connectivity positions Offsets(first):Offsets(last+1)-1, so adjacent
  ! nonempty ranks share one offset entry but no connectivity. Ranks
  ! without cells take part in every collective read with empty
  ! selections. Vertex rows are read in the same block pattern.
  !
  use Option_module
  use String_module
  use HDF5_Aux_module

  implicit none

  type(grid_unstructured_type) :: ugrid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  type(unstructured_polyhedra_type), pointer :: pgrid
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: file_id, grp_id, dset_offsets, dset_conn, dset_vert
  integer(HID_T) :: xfer_id
  integer(HSIZE_T) :: dims(2)
  PetscInt64 :: num_cells_global, num_vertices_global, conn_length
  PetscInt64 :: first_cell, ncell_local, first_vertex, nvert_local
  PetscInt64 :: conn_start, conn_count, conn_sum
  PetscInt64, allocatable, target :: offsets8(:), conn8(:)
  PetscInt :: icell, ivert, nerr, status, bad_position, nfaces, nface_nodes
  PetscInt :: ibeg, iend, nranks, myrank
  PetscBool :: found
  PetscErrorCode :: ierr
  integer :: hdf5_err

  pgrid => ugrid%polyhedra_grid
  pgrid%representation = POLYHEDRA_REPR_XDMF_HDF5
  nranks = option%comm%size
  myrank = option%comm%rank

#ifdef SERIAL_HDF5
  if (option%comm%size > 1) then
    option%io_buffer = 'UNSTRUCTURED_POLYHEDRA HDF5 meshes require a &
      &parallel HDF5 build when run on more than one process.'
    call PrintErrMsg(option)
  endif
#endif

  option%io_buffer = 'Reading polyhedra HDF5 mesh "' // trim(filename) // '".'
  call PrintMsg(option)

  call HDF5FileOpenReadOnly(filename,file_id,PETSC_TRUE,'',option)

  ! schema declarations
  if (.not.PolyH5LinkExists(file_id,'Domain') .or. &
      .not.PolyH5LinkExists(file_id,'Domain/CellsPolyhedra')) then
    option%io_buffer = 'Group /Domain/CellsPolyhedra not found in "' // &
      trim(filename) // '". UNSTRUCTURED_POLYHEDRA HDF5 meshes use &
      &/Domain/Vertices, /Domain/CellsPolyhedra/ConnectivityWithXDMFCodes, &
      &and /Domain/CellsPolyhedra/Offsets; the legacy /Domain/Cells table &
      &is read only by TYPE UNSTRUCTURED.'
    call PrintErrMsg(option)
  endif
  call h5gopen_f(file_id,'Domain/CellsPolyhedra',grp_id,hdf5_err)
  call PolyH5CheckIntAttribute(grp_id,'SchemaVersion',1_8,filename,option)
  call PolyH5CheckIntAttribute(grp_id,'NodeIndexBase',1_8,filename,option)
  call PolyH5CheckIntAttribute(grp_id,'OffsetIndexBase',1_8,filename,option)
  call PolyH5ReadStringAttribute(grp_id,'CellTypeConvention',string,found)
  if (.not.found .or. trim(string) /= 'XDMF') then
    option%io_buffer = 'Attribute /Domain/CellsPolyhedra/&
      &CellTypeConvention must be "XDMF" in "' // trim(filename) // '"'
    if (found) option%io_buffer = trim(option%io_buffer) // &
      ' (found "' // trim(string) // '")'
    option%io_buffer = trim(option%io_buffer) // '.'
    call PrintErrMsg(option)
  endif
  call h5gclose_f(grp_id,hdf5_err)

  ! dataset extents
  call PolyH5OpenDataset(file_id,'Domain/CellsPolyhedra/Offsets',1, &
                         H5T_INTEGER_F,8,filename,dset_offsets,dims,option)
  num_cells_global = int(dims(1),kind=8) - 1
  call PolyH5OpenDataset(file_id, &
                         'Domain/CellsPolyhedra/ConnectivityWithXDMFCodes', &
                         1,H5T_INTEGER_F,8,filename,dset_conn,dims,option)
  conn_length = int(dims(1),kind=8)
  call PolyH5OpenDataset(file_id,'Domain/Vertices',2,H5T_FLOAT_F,-1, &
                         filename,dset_vert,dims,option)
  ! the on-disk [Nv,3] layout appears as (3,Nv) through the Fortran API
  if (dims(1) /= 3) then
    option%io_buffer = 'Dataset /Domain/Vertices in "' // trim(filename) // &
      '" must have shape [number of vertices, 3].'
    call PrintErrMsg(option)
  endif
  num_vertices_global = int(dims(2),kind=8)

  if (num_cells_global < 1 .or. conn_length < 1) then
    option%io_buffer = 'Polyhedra mesh "' // trim(filename) // &
      '" contains no cells.'
    call PrintErrMsg(option)
  endif
  if (num_vertices_global < 1) then
    option%io_buffer = 'Polyhedra mesh "' // trim(filename) // &
      '" contains no vertices.'
    call PrintErrMsg(option)
  endif
  if (num_cells_global > MAX_PETSC_INT .or. &
      num_vertices_global > MAX_PETSC_INT) then
    option%io_buffer = 'The number of cells or vertices in "' // &
      trim(filename) // '" exceeds the PETSc integer range. Configure &
      &PETSc with --with-64-bit-indices.'
    call PrintErrMsg(option)
  endif
  pgrid%num_cells_global = int(num_cells_global)
  pgrid%num_vertices_global = int(num_vertices_global)
  ugrid%nmax = pgrid%num_cells_global
  ugrid%num_vertices_global = pgrid%num_vertices_global

  call h5pcreate_f(H5P_DATASET_XFER_F,xfer_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(xfer_id,H5FD_MPIO_COLLECTIVE_F,hdf5_err)
#endif

  ! offsets of this rank's cells
  call PolyBlockRange(num_cells_global,nranks,myrank,first_cell, &
                      ncell_local)
  allocate(offsets8(ncell_local+1))
  offsets8 = 0
  if (ncell_local > 0) then
    call PolyH5ReadInt64(dset_offsets,first_cell-1,ncell_local+1, &
                         offsets8,xfer_id,hdf5_err)
  else
    call PolyH5ReadInt64(dset_offsets,0_8,0_8,offsets8,xfer_id,hdf5_err)
  endif
  call PolyCheckHDF5Error(hdf5_err,'Offsets',filename,option)

  nerr = 0
  if (ncell_local > 0) then
    if (first_cell == 1 .and. offsets8(1) /= 1) then
      nerr = nerr + 1
      option%io_buffer = 'Offsets(1) must be 1 (1-based offsets).'
      call PrintMsgByRank(option)
    endif
    if (first_cell + ncell_local - 1 == num_cells_global .and. &
        offsets8(ncell_local+1) /= conn_length + 1) then
      nerr = nerr + 1
      option%io_buffer = 'Offsets(Nc+1) must equal the connectivity &
        &length + 1 (' // trim(PolyIntString(conn_length+1)) // '), not ' // &
        trim(PolyIntString(offsets8(ncell_local+1))) // '.'
      call PrintMsgByRank(option)
    endif
    do icell = 1, int(ncell_local)
      if (offsets8(icell) < 1 .or. offsets8(icell+1) > conn_length + 1 .or. &
          offsets8(icell+1) <= offsets8(icell)) then
        nerr = nerr + 1
        if (nerr <= MAX_ERRORS_PRINTED) then
          option%io_buffer = 'Offsets of natural cell ' // &
            trim(PolyIntString(first_cell+icell-1)) // ' are not strictly &
            &increasing within [1,L+1]: ' // &
            trim(PolyIntString(offsets8(icell))) // ', ' // &
            trim(PolyIntString(offsets8(icell+1))) // '.'
          call PrintMsgByRank(option)
        endif
      endif
    enddo
  endif
  call PolyCheckErrors(nerr,'Invalid /Domain/CellsPolyhedra/Offsets in "' // &
                       trim(filename) // '".',option)

  ! connectivity slice of this rank's cells
  conn_start = 0
  conn_count = 0
  if (ncell_local > 0) then
    conn_start = offsets8(1) - 1
    conn_count = offsets8(ncell_local+1) - offsets8(1)
  endif
  nerr = 0
  if (conn_count > MAX_PETSC_INT) then
    nerr = 1
    option%io_buffer = 'Local connectivity slice exceeds the PETSc integer &
      &range; use more processes or 64-bit PETSc indices.'
    call PrintMsgByRank(option)
  endif
  call MPI_Allreduce(conn_count,conn_sum,ONE_INTEGER_MPI,MPI_INTEGER8, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  if (conn_sum /= conn_length) then
    nerr = nerr + 1
    option%io_buffer = 'Connectivity read counts sum to ' // &
      trim(PolyIntString(conn_sum)) // ' instead of ' // &
      trim(PolyIntString(conn_length)) // '.'
    call PrintMsg(option)
  endif
  call PolyCheckErrors(nerr,'Invalid connectivity layout in "' // &
                       trim(filename) // '".',option)

  allocate(conn8(max(conn_count,1_8)))
  conn8 = 0
  call PolyH5ReadInt64(dset_conn,conn_start,conn_count,conn8,xfer_id, &
                       hdf5_err)
  call PolyCheckHDF5Error(hdf5_err,'ConnectivityWithXDMFCodes',filename, &
                          option)

  ! narrow to PetscInt with range checks and validate each record
  pgrid%io_num_cells = int(ncell_local)
  pgrid%io_first_cell = int(first_cell)
  allocate(pgrid%io_record_offset(pgrid%io_num_cells+1))
  allocate(pgrid%io_records(max(conn_count,1_8)))
  pgrid%io_records = 0
  do icell = 1, pgrid%io_num_cells + 1
    pgrid%io_record_offset(icell) = int(offsets8(icell) - offsets8(1)) + 1
  enddo
  if (pgrid%io_num_cells == 0) pgrid%io_record_offset(1) = 1
  nerr = 0
  do ivert = 1, int(conn_count)
    if (abs(conn8(ivert)) > MAX_PETSC_INT) then
      nerr = nerr + 1
      if (nerr <= MAX_ERRORS_PRINTED) then
        option%io_buffer = 'Connectivity value at position ' // &
          trim(PolyIntString(conn_start+ivert)) // &
          ' exceeds the PETSc integer range.'
        call PrintMsgByRank(option)
      endif
    else
      pgrid%io_records(ivert) = int(conn8(ivert))
    endif
  enddo
  deallocate(conn8)
  if (nerr == 0) then
    do icell = 1, pgrid%io_num_cells
      ibeg = pgrid%io_record_offset(icell)
      iend = pgrid%io_record_offset(icell+1) - 1
      call XDMFRecordCheck(pgrid%io_records(ibeg:iend), &
                           pgrid%num_vertices_global,nfaces,nface_nodes, &
                           status,bad_position)
      if (status /= XDMF_RECORD_OK) then
        nerr = nerr + 1
        if (nerr <= MAX_ERRORS_PRINTED) then
          option%io_buffer = 'Natural cell ' // &
            trim(StringWrite(pgrid%io_first_cell+icell-1)) // ': ' // &
            trim(XDMFRecordStatusString(status))
          if (bad_position > 0) then
            option%io_buffer = trim(option%io_buffer) // &
              ' (connectivity position ' // &
              trim(PolyIntString(offsets8(icell)+bad_position-1)) // ')'
          endif
          option%io_buffer = trim(option%io_buffer) // '.'
          call PrintMsgByRank(option)
        endif
      endif
    enddo
  endif
  deallocate(offsets8)
  call PolyCheckErrors(nerr,'Invalid cell records in &
    &/Domain/CellsPolyhedra/ConnectivityWithXDMFCodes of "' // &
    trim(filename) // '".',option)

  ! this rank's block of vertex rows
  call PolyBlockRange(num_vertices_global,nranks,myrank,first_vertex, &
                      nvert_local)
  allocate(pgrid%io_vertex_xyz(3,max(nvert_local,1_8)))
  pgrid%io_vertex_xyz = 0.d0
  call PolyH5ReadVertices(dset_vert,first_vertex-1,nvert_local, &
                          pgrid%io_vertex_xyz,xfer_id,hdf5_err)
  call PolyCheckHDF5Error(hdf5_err,'Vertices',filename,option)
  pgrid%io_first_vertex = int(first_vertex)
  pgrid%io_num_vertices = int(nvert_local)

  call h5pclose_f(xfer_id,hdf5_err)
  call HDF5DatasetClose(dset_offsets,option)
  call HDF5DatasetClose(dset_conn,option)
  call HDF5DatasetClose(dset_vert,option)
  call HDF5FileClose(file_id,option)

  option%io_buffer = '  ' // trim(PolyIntString(num_cells_global)) // &
    ' cells, ' // trim(PolyIntString(num_vertices_global)) // &
    ' vertices, ' // trim(PolyIntString(conn_length)) // &
    ' connectivity integers.'
  call PrintMsg(option)
  call PolyReportMemory('I/O block records (bytes)', &
                        int(size(pgrid%io_records),kind=8)*8,option)

end subroutine UGridPolyhedraReadHDF5

! ************************************************************************** !

subroutine UGridPolyhedraDecomposeDistributed(ugrid,option)
  !
  ! Decomposes a polyhedra mesh read by UGridPolyhedraReadHDF5() and
  ! builds the local/ghost topology and geometry. See the module header.
  !
  use Option_module
  use String_module
  use Utility_module, only : DeallocateArray

  implicit none

  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option

  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscInt :: nranks, nio, icell, iface, nerr, n, i, count
  PetscInt :: nfaces, nface_nodes, ibeg, iend, pos, dest
  PetscInt :: nocc_io, nfn_io, nlocal_new, global_offset_new, ngmax
  PetscInt :: nghost, ghosted_id, local_id, occ
  PetscInt :: status, bad_position
  ! I/O-block cell faces (natural node ids) and their neighbors
  PetscInt, allocatable :: io_face_offset(:), io_fv_offset(:), io_fv(:)
  PetscInt, allocatable :: io_nbr_nat(:), io_nbr_face(:)
  ! partitioning
  PetscInt, allocatable :: adj_ia(:), adj_ja(:)
  PetscInt, pointer :: dest_rank(:), new_petsc(:)
  Mat :: adj_mat
  IS :: is_new, is_num
  ! migrated/ghost cell messages
  PetscInt, allocatable :: owned_buf(:), owned_pos(:)
  PetscInt, allocatable :: ghost_buf(:), ghost_pos(:)
  PetscInt, allocatable :: ghost_petsc(:), ghost_natural(:)
  PetscInt, allocatable :: send_buf(:), recv_buf(:)
  PetscInt, allocatable :: send_count(:), recv_count(:), cursor(:)
  PetscInt, allocatable :: rank_offset(:)
  PetscInt, allocatable :: int_array(:), int_array2(:)
  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid
  nranks = option%comm%size
  nio = pgrid%io_num_cells

  call PrintMsg(option,'Decomposing polyhedra mesh (distributed).')

  if (pgrid%num_cells_global < nranks) then
    option%io_buffer = 'The polyhedra mesh has ' // &
      trim(StringWrite(pgrid%num_cells_global)) // ' cells, fewer than &
      &the ' // trim(StringWrite(nranks)) // ' MPI processes; every &
      &process must own at least one cell.'
    call PrintErrMsg(option)
  endif

  ! ---- expand I/O-block records into faces ---------------------------------
  allocate(io_face_offset(nio+1))
  io_face_offset(1) = 1
  nfn_io = 0
  do icell = 1, nio
    ibeg = pgrid%io_record_offset(icell)
    iend = pgrid%io_record_offset(icell+1) - 1
    call XDMFRecordCheck(pgrid%io_records(ibeg:iend), &
                         pgrid%num_vertices_global,nfaces,nface_nodes, &
                         status,bad_position)
    io_face_offset(icell+1) = io_face_offset(icell) + nfaces
    nfn_io = nfn_io + nface_nodes
  enddo
  nocc_io = io_face_offset(nio+1) - 1
  allocate(io_fv_offset(nocc_io+1))
  allocate(io_fv(max(nfn_io,1)))
  io_fv_offset(1) = 1
  pos = 0
  do icell = 1, nio
    ibeg = pgrid%io_record_offset(icell)
    iend = pgrid%io_record_offset(icell+1) - 1
    nfaces = io_face_offset(icell+1) - io_face_offset(icell)
    allocate(int_array(nfaces+1))
    call XDMFRecordCheck(pgrid%io_records(ibeg:iend), &
                         pgrid%num_vertices_global,i,nface_nodes, &
                         status,bad_position)
    call XDMFRecordExpand(pgrid%io_records(ibeg:iend),int_array, &
                          io_fv(pos+1:pos+nface_nodes))
    do iface = 1, nfaces
      occ = io_face_offset(icell) + iface - 1
      io_fv_offset(occ+1) = pos + int_array(iface+1)
    enddo
    pos = pos + nface_nodes
    deallocate(int_array)
  enddo

  ! ---- match faces through the distributed directory -----------------------
  allocate(io_nbr_nat(max(nocc_io,1)))
  allocate(io_nbr_face(max(nocc_io,1)))
  call PolyMatchFaces(nio,pgrid%io_first_cell,io_face_offset,io_fv_offset, &
                      io_fv,io_nbr_nat,io_nbr_face,option)

  ! ---- dual graph and partition ---------------------------------------------
  allocate(adj_ia(nio+1))
  allocate(adj_ja(max(nocc_io,1)))
  adj_ia(1) = 0
  nerr = 0
  count = 0
  do icell = 1, nio
    ibeg = count
    do occ = io_face_offset(icell), io_face_offset(icell+1)-1
      if (io_nbr_nat(occ) > 0) then
        count = count + 1
        adj_ja(count) = io_nbr_nat(occ) - 1
      endif
    enddo
    n = count - ibeg
    if (n > 1) then
      call PolySortInt(n,adj_ja(ibeg+1:count),ierr);CHKERRQ(ierr)
      do i = ibeg+2, count
        if (adj_ja(i) == adj_ja(i-1)) then
          nerr = nerr + 1
          if (nerr <= MAX_ERRORS_PRINTED) then
            option%io_buffer = 'Natural cells ' // &
              trim(StringWrite(pgrid%io_first_cell+icell-1)) // ' and ' // &
              trim(StringWrite(adj_ja(i)+1)) // ' share more than one &
              &face; this topology is not supported.'
            call PrintMsgByRank(option)
          endif
        endif
      enddo
    endif
    adj_ia(icell+1) = count
  enddo
  call PolyCheckErrors(nerr,'Unsupported polyhedra mesh topology.',option)

  call MatCreateMPIAdj(option%mycomm,nio,pgrid%num_cells_global,adj_ia, &
                       adj_ja,PETSC_NULL_INTEGER_ARRAY,adj_mat, &
                       ierr);CHKERRQ(ierr)
  call UGridPartition(ugrid,option,adj_mat,is_new,nlocal_new)
  ! the Fortran MatCreateMPIAdj() does not take ownership of adj_ia/adj_ja
  call MatDestroy(adj_mat,ierr);CHKERRQ(ierr)
  deallocate(adj_ia)
  deallocate(adj_ja)
  call ISPartitioningToNumbering(is_new,is_num,ierr);CHKERRQ(ierr)
  call ISGetIndices(is_new,dest_rank,ierr);CHKERRQ(ierr)
  call ISGetIndices(is_num,new_petsc,ierr);CHKERRQ(ierr)

  global_offset_new = 0
  call MPI_Exscan(nlocal_new,global_offset_new,ONE_INTEGER_MPI, &
                  MPIU_INTEGER,MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)

  ! ---- migrate cells to their new owners -----------------------------------
  ! message: petsc id (0-based), natural id, record length, record,
  !          number of faces, (neighbor natural id, neighbor face) per face
  allocate(send_count(nranks))
  allocate(recv_count(nranks))
  allocate(cursor(nranks))
  send_count = 0
  do icell = 1, nio
    dest = dest_rank(icell) + 1
    send_count(dest) = send_count(dest) + PolyMessageLength( &
      pgrid%io_record_offset(icell+1)-pgrid%io_record_offset(icell), &
      io_face_offset(icell+1)-io_face_offset(icell))
  enddo
  allocate(send_buf(max(sum(send_count),1)))
  cursor(1) = 0
  do i = 2, nranks
    cursor(i) = cursor(i-1) + send_count(i-1)
  enddo
  do icell = 1, nio
    dest = dest_rank(icell) + 1
    ibeg = pgrid%io_record_offset(icell)
    iend = pgrid%io_record_offset(icell+1) - 1
    call PolyPackMessage(new_petsc(icell),pgrid%io_first_cell+icell-1, &
                         pgrid%io_records(ibeg:iend), &
                         io_nbr_nat(io_face_offset(icell): &
                                    io_face_offset(icell+1)-1), &
                         io_nbr_face(io_face_offset(icell): &
                                     io_face_offset(icell+1)-1), &
                         send_buf,cursor(dest))
  enddo
  call ISRestoreIndices(is_new,dest_rank,ierr);CHKERRQ(ierr)
  call ISRestoreIndices(is_num,new_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(is_new,ierr);CHKERRQ(ierr)
  call ISDestroy(is_num,ierr);CHKERRQ(ierr)
  deallocate(io_face_offset,io_fv_offset,io_fv,io_nbr_nat,io_nbr_face)
  call DeallocateArray(pgrid%io_record_offset)
  call DeallocateArray(pgrid%io_records)

  call PolyExchangeInt(send_buf,send_count,owned_buf,recv_count, &
                       'cell migration',option)
  deallocate(send_buf)

  ! index the received messages by local id
  allocate(owned_pos(nlocal_new))
  owned_pos = 0
  nerr = 0
  pos = 0
  do while (pos < sum(recv_count))
    local_id = owned_buf(pos+1) - global_offset_new + 1
    if (local_id < 1 .or. local_id > nlocal_new) then
      nerr = nerr + 1
    else if (owned_pos(local_id) /= 0) then
      nerr = nerr + 1
    else
      owned_pos(local_id) = pos + 1
    endif
    pos = pos + PolyMessageLengthAt(owned_buf,pos+1)
  enddo
  if (any(owned_pos == 0)) nerr = nerr + 1
  call PolyCheckErrors(nerr,'Inconsistent cell migration in polyhedra &
                       &decomposition.',option)

  ! ---- numbering: natural <-> PETSc -----------------------------------------
  ugrid%nlmax = nlocal_new
  ugrid%global_offset = global_offset_new
  allocate(ugrid%cell_ids_petsc(nlocal_new))
  allocate(int_array(nlocal_new))
  allocate(int_array2(nlocal_new))
  do local_id = 1, nlocal_new
    ugrid%cell_ids_petsc(local_id) = global_offset_new + local_id
    int_array(local_id) = owned_buf(owned_pos(local_id)+1) - 1
    int_array2(local_id) = global_offset_new + local_id - 1
  enddo
  call AOCreateMemoryScalable(option%mycomm,nlocal_new,int_array,int_array2, &
                              ugrid%ao_natural_to_petsc,ierr);CHKERRQ(ierr)
  deallocate(int_array,int_array2)

  ! neighbors of owned cells that are owned elsewhere become ghosts
  count = 0
  do local_id = 1, nlocal_new
    count = count + PolyMessageNumNeighbors(owned_buf,owned_pos(local_id))
  enddo
  allocate(int_array(max(count,1)))
  count = 0
  do local_id = 1, nlocal_new
    call PolyMessageAppendNeighbors(owned_buf,owned_pos(local_id), &
                                    int_array,count)
  enddo
  if (count > 0) then
    call PetscSortRemoveDupsInt(count,int_array,ierr);CHKERRQ(ierr)
  endif
  allocate(int_array2(max(count,1)))
  int_array2(1:count) = int_array(1:count) - 1
  call AOApplicationToPetsc(ugrid%ao_natural_to_petsc,count,int_array2, &
                            ierr);CHKERRQ(ierr)
  nghost = 0
  do i = 1, count
    if (int_array2(i) < global_offset_new .or. &
        int_array2(i) >= global_offset_new + nlocal_new) then
      nghost = nghost + 1
      int_array2(nghost) = int_array2(i)
      int_array(nghost) = int_array(i)
    endif
  enddo
  allocate(ghost_petsc(max(nghost,1)))
  allocate(ghost_natural(max(nghost,1)))
  ghost_petsc(1:nghost) = int_array2(1:nghost)
  ghost_natural(1:nghost) = int_array(1:nghost)
  deallocate(int_array,int_array2)
  if (nghost > 1) then
    call PolySortIntWithArray(nghost,ghost_petsc,ghost_natural, &
                               ierr);CHKERRQ(ierr)
  endif
  ngmax = nlocal_new + nghost
  ugrid%num_ghost_cells = nghost
  ugrid%ngmax = ngmax
  allocate(ugrid%ghost_cell_ids_petsc(max(nghost,1)))
  ugrid%ghost_cell_ids_petsc(1:nghost) = ghost_petsc(1:nghost) + 1
  if (nghost == 0) then
    deallocate(ugrid%ghost_cell_ids_petsc)
    nullify(ugrid%ghost_cell_ids_petsc)
  endif
  allocate(ugrid%cell_ids_natural(ngmax))
  do local_id = 1, nlocal_new
    ugrid%cell_ids_natural(local_id) = owned_buf(owned_pos(local_id)+1)
  enddo
  ugrid%cell_ids_natural(nlocal_new+1:ngmax) = ghost_natural(1:nghost)

  ! ---- fetch ghost cell records from their owners --------------------------
  allocate(rank_offset(nranks+1))
  call PolyRankOffsets(nlocal_new,rank_offset,option)
  send_count = 0
  do i = 1, nghost
    dest = PolyOwnerFromOffsets(ghost_petsc(i),rank_offset)
    send_count(dest) = send_count(dest) + 1
  enddo
  ! ghost_petsc is sorted and rank ranges ascend, so requests are grouped
  call PolyExchangeInt(ghost_petsc,send_count,recv_buf,recv_count, &
                       'ghost requests',option)
  send_count = 0
  do dest = 1, nranks
    ibeg = sum(recv_count(1:dest-1))
    do i = ibeg+1, ibeg+recv_count(dest)
      local_id = recv_buf(i) - global_offset_new + 1
      send_count(dest) = send_count(dest) + &
        PolyMessageLengthAt(owned_buf,owned_pos(local_id))
    enddo
  enddo
  allocate(send_buf(max(sum(send_count),1)))
  pos = 0
  do i = 1, sum(recv_count)
    local_id = recv_buf(i) - global_offset_new + 1
    n = PolyMessageLengthAt(owned_buf,owned_pos(local_id))
    send_buf(pos+1:pos+n) = &
      owned_buf(owned_pos(local_id):owned_pos(local_id)+n-1)
    pos = pos + n
  enddo
  deallocate(recv_buf)
  call PolyExchangeInt(send_buf,send_count,ghost_buf,recv_count, &
                       'ghost records',option)
  deallocate(send_buf)
  allocate(ghost_pos(max(nghost,1)))
  ghost_pos = 0
  nerr = 0
  pos = 0
  do while (pos < sum(recv_count))
    call PolyFindInt(ghost_buf(pos+1),nghost,ghost_petsc,i, &
                      ierr);CHKERRQ(ierr)
    if (i < 0) then
      nerr = nerr + 1
    else
      ghost_pos(i+1) = pos + 1
    endif
    pos = pos + PolyMessageLengthAt(ghost_buf,pos+1)
  enddo
  if (nghost > 0) then
    if (any(ghost_pos(1:nghost) == 0)) nerr = nerr + 1
  endif
  call PolyCheckErrors(nerr,'Inconsistent ghost cell exchange in polyhedra &
                       &decomposition.',option)
  deallocate(rank_offset,ghost_petsc,ghost_natural)

  ! ---- local topology, vertices, geometry ----------------------------------
  call PolyBuildLocalGrid(ugrid,owned_buf,owned_pos,ghost_buf,ghost_pos, &
                          option)
  deallocate(owned_buf,owned_pos,ghost_buf,ghost_pos)
  deallocate(send_count,recv_count,cursor)

  ! all cells use the polyhedra branch of shared routines, including
  ! conventional ones; cell_xdmf_code keeps the original type
  allocate(ugrid%cell_type(ngmax))
  ugrid%cell_type = POLY_TYPE

  ! ---- neighbor lists in the legacy layout (padded to the local maximum) ---
  n = 0
  do local_id = 1, ugrid%nlmax
    count = 0
    do occ = pgrid%cell_face_offset(local_id), &
             pgrid%cell_face_offset(local_id+1)-1
      if (pgrid%face_nbr_ghosted(occ) > 0) count = count + 1
    enddo
    n = max(n,count)
  enddo
  ugrid%max_ndual_per_cell = n
  allocate(ugrid%cell_neighbors_local_ghosted(0:n,ugrid%nlmax))
  ugrid%cell_neighbors_local_ghosted = 0
  do local_id = 1, ugrid%nlmax
    count = 0
    do occ = pgrid%cell_face_offset(local_id), &
             pgrid%cell_face_offset(local_id+1)-1
      ghosted_id = pgrid%face_nbr_ghosted(occ)
      if (ghosted_id > 0) then
        count = count + 1
        ! ghosted neighbors are flagged negative (see UGridNaturalToPetsc)
        if (ghosted_id > ugrid%nlmax) ghosted_id = -ghosted_id
        ugrid%cell_neighbors_local_ghosted(count,local_id) = ghosted_id
      endif
    enddo
    ugrid%cell_neighbors_local_ghosted(0,local_id) = count
  enddo

  call PrintMsg(option,'Polyhedra mesh decomposition complete.')

end subroutine UGridPolyhedraDecomposeDistributed

! ************************************************************************** !

function UGridPolyhedraXDMFStreamLength(ugrid,option)
  !
  ! Length of the global mixed-topology connectivity stream written by
  ! UGridPolyhedraWriteXDMFMesh() (sum of the record lengths of all cells).
  !
  use Option_module

  implicit none

  type(grid_unstructured_type) :: ugrid
  type(option_type) :: option
  PetscInt64 :: UGridPolyhedraXDMFStreamLength

  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscInt64 :: local_length
  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid
  local_length = int(pgrid%cell_record_offset(ugrid%nlmax+1) - &
                     pgrid%cell_record_offset(1),kind=8)
  call MPI_Allreduce(local_length,UGridPolyhedraXDMFStreamLength, &
                     ONE_INTEGER_MPI,MPI_INTEGER8,MPI_SUM,option%mycomm, &
                     ierr);CHKERRQ(ierr)

end function UGridPolyhedraXDMFStreamLength

! ************************************************************************** !

subroutine UGridPolyhedraWriteXDMFMesh(ugrid,loc_id,option,stream_length)
  !
  ! Writes the mesh of a polyhedra grid read from the HDF5 XDMF-code schema
  ! for XDMF3 visualization, under loc_id:
  !
  !   Vertices  float64 [Nv, 3]  rows in natural vertex order
  !   Cells     int64   [L]      mixed-topology stream in natural cell order
  !
  ! Cells holds each cell's original record (codes 6-9 stay compact; code
  ! 16 keeps its face records) with node ids converted to the 0-based
  ! vertex row indices required by XDMF; codes and counts are unchanged.
  ! Natural cell order matches cell fields written from natural vectors.
  !
  ! Owned records are routed to the ranks holding contiguous blocks of
  ! natural ids, and every rank writes its block (and its I/O block of
  ! vertex rows) as a hyperslab; no rank assembles the output mesh.
  !
  use Option_module

  implicit none

  type(grid_unstructured_type) :: ugrid
  integer(HID_T) :: loc_id
  type(option_type) :: option
  PetscInt64, intent(out) :: stream_length

  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscInt :: nranks, local_id, dest, reclen, ibeg, pos, i, nerr, nblock
  PetscInt :: count
  PetscInt64 :: first_cell, nblock8, local_length, offset
  PetscInt, allocatable :: send_buf(:), recv_buf(:), send_count(:)
  PetscInt, allocatable :: recv_count(:), cursor(:), block_pos(:)
  PetscInt, allocatable :: record(:)
  PetscInt64, allocatable, target :: stream(:)
  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid
  nranks = option%comm%size

  ! route owned records to the natural-block owners
  allocate(send_count(nranks))
  allocate(cursor(nranks))
  send_count = 0
  do local_id = 1, ugrid%nlmax
    dest = PolyBlockOwner(int(ugrid%cell_ids_natural(local_id),kind=8), &
                          int(pgrid%num_cells_global,kind=8),nranks) + 1
    reclen = pgrid%cell_record_offset(local_id+1) - &
             pgrid%cell_record_offset(local_id)
    send_count(dest) = send_count(dest) + 2 + reclen
  enddo
  allocate(send_buf(max(sum(send_count),1)))
  cursor(1) = 0
  do i = 2, nranks
    cursor(i) = cursor(i-1) + send_count(i-1)
  enddo
  do local_id = 1, ugrid%nlmax
    dest = PolyBlockOwner(int(ugrid%cell_ids_natural(local_id),kind=8), &
                          int(pgrid%num_cells_global,kind=8),nranks) + 1
    ibeg = pgrid%cell_record_offset(local_id)
    reclen = pgrid%cell_record_offset(local_id+1) - ibeg
    pos = cursor(dest)
    send_buf(pos+1) = ugrid%cell_ids_natural(local_id)
    send_buf(pos+2) = reclen
    send_buf(pos+3:pos+2+reclen) = pgrid%cell_records(ibeg:ibeg+reclen-1)
    cursor(dest) = pos + 2 + reclen
  enddo
  call PolyExchangeInt(send_buf,send_count,recv_buf,recv_count, &
                       'output cell records',option,PETSC_FALSE)
  deallocate(send_buf,send_count,cursor)

  ! order this block's records by natural id
  call PolyBlockRange(int(pgrid%num_cells_global,kind=8),nranks, &
                      int(option%comm%rank),first_cell,nblock8)
  nblock = int(nblock8)
  allocate(block_pos(max(nblock,1)))
  block_pos = 0
  nerr = 0
  local_length = 0
  pos = 0
  do while (pos < sum(recv_count))
    i = recv_buf(pos+1) - int(first_cell) + 1
    if (i < 1 .or. i > nblock) then
      nerr = nerr + 1
    else if (block_pos(i) /= 0) then
      nerr = nerr + 1
    else
      block_pos(i) = pos + 1
      local_length = local_length + recv_buf(pos+2)
    endif
    pos = pos + 2 + recv_buf(pos+2)
  enddo
  if (nblock > 0) then
    if (any(block_pos(1:nblock) == 0)) nerr = nerr + 1
  endif
  call PolyCheckErrors(nerr,'Inconsistent cell records in polyhedra &
                       &output.',option)

  allocate(stream(max(local_length,1_8)))
  count = 0
  do i = 1, nblock
    pos = block_pos(i)
    reclen = recv_buf(pos+1)
    allocate(record(reclen))
    record = recv_buf(pos+2:pos+1+reclen)
    call XDMFRecordShiftNodes(record,-1)
    stream(count+1:count+reclen) = int(record,kind=8)
    count = count + reclen
    deallocate(record)
  enddo
  deallocate(recv_buf,recv_count,block_pos)

  offset = 0
  call MPI_Exscan(local_length,offset,ONE_INTEGER_MPI,MPI_INTEGER8, &
                  MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  if (option%comm%rank == 0) offset = 0
  call MPI_Allreduce(local_length,stream_length,ONE_INTEGER_MPI, &
                     MPI_INTEGER8,MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)

  call PolyH5WriteInt64(loc_id,'Cells',stream_length,offset,local_length, &
                        stream,option)
  deallocate(stream)

  call PolyH5WriteVertices(loc_id,'Vertices', &
                           int(pgrid%num_vertices_global,kind=8), &
                           int(pgrid%io_first_vertex-1,kind=8), &
                           int(pgrid%io_num_vertices,kind=8), &
                           pgrid%io_vertex_xyz,option)

end subroutine UGridPolyhedraWriteXDMFMesh

! ************************************************************************** !

subroutine PolyMatchFaces(nio,first_cell,face_offset,fv_offset,fv, &
                          nbr_natural,nbr_face,option)
  !
  ! Matches the cell-face occurrences of this rank's I/O block with those
  ! of all other blocks. Each occurrence is routed to a directory rank
  ! chosen by hashing its sorted node ids (the full key). The directory
  ! sorts occurrences by full key, so hash collisions only affect load;
  ! equal keys identify a shared face. Exactly two occurrences with
  ! reversed cyclic order make an interior face; one makes a boundary face;
  ! more than two (non-manifold), equal winding, or a cell sharing a face
  ! with itself is an error. Returns, for each occurrence, the natural id
  ! and cell-local face index of the matching occurrence (0 on the
  ! boundary).
  !
  use Option_module
  use String_module

  implicit none

  PetscInt, intent(in) :: nio
  PetscInt, intent(in) :: first_cell
  PetscInt, intent(in) :: face_offset(:)
  PetscInt, intent(in) :: fv_offset(:)
  PetscInt, intent(in) :: fv(:)
  PetscInt, intent(out) :: nbr_natural(:)
  PetscInt, intent(out) :: nbr_face(:)
  type(option_type) :: option

  PetscInt :: nranks, nocc, icell, iface, occ, nv, dest, pos, i, j, g
  PetscInt :: nmsg, gbeg, gend, nerr, m1, m2, isrc
  PetscInt, allocatable :: key(:), occ_dest(:), occ_slot(:)
  PetscInt, allocatable :: send_count(:), recv_count(:), cursor(:)
  PetscInt, allocatable :: occ_count(:), occ_cursor(:)
  PetscInt, allocatable :: send_buf(:), recv_buf(:)
  PetscInt, allocatable :: msg_pos(:), key_offset(:), keys(:), perm(:)
  PetscInt, allocatable :: reply(:), reply_recv(:), reply_count(:)
  PetscInt, allocatable :: msgs_from(:)
  PetscErrorCode :: ierr

  nranks = option%comm%size
  nocc = face_offset(nio+1) - 1

  ! route each occurrence: nv, cyclic nodes, natural cell id, local face
  allocate(occ_dest(max(nocc,1)))
  allocate(occ_slot(max(nocc,1)))
  allocate(send_count(nranks))
  allocate(occ_count(nranks))
  send_count = 0
  occ_count = 0
  do occ = 1, nocc
    nv = fv_offset(occ+1) - fv_offset(occ)
    allocate(key(nv))
    key = fv(fv_offset(occ):fv_offset(occ+1)-1)
    call PolySortInt(nv,key,ierr);CHKERRQ(ierr)
    dest = PolyKeyOwner(key,nranks) + 1
    deallocate(key)
    occ_dest(occ) = dest
    send_count(dest) = send_count(dest) + nv + 3
    occ_count(dest) = occ_count(dest) + 1
  enddo
  allocate(send_buf(max(sum(send_count),1)))
  allocate(cursor(nranks))
  allocate(occ_cursor(nranks))
  cursor(1) = 0
  occ_cursor(1) = 0
  do i = 2, nranks
    cursor(i) = cursor(i-1) + send_count(i-1)
    occ_cursor(i) = occ_cursor(i-1) + occ_count(i-1)
  enddo
  do icell = 1, nio
    do iface = 1, face_offset(icell+1) - face_offset(icell)
      occ = face_offset(icell) + iface - 1
      dest = occ_dest(occ)
      nv = fv_offset(occ+1) - fv_offset(occ)
      pos = cursor(dest)
      send_buf(pos+1) = nv
      send_buf(pos+2:pos+nv+1) = fv(fv_offset(occ):fv_offset(occ+1)-1)
      send_buf(pos+nv+2) = first_cell + icell - 1
      send_buf(pos+nv+3) = iface
      cursor(dest) = pos + nv + 3
      occ_cursor(dest) = occ_cursor(dest) + 1
      occ_slot(occ) = occ_cursor(dest)
    enddo
  enddo
  call PolyExchangeInt(send_buf,send_count,recv_buf,recv_count, &
                       'face directory',option)
  deallocate(send_buf)

  ! directory: index received occurrences and their sorted keys
  nmsg = 0
  pos = 0
  do while (pos < sum(recv_count))
    nmsg = nmsg + 1
    pos = pos + recv_buf(pos+1) + 3
  enddo
  allocate(msg_pos(max(nmsg,1)))
  allocate(msgs_from(nranks))
  allocate(key_offset(nmsg+1))
  allocate(keys(max(sum(recv_count),1)))
  msgs_from = 0
  key_offset(1) = 1
  pos = 0
  isrc = 1
  do i = 1, nmsg
    do while (isrc < nranks .and. pos >= sum(recv_count(1:isrc)))
      isrc = isrc + 1
    enddo
    msgs_from(isrc) = msgs_from(isrc) + 1
    msg_pos(i) = pos + 1
    nv = recv_buf(pos+1)
    keys(key_offset(i):key_offset(i)+nv-1) = recv_buf(pos+2:pos+nv+1)
    call PolySortInt(nv,keys(key_offset(i):key_offset(i)+nv-1), &
                      ierr);CHKERRQ(ierr)
    key_offset(i+1) = key_offset(i) + nv
    pos = pos + nv + 3
  enddo
  allocate(perm(max(nmsg,1)))
  call XDMFSortKeys(nmsg,key_offset,keys,perm)

  ! reply per received occurrence: matching natural cell id and local face
  allocate(reply(max(2*nmsg,1)))
  reply = 0
  nerr = 0
  gbeg = 1
  do while (gbeg <= nmsg)
    gend = gbeg
    do while (gend < nmsg)
      if (XDMFCompareKeys( &
            keys(key_offset(perm(gbeg)):key_offset(perm(gbeg)+1)-1), &
            keys(key_offset(perm(gend+1)):key_offset(perm(gend+1)+1)-1)) &
          /= 0) exit
      gend = gend + 1
    enddo
    g = gend - gbeg + 1
    if (g == 2) then
      m1 = perm(gbeg)
      m2 = perm(gend)
      nv = recv_buf(msg_pos(m1))
      i = msg_pos(m1) + nv + 1 ! natural cell id of m1
      j = msg_pos(m2) + nv + 1 ! natural cell id of m2
      if (recv_buf(i) == recv_buf(j)) then
        nerr = nerr + 1
        if (nerr <= MAX_ERRORS_PRINTED) then
          option%io_buffer = 'Natural cell ' // &
            trim(StringWrite(recv_buf(i))) // ' lists the same face twice.'
          call PrintMsgByRank(option)
        endif
      else if (.not.XDMFFaceIsReversed( &
                 recv_buf(msg_pos(m1)+1:msg_pos(m1)+nv), &
                 recv_buf(msg_pos(m2)+1:msg_pos(m2)+nv))) then
        nerr = nerr + 1
        if (nerr <= MAX_ERRORS_PRINTED) then
          option%io_buffer = 'Natural cells ' // &
            trim(StringWrite(recv_buf(i))) // ' and ' // &
            trim(StringWrite(recv_buf(j))) // ' list their shared face &
            &with the same winding; each cell must wind its faces outward &
            &(for wedges, nodes 1-2-3 must wind toward nodes 4-5-6).'
          call PrintMsgByRank(option)
        endif
      else
        reply(2*m1-1) = recv_buf(j)
        reply(2*m1) = recv_buf(j+1)
        reply(2*m2-1) = recv_buf(i)
        reply(2*m2) = recv_buf(i+1)
      endif
    else if (g > 2) then
      nerr = nerr + 1
      if (nerr <= MAX_ERRORS_PRINTED) then
        m1 = perm(gbeg)
        nv = recv_buf(msg_pos(m1))
        option%io_buffer = 'A face is shared by ' // trim(StringWrite(g)) // &
          ' cells (non-manifold), including natural cell ' // &
          trim(StringWrite(recv_buf(msg_pos(m1)+nv+1))) // '.'
        call PrintMsgByRank(option)
      endif
    endif
    gbeg = gend + 1
  enddo
  deallocate(recv_buf,msg_pos,key_offset,keys,perm)
  call PolyCheckErrors(nerr,'Invalid face sharing in polyhedra mesh.',option)

  ! return replies in the order the occurrences were received
  allocate(reply_count(nranks))
  reply_count = 2*msgs_from
  call PolyExchangeInt(reply,reply_count,reply_recv,recv_count, &
                       'face directory replies',option)
  deallocate(reply)
  do occ = 1, nocc
    nbr_natural(occ) = reply_recv(2*occ_slot(occ)-1)
    nbr_face(occ) = reply_recv(2*occ_slot(occ))
  enddo
  deallocate(reply_recv,reply_count,msgs_from)
  deallocate(occ_dest,occ_slot,send_count,occ_count,cursor,occ_cursor)
  if (allocated(recv_count)) deallocate(recv_count)

end subroutine PolyMatchFaces

! ************************************************************************** !

subroutine PolyBuildLocalGrid(ugrid,owned_buf,owned_pos,ghost_buf,ghost_pos, &
                              option)
  !
  ! Builds the local (owned + ghost) topology from migrated cell messages,
  ! fetches the referenced vertex coordinates from their I/O-block owners,
  ! computes geometry, and fills the polyhedra grid fields used downstream.
  ! Padded 2D arrays kept for existing consumers (cell_faceids,
  ! cell_vertids, face_vertids, ugrid%cell_vertices) are sized by the
  ! local maxima, not by global maxima.
  !
  use Option_module
  use String_module

  implicit none

  type(grid_unstructured_type) :: ugrid
  PetscInt, intent(in) :: owned_buf(:)
  PetscInt, intent(in) :: owned_pos(:)
  PetscInt, intent(in) :: ghost_buf(:)
  PetscInt, intent(in) :: ghost_pos(:)
  type(option_type) :: option

  type(unstructured_polyhedra_type), pointer :: pgrid
  PetscInt :: nlmax, ngmax, ghosted_id, pos, reclen, nfaces, nface_nodes
  PetscInt :: nocc, nfn, nrec, occ, iface, i, j, n, nv, nvert_local, count
  PetscInt :: status, bad_position, nerr, max_nface, max_nvert_face
  PetscInt :: max_nvert_cell, nbr, nbr_occ, nranks, dest, vid, ibeg
  PetscInt, allocatable :: msg_nbr_nat(:), msg_nbr_face(:)
  PetscInt, allocatable :: fv_natural(:), nat_sorted(:), nat_ghosted(:)
  PetscInt, allocatable :: face_offset(:), cell_vert(:)
  PetscInt, allocatable :: send_count(:), recv_count(:), recv_ids(:)
  PetscReal, allocatable :: send_xyz(:), recv_xyz(:), xyz(:,:)
  PetscInt, allocatable :: local_nodes(:)
  PetscReal :: area, centroid(3), unit_normal(3), nonplanarity
  PetscReal :: max_nonplanarity, volume, extent, max_rel_mismatch
  PetscErrorCode :: ierr

  pgrid => ugrid%polyhedra_grid
  nlmax = ugrid%nlmax
  ngmax = ugrid%ngmax
  nranks = option%comm%size

  ! ---- sizes -----------------------------------------------------------------
  nocc = 0
  nfn = 0
  nrec = 0
  do ghosted_id = 1, ngmax
    call PolyMessageAt(ghosted_id,nlmax,owned_pos,ghost_pos,pos)
    reclen = PolyBuf(ghosted_id,nlmax,owned_buf,ghost_buf,pos+3)
    nrec = nrec + reclen
    call PolyMessageCheck(ghosted_id,nlmax,owned_buf,ghost_buf,pos, &
                          pgrid%num_vertices_global,nfaces,nface_nodes)
    nocc = nocc + nfaces
    nfn = nfn + nface_nodes
  enddo

  allocate(pgrid%cell_xdmf_code(ngmax))
  allocate(pgrid%cell_record_offset(ngmax+1))
  allocate(pgrid%cell_records(max(nrec,1)))
  allocate(pgrid%cell_face_offset(ngmax+1))
  allocate(pgrid%face_vert_offset(nocc+1))
  allocate(pgrid%face_vert_ids(max(nfn,1)))
  allocate(pgrid%face_nbr_ghosted(max(nocc,1)))
  allocate(pgrid%face_nbr_face(max(nocc,1)))
  allocate(msg_nbr_nat(max(nocc,1)))
  allocate(msg_nbr_face(max(nocc,1)))
  allocate(fv_natural(max(nfn,1)))

  ! ---- unpack records and faces (natural node ids) -------------------------
  pgrid%cell_record_offset(1) = 1
  pgrid%cell_face_offset(1) = 1
  pgrid%face_vert_offset(1) = 1
  occ = 0
  nfn = 0
  do ghosted_id = 1, ngmax
    call PolyMessageAt(ghosted_id,nlmax,owned_pos,ghost_pos,pos)
    reclen = PolyBuf(ghosted_id,nlmax,owned_buf,ghost_buf,pos+3)
    i = pgrid%cell_record_offset(ghosted_id)
    do j = 1, reclen
      pgrid%cell_records(i+j-1) = &
        PolyBuf(ghosted_id,nlmax,owned_buf,ghost_buf,pos+3+j)
    enddo
    pgrid%cell_record_offset(ghosted_id+1) = i + reclen
    pgrid%cell_xdmf_code(ghosted_id) = pgrid%cell_records(i)
    call XDMFRecordCheck(pgrid%cell_records(i:i+reclen-1), &
                         pgrid%num_vertices_global,nfaces,nface_nodes, &
                         status,bad_position)
    allocate(face_offset(nfaces+1))
    call XDMFRecordExpand(pgrid%cell_records(i:i+reclen-1),face_offset, &
                          fv_natural(nfn+1:nfn+nface_nodes))
    ! neighbor data follow the record: nfaces, (natural, face) pairs
    j = pos + 4 + reclen
    do iface = 1, nfaces
      occ = occ + 1
      pgrid%face_vert_offset(occ+1) = nfn + face_offset(iface+1)
      msg_nbr_nat(occ) = PolyBuf(ghosted_id,nlmax,owned_buf,ghost_buf, &
                                 j+2*iface-1)
      msg_nbr_face(occ) = PolyBuf(ghosted_id,nlmax,owned_buf,ghost_buf, &
                                  j+2*iface)
    enddo
    nfn = nfn + nface_nodes
    pgrid%cell_face_offset(ghosted_id+1) = occ + 1
    deallocate(face_offset)
  enddo

  ! ---- neighbors: natural id -> ghosted id ---------------------------------
  allocate(nat_sorted(ngmax))
  allocate(nat_ghosted(ngmax))
  do ghosted_id = 1, ngmax
    nat_sorted(ghosted_id) = ugrid%cell_ids_natural(ghosted_id)
    nat_ghosted(ghosted_id) = ghosted_id
  enddo
  call PolySortIntWithArray(ngmax,nat_sorted,nat_ghosted,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, ngmax
    do occ = pgrid%cell_face_offset(ghosted_id), &
             pgrid%cell_face_offset(ghosted_id+1)-1
      pgrid%face_nbr_face(occ) = 0
      if (msg_nbr_nat(occ) == 0) then
        pgrid%face_nbr_ghosted(occ) = 0
      else
        call PolyFindInt(msg_nbr_nat(occ),ngmax,nat_sorted,i, &
                          ierr);CHKERRQ(ierr)
        if (i < 0) then
          pgrid%face_nbr_ghosted(occ) = -1
        else
          nbr = nat_ghosted(i+1)
          pgrid%face_nbr_ghosted(occ) = nbr
          pgrid%face_nbr_face(occ) = pgrid%cell_face_offset(nbr) + &
                                     msg_nbr_face(occ) - 1
        endif
      endif
    enddo
  enddo
  deallocate(nat_sorted,nat_ghosted,msg_nbr_nat,msg_nbr_face)
  nerr = 0
  do ghosted_id = 1, ngmax
    do occ = pgrid%cell_face_offset(ghosted_id), &
             pgrid%cell_face_offset(ghosted_id+1)-1
      nbr = pgrid%face_nbr_ghosted(occ)
      if (nbr > 0) then
        nbr_occ = pgrid%face_nbr_face(occ)
        if (nbr_occ < pgrid%cell_face_offset(nbr) .or. &
            nbr_occ >= pgrid%cell_face_offset(nbr+1)) then
          nerr = nerr + 1
        else if (pgrid%face_nbr_ghosted(nbr_occ) /= ghosted_id .or. &
                 pgrid%face_nbr_face(nbr_occ) /= occ) then
          nerr = nerr + 1
        endif
      endif
    enddo
    if (ghosted_id <= nlmax) then
      ! every neighbor of an owned cell must be in the halo
      do occ = pgrid%cell_face_offset(ghosted_id), &
               pgrid%cell_face_offset(ghosted_id+1)-1
        if (pgrid%face_nbr_ghosted(occ) < 0) nerr = nerr + 1
      enddo
    endif
  enddo
  call PolyCheckErrors(nerr,'Inconsistent face neighbors after polyhedra &
                       &decomposition.',option)

  ! ---- local vertices ----------------------------------------------------------
  allocate(local_nodes(max(nfn,1)))
  local_nodes(1:nfn) = fv_natural(1:nfn)
  nvert_local = nfn
  if (nvert_local > 0) then
    call PetscSortRemoveDupsInt(nvert_local,local_nodes,ierr);CHKERRQ(ierr)
  endif
  allocate(ugrid%vertex_ids_natural(max(nvert_local,1)))
  ugrid%vertex_ids_natural(1:nvert_local) = local_nodes(1:nvert_local)
  deallocate(local_nodes)
  do i = 1, nfn
    call PolyFindInt(fv_natural(i),nvert_local,ugrid%vertex_ids_natural, &
                      j,ierr);CHKERRQ(ierr)
    pgrid%face_vert_ids(i) = j + 1
  enddo
  deallocate(fv_natural)

  ! request coordinates from the owners of the vertex I/O blocks
  allocate(send_count(nranks))
  send_count = 0
  do i = 1, nvert_local
    dest = PolyBlockOwner(int(ugrid%vertex_ids_natural(i),kind=8), &
                          int(pgrid%num_vertices_global,kind=8),nranks) + 1
    send_count(dest) = send_count(dest) + 1
  enddo
  ! vertex_ids_natural is sorted, so requests are grouped by owner
  call PolyExchangeInt(ugrid%vertex_ids_natural,send_count,recv_ids, &
                       recv_count,'vertex requests',option)
  n = sum(recv_count)
  allocate(send_xyz(max(3*n,1)))
  nerr = 0
  do i = 1, n
    vid = recv_ids(i) - pgrid%io_first_vertex + 1
    if (vid < 1 .or. vid > pgrid%io_num_vertices) then
      nerr = nerr + 1
      cycle
    endif
    send_xyz(3*i-2:3*i) = pgrid%io_vertex_xyz(:,vid)
  enddo
  call PolyCheckErrors(nerr,'Inconsistent vertex request in polyhedra &
                       &decomposition.',option)
  call PolyExchangeReal(send_xyz,3*recv_count,recv_xyz,3*send_count, &
                        'vertex coordinates',option)
  deallocate(send_xyz,recv_ids,recv_count)
  ugrid%num_vertices_local = nvert_local
  pgrid%num_vertices_local = nvert_local
  allocate(ugrid%vertices(max(nvert_local,1)))
  allocate(pgrid%vertex_coordinates(max(nvert_local,1)))
  do i = 1, nvert_local
    ugrid%vertices(i)%id = ugrid%vertex_ids_natural(i)
    ugrid%vertices(i)%x = recv_xyz(3*i-2)
    ugrid%vertices(i)%y = recv_xyz(3*i-1)
    ugrid%vertices(i)%z = recv_xyz(3*i)
    pgrid%vertex_coordinates(i) = ugrid%vertices(i)
  enddo
  deallocate(recv_xyz,send_count)

  ! ---- per-cell counts and padded legacy arrays ----------------------------
  pgrid%num_cells_local = nlmax
  pgrid%num_faces_local = nocc
  max_nface = 0
  max_nvert_face = 0
  max_nvert_cell = 0
  allocate(pgrid%cell_ids(ngmax))
  allocate(pgrid%cell_nfaces(ngmax))
  allocate(pgrid%cell_nverts(ngmax))
  allocate(cell_vert(max(nfn,1)))
  do ghosted_id = 1, ngmax
    pgrid%cell_ids(ghosted_id) = ugrid%cell_ids_natural(ghosted_id)
    nfaces = pgrid%cell_face_offset(ghosted_id+1) - &
             pgrid%cell_face_offset(ghosted_id)
    pgrid%cell_nfaces(ghosted_id) = nfaces
    max_nface = max(max_nface,nfaces)
    ibeg = pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id))
    n = pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id+1)) - ibeg
    cell_vert(1:n) = pgrid%face_vert_ids(ibeg:ibeg+n-1)
    call PetscSortRemoveDupsInt(n,cell_vert,ierr);CHKERRQ(ierr)
    pgrid%cell_nverts(ghosted_id) = n
    max_nvert_cell = max(max_nvert_cell,n)
  enddo
  do occ = 1, nocc
    max_nvert_face = max(max_nvert_face,pgrid%face_vert_offset(occ+1) - &
                                        pgrid%face_vert_offset(occ))
  enddo
  pgrid%max_nface_per_cell = max_nface
  pgrid%max_nvert_per_face = max_nvert_face
  pgrid%max_nvert_per_cell = max_nvert_cell
  ugrid%max_nvert_per_cell = max_nvert_cell

  allocate(pgrid%cell_faceids(max_nface,ngmax))
  allocate(pgrid%cell_vertids(max_nvert_cell,ngmax))
  allocate(ugrid%cell_vertices(0:max_nvert_cell,ngmax))
  pgrid%cell_faceids = 0
  pgrid%cell_vertids = 0
  ugrid%cell_vertices = 0
  do ghosted_id = 1, ngmax
    do iface = 1, pgrid%cell_nfaces(ghosted_id)
      pgrid%cell_faceids(iface,ghosted_id) = &
        pgrid%cell_face_offset(ghosted_id) + iface - 1
    enddo
    ibeg = pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id))
    n = pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id+1)) - ibeg
    cell_vert(1:n) = pgrid%face_vert_ids(ibeg:ibeg+n-1)
    call PetscSortRemoveDupsInt(n,cell_vert,ierr);CHKERRQ(ierr)
    pgrid%cell_vertids(1:n,ghosted_id) = cell_vert(1:n)
    ugrid%cell_vertices(0,ghosted_id) = n
    ugrid%cell_vertices(1:n,ghosted_id) = cell_vert(1:n)
  enddo
  deallocate(cell_vert)

  allocate(pgrid%face_ids(max(nocc,1)))
  allocate(pgrid%face_cellids(max(nocc,1)))
  allocate(pgrid%face_nverts(max(nocc,1)))
  allocate(pgrid%face_vertids(max_nvert_face,max(nocc,1)))
  pgrid%face_vertids = 0
  do ghosted_id = 1, ngmax
    do occ = pgrid%cell_face_offset(ghosted_id), &
             pgrid%cell_face_offset(ghosted_id+1)-1
      pgrid%face_ids(occ) = occ
      pgrid%face_cellids(occ) = ghosted_id
      ibeg = pgrid%face_vert_offset(occ)
      nv = pgrid%face_vert_offset(occ+1) - ibeg
      pgrid%face_nverts(occ) = nv
      pgrid%face_vertids(1:nv,occ) = pgrid%face_vert_ids(ibeg:ibeg+nv-1)
    enddo
  enddo

  ! ---- geometry -------------------------------------------------------------
  allocate(pgrid%face_areas(max(nocc,1)))
  allocate(pgrid%face_centroids(max(nocc,1)))
  allocate(pgrid%face_normals(max(nocc,1)))
  allocate(pgrid%cell_volumes(ngmax))
  allocate(pgrid%cell_centroids(ngmax))
  nerr = 0
  max_nonplanarity = 0.d0
  do occ = 1, nocc
    ibeg = pgrid%face_vert_offset(occ)
    nv = pgrid%face_vert_offset(occ+1) - ibeg
    allocate(xyz(3,nv))
    do i = 1, nv
      j = pgrid%face_vert_ids(ibeg+i-1)
      xyz(:,i) = [ugrid%vertices(j)%x,ugrid%vertices(j)%y,ugrid%vertices(j)%z]
    enddo
    call XDMFPolygonGeometry(xyz,area,centroid,unit_normal,nonplanarity)
    extent = 0.d0
    do i = 1, nv
      extent = max(extent,sqrt(sum((xyz(:,i)-xyz(:,1))**2)))
    enddo
    if (area <= ZERO_AREA_TOL*extent*extent) then
      nerr = nerr + 1
      if (nerr <= MAX_ERRORS_PRINTED) then
        ghosted_id = pgrid%face_cellids(occ)
        option%io_buffer = 'Face ' // &
          trim(StringWrite(occ-pgrid%cell_face_offset(ghosted_id)+1)) // &
          ' of natural cell ' // &
          trim(StringWrite(ugrid%cell_ids_natural(ghosted_id))) // &
          ' has zero area.'
        call PrintMsgByRank(option)
      endif
    endif
    max_nonplanarity = max(max_nonplanarity,nonplanarity)
    pgrid%face_areas(occ) = area
    pgrid%face_centroids(occ)%x = centroid(1)
    pgrid%face_centroids(occ)%y = centroid(2)
    pgrid%face_centroids(occ)%z = centroid(3)
    pgrid%face_normals(occ)%x = unit_normal(1)
    pgrid%face_normals(occ)%y = unit_normal(2)
    pgrid%face_normals(occ)%z = unit_normal(3)
    deallocate(xyz)
  enddo
  call PolyCheckErrors(nerr,'Degenerate faces in polyhedra mesh.',option)

  do ghosted_id = 1, ngmax
    n = pgrid%cell_nverts(ghosted_id)
    allocate(xyz(3,n))
    do i = 1, n
      j = pgrid%cell_vertids(i,ghosted_id)
      xyz(:,i) = [ugrid%vertices(j)%x,ugrid%vertices(j)%y,ugrid%vertices(j)%z]
    enddo
    nfaces = pgrid%cell_nfaces(ghosted_id)
    ibeg = pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id))
    count = pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id+1)) - &
            ibeg
    allocate(face_offset(nfaces+1))
    allocate(local_nodes(count))
    do iface = 1, nfaces+1
      face_offset(iface) = &
        pgrid%face_vert_offset(pgrid%cell_face_offset(ghosted_id)+iface-1) - &
        ibeg + 1
    enddo
    ! cell-local node index of each face node
    do i = 1, count
      call PolyFindInt(pgrid%face_vert_ids(ibeg+i-1),n, &
                        pgrid%cell_vertids(1:n,ghosted_id),j, &
                        ierr);CHKERRQ(ierr)
      local_nodes(i) = j + 1
    enddo
    call XDMFPolyhedronGeometry(xyz,face_offset,local_nodes,volume,centroid)
    extent = 0.d0
    do i = 1, 3
      extent = max(extent,maxval(xyz(i,:))-minval(xyz(i,:)))
    enddo
    if (volume <= ZERO_VOLUME_TOL*extent**3) then
      nerr = nerr + 1
      if (nerr <= MAX_ERRORS_PRINTED) then
        option%io_buffer = 'Natural cell ' // &
          trim(StringWrite(ugrid%cell_ids_natural(ghosted_id))) // &
          ' (XDMF code ' // &
          trim(StringWrite(pgrid%cell_xdmf_code(ghosted_id))) // &
          ') has non-positive volume ' // trim(StringWrite(volume)) // &
          '; its faces must wind outward (right-hand rule).'
        if (pgrid%cell_xdmf_code(ghosted_id) == XDMF_CODE_WEDGE) then
          option%io_buffer = trim(option%io_buffer) // ' Wedge nodes &
            &1-2-3 must wind toward nodes 4-5-6.'
        endif
        call PrintMsgByRank(option)
      endif
    endif
    pgrid%cell_volumes(ghosted_id) = volume
    pgrid%cell_centroids(ghosted_id)%x = centroid(1)
    pgrid%cell_centroids(ghosted_id)%y = centroid(2)
    pgrid%cell_centroids(ghosted_id)%z = centroid(3)
    deallocate(xyz,face_offset,local_nodes)
  enddo
  call PolyCheckErrors(nerr,'Inverted or degenerate cells in polyhedra &
                       &mesh.',option)

  ! both copies of a shared face must agree
  max_rel_mismatch = 0.d0
  do occ = 1, nocc
    nbr_occ = pgrid%face_nbr_face(occ)
    if (nbr_occ <= 0) cycle
    max_rel_mismatch = max(max_rel_mismatch, &
      abs(pgrid%face_areas(occ)-pgrid%face_areas(nbr_occ))/ &
        pgrid%face_areas(occ), &
      abs(pgrid%face_normals(occ)%x+pgrid%face_normals(nbr_occ)%x) + &
      abs(pgrid%face_normals(occ)%y+pgrid%face_normals(nbr_occ)%y) + &
      abs(pgrid%face_normals(occ)%z+pgrid%face_normals(nbr_occ)%z))
  enddo
  nerr = 0
  if (max_rel_mismatch > SHARED_FACE_TOL) then
    nerr = 1
    option%io_buffer = 'Shared face geometry mismatch: ' // &
      trim(StringWrite(max_rel_mismatch))
    call PrintMsgByRank(option)
  endif
  call PolyCheckErrors(nerr,'Inconsistent shared-face geometry in polyhedra &
                       &mesh.',option)

  call MPI_Allreduce(MPI_IN_PLACE,max_nonplanarity,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)
  option%io_buffer = '  Maximum face nonplanarity (distance/sqrt(area)): ' // &
    trim(StringWrite(max_nonplanarity))
  call PrintMsg(option)

  call PolyReportMemory('local topology (bytes)', &
    int(size(pgrid%cell_records)+size(pgrid%face_vert_ids)+ &
        size(pgrid%face_vertids)+size(pgrid%cell_faceids)+ &
        size(pgrid%cell_vertids)+size(ugrid%cell_vertices)+ &
        4*nocc+4*ngmax,kind=8)*8 + &
    int(nocc*8+ngmax*4+nvert_local*4,kind=8)*8,option)
  option%io_buffer = '  Local cells (owned/ghost) max over ranks: '
  count = ngmax - nlmax
  call MPI_Allreduce(MPI_IN_PLACE,count,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  n = nlmax
  call MPI_Allreduce(MPI_IN_PLACE,n,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  option%io_buffer = trim(option%io_buffer) // ' ' // &
    trim(StringWrite(n)) // '/' // trim(StringWrite(count))
  call PrintMsg(option)

end subroutine PolyBuildLocalGrid

! ************************************************************************** !
! Cell messages: petsc id (0-based), natural id, record length, record,
! number of faces, then (neighbor natural id, neighbor face) per face.
! ************************************************************************** !

function PolyMessageLength(reclen,nfaces)

  implicit none

  PetscInt, intent(in) :: reclen
  PetscInt, intent(in) :: nfaces
  PetscInt :: PolyMessageLength

  PolyMessageLength = 3 + reclen + 1 + 2*nfaces

end function PolyMessageLength

! ************************************************************************** !

function PolyMessageLengthAt(buf,pos)
  !
  ! Length of the message starting at buf(pos).
  !
  implicit none

  PetscInt, intent(in) :: buf(:)
  PetscInt, intent(in) :: pos
  PetscInt :: PolyMessageLengthAt

  PetscInt :: reclen

  reclen = buf(pos+2)
  PolyMessageLengthAt = PolyMessageLength(reclen,buf(pos+3+reclen))

end function PolyMessageLengthAt

! ************************************************************************** !

subroutine PolyPackMessage(petsc_id,natural_id,record,nbr_natural,nbr_face, &
                           buf,pos)
  !
  ! Appends a cell message to buf after position pos (updated).
  !
  implicit none

  PetscInt, intent(in) :: petsc_id
  PetscInt, intent(in) :: natural_id
  PetscInt, intent(in) :: record(:)
  PetscInt, intent(in) :: nbr_natural(:)
  PetscInt, intent(in) :: nbr_face(:)
  PetscInt :: buf(:)
  PetscInt :: pos

  PetscInt :: reclen, nfaces, iface

  reclen = size(record)
  nfaces = size(nbr_natural)
  buf(pos+1) = petsc_id
  buf(pos+2) = natural_id
  buf(pos+3) = reclen
  buf(pos+4:pos+3+reclen) = record
  pos = pos + 3 + reclen + 1
  buf(pos) = nfaces
  do iface = 1, nfaces
    buf(pos+2*iface-1) = nbr_natural(iface)
    buf(pos+2*iface) = nbr_face(iface)
  enddo
  pos = pos + 2*nfaces

end subroutine PolyPackMessage

! ************************************************************************** !

function PolyMessageNumNeighbors(buf,pos)
  !
  ! Number of faces with a neighbor in the message starting at buf(pos).
  !
  implicit none

  PetscInt, intent(in) :: buf(:)
  PetscInt, intent(in) :: pos
  PetscInt :: PolyMessageNumNeighbors

  PetscInt :: reclen, nfaces, iface, j

  reclen = buf(pos+2)
  j = pos + 3 + reclen
  nfaces = buf(j)
  PolyMessageNumNeighbors = 0
  do iface = 1, nfaces
    if (buf(j+2*iface-1) > 0) &
      PolyMessageNumNeighbors = PolyMessageNumNeighbors + 1
  enddo

end function PolyMessageNumNeighbors

! ************************************************************************** !

subroutine PolyMessageAppendNeighbors(buf,pos,list,count)
  !
  ! Appends the neighbor natural ids of the message at buf(pos) to list.
  !
  implicit none

  PetscInt, intent(in) :: buf(:)
  PetscInt, intent(in) :: pos
  PetscInt :: list(:)
  PetscInt :: count

  PetscInt :: reclen, nfaces, iface, j

  reclen = buf(pos+2)
  j = pos + 3 + reclen
  nfaces = buf(j)
  do iface = 1, nfaces
    if (buf(j+2*iface-1) > 0) then
      count = count + 1
      list(count) = buf(j+2*iface-1)
    endif
  enddo

end subroutine PolyMessageAppendNeighbors

! ************************************************************************** !

subroutine PolyMessageAt(ghosted_id,nlmax,owned_pos,ghost_pos,pos)
  !
  ! Position, in the owned or ghost buffer, of a ghosted cell's message
  ! minus one (so fields are at pos+1, pos+2, ...).
  !
  implicit none

  PetscInt, intent(in) :: ghosted_id
  PetscInt, intent(in) :: nlmax
  PetscInt, intent(in) :: owned_pos(:)
  PetscInt, intent(in) :: ghost_pos(:)
  PetscInt, intent(out) :: pos

  if (ghosted_id <= nlmax) then
    pos = owned_pos(ghosted_id) - 1
  else
    pos = ghost_pos(ghosted_id-nlmax) - 1
  endif

end subroutine PolyMessageAt

! ************************************************************************** !

function PolyBuf(ghosted_id,nlmax,owned_buf,ghost_buf,i)
  !
  ! Entry i of the owned (ghosted_id <= nlmax) or ghost message buffer.
  !
  implicit none

  PetscInt, intent(in) :: ghosted_id
  PetscInt, intent(in) :: nlmax
  PetscInt, intent(in) :: owned_buf(:)
  PetscInt, intent(in) :: ghost_buf(:)
  PetscInt, intent(in) :: i
  PetscInt :: PolyBuf

  if (ghosted_id <= nlmax) then
    PolyBuf = owned_buf(i)
  else
    PolyBuf = ghost_buf(i)
  endif

end function PolyBuf

! ************************************************************************** !

subroutine PolyMessageCheck(ghosted_id,nlmax,owned_buf,ghost_buf,pos, &
                            num_vertices_global,nfaces,nface_nodes)
  !
  ! Face and face-node counts of the record in a cell message.
  !
  implicit none

  PetscInt, intent(in) :: ghosted_id
  PetscInt, intent(in) :: nlmax
  PetscInt, intent(in) :: owned_buf(:)
  PetscInt, intent(in) :: ghost_buf(:)
  PetscInt, intent(in) :: pos
  PetscInt, intent(in) :: num_vertices_global
  PetscInt, intent(out) :: nfaces
  PetscInt, intent(out) :: nface_nodes

  PetscInt :: reclen, status, bad_position

  if (ghosted_id <= nlmax) then
    reclen = owned_buf(pos+3)
    call XDMFRecordCheck(owned_buf(pos+4:pos+3+reclen),num_vertices_global, &
                         nfaces,nface_nodes,status,bad_position)
  else
    reclen = ghost_buf(pos+3)
    call XDMFRecordCheck(ghost_buf(pos+4:pos+3+reclen),num_vertices_global, &
                         nfaces,nface_nodes,status,bad_position)
  endif

end subroutine PolyMessageCheck

! ************************************************************************** !
! Distribution helpers
! ************************************************************************** !

subroutine PolyBlockRange(n_global,nranks,rank,first,count)
  !
  ! Contiguous block of the 1-based ids 1..n_global assigned to a rank:
  ! count = q + (1 if rank < s), first = 1 + rank*q + min(rank,s), with
  ! q = n_global/nranks and s = mod(n_global,nranks).
  !
  implicit none

  PetscInt64, intent(in) :: n_global
  PetscInt, intent(in) :: nranks
  PetscInt, intent(in) :: rank
  PetscInt64, intent(out) :: first
  PetscInt64, intent(out) :: count

  PetscInt64 :: q, s, r

  q = n_global/int(nranks,kind=8)
  s = mod(n_global,int(nranks,kind=8))
  r = int(rank,kind=8)
  count = q
  if (r < s) count = count + 1
  first = 1 + r*q + min(r,s)

end subroutine PolyBlockRange

! ************************************************************************** !

function PolyBlockOwner(id,n_global,nranks)
  !
  ! Zero-based rank whose PolyBlockRange() block contains the 1-based id.
  !
  implicit none

  PetscInt64, intent(in) :: id
  PetscInt64, intent(in) :: n_global
  PetscInt, intent(in) :: nranks
  PetscInt :: PolyBlockOwner

  PetscInt64 :: q, s

  q = n_global/int(nranks,kind=8)
  s = mod(n_global,int(nranks,kind=8))
  if (id <= s*(q+1)) then
    PolyBlockOwner = int((id-1)/(q+1))
  else
    PolyBlockOwner = int(s + (id-1-s*(q+1))/q)
  endif

end function PolyBlockOwner

! ************************************************************************** !

subroutine PolyRankOffsets(nlocal,rank_offset,option)
  !
  ! rank_offset(r+1) is the 0-based global offset of rank r's cells;
  ! rank_offset(nranks+1) is the global count. O(number of ranks).
  !
  use Option_module

  implicit none

  PetscInt, intent(in) :: nlocal
  PetscInt, intent(out) :: rank_offset(:)
  type(option_type) :: option

  PetscInt :: i, nranks
  PetscErrorCode :: ierr

  nranks = option%comm%size
  call MPI_Allgather(nlocal,ONE_INTEGER_MPI,MPIU_INTEGER,rank_offset(2:), &
                     ONE_INTEGER_MPI,MPIU_INTEGER,option%mycomm, &
                     ierr);CHKERRQ(ierr)
  rank_offset(1) = 0
  do i = 2, nranks+1
    rank_offset(i) = rank_offset(i-1) + rank_offset(i)
  enddo

end subroutine PolyRankOffsets

! ************************************************************************** !

function PolyOwnerFromOffsets(petsc_id,rank_offset)
  !
  ! 1-based rank index owning a 0-based PETSc id (binary search).
  !
  implicit none

  PetscInt, intent(in) :: petsc_id
  PetscInt, intent(in) :: rank_offset(:)
  PetscInt :: PolyOwnerFromOffsets

  PetscInt :: lo, hi, mid

  lo = 1
  hi = size(rank_offset) - 1
  do while (lo < hi)
    mid = (lo + hi + 1)/2
    if (rank_offset(mid) <= petsc_id) then
      lo = mid
    else
      hi = mid - 1
    endif
  enddo
  PolyOwnerFromOffsets = lo

end function PolyOwnerFromOffsets

! ************************************************************************** !

function PolyKeyOwner(key,nranks)
  !
  ! Zero-based directory rank for a sorted face key.
  !
  implicit none

  PetscInt, intent(in) :: key(:)
  PetscInt, intent(in) :: nranks
  PetscInt :: PolyKeyOwner

  PetscInt64, parameter :: PRIME = 2147483647_8
  PetscInt64 :: h
  PetscInt :: i

  h = int(size(key),kind=8)
  do i = 1, size(key)
    h = mod(h*31_8 + mod(int(key(i),kind=8),PRIME),PRIME)
  enddo
  PolyKeyOwner = int(mod(h,int(nranks,kind=8)))

end function PolyKeyOwner

! ************************************************************************** !

subroutine PolyExchangeInt(send_buf,send_count,recv_buf,recv_count,label, &
                           option,report)
  !
  ! Personalized all-to-all exchange of PetscInt data. send_buf holds the
  ! data for rank 0, then rank 1, etc.; send_count(r+1) entries go to rank
  ! r. recv_buf is allocated (at least one entry) and ordered by source.
  !
  use Option_module

  implicit none

  PetscInt, intent(in) :: send_buf(:)
  PetscInt, intent(in) :: send_count(:)
  PetscInt, allocatable, intent(out) :: recv_buf(:)
  PetscInt, allocatable, intent(out) :: recv_count(:)
  character(len=*) :: label
  type(option_type) :: option
  PetscBool, optional :: report

  PetscMPIInt, allocatable :: scount(:), sdispl(:), rcount(:), rdispl(:)
  PetscInt :: nranks, ntotal
  PetscBool :: report_
  PetscErrorCode :: ierr

  report_ = PETSC_TRUE
  if (present(report)) report_ = report

  nranks = option%comm%size
  allocate(scount(nranks),sdispl(nranks),rcount(nranks),rdispl(nranks))
  call PolyExchangeCounts(send_count,scount,sdispl,rcount,rdispl, &
                          ntotal,label,option)
  allocate(recv_count(nranks))
  recv_count = rcount
  allocate(recv_buf(max(ntotal,1)))
  call MPI_Alltoallv(send_buf,scount,sdispl,MPIU_INTEGER, &
                     recv_buf,rcount,rdispl,MPIU_INTEGER, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  if (report_) then
    call PolyReportTransfer(label,int(sum(send_count),kind=8)*8, &
                            int(ntotal,kind=8)*8,option)
  endif
  deallocate(scount,sdispl,rcount,rdispl)

end subroutine PolyExchangeInt

! ************************************************************************** !

subroutine PolyExchangeReal(send_buf,send_count,recv_buf,recv_count,label, &
                            option)
  !
  ! As PolyExchangeInt() for PetscReal data. recv_count is input: the
  ! number of reals expected from each rank.
  !
  use Option_module

  implicit none

  PetscReal, intent(in) :: send_buf(:)
  PetscInt, intent(in) :: send_count(:)
  PetscReal, allocatable, intent(out) :: recv_buf(:)
  PetscInt, intent(in) :: recv_count(:)
  character(len=*) :: label
  type(option_type) :: option

  PetscMPIInt, allocatable :: scount(:), sdispl(:), rcount(:), rdispl(:)
  PetscInt :: nranks, ntotal
  PetscErrorCode :: ierr

  nranks = option%comm%size
  allocate(scount(nranks),sdispl(nranks),rcount(nranks),rdispl(nranks))
  call PolyExchangeCounts(send_count,scount,sdispl,rcount,rdispl, &
                          ntotal,label,option)
  if (any(rcount /= recv_count)) then
    option%io_buffer = 'Mismatched counts in ' // trim(label) // '.'
    call PrintErrMsgByRank(option)
  endif
  allocate(recv_buf(max(ntotal,1)))
  call MPI_Alltoallv(send_buf,scount,sdispl,MPI_DOUBLE_PRECISION, &
                     recv_buf,rcount,rdispl,MPI_DOUBLE_PRECISION, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  call PolyReportTransfer(label,int(sum(send_count),kind=8)*8, &
                          int(ntotal,kind=8)*8,option)
  deallocate(scount,sdispl,rcount,rdispl)

end subroutine PolyExchangeReal

! ************************************************************************** !

subroutine PolyExchangeCounts(send_count,scount,sdispl,rcount,rdispl, &
                              ntotal,label,option)
  !
  ! Exchanges per-rank counts and computes MPI displacements, checking that
  ! the totals fit in MPI's 32-bit counts on every rank (collectively).
  !
  use Option_module

  implicit none

  PetscInt, intent(in) :: send_count(:)
  PetscMPIInt, intent(out) :: scount(:), sdispl(:), rcount(:), rdispl(:)
  PetscInt, intent(out) :: ntotal
  character(len=*) :: label
  type(option_type) :: option

  PetscInt :: i, nranks
  PetscInt64 :: stotal, rtotal
  PetscMPIInt :: iflag
  PetscErrorCode :: ierr

  nranks = option%comm%size
  stotal = 0
  iflag = 0
  do i = 1, nranks
    stotal = stotal + send_count(i)
  enddo
  if (stotal > MAX_MPI_COUNT) iflag = 1
  call MPI_Allreduce(MPI_IN_PLACE,iflag,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  if (iflag > 0) then
    option%io_buffer = 'Polyhedra exchange "' // trim(label) // '" exceeds &
      &2^31-1 entries on a process; use more processes.'
    call PrintErrMsg(option)
  endif
  scount = int(send_count)
  call MPI_Alltoall(scount,ONE_INTEGER_MPI,MPI_INTEGER,rcount, &
                    ONE_INTEGER_MPI,MPI_INTEGER,option%mycomm, &
                    ierr);CHKERRQ(ierr)
  rtotal = 0
  do i = 1, nranks
    rtotal = rtotal + rcount(i)
  enddo
  iflag = 0
  if (rtotal > MAX_MPI_COUNT) iflag = 1
  call MPI_Allreduce(MPI_IN_PLACE,iflag,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  if (iflag > 0) then
    option%io_buffer = 'Polyhedra exchange "' // trim(label) // '" exceeds &
      &2^31-1 entries on a process; use more processes.'
    call PrintErrMsg(option)
  endif
  sdispl(1) = 0
  rdispl(1) = 0
  do i = 2, nranks
    sdispl(i) = sdispl(i-1) + scount(i-1)
    rdispl(i) = rdispl(i-1) + rcount(i-1)
  enddo
  ntotal = int(rtotal)

end subroutine PolyExchangeCounts

! ************************************************************************** !

subroutine PolyCheckErrors(nerr,message,option)
  !
  ! Collective error checkpoint: stops on all ranks if any rank reports an
  ! error. Rank-specific details must be printed beforehand.
  !
  use Option_module

  implicit none

  PetscInt, intent(in) :: nerr
  character(len=*) :: message
  type(option_type) :: option

  PetscInt :: nerr_global
  PetscErrorCode :: ierr

  call MPI_Allreduce(nerr,nerr_global,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                     option%mycomm,ierr);CHKERRQ(ierr)
  if (nerr_global > 0) then
    option%io_buffer = message
    call PrintErrMsg(option)
  endif

end subroutine PolyCheckErrors

! ************************************************************************** !

subroutine PolyReportTransfer(label,bytes_sent,bytes_received,option)
  !
  ! Prints the maximum per-rank volume of an exchange.
  !
  use Option_module
  use String_module

  implicit none

  character(len=*) :: label
  PetscInt64, intent(in) :: bytes_sent
  PetscInt64, intent(in) :: bytes_received
  type(option_type) :: option

  PetscInt64 :: values(2)
  PetscErrorCode :: ierr

  values = [bytes_sent,bytes_received]
  call MPI_Allreduce(MPI_IN_PLACE,values,TWO_INTEGER_MPI,MPI_INTEGER8, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  option%io_buffer = '  ' // trim(label) // ': max bytes sent/received &
    &per process: ' // trim(PolyIntString(values(1))) // '/' // &
    trim(PolyIntString(values(2)))
  call PrintMsg(option)

end subroutine PolyReportTransfer

! ************************************************************************** !

subroutine PolyReportMemory(label,bytes,option)
  !
  ! Prints the maximum and total over ranks of a per-rank storage size.
  !
  use Option_module
  use String_module

  implicit none

  character(len=*) :: label
  PetscInt64, intent(in) :: bytes
  type(option_type) :: option

  PetscInt64 :: max_bytes, sum_bytes
  PetscErrorCode :: ierr

  call MPI_Allreduce(bytes,max_bytes,ONE_INTEGER_MPI,MPI_INTEGER8, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(bytes,sum_bytes,ONE_INTEGER_MPI,MPI_INTEGER8, &
                     MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  option%io_buffer = '  ' // trim(label) // ': max per process ' // &
    trim(PolyIntString(max_bytes)) // ', total ' // &
    trim(PolyIntString(sum_bytes))
  call PrintMsg(option)

end subroutine PolyReportMemory

! ************************************************************************** !

function PolyIntString(i)
  !
  ! Formats a 64-bit integer.
  !
  implicit none

  PetscInt64, intent(in) :: i
  character(len=32) :: PolyIntString

  write(PolyIntString,'(i0)') i

end function PolyIntString

! ************************************************************************** !
! Wrappers for PETSc sorting/searching, whose lengths are PetscCount
! ************************************************************************** !

subroutine PolySortInt(n,a,ierr)

  implicit none

  PetscInt, intent(in) :: n
  PetscInt :: a(:)
  PetscErrorCode :: ierr

  PetscCount :: count

  count = n
  call PetscSortInt(count,a,ierr);CHKERRQ(ierr)

end subroutine PolySortInt

! ************************************************************************** !

subroutine PolySortIntWithArray(n,a,b,ierr)

  implicit none

  PetscInt, intent(in) :: n
  PetscInt :: a(:)
  PetscInt :: b(:)
  PetscErrorCode :: ierr

  PetscCount :: count

  count = n
  call PetscSortIntWithArray(count,a,b,ierr);CHKERRQ(ierr)

end subroutine PolySortIntWithArray

! ************************************************************************** !

subroutine PolyFindInt(key,n,a,loc,ierr)
  !
  ! Binary search of a sorted array; loc is 0-based, negative if absent.
  !
  implicit none

  PetscInt, intent(in) :: key
  PetscInt, intent(in) :: n
  PetscInt, intent(in) :: a(:)
  PetscInt, intent(out) :: loc
  PetscErrorCode :: ierr

  PetscCount :: count

  count = n
  call PetscFindInt(key,count,a,loc,ierr);CHKERRQ(ierr)

end subroutine PolyFindInt

! ************************************************************************** !
! HDF5 helpers
! ************************************************************************** !

function PolyH5LinkExists(loc_id,name)

  implicit none

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  PetscBool :: PolyH5LinkExists

  logical :: exists
  integer :: hdf5_err

  call h5lexists_f(loc_id,name,exists,hdf5_err)
  PolyH5LinkExists = (hdf5_err >= 0 .and. exists)

end function PolyH5LinkExists

! ************************************************************************** !

subroutine PolyH5OpenDataset(file_id,path,expected_rank,expected_class, &
                             expected_size,filename,dset_id,dims,option)
  !
  ! Opens a dataset and checks its rank, type class, and (if
  ! expected_size > 0) its type size in bytes. dims are in Fortran order.
  !
  use Option_module
  use String_module

  implicit none

  integer(HID_T), intent(in) :: file_id
  character(len=*), intent(in) :: path
  integer, intent(in) :: expected_rank
  integer, intent(in) :: expected_class
  integer, intent(in) :: expected_size
  character(len=*), intent(in) :: filename
  integer(HID_T), intent(out) :: dset_id
  integer(HSIZE_T), intent(out) :: dims(2)
  type(option_type) :: option

  integer(HID_T) :: space_id, type_id
  integer(HSIZE_T) :: maxdims(2)
  integer(SIZE_T) :: type_size
  integer :: ndims, type_class, hdf5_err

  if (.not.PolyH5LinkExists(file_id,path)) then
    option%io_buffer = 'Dataset /' // trim(path) // ' not found in "' // &
      trim(filename) // '".'
    call PrintErrMsg(option)
  endif
  call h5dopen_f(file_id,path,dset_id,hdf5_err)
  call h5dget_space_f(dset_id,space_id,hdf5_err)
  call h5sget_simple_extent_ndims_f(space_id,ndims,hdf5_err)
  if (ndims /= expected_rank) then
    option%io_buffer = 'Dataset /' // trim(path) // ' in "' // &
      trim(filename) // '" must have rank ' // &
      trim(PolyIntString(int(expected_rank,kind=8))) // '.'
    call PrintErrMsg(option)
  endif
  dims = 0
  call h5sget_simple_extent_dims_f(space_id,dims(1:ndims),maxdims(1:ndims), &
                                   hdf5_err)
  call h5sclose_f(space_id,hdf5_err)
  call h5dget_type_f(dset_id,type_id,hdf5_err)
  call h5tget_class_f(type_id,type_class,hdf5_err)
  call h5tget_size_f(type_id,type_size,hdf5_err)
  call h5tclose_f(type_id,hdf5_err)
  if (type_class /= expected_class .or. &
      (expected_size > 0 .and. type_size /= expected_size)) then
    if (expected_class == H5T_INTEGER_F) then
      option%io_buffer = 'Dataset /' // trim(path) // ' in "' // &
        trim(filename) // '" must hold 64-bit integers.'
    else
      option%io_buffer = 'Dataset /' // trim(path) // ' in "' // &
        trim(filename) // '" must hold floating point values.'
    endif
    call PrintErrMsg(option)
  endif

end subroutine PolyH5OpenDataset

! ************************************************************************** !

subroutine PolyH5ReadInt64(dset_id,start,count,buf,xfer_id,hdf5_err)
  !
  ! Collective read of count entries of a rank-1 integer dataset starting
  ! at the 0-based position start. count may be zero (empty selection).
  !
  implicit none

  integer(HID_T), intent(in) :: dset_id
  PetscInt64, intent(in) :: start
  PetscInt64, intent(in) :: count
  PetscInt64, target :: buf(:)
  integer(HID_T), intent(in) :: xfer_id
  integer, intent(out) :: hdf5_err

  integer(HID_T) :: file_space_id, mem_space_id
  integer(HSIZE_T) :: hstart(1), hcount(1)
  type(c_ptr) :: f_ptr
  integer :: err

  call h5dget_space_f(dset_id,file_space_id,err)
  if (count > 0) then
    hstart(1) = int(start,kind=HSIZE_T)
    hcount(1) = int(count,kind=HSIZE_T)
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,hstart, &
                               hcount,err)
    call h5screate_simple_f(1,hcount,mem_space_id,err)
  else
    call h5sselect_none_f(file_space_id,err)
    hcount(1) = 1
    call h5screate_simple_f(1,hcount,mem_space_id,err)
    call h5sselect_none_f(mem_space_id,err)
  endif
  f_ptr = c_loc(buf(1))
  call h5dread_f(dset_id,h5kind_to_type(8,H5_INTEGER_KIND),f_ptr,hdf5_err, &
                 mem_space_id,file_space_id,xfer_id)
  call h5sclose_f(mem_space_id,err)
  call h5sclose_f(file_space_id,err)

end subroutine PolyH5ReadInt64

! ************************************************************************** !

subroutine PolyH5ReadVertices(dset_id,start,count,buf,xfer_id,hdf5_err)
  !
  ! Collective read of count vertex rows (0-based start row) of the [Nv,3]
  ! vertex dataset into buf(3,count). count may be zero.
  !
  implicit none

  integer(HID_T), intent(in) :: dset_id
  PetscInt64, intent(in) :: start
  PetscInt64, intent(in) :: count
  PetscReal, target :: buf(:,:)
  integer(HID_T), intent(in) :: xfer_id
  integer, intent(out) :: hdf5_err

  integer(HID_T) :: file_space_id, mem_space_id
  integer(HSIZE_T) :: hstart(2), hcount(2)
  type(c_ptr) :: f_ptr
  integer :: err

  call h5dget_space_f(dset_id,file_space_id,err)
  if (count > 0) then
    hstart = [0_HSIZE_T,int(start,kind=HSIZE_T)]
    hcount = [3_HSIZE_T,int(count,kind=HSIZE_T)]
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,hstart, &
                               hcount,err)
    call h5screate_simple_f(2,hcount,mem_space_id,err)
  else
    call h5sselect_none_f(file_space_id,err)
    hcount = [3_HSIZE_T,1_HSIZE_T]
    call h5screate_simple_f(2,hcount,mem_space_id,err)
    call h5sselect_none_f(mem_space_id,err)
  endif
  f_ptr = c_loc(buf(1,1))
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,f_ptr,hdf5_err,mem_space_id, &
                 file_space_id,xfer_id)
  call h5sclose_f(mem_space_id,err)
  call h5sclose_f(file_space_id,err)

end subroutine PolyH5ReadVertices

! ************************************************************************** !

subroutine PolyH5WriteInt64(loc_id,name,total,start,count,buf,option)
  !
  ! Collectively creates (or opens) a rank-1 int64 dataset of length total
  ! and writes count entries at the 0-based position start.
  !
  use Option_module

  implicit none

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  PetscInt64, intent(in) :: total
  PetscInt64, intent(in) :: start
  PetscInt64, intent(in) :: count
  PetscInt64, target :: buf(:)
  type(option_type) :: option

  integer(HID_T) :: dset_id, file_space_id, mem_space_id, xfer_id
  integer(HSIZE_T) :: dims(1), hstart(1), hcount(1)
  type(c_ptr) :: f_ptr
  integer :: hdf5_err, err

  dims(1) = int(total,kind=HSIZE_T)
  if (PolyH5LinkExists(loc_id,name)) then
    call h5dopen_f(loc_id,name,dset_id,err)
  else
    call h5screate_simple_f(1,dims,file_space_id,err)
    call h5dcreate_f(loc_id,name,H5T_STD_I64LE,file_space_id,dset_id,err)
    call h5sclose_f(file_space_id,err)
  endif
  call h5dget_space_f(dset_id,file_space_id,err)
  if (count > 0) then
    hstart(1) = int(start,kind=HSIZE_T)
    hcount(1) = int(count,kind=HSIZE_T)
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,hstart, &
                               hcount,err)
    call h5screate_simple_f(1,hcount,mem_space_id,err)
  else
    call h5sselect_none_f(file_space_id,err)
    hcount(1) = 1
    call h5screate_simple_f(1,hcount,mem_space_id,err)
    call h5sselect_none_f(mem_space_id,err)
  endif
  call h5pcreate_f(H5P_DATASET_XFER_F,xfer_id,err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(xfer_id,H5FD_MPIO_COLLECTIVE_F,err)
#endif
  f_ptr = c_loc(buf(1))
  call h5dwrite_f(dset_id,h5kind_to_type(8,H5_INTEGER_KIND),f_ptr, &
                  hdf5_err,mem_space_id,file_space_id,xfer_id)
  call h5pclose_f(xfer_id,err)
  call h5sclose_f(mem_space_id,err)
  call h5sclose_f(file_space_id,err)
  call h5dclose_f(dset_id,err)
  call PolyCheckHDF5Error(hdf5_err,name,'output file',option)

end subroutine PolyH5WriteInt64

! ************************************************************************** !

subroutine PolyH5WriteVertices(loc_id,name,total,start,count,buf,option)
  !
  ! Collectively creates (or opens) a float64 [total,3] dataset and writes
  ! count rows (buf(3,count)) at the 0-based row start.
  !
  use Option_module

  implicit none

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  PetscInt64, intent(in) :: total
  PetscInt64, intent(in) :: start
  PetscInt64, intent(in) :: count
  PetscReal, target :: buf(:,:)
  type(option_type) :: option

  integer(HID_T) :: dset_id, file_space_id, mem_space_id, xfer_id
  integer(HSIZE_T) :: dims(2), hstart(2), hcount(2)
  type(c_ptr) :: f_ptr
  integer :: hdf5_err, err

  dims = [3_HSIZE_T,int(total,kind=HSIZE_T)]
  if (PolyH5LinkExists(loc_id,name)) then
    call h5dopen_f(loc_id,name,dset_id,err)
  else
    call h5screate_simple_f(2,dims,file_space_id,err)
    call h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,file_space_id,dset_id,err)
    call h5sclose_f(file_space_id,err)
  endif
  call h5dget_space_f(dset_id,file_space_id,err)
  if (count > 0) then
    hstart = [0_HSIZE_T,int(start,kind=HSIZE_T)]
    hcount = [3_HSIZE_T,int(count,kind=HSIZE_T)]
    call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,hstart, &
                               hcount,err)
    call h5screate_simple_f(2,hcount,mem_space_id,err)
  else
    call h5sselect_none_f(file_space_id,err)
    hcount = [3_HSIZE_T,1_HSIZE_T]
    call h5screate_simple_f(2,hcount,mem_space_id,err)
    call h5sselect_none_f(mem_space_id,err)
  endif
  call h5pcreate_f(H5P_DATASET_XFER_F,xfer_id,err)
#ifndef SERIAL_HDF5
  call h5pset_dxpl_mpio_f(xfer_id,H5FD_MPIO_COLLECTIVE_F,err)
#endif
  f_ptr = c_loc(buf(1,1))
  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,f_ptr,hdf5_err,mem_space_id, &
                  file_space_id,xfer_id)
  call h5pclose_f(xfer_id,err)
  call h5sclose_f(mem_space_id,err)
  call h5sclose_f(file_space_id,err)
  call h5dclose_f(dset_id,err)
  call PolyCheckHDF5Error(hdf5_err,name,'output file',option)

end subroutine PolyH5WriteVertices

! ************************************************************************** !

subroutine PolyCheckHDF5Error(hdf5_err,dataset_name,filename,option)
  !
  ! Collective check of a dataset read.
  !
  use Option_module

  implicit none

  integer, intent(in) :: hdf5_err
  character(len=*) :: dataset_name
  character(len=*) :: filename
  type(option_type) :: option

  PetscInt :: nerr

  nerr = 0
  if (hdf5_err < 0) nerr = 1
  call PolyCheckErrors(nerr,'Error reading or writing dataset "' // &
                       trim(dataset_name) // '" from "' // trim(filename) // &
                       '".',option)

end subroutine PolyCheckHDF5Error

! ************************************************************************** !

subroutine PolyH5CheckIntAttribute(loc_id,name,expected,filename,option)
  !
  ! Requires a scalar integer attribute with the expected value.
  !
  use Option_module
  use String_module

  implicit none

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  PetscInt64, intent(in) :: expected
  character(len=*), intent(in) :: filename
  type(option_type) :: option

  integer(HID_T) :: attr_id, type_id, space_id
  integer(HSSIZE_T) :: npoints
  PetscInt64, target :: value(1)
  type(c_ptr) :: f_ptr
  logical :: exists
  integer :: type_class, hdf5_err
  PetscBool :: valid

  valid = PETSC_FALSE
  value = -huge(value(1))
  call h5aexists_f(loc_id,name,exists,hdf5_err)
  if (exists) then
    call h5aopen_f(loc_id,name,attr_id,hdf5_err)
    call h5aget_type_f(attr_id,type_id,hdf5_err)
    call h5tget_class_f(type_id,type_class,hdf5_err)
    call h5tclose_f(type_id,hdf5_err)
    call h5aget_space_f(attr_id,space_id,hdf5_err)
    call h5sget_simple_extent_npoints_f(space_id,npoints,hdf5_err)
    call h5sclose_f(space_id,hdf5_err)
    if (type_class == H5T_INTEGER_F .and. npoints == 1) then
      f_ptr = c_loc(value(1))
      call h5aread_f(attr_id,h5kind_to_type(8,H5_INTEGER_KIND),f_ptr, &
                     hdf5_err)
      valid = (hdf5_err >= 0 .and. value(1) == expected)
    endif
    call h5aclose_f(attr_id,hdf5_err)
  endif
  if (.not.valid) then
    option%io_buffer = 'Attribute /Domain/CellsPolyhedra/' // trim(name) // &
      ' must be the integer ' // trim(PolyIntString(expected)) // ' in "' // &
      trim(filename) // '"'
    if (exists .and. value(1) /= -huge(value(1))) then
      option%io_buffer = trim(option%io_buffer) // ' (found ' // &
        trim(PolyIntString(value(1))) // ')'
    endif
    option%io_buffer = trim(option%io_buffer) // '.'
    call PrintErrMsg(option)
  endif

end subroutine PolyH5CheckIntAttribute

! ************************************************************************** !

subroutine PolyH5ReadStringAttribute(loc_id,name,value,found)
  !
  ! Reads a scalar string attribute stored as either a fixed-length or a
  ! variable-length string (h5py writes the latter by default).
  !
  implicit none

  integer(HID_T), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: value
  PetscBool, intent(out) :: found

  integer(HID_T) :: attr_id, type_id, mem_type_id, space_id
  integer(SIZE_T) :: type_size
  type(c_ptr), target :: vl_ptr(1)
  type(c_ptr) :: f_ptr
  character(kind=c_char), pointer :: chars(:)
  character(kind=c_char), allocatable, target :: fixed(:)
  logical :: exists, is_variable
  integer :: type_class, hdf5_err, i

  value = ''
  found = PETSC_FALSE
  call h5aexists_f(loc_id,name,exists,hdf5_err)
  if (.not.exists) return
  call h5aopen_f(loc_id,name,attr_id,hdf5_err)
  call h5aget_type_f(attr_id,type_id,hdf5_err)
  call h5tget_class_f(type_id,type_class,hdf5_err)
  if (type_class == H5T_STRING_F) then
    call h5tis_variable_str_f(type_id,is_variable,hdf5_err)
    call h5tcopy_f(type_id,mem_type_id,hdf5_err)
    if (is_variable) then
      vl_ptr(1) = c_null_ptr
      f_ptr = c_loc(vl_ptr(1))
      call h5aread_f(attr_id,mem_type_id,f_ptr,hdf5_err)
      if (hdf5_err >= 0 .and. c_associated(vl_ptr(1))) then
        call c_f_pointer(vl_ptr(1),chars,[len(value)])
        do i = 1, len(value)
          if (chars(i) == c_null_char) exit
          value(i:i) = chars(i)
        enddo
        found = PETSC_TRUE
        call h5aget_space_f(attr_id,space_id,hdf5_err)
        call h5dvlen_reclaim_f(mem_type_id,space_id,H5P_DEFAULT_F, &
                               f_ptr,hdf5_err)
        call h5sclose_f(space_id,hdf5_err)
      endif
    else
      call h5tget_size_f(type_id,type_size,hdf5_err)
      allocate(fixed(type_size))
      fixed = c_null_char
      f_ptr = c_loc(fixed(1))
      call h5aread_f(attr_id,mem_type_id,f_ptr,hdf5_err)
      if (hdf5_err >= 0) then
        do i = 1, min(int(type_size),len(value))
          if (fixed(i) == c_null_char) exit
          value(i:i) = fixed(i)
        enddo
        found = PETSC_TRUE
      endif
      deallocate(fixed)
    endif
    call h5tclose_f(mem_type_id,hdf5_err)
  endif
  call h5tclose_f(type_id,hdf5_err)
  call h5aclose_f(attr_id,hdf5_err)

end subroutine PolyH5ReadStringAttribute

end module Grid_Unstructured_Poly_HDF5_module
