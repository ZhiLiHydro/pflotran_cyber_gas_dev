module Grid_Unstructured_XDMF_module
  !
  ! XDMF mixed-topology cell records shared by polyhedra mesh input and
  ! output. A record is an XDMF cell code followed by its payload:
  !
  !   6  (tetrahedron):  4 node ids
  !   7  (pyramid):      5 node ids
  !   8  (wedge):        6 node ids
  !   9  (hexahedron):   8 node ids
  !   16 (polyhedron):   nfaces, then for each face its node count and its
  !                      node ids in cyclic order (first node not repeated)
  !
  ! These codes are neither Exodus element codes nor PFLOTRAN's cell-type
  ! enumerations (HEX_TYPE, etc.). Node ids are 1-based.
  !
  ! Conventional cells are expanded into faces with the templates of
  ! UCellGetFaceVertices() (grid_unstructured_cell.F90), so a cell's local
  ! face order is:
  !   tetrahedron: (1,2,4) (2,3,4) (1,4,3) (1,3,2)
  !   pyramid:     (1,2,5) (2,3,5) (3,4,5) (4,1,5) (1,4,3,2)
  !   wedge:       (1,2,5,4) (2,3,6,5) (3,1,4,6) (1,3,2) (4,5,6)
  !   hexahedron:  (1,2,6,5) (2,3,7,6) (3,4,8,7) (4,1,5,8) (1,4,3,2) (5,6,7,8)
  ! All template faces wind outward. For the wedge, nodes 1-2-3 wind toward
  ! nodes 4-5-6 (the Exodus convention used by PFLOTRAN's XDMF writer); a
  ! wedge in the opposite (VTK) convention is inverted and is rejected.
  !
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: XDMF_CODE_TETRAHEDRON = 6
  PetscInt, parameter, public :: XDMF_CODE_PYRAMID = 7
  PetscInt, parameter, public :: XDMF_CODE_WEDGE = 8
  PetscInt, parameter, public :: XDMF_CODE_HEXAHEDRON = 9
  PetscInt, parameter, public :: XDMF_CODE_POLYHEDRON = 16

  PetscInt, parameter, public :: XDMF_POLYHEDRON_MIN_FACES = 4
  PetscInt, parameter, public :: XDMF_POLYGON_MIN_NODES = 3

  ! status returned by XDMFRecordCheck()
  PetscInt, parameter, public :: XDMF_RECORD_OK = 0
  PetscInt, parameter, public :: XDMF_RECORD_EMPTY = 1
  PetscInt, parameter, public :: XDMF_RECORD_UNSUPPORTED_CODE = 2
  PetscInt, parameter, public :: XDMF_RECORD_TRUNCATED = 3
  PetscInt, parameter, public :: XDMF_RECORD_TRAILING_DATA = 4
  PetscInt, parameter, public :: XDMF_RECORD_TOO_FEW_FACES = 5
  PetscInt, parameter, public :: XDMF_RECORD_TOO_FEW_NODES = 6
  PetscInt, parameter, public :: XDMF_RECORD_NODE_OUT_OF_RANGE = 7
  PetscInt, parameter, public :: XDMF_RECORD_REPEATED_NODE = 8
  PetscInt, parameter, public :: XDMF_RECORD_NOT_CLOSED = 9

  public :: XDMFIsConventionalCode, &
            XDMFConventionalNumNodes, &
            XDMFRecordCheck, &
            XDMFRecordExpand, &
            XDMFRecordStatusString, &
            XDMFFaceIsReversed, &
            XDMFCompareKeys, &
            XDMFSortKeys, &
            XDMFPolygonGeometry, &
            XDMFPolyhedronGeometry

contains

! ************************************************************************** !

function XDMFIsConventionalCode(code)
  !
  ! Returns true for the fixed-template XDMF codes 6 through 9.
  !
  implicit none

  PetscInt, intent(in) :: code
  PetscBool :: XDMFIsConventionalCode

  select case(code)
    case(XDMF_CODE_TETRAHEDRON,XDMF_CODE_PYRAMID,XDMF_CODE_WEDGE, &
         XDMF_CODE_HEXAHEDRON)
      XDMFIsConventionalCode = PETSC_TRUE
    case default
      XDMFIsConventionalCode = PETSC_FALSE
  end select

end function XDMFIsConventionalCode

! ************************************************************************** !

function XDMFConventionalNumNodes(code)
  !
  ! Number of node ids following a conventional code; 0 otherwise.
  !
  implicit none

  PetscInt, intent(in) :: code
  PetscInt :: XDMFConventionalNumNodes

  select case(code)
    case(XDMF_CODE_TETRAHEDRON)
      XDMFConventionalNumNodes = 4
    case(XDMF_CODE_PYRAMID)
      XDMFConventionalNumNodes = 5
    case(XDMF_CODE_WEDGE)
      XDMFConventionalNumNodes = 6
    case(XDMF_CODE_HEXAHEDRON)
      XDMFConventionalNumNodes = 8
    case default
      XDMFConventionalNumNodes = 0
  end select

end function XDMFConventionalNumNodes

! ************************************************************************** !

subroutine XDMFConventionalTemplate(code,nfaces,face_offset,face_nodes)
  !
  ! Returns the face template of a conventional cell: face iface uses the
  ! cell-local node positions face_nodes(face_offset(iface):
  ! face_offset(iface+1)-1). See the module header for the ordering.
  !
  implicit none

  PetscInt, intent(in) :: code
  PetscInt, intent(out) :: nfaces
  PetscInt, intent(out) :: face_offset(7)
  PetscInt, intent(out) :: face_nodes(24)

  face_offset = 0
  face_nodes = 0
  select case(code)
    case(XDMF_CODE_TETRAHEDRON)
      nfaces = 4
      face_offset(1:5) = [1,4,7,10,13]
      face_nodes(1:12) = [1,2,4, 2,3,4, 1,4,3, 1,3,2]
    case(XDMF_CODE_PYRAMID)
      nfaces = 5
      face_offset(1:6) = [1,4,7,10,13,17]
      face_nodes(1:16) = [1,2,5, 2,3,5, 3,4,5, 4,1,5, 1,4,3,2]
    case(XDMF_CODE_WEDGE)
      nfaces = 5
      face_offset(1:6) = [1,5,9,13,16,19]
      face_nodes(1:18) = [1,2,5,4, 2,3,6,5, 3,1,4,6, 1,3,2, 4,5,6]
    case(XDMF_CODE_HEXAHEDRON)
      nfaces = 6
      face_offset(1:7) = [1,5,9,13,17,21,25]
      face_nodes(1:24) = [1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,4,3,2, &
                          5,6,7,8]
    case default
      nfaces = 0
  end select

end subroutine XDMFConventionalTemplate

! ************************************************************************** !

subroutine XDMFRecordCheck(record,num_vertices_global,nfaces,nface_nodes, &
                           status,bad_position)
  !
  ! Validates one XDMF cell record without expanding it. The record must be
  ! consumed exactly. On success, returns the number of faces and the total
  ! number of face-node entries (sum over faces of the face node count).
  ! Checks node ids in [1,num_vertices_global], at least 4 faces per code 16
  ! record, at least 3 nodes per face, no repeated node within a face, and
  ! a closed, consistently oriented face shell (every directed edge appears
  ! once and its reverse appears once).
  !
  ! bad_position is the 1-based position within the record at which the
  ! problem was detected (0 if not applicable).
  !
  implicit none

  PetscInt, intent(in) :: record(:)
  PetscInt, intent(in) :: num_vertices_global
  PetscInt, intent(out) :: nfaces
  PetscInt, intent(out) :: nface_nodes
  PetscInt, intent(out) :: status
  PetscInt, intent(out) :: bad_position

  PetscInt :: n, code, nnode, iface, pos, nv, i, j
  PetscInt :: tnfaces, toffset(7), tnodes(24)

  n = size(record)
  nfaces = 0
  nface_nodes = 0
  status = XDMF_RECORD_OK
  bad_position = 0

  if (n < 1) then
    status = XDMF_RECORD_EMPTY
    return
  endif
  code = record(1)
  if (XDMFIsConventionalCode(code)) then
    nnode = XDMFConventionalNumNodes(code)
    if (n < nnode + 1) then
      status = XDMF_RECORD_TRUNCATED
      bad_position = n
      return
    else if (n > nnode + 1) then
      status = XDMF_RECORD_TRAILING_DATA
      bad_position = nnode + 2
      return
    endif
    do i = 2, n
      if (record(i) < 1 .or. record(i) > num_vertices_global) then
        status = XDMF_RECORD_NODE_OUT_OF_RANGE
        bad_position = i
        return
      endif
      do j = 2, i-1
        if (record(j) == record(i)) then
          status = XDMF_RECORD_REPEATED_NODE
          bad_position = i
          return
        endif
      enddo
    enddo
    call XDMFConventionalTemplate(code,tnfaces,toffset,tnodes)
    nfaces = tnfaces
    nface_nodes = toffset(tnfaces+1) - 1
  else if (code == XDMF_CODE_POLYHEDRON) then
    if (n < 2) then
      status = XDMF_RECORD_TRUNCATED
      bad_position = n
      return
    endif
    nfaces = record(2)
    if (nfaces < XDMF_POLYHEDRON_MIN_FACES) then
      status = XDMF_RECORD_TOO_FEW_FACES
      bad_position = 2
      return
    endif
    pos = 2
    do iface = 1, nfaces
      pos = pos + 1
      if (pos > n) then
        status = XDMF_RECORD_TRUNCATED
        bad_position = n
        return
      endif
      nv = record(pos)
      if (nv < XDMF_POLYGON_MIN_NODES) then
        status = XDMF_RECORD_TOO_FEW_NODES
        bad_position = pos
        return
      endif
      if (nv > n - pos) then
        status = XDMF_RECORD_TRUNCATED
        bad_position = n
        return
      endif
      do i = pos+1, pos+nv
        if (record(i) < 1 .or. record(i) > num_vertices_global) then
          status = XDMF_RECORD_NODE_OUT_OF_RANGE
          bad_position = i
          return
        endif
        do j = pos+1, i-1
          if (record(j) == record(i)) then
            status = XDMF_RECORD_REPEATED_NODE
            bad_position = i
            return
          endif
        enddo
      enddo
      nface_nodes = nface_nodes + nv
      pos = pos + nv
    enddo
    if (pos /= n) then
      status = XDMF_RECORD_TRAILING_DATA
      bad_position = pos + 1
      return
    endif
  else
    status = XDMF_RECORD_UNSUPPORTED_CODE
    bad_position = 1
    return
  endif

  if (.not.XDMFRecordIsClosed(record,nfaces,nface_nodes)) then
    status = XDMF_RECORD_NOT_CLOSED
    bad_position = 0
    nfaces = 0
    nface_nodes = 0
  endif

end subroutine XDMFRecordCheck

! ************************************************************************** !

function XDMFRecordIsClosed(record,nfaces,nface_nodes)
  !
  ! Returns true if every directed face edge of a (valid) record occurs once
  ! and its reverse occurs once, i.e. the faces form a closed shell with
  ! consistent winding in which each edge is shared by exactly two faces.
  !
  implicit none

  PetscInt, intent(in) :: record(:)
  PetscInt, intent(in) :: nfaces
  PetscInt, intent(in) :: nface_nodes
  PetscBool :: XDMFRecordIsClosed

  PetscInt :: face_offset(nfaces+1)
  PetscInt :: face_nodes(nface_nodes)
  PetscInt :: edge_a(nface_nodes), edge_b(nface_nodes)
  PetscInt :: iface, i, j, nedge, nv, ibeg
  PetscBool :: found

  call XDMFRecordExpand(record,face_offset,face_nodes)
  nedge = 0
  do iface = 1, nfaces
    ibeg = face_offset(iface)
    nv = face_offset(iface+1) - ibeg
    do i = 0, nv-1
      nedge = nedge + 1
      edge_a(nedge) = face_nodes(ibeg+i)
      edge_b(nedge) = face_nodes(ibeg+mod(i+1,nv))
    enddo
  enddo

  XDMFRecordIsClosed = PETSC_FALSE
  do i = 1, nedge
    found = PETSC_FALSE
    do j = 1, nedge
      if (j /= i .and. edge_a(j) == edge_a(i) .and. &
          edge_b(j) == edge_b(i)) return
      if (edge_a(j) == edge_b(i) .and. edge_b(j) == edge_a(i)) then
        if (found) return
        found = PETSC_TRUE
      endif
    enddo
    if (.not.found) return
  enddo
  XDMFRecordIsClosed = PETSC_TRUE

end function XDMFRecordIsClosed

! ************************************************************************** !

subroutine XDMFRecordExpand(record,face_offset,face_nodes)
  !
  ! Expands a record validated by XDMFRecordCheck() into faces: face iface
  ! has nodes face_nodes(face_offset(iface):face_offset(iface+1)-1), in the
  ! cyclic order given by the record (code 16) or the template (codes 6-9).
  ! face_offset must hold nfaces+1 entries and face_nodes nface_nodes.
  !
  implicit none

  PetscInt, intent(in) :: record(:)
  PetscInt, intent(out) :: face_offset(:)
  PetscInt, intent(out) :: face_nodes(:)

  PetscInt :: code, nfaces, iface, pos, nv, count, i
  PetscInt :: toffset(7), tnodes(24)

  code = record(1)
  if (code == XDMF_CODE_POLYHEDRON) then
    nfaces = record(2)
    pos = 2
    count = 0
    face_offset(1) = 1
    do iface = 1, nfaces
      pos = pos + 1
      nv = record(pos)
      face_nodes(count+1:count+nv) = record(pos+1:pos+nv)
      count = count + nv
      face_offset(iface+1) = count + 1
      pos = pos + nv
    enddo
  else
    call XDMFConventionalTemplate(code,nfaces,toffset,tnodes)
    face_offset(1:nfaces+1) = toffset(1:nfaces+1)
    do i = 1, toffset(nfaces+1)-1
      ! record(1) is the code, so cell-local node k is record(k+1)
      face_nodes(i) = record(tnodes(i)+1)
    enddo
  endif

end subroutine XDMFRecordExpand

! ************************************************************************** !

function XDMFRecordStatusString(status)
  !
  ! Describes a status returned by XDMFRecordCheck().
  !
  implicit none

  PetscInt, intent(in) :: status
  character(len=MAXWORDLENGTH*2) :: XDMFRecordStatusString

  select case(status)
    case(XDMF_RECORD_OK)
      XDMFRecordStatusString = 'valid record'
    case(XDMF_RECORD_EMPTY)
      XDMFRecordStatusString = 'empty record'
    case(XDMF_RECORD_UNSUPPORTED_CODE)
      XDMFRecordStatusString = 'unsupported XDMF cell code (supported: &
        &6, 7, 8, 9, 16)'
    case(XDMF_RECORD_TRUNCATED)
      XDMFRecordStatusString = 'record is truncated'
    case(XDMF_RECORD_TRAILING_DATA)
      XDMFRecordStatusString = 'record has trailing data'
    case(XDMF_RECORD_TOO_FEW_FACES)
      XDMFRecordStatusString = 'code 16 record has fewer than 4 faces'
    case(XDMF_RECORD_TOO_FEW_NODES)
      XDMFRecordStatusString = 'face has fewer than 3 nodes'
    case(XDMF_RECORD_NODE_OUT_OF_RANGE)
      XDMFRecordStatusString = 'node id outside [1,number of vertices]'
    case(XDMF_RECORD_REPEATED_NODE)
      XDMFRecordStatusString = 'node repeated within a face or cell'
    case(XDMF_RECORD_NOT_CLOSED)
      XDMFRecordStatusString = 'faces do not form a closed, consistently &
        &wound shell'
    case default
      XDMFRecordStatusString = 'unknown record status'
  end select

end function XDMFRecordStatusString

! ************************************************************************** !

function XDMFFaceIsReversed(a,b)
  !
  ! Returns true if face b lists the nodes of face a in reverse cyclic
  ! order (as the two cells sharing a face must). Both faces must have the
  ! same length and node set.
  !
  implicit none

  PetscInt, intent(in) :: a(:)
  PetscInt, intent(in) :: b(:)
  PetscBool :: XDMFFaceIsReversed

  PetscInt :: n, i, j, k

  XDMFFaceIsReversed = PETSC_FALSE
  n = size(a)
  if (size(b) /= n) return
  j = 0
  do i = 1, n
    if (b(i) == a(1)) then
      j = i
      exit
    endif
  enddo
  if (j == 0) return
  do k = 0, n-1
    if (a(1+k) /= b(modulo(j-1-k,n)+1)) return
  enddo
  XDMFFaceIsReversed = PETSC_TRUE

end function XDMFFaceIsReversed

! ************************************************************************** !

function XDMFCompareKeys(a,b)
  !
  ! Orders sorted face keys: first by length, then lexicographically.
  ! Returns -1, 0, or 1.
  !
  implicit none

  PetscInt, intent(in) :: a(:)
  PetscInt, intent(in) :: b(:)
  PetscInt :: XDMFCompareKeys

  PetscInt :: i

  XDMFCompareKeys = 0
  if (size(a) /= size(b)) then
    XDMFCompareKeys = merge(-1,1,size(a) < size(b))
    return
  endif
  do i = 1, size(a)
    if (a(i) /= b(i)) then
      XDMFCompareKeys = merge(-1,1,a(i) < b(i))
      return
    endif
  enddo

end function XDMFCompareKeys

! ************************************************************************** !

subroutine XDMFSortKeys(nkeys,key_offset,keys,perm)
  !
  ! Stable merge sort of variable-length keys stored in CSR form (key i is
  ! keys(key_offset(i):key_offset(i+1)-1)). On return, perm lists the key
  ! indices in ascending XDMFCompareKeys() order.
  !
  implicit none

  PetscInt, intent(in) :: nkeys
  PetscInt, intent(in) :: key_offset(:)
  PetscInt, intent(in) :: keys(:)
  PetscInt, intent(out) :: perm(:)

  PetscInt, allocatable :: work(:)
  PetscInt :: width, lo, mid, hi, i, j, k

  do i = 1, nkeys
    perm(i) = i
  enddo
  if (nkeys < 2) return
  allocate(work(nkeys))
  width = 1
  do while (width < nkeys)
    lo = 1
    do while (lo <= nkeys)
      mid = min(lo+width,nkeys+1)
      hi = min(lo+2*width,nkeys+1)
      i = lo
      j = mid
      k = lo
      do while (i < mid .and. j < hi)
        if (XDMFCompareKeys(keys(key_offset(perm(j)): &
                                 key_offset(perm(j)+1)-1), &
                            keys(key_offset(perm(i)): &
                                 key_offset(perm(i)+1)-1)) < 0) then
          work(k) = perm(j)
          j = j + 1
        else
          work(k) = perm(i)
          i = i + 1
        endif
        k = k + 1
      enddo
      do while (i < mid)
        work(k) = perm(i)
        i = i + 1
        k = k + 1
      enddo
      do while (j < hi)
        work(k) = perm(j)
        j = j + 1
        k = k + 1
      enddo
      lo = hi
    enddo
    perm(1:nkeys) = work(1:nkeys)
    width = 2*width
  enddo
  deallocate(work)

end subroutine XDMFSortKeys

! ************************************************************************** !

subroutine XDMFPolygonGeometry(xyz,area,centroid,unit_normal,nonplanarity)
  !
  ! Geometry of a polygon face with nodes xyz(3,n) in cyclic order.
  !
  ! The face is triangulated as a fan about its node average c:
  ! triangles (c,p_i,p_i+1). The fan depends only on the cyclic node order,
  ! not on the starting node or the winding direction, so the two copies of
  ! a shared face produce the same triangles. The area vector is the sum of
  ! the triangle area vectors, which is exact for any planar simple polygon
  ! (convex or not) and, for a warped face, is the area vector of the
  ! triangulated surface. area is its magnitude (the projected area that
  ! carries a uniform flux), unit_normal its direction (outward for an
  ! outward-wound face), and centroid the area-weighted centroid using
  ! triangle areas signed by unit_normal. nonplanarity is the maximum node
  ! distance from the plane through c normal to unit_normal, divided by
  ! sqrt(area). Coordinates are shifted by c to limit cancellation for
  ! large map coordinates.
  !
  ! Returns area = 0 for a degenerate face.
  !
  implicit none

  PetscReal, intent(in) :: xyz(:,:)
  PetscReal, intent(out) :: area
  PetscReal, intent(out) :: centroid(3)
  PetscReal, intent(out) :: unit_normal(3)
  PetscReal, intent(out) :: nonplanarity

  PetscReal :: c(3), p(3), q(3), tri(3,size(xyz,2)), area_vec(3)
  PetscReal :: w, wsum
  PetscInt :: n, i

  n = size(xyz,2)
  c = sum(xyz,dim=2)/dble(n)
  area_vec = 0.d0
  do i = 1, n
    p = xyz(:,i) - c
    q = xyz(:,mod(i,n)+1) - c
    tri(:,i) = 0.5d0*Cross(p,q)
    area_vec = area_vec + tri(:,i)
  enddo
  area = sqrt(dot_product(area_vec,area_vec))
  centroid = c
  unit_normal = 0.d0
  nonplanarity = 0.d0
  if (area <= 0.d0) return
  unit_normal = area_vec/area

  centroid = 0.d0
  wsum = 0.d0
  do i = 1, n
    p = xyz(:,i) - c
    q = xyz(:,mod(i,n)+1) - c
    w = dot_product(tri(:,i),unit_normal)
    centroid = centroid + w*(p+q)/3.d0
    wsum = wsum + w
  enddo
  centroid = c + centroid/wsum

  do i = 1, n
    nonplanarity = max(nonplanarity, &
                       abs(dot_product(xyz(:,i)-c,unit_normal)))
  enddo
  nonplanarity = nonplanarity/sqrt(area)

end subroutine XDMFPolygonGeometry

! ************************************************************************** !

subroutine XDMFPolyhedronGeometry(xyz,face_offset,face_nodes,volume, &
                                  centroid)
  !
  ! Signed volume and volume-weighted centroid of a closed polyhedron whose
  ! faces reference the node coordinates xyz(3,:) through cell-local node
  ! indices: face iface uses xyz(:,face_nodes(face_offset(iface):
  ! face_offset(iface+1)-1)). Faces are triangulated exactly as in
  ! XDMFPolygonGeometry(), and the volume sums the signed tetrahedra formed
  ! by each triangle and the node average r of the cell (coordinates are
  ! shifted by r). The volume is positive for outward-wound faces; it is
  ! exact for planar faces whether or not the cell is convex.
  !
  implicit none

  PetscReal, intent(in) :: xyz(:,:)
  PetscInt, intent(in) :: face_offset(:)
  PetscInt, intent(in) :: face_nodes(:)
  PetscReal, intent(out) :: volume
  PetscReal, intent(out) :: centroid(3)

  PetscReal :: r(3), c(3), p(3), q(3), v
  PetscInt :: nfaces, iface, ibeg, nv, i

  nfaces = size(face_offset) - 1
  r = sum(xyz,dim=2)/dble(size(xyz,2))
  volume = 0.d0
  centroid = 0.d0
  do iface = 1, nfaces
    ibeg = face_offset(iface)
    nv = face_offset(iface+1) - ibeg
    c = 0.d0
    do i = 0, nv-1
      c = c + xyz(:,face_nodes(ibeg+i))
    enddo
    c = c/dble(nv) - r
    do i = 0, nv-1
      p = xyz(:,face_nodes(ibeg+i)) - r
      q = xyz(:,face_nodes(ibeg+mod(i+1,nv))) - r
      v = dot_product(c,Cross(p,q))/6.d0
      volume = volume + v
      ! tetrahedron (0,c,p,q) has centroid (c+p+q)/4 relative to r
      centroid = centroid + v*(c+p+q)*0.25d0
    enddo
  enddo
  if (volume /= 0.d0) then
    centroid = r + centroid/volume
  else
    centroid = r
  endif

end subroutine XDMFPolyhedronGeometry

! ************************************************************************** !

function Cross(a,b)

  implicit none

  PetscReal, intent(in) :: a(3), b(3)
  PetscReal :: Cross(3)

  Cross(1) = a(2)*b(3) - a(3)*b(2)
  Cross(2) = a(3)*b(1) - a(1)*b(3)
  Cross(3) = a(1)*b(2) - a(2)*b(1)

end function Cross

end module Grid_Unstructured_XDMF_module
