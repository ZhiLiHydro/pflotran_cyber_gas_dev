#!/usr/bin/env python3
"""
Generates HDF5 polyhedra meshes in the PFLOTRAN mixed XDMF-code schema:

  /Domain/Vertices                               float64 [Nv, 3]
  /Domain/CellsPolyhedra/ConnectivityWithXDMFCodes int64 [L]
  /Domain/CellsPolyhedra/Offsets                   int64 [Nc + 1]

with attributes on /Domain/CellsPolyhedra:

  SchemaVersion = 1, CellTypeConvention = "XDMF",
  NodeIndexBase = 1, OffsetIndexBase = 1

Node IDs and offsets are 1-based. Cell records use XDMF codes 6 (tet),
7 (pyramid), 8 (wedge), 9 (hex), and 16 (arbitrary polyhedron:
16, nfaces, n1, ids..., n2, ids..., ...). Conventional cells use the node
order of PFLOTRAN's cell templates (grid_unstructured_cell.F90); for the
wedge, nodes 1-2-3 wind toward nodes 4-5-6 (Exodus convention). No
geometry (volumes, centroids, areas) is stored; PFLOTRAN computes it.

Each file also holds a boundary region /Regions/top defined by
"Cell Ids" (1-based natural cell ids) and "Face Ids" (1-based cell-local
face indices): the faces whose vertices all lie on the top of the domain.
Cell-local faces are numbered in template order for codes 6-9 (TEMPLATES
below, matching UCellGetFaceVertices) and in record order for code 16.
The same faces are also written as the face set /Regions/top_vertices:
"FaceVertexIDs" (1-based node ids, faces concatenated) and "FaceOffsets"
(1-based, Nf + 1 entries), with SchemaVersion, NodeIndexBase, and
OffsetIndexBase attributes. Winding and starting node are not significant;
every other face is written reversed and rotated to exercise that.

Files written (all describe outward-wound, closed cells):

  mixed_polyhedra.h5        mixed.ugp with conventional cells stored
                            compactly (codes 6-9) where the faces match a
                            template, and code 16 otherwise
  mixed_polyhedra_code16.h5 mixed.ugp with every cell stored as code 16
  pentagonal_prisms.h5      two layers of two pentagonal prisms; cells share
                            pentagonal (between layers) and quadrilateral
                            (within a layer) faces

Usage: python generate_polyhedra_h5_grid.py
"""
import itertools
import sys

import h5py
import numpy as np

XDMF_TET = 6
XDMF_PYRAMID = 7
XDMF_WEDGE = 8
XDMF_HEX = 9
XDMF_POLYHEDRON = 16

# Cell-local faces (1-based node positions, outward winding) matching
# UCellGetFaceVertices() in grid_unstructured_cell.F90.
TEMPLATES = {
    XDMF_TET: (4, [(1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2)]),
    XDMF_PYRAMID: (5, [(1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5),
                       (1, 4, 3, 2)]),
    XDMF_WEDGE: (6, [(1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6), (1, 3, 2),
                     (4, 5, 6)]),
    XDMF_HEX: (8, [(1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8),
                   (1, 4, 3, 2), (5, 6, 7, 8)]),
}


def canonical_cycle(face):
    """Rotation of a cyclic face list that starts at its smallest node."""
    i = face.index(min(face))
    return tuple(face[i:] + face[:i])


def face_area_vector(xyz):
    c = xyz.mean(axis=0)
    a = np.zeros(3)
    n = len(xyz)
    for i in range(n):
        a += 0.5 * np.cross(xyz[i] - c, xyz[(i + 1) % n] - c)
    return a


def cell_volume(faces, vertices):
    """Signed volume from 1-based face node lists; positive if outward."""
    ids = sorted({v for f in faces for v in f})
    r = vertices[[i - 1 for i in ids]].mean(axis=0)
    vol = 0.0
    for f in faces:
        xyz = vertices[[i - 1 for i in f]] - r
        c = xyz.mean(axis=0)
        n = len(f)
        for i in range(n):
            vol += np.dot(c, np.cross(xyz[i], xyz[(i + 1) % n])) / 6.0
    return vol


def check_closed_and_oriented(faces):
    """Each directed edge appears once and its reverse appears once."""
    edges = set()
    for f in faces:
        n = len(f)
        for i in range(n):
            e = (f[i], f[(i + 1) % n])
            if e in edges:
                return False
            edges.add(e)
    return all((b, a) in edges for (a, b) in edges)


def match_template(faces):
    """Return (code, ordered nodes) if faces match a conventional template."""
    nodes = sorted({v for f in faces for v in f})
    target = sorted(canonical_cycle(list(f)) for f in faces)
    for code, (nnode, tfaces) in TEMPLATES.items():
        if nnode != len(nodes) or len(tfaces) != len(faces):
            continue
        for perm in itertools.permutations(nodes):
            mapped = sorted(canonical_cycle([perm[i - 1] for i in tf])
                            for tf in tfaces)
            if mapped == target:
                return code, list(perm)
    return None


def polyhedron_record(faces):
    rec = [XDMF_POLYHEDRON, len(faces)]
    for f in faces:
        rec.append(len(f))
        rec.extend(f)
    return rec


def record_faces(rec):
    """Faces of a record (node ids) in PFLOTRAN's cell-local face order."""
    code = rec[0]
    if code == XDMF_POLYHEDRON:
        faces = []
        p = 2
        for _ in range(rec[1]):
            n = rec[p]
            faces.append(list(rec[p + 1:p + 1 + n]))
            p += 1 + n
        return faces
    return [[rec[i] for i in f] for f in TEMPLATES[code][1]]


def top_faces(records, vertices):
    """(cell id, face index) pairs of faces lying on the domain top."""
    ztop = vertices[:, 2].max()
    cell_ids = []
    face_ids = []
    for icell, rec in enumerate(records, 1):
        for iface, face in enumerate(record_faces(rec), 1):
            if np.all(np.isclose(vertices[[v - 1 for v in face], 2], ztop)):
                cell_ids.append(icell)
                face_ids.append(iface)
    return cell_ids, face_ids


def write_mesh(filename, vertices, records):
    offsets = [1]
    conn = []
    for rec in records:
        conn.extend(rec)
        offsets.append(offsets[-1] + len(rec))
    with h5py.File(filename, 'w') as h5:
        dom = h5.create_group('Domain')
        dom.create_dataset('Vertices', data=np.asarray(vertices, 'f8'))
        grp = dom.create_group('CellsPolyhedra')
        grp.create_dataset('ConnectivityWithXDMFCodes',
                           data=np.asarray(conn, 'i8'))
        grp.create_dataset('Offsets', data=np.asarray(offsets, 'i8'))
        grp.attrs['SchemaVersion'] = np.int64(1)
        grp.attrs['CellTypeConvention'] = 'XDMF'
        grp.attrs['NodeIndexBase'] = np.int64(1)
        grp.attrs['OffsetIndexBase'] = np.int64(1)
        cell_ids, face_ids = top_faces(records, np.asarray(vertices))
        top = h5.create_group('Regions').create_group('top')
        top.create_dataset('Cell Ids', data=np.asarray(cell_ids, 'i8'))
        top.create_dataset('Face Ids', data=np.asarray(face_ids, 'i8'))
        ids = []
        offsets = [1]
        for n, (icell, iface) in enumerate(zip(cell_ids, face_ids)):
            face = record_faces(records[icell - 1])[iface - 1]
            if n % 2 == 1:
                face = list(reversed(face))
                face = face[1:] + face[:1]
            ids.extend(face)
            offsets.append(offsets[-1] + len(face))
        fset = h5['Regions'].create_group('top_vertices')
        fset.create_dataset('FaceVertexIDs', data=np.asarray(ids, 'i8'))
        fset.create_dataset('FaceOffsets', data=np.asarray(offsets, 'i8'))
        fset.attrs['SchemaVersion'] = np.int64(1)
        fset.attrs['NodeIndexBase'] = np.int64(1)
        fset.attrs['OffsetIndexBase'] = np.int64(1)
    print('Wrote {}: {} cells, {} vertices, {} connectivity integers, '
          '{} top faces'.format(filename, len(records), len(vertices),
                                len(conn), len(cell_ids)))


def validate(cells, vertices):
    for icell, faces in enumerate(cells, 1):
        if not check_closed_and_oriented(faces):
            sys.exit('Cell {} is not closed/consistently oriented'.format(icell))
        vol = cell_volume(faces, vertices)
        if vol <= 0.0:
            sys.exit('Cell {} has non-positive volume {}'.format(icell, vol))


def read_ugp(filename):
    with open(filename) as f:
        lines = [l.split() for l in f if l.strip()]
    pos = 0
    ncell = int(lines[pos][1])
    pos += 1
    nfaces = []
    for i in range(ncell):
        nfaces.append(int(lines[pos + i][1]))
    pos += ncell
    nface_total = int(lines[pos][1])
    pos += 1
    cells = [[] for _ in range(ncell)]
    for i in range(nface_total):
        tok = lines[pos + i]
        cell = int(tok[1])
        nv = int(tok[2])
        cells[cell - 1].append([int(v) for v in tok[3:3 + nv]])
    pos += nface_total
    nvert = int(lines[pos][1])
    pos += 1
    vertices = np.array([[float(x) for x in lines[pos + i][:3]]
                         for i in range(nvert)])
    for i, faces in enumerate(cells):
        assert len(faces) == nfaces[i]
    return cells, vertices


def mixed_meshes():
    cells, vertices = read_ugp('mixed.ugp')
    validate(cells, vertices)
    compact = []
    counts = {}
    for faces in cells:
        match = match_template(faces)
        if match is None:
            compact.append(polyhedron_record(faces))
            counts[XDMF_POLYHEDRON] = counts.get(XDMF_POLYHEDRON, 0) + 1
        else:
            code, nodes = match
            compact.append([code] + nodes)
            counts[code] = counts.get(code, 0) + 1
    print('mixed.ugp cell codes: {}'.format(counts))
    write_mesh('mixed_polyhedra.h5', vertices, compact)
    write_mesh('mixed_polyhedra_code16.h5', vertices,
               [polyhedron_record(faces) for faces in cells])


def pentagonal_prisms():
    """
    Pentagonal prisms tiled in x and stacked in z. Each layer holds two
    prisms whose pentagons are a unit square with one corner cut, mirrored
    about x = 1 so that the pair shares a quadrilateral side face. Layers
    share pentagonal faces.
    """
    # pentagon of the left prism (counterclockwise from above)
    left = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.3, 1.0), (0.0, 0.7)]
    # right prism mirrored about x = 1 (counterclockwise from above)
    right = [(1.0, 0.0), (2.0, 0.0), (2.0, 0.7), (1.7, 1.0), (1.0, 1.0)]
    nlayer = 2
    dz = 0.5
    xy = []
    index = {}

    def node(p, k):
        key = (p[0], p[1], k)
        if key not in index:
            xy.append((p[0], p[1], k * dz))
            index[key] = len(xy)
        return index[key]

    cells = []
    for k in range(nlayer):
        for pent in (left, right):
            bottom = [node(p, k) for p in pent]
            top = [node(p, k + 1) for p in pent]
            faces = [list(reversed(bottom)), list(top)]
            for i in range(5):
                j = (i + 1) % 5
                faces.append([bottom[i], bottom[j], top[j], top[i]])
            cells.append(faces)
    vertices = np.array(xy)
    validate(cells, vertices)
    write_mesh('pentagonal_prisms.h5', vertices,
               [polyhedron_record(faces) for faces in cells])


if __name__ == '__main__':
    mixed_meshes()
    pentagonal_prisms()
