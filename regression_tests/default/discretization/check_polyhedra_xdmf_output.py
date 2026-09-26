#!/usr/bin/env python3
"""
Checks PFLOTRAN HDF5/XDMF visualization output of polyhedra meshes read from
the HDF5 XDMF-code schema (see generate_polyhedra_h5_grid.py).

For each (output, input mesh) pair below whose output exists, verifies that

  - /Domain/Cells is int64 and equals the input connectivity with only the
    node ids converted to 0-based indices (codes, face counts, and face node
    counts unchanged; conventional records stay compact),
  - /Domain/Vertices equals the input vertices,
  - the "Natural ID" field of every snapshot is 1..Nc (cell fields are in
    the same order as the connectivity),
  - /Domain/XC, YC, ZC equal centroids computed from the written stream,
  - each .xmf file declares XDMF 3.0, a Mixed topology with Nc elements,
    and an 8-byte integer stream of the written length.

If face flow rates / face velocities were written, it also verifies that
each has one row per cell holding the face count followed by at least
that many per-face values (cell-local face order, zero padding), that the
flow rates of the two copies of every shared face are equal and opposite
and the face velocities equal, and that boundary faces carry flow only if
they belong to the boundary region /Regions/top.

Run as the python_post_process_script of the polyhedra regression tests.
Exits nonzero on the first failure.
"""
import glob
import re
import sys

import h5py
import numpy as np

CHECKS = [
    ('mixed_polyhedra_h5-np4', 'mixed_polyhedra.h5'),
    ('pentagonal_prisms_h5-np2', 'pentagonal_prisms.h5'),
    ('pentagonal_prisms_h5_flowrate-np2', 'pentagonal_prisms.h5'),
]

TEMPLATES = {
    6: [(1, 2, 4), (2, 3, 4), (1, 4, 3), (1, 3, 2)],
    7: [(1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 1, 5), (1, 4, 3, 2)],
    8: [(1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6), (1, 3, 2), (4, 5, 6)],
    9: [(1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8), (1, 4, 3, 2),
        (5, 6, 7, 8)],
}
NUM_NODES = {6: 4, 7: 5, 8: 6, 9: 8}


def parse(stream):
    """Split a mixed stream into (record, faces) per cell."""
    cells = []
    p = 0
    while p < len(stream):
        code = int(stream[p])
        if code == 16:
            q = p + 2
            faces = []
            for _ in range(int(stream[p + 1])):
                n = int(stream[q])
                faces.append([int(v) for v in stream[q + 1:q + 1 + n]])
                q += 1 + n
        elif code in NUM_NODES:
            q = p + 1 + NUM_NODES[code]
            nodes = [int(v) for v in stream[p + 1:q]]
            faces = [[nodes[i - 1] for i in f] for f in TEMPLATES[code]]
        else:
            raise ValueError('unsupported code {} at {}'.format(code, p))
        cells.append(([int(v) for v in stream[p:q]], faces))
        p = q
    return cells


def shift_nodes(record, shift):
    record = list(record)
    if record[0] == 16:
        q = 2
        for _ in range(record[1]):
            n = record[q]
            record[q + 1:q + 1 + n] = [v + shift for v in record[q + 1:q + 1 + n]]
            q += 1 + n
    else:
        record[1:] = [v + shift for v in record[1:]]
    return record


def centroid(faces, vertices):
    ids = sorted({v for f in faces for v in f})
    r = vertices[ids].mean(axis=0)
    vol = 0.0
    cen = np.zeros(3)
    for f in faces:
        x = vertices[f] - r
        c = x.mean(axis=0)
        for i in range(len(f)):
            p, q = x[i], x[(i + 1) % len(f)]
            v = np.dot(c, np.cross(p, q)) / 6.0
            vol += v
            cen += v * (c + p + q) / 4.0
    return r + cen / vol


def check_face_output(h5, cells, top_faces):
    """Checks face flow rates and face velocities, if present."""
    owners = {}
    for icell, (_, faces) in enumerate(cells):
        for iface, face in enumerate(faces):
            owners.setdefault(tuple(sorted(face)), []).append((icell, iface))
    shared = [v for v in owners.values() if len(v) == 2]
    boundary = [v[0] for v in owners.values() if len(v) == 1]
    nfaces = np.array([len(f) for _, f in cells])
    nchecked = 0
    for g in [g for g in h5 if g != 'Domain']:
        for name in h5[g]:
            if 'Flowrate' not in name and 'Velocity at cell face' not in name:
                continue
            d = h5[g][name][:]
            assert d.shape[0] == len(cells) and d.shape[1] >= nfaces.max() + 1, \
                '{}: shape {}'.format(name, d.shape)
            assert np.array_equal(d[:, 0], nfaces), name + ': face counts'
            q = d[:, 1:]
            scale = max(np.abs(q).max(), 1e-300)
            for i, n in enumerate(nfaces):
                assert not q[i, n:].any(), name + ': nonzero padding'
            sign = 1.0 if 'Flowrate' in name else -1.0
            err = max(abs(q[a] + sign * q[b]) for a, b in shared)
            assert err <= 1e-12 * scale, '{}: shared faces differ by {}' \
                .format(name, err)
            if 'Flowrate' in name:
                for c, f in boundary:
                    if (c + 1, f + 1) not in top_faces:
                        assert abs(q[c, f]) <= 1e-12 * scale, \
                            '{}: flow on boundary face {} of cell {}' \
                            .format(name, f + 1, c + 1)
            nchecked += 1
    return nchecked


def check(prefix, mesh):
    with h5py.File(mesh, 'r') as h5:
        vin = h5['Domain/Vertices'][:]
        cin = h5['Domain/CellsPolyhedra/ConnectivityWithXDMFCodes'][:]
        top_faces = set(zip(h5['Regions/top/Cell Ids'][:].tolist(),
                            h5['Regions/top/Face Ids'][:].tolist()))
    expected = [x for rec, _ in parse(cin) for x in shift_nodes(rec, -1)]

    with h5py.File(prefix + '.h5', 'r') as h5:
        dom = h5['Domain']
        stream = dom['Cells'][:]
        assert stream.dtype == np.int64, 'Cells is {}'.format(stream.dtype)
        assert np.array_equal(stream, expected), 'Cells stream differs'
        assert np.array_equal(dom['Vertices'][:], vin), 'Vertices differ'
        cells = parse(stream)
        ncell = len(cells)
        xyz = np.c_[dom['XC'][:], dom['YC'][:], dom['ZC'][:]]
        cen = np.array([centroid(f, dom['Vertices'][:]) for _, f in cells])
        err = np.abs(cen - xyz).max()
        assert err < 1e-10, 'cell centers differ by {}'.format(err)
        groups = [g for g in h5 if g != 'Domain']
        assert groups, 'no snapshots'
        for g in groups:
            names = [k for k in h5[g] if k.startswith('Natural ID')]
            assert names, 'no Natural ID in ' + g
            nat = np.asarray(h5[g][names[0]][:]).ravel()
            assert np.array_equal(nat, np.arange(1, ncell + 1)), \
                'field order differs in ' + g
        nface_checked = check_face_output(h5, cells, top_faces)

    xmfs = sorted(glob.glob(prefix + '-[0-9]*.xmf'))
    assert xmfs, 'no .xmf files'
    for xmf in xmfs:
        text = open(xmf).read()
        assert '<Xdmf Version="3.0">' in text, xmf + ': not XDMF 3.0'
        assert re.search(r'Topology Type="Mixed" NumberOfElements="{}"'
                         .format(ncell), text), xmf + ': topology'
        assert re.search(r'DataType="Int" Precision="8" Dimensions="{}"'
                         .format(len(stream)), text), xmf + ': stream item'
    print('{}: OK ({} cells, {} stream integers, {} snapshots, {} xmf, '
          '{} face datasets)'.format(prefix, ncell, len(stream), len(groups),
                                     len(xmfs), nface_checked))


if __name__ == '__main__':
    nchecked = 0
    for prefix, mesh in CHECKS:
        if glob.glob(prefix + '.h5'):
            check(prefix, mesh)
            nchecked += 1
    if nchecked == 0:
        sys.exit('No polyhedra XDMF output found to check.')
