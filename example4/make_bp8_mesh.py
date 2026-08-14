#!/usr/bin/env python3
"""
Generates a flat, uniform triangular mesh for the SEAS BP8 whole-space
benchmark and writes it in the (non-standard, project-local) GTS-like
format consumed by NikkhooWalter2015/calc_nikkhoo_fs.f90 and
NikkhooWalter2015/calc_nikkhoo.f90: a header line "NV NE NC" (NE is
unused/always 0 in this codebase's mesh files -- confirmed against
example2/triangular_mesh.gts, whose face section lists vertex-index
triples directly rather than the edge-index triples a real GTS file
would use), NV lines of "x y z", then NC lines of 1-based vertex-index
triples.

The mesh is flat at z=0 by construction: calc_nikkhoo_fs.f90's local-frame
convention (see calc_local_coordinate's "horizontal triangle" branch)
gives an unambiguous strike/dip/normal frame a1=(1,0,0), a2=(0,1,0),
a3=(0,0,1) for z=0 triangles, so mesh x maps to BP8's fault-plane
coordinate x2 (strike) and mesh y maps to x3 (dip); the whole-space fault
normal is BP8's x1.

Triangle winding is kept consistent (normal = +z for every cell, via
cross(p2-p1, p3-p1)) so calc_nikkhoo_fs.f90 does not need a winding-order
correction pass.

Usage:
    python3 make_bp8_mesh.py --domain 1200 --cell-size 10 --out triangular_mesh.gts
"""
import argparse
import numpy as np


def build_mesh(domain_size, cell_size):
    half = domain_size / 2.0
    n = int(round(domain_size / cell_size))
    if abs(n * cell_size - domain_size) > 1e-6:
        raise ValueError(f"domain_size ({domain_size}) must be an integer multiple of cell_size ({cell_size})")

    coords = -half + cell_size * np.arange(n + 1)  # n+1 grid lines, exact multiples of cell_size
    nv_per_side = n + 1

    vertices = np.zeros(((n + 1) * (n + 1), 3))
    for j in range(n + 1):        # y index
        for i in range(n + 1):    # x index
            vid = j * nv_per_side + i
            vertices[vid, 0] = coords[i]
            vertices[vid, 1] = coords[j]
            vertices[vid, 2] = 0.0

    def vid(i, j):
        return j * nv_per_side + i + 1  # 1-based

    faces = []
    for j in range(n):
        for i in range(n):
            v00 = vid(i, j)
            v10 = vid(i + 1, j)
            v01 = vid(i, j + 1)
            v11 = vid(i + 1, j + 1)
            # Both triangles wound counter-clockwise when viewed from +z,
            # i.e. cross(p2-p1, p3-p1) points +z for both.
            faces.append((v00, v10, v11))
            faces.append((v00, v11, v01))

    return vertices, faces


def write_gts(path, vertices, faces):
    with open(path, 'w') as f:
        f.write(f"{len(vertices)} 0 {len(faces)}\n")
        for x, y, z in vertices:
            f.write(f"{x:15.7e} {y:15.7e} {z:15.7e}\n")
        for a, b, c in faces:
            f.write(f"{a} {b} {c}\n")


def check_winding(vertices, faces, n_check=50):
    import random
    bad = 0
    sample = faces if len(faces) <= n_check else random.sample(faces, n_check)
    for a, b, c in sample:
        p1, p2, p3 = vertices[a-1], vertices[b-1], vertices[c-1]
        nv = np.cross(p2 - p1, p3 - p1)
        if nv[2] <= 0:
            bad += 1
    return bad


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--domain', type=float, default=1200.0, help='square domain side length (m)')
    ap.add_argument('--cell-size', type=float, default=10.0, help='uniform cell size (m)')
    ap.add_argument('--out', type=str, default='triangular_mesh.gts')
    args = ap.parse_args()

    vertices, faces = build_mesh(args.domain, args.cell_size)
    write_gts(args.out, vertices, faces)

    bad = check_winding(vertices, faces)
    n_side = int(round(args.domain / args.cell_size))
    print(f"Wrote {args.out}: {len(vertices)} vertices, {len(faces)} faces "
          f"({n_side}x{n_side} cells, {args.cell_size} m each, {args.domain} m domain)")
    print(f"Winding check: {bad}/50 sampled faces have non-positive z-normal (expect 0)")
