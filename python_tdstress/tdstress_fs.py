"""
Triangular Dislocation Stress calculation in elastic full-space.

Translation from MATLAB code by Nikkhoo & Walter (2015).

Reference:
Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
artefact-free solution. Geophysical Journal International.
"""

import numpy as np
from .td_utils import coord_trans, tens_trans, trimodefinder
from .ang_dislocation import td_setup_s


def tdstress_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam):
    """
    Calculate stresses and strains for a triangular dislocation in elastic full-space.

    Parameters
    ----------
    X, Y, Z : array_like
        Coordinates of calculation points in EFCS (East, North, Up).
        Must have the same size.
    P1, P2, P3 : array_like, shape (3,)
        Coordinates of TD vertices in EFCS
    Ss, Ds, Ts : float
        TD slip vector components (Strike-slip, Dip-slip, Tensile-slip)
    mu, lam : float
        Lame constants

    Returns
    -------
    Stress : ndarray, shape (n, 6)
        Stress tensor components [Sxx, Syy, Szz, Sxy, Sxz, Syz]
    Strain : ndarray, shape (n, 6)
        Strain tensor components [Exx, Eyy, Ezz, Exy, Exz, Eyz]
    """
    # Convert to arrays and flatten
    X = np.atleast_1d(X).flatten()
    Y = np.atleast_1d(Y).flatten()
    Z = np.atleast_1d(Z).flatten()

    P1 = np.atleast_1d(P1).flatten()
    P2 = np.atleast_1d(P2).flatten()
    P3 = np.atleast_1d(P3).flatten()

    n_points = len(X)

    # Poisson's ratio
    nu = 1 / (1 + lam / mu) / 2

    # Burgers vector components
    bx = Ts  # Tensile-slip
    by = Ss  # Strike-slip
    bz = Ds  # Dip-slip

    # Calculate unit normal, strike, and dip vectors
    eY = np.array([0, 1, 0])
    eZ = np.array([0, 0, 1])

    Vnorm = np.cross(P2 - P1, P3 - P1)
    Vnorm = Vnorm / np.linalg.norm(Vnorm)

    Vstrike = np.cross(eZ, Vnorm)
    if np.linalg.norm(Vstrike) == 0:
        Vstrike = eY * Vnorm[2]
        # For horizontal elements (image dislocation case)
        if P1[2] > 0:
            Vstrike = -Vstrike
    Vstrike = Vstrike / np.linalg.norm(Vstrike)

    Vdip = np.cross(Vnorm, Vstrike)

    # Transformation matrix (rows are unit vectors)
    A = np.array([Vnorm, Vstrike, Vdip])

    # Transform coordinates from EFCS into TDCS
    p1 = np.zeros(3)
    p2 = np.zeros(3)
    p3 = np.zeros(3)

    x, y, z = coord_trans(X - P2[0], Y - P2[1], Z - P2[2], A)
    p1[0], p1[1], p1[2] = coord_trans(P1[0] - P2[0], P1[1] - P2[1], P1[2] - P2[2], A)
    p3[0], p3[1], p3[2] = coord_trans(P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2], A)

    # Calculate unit vectors along TD sides in TDCS
    e12 = (p2 - p1) / np.linalg.norm(p2 - p1)
    e13 = (p3 - p1) / np.linalg.norm(p3 - p1)
    e23 = (p3 - p2) / np.linalg.norm(p3 - p2)

    # Calculate TD angles
    A_angle = np.arccos(np.dot(e12, e13))
    B_angle = np.arccos(-np.dot(e12, e23))
    C_angle = np.arccos(np.dot(e23, e13))

    # Determine configuration
    Trimode = trimodefinder(y, z, x, p1, p2, p3)

    casepLog = (Trimode == 1)
    casenLog = (Trimode == -1)
    casezLog = (Trimode == 0)

    # Initialize strain arrays
    exx = np.zeros(n_points)
    eyy = np.zeros(n_points)
    ezz = np.zeros(n_points)
    exy = np.zeros(n_points)
    exz = np.zeros(n_points)
    eyz = np.zeros(n_points)

    # Configuration I (casepLog)
    if np.any(casepLog):
        # First angular dislocation
        exx1, eyy1, ezz1, exy1, exz1, eyz1 = td_setup_s(
            x[casepLog], y[casepLog], z[casepLog], A_angle,
            bx, by, bz, nu, p1, -e13)

        # Second angular dislocation
        exx2, eyy2, ezz2, exy2, exz2, eyz2 = td_setup_s(
            x[casepLog], y[casepLog], z[casepLog], B_angle,
            bx, by, bz, nu, p2, e12)

        # Third angular dislocation
        exx3, eyy3, ezz3, exy3, exz3, eyz3 = td_setup_s(
            x[casepLog], y[casepLog], z[casepLog], C_angle,
            bx, by, bz, nu, p3, e23)

        exx[casepLog] = exx1 + exx2 + exx3
        eyy[casepLog] = eyy1 + eyy2 + eyy3
        ezz[casepLog] = ezz1 + ezz2 + ezz3
        exy[casepLog] = exy1 + exy2 + exy3
        exz[casepLog] = exz1 + exz2 + exz3
        eyz[casepLog] = eyz1 + eyz2 + eyz3

    # Configuration II (casenLog)
    if np.any(casenLog):
        # First angular dislocation
        exx1, eyy1, ezz1, exy1, exz1, eyz1 = td_setup_s(
            x[casenLog], y[casenLog], z[casenLog], -A_angle,
            bx, by, bz, nu, p1, e13)

        # Second angular dislocation
        exx2, eyy2, ezz2, exy2, exz2, eyz2 = td_setup_s(
            x[casenLog], y[casenLog], z[casenLog], -B_angle,
            bx, by, bz, nu, p2, -e12)

        # Third angular dislocation
        exx3, eyy3, ezz3, exy3, exz3, eyz3 = td_setup_s(
            x[casenLog], y[casenLog], z[casenLog], -C_angle,
            bx, by, bz, nu, p3, -e23)

        exx[casenLog] = exx1 + exx2 + exx3
        eyy[casenLog] = eyy1 + eyy2 + eyy3
        ezz[casenLog] = ezz1 + ezz2 + ezz3
        exy[casenLog] = exy1 + exy2 + exy3
        exz[casenLog] = exz1 + exz2 + exz3
        eyz[casenLog] = eyz1 + eyz2 + eyz3

    # Configuration III (casezLog) - singular points
    if np.any(casezLog):
        exx[casezLog] = np.nan
        eyy[casezLog] = np.nan
        ezz[casezLog] = np.nan
        exy[casezLog] = np.nan
        exz[casezLog] = np.nan
        eyz[casezLog] = np.nan

    # Transform strain tensor from TDCS to EFCS
    Exx, Eyy, Ezz, Exy, Exz, Eyz = tens_trans(
        exx, eyy, ezz, exy, exz, eyz, A.T)

    # Calculate stress tensor
    Sxx = 2 * mu * Exx + lam * (Exx + Eyy + Ezz)
    Syy = 2 * mu * Eyy + lam * (Exx + Eyy + Ezz)
    Szz = 2 * mu * Ezz + lam * (Exx + Eyy + Ezz)
    Sxy = 2 * mu * Exy
    Sxz = 2 * mu * Exz
    Syz = 2 * mu * Eyz

    # Stack results
    Stress = np.column_stack([Sxx, Syy, Szz, Sxy, Sxz, Syz])
    Strain = np.column_stack([Exx, Eyy, Ezz, Exy, Exz, Eyz])

    return Stress, Strain
