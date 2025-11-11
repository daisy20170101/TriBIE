"""
Angular dislocation free surface correction setup functions.

These functions handle the setup and calculation of free surface corrections
for angular dislocation pairs on triangular dislocation sides.

Translation from MATLAB code by Nikkhoo & Walter (2015).
"""

import numpy as np
from .td_utils import coord_trans, tens_trans
from .ang_dislocation_fsc import ang_dis_strain_fsc


def ang_setup_fsc_s(X, Y, Z, bX, bY, bZ, PA, PB, mu, lam):
    """
    Calculate Free Surface Correction for angular dislocation pair.

    This function calculates the free surface correction to strains and
    stresses associated with an angular dislocation pair on each TD side.

    Parameters
    ----------
    X, Y, Z : array_like
        Coordinates of calculation points in EFCS
    bX, bY, bZ : float
        Burgers vector components in EFCS
    PA, PB : array_like, shape (3,)
        Coordinates of angular dislocation endpoints (side vertices)
    mu, lam : float
        Lame constants

    Returns
    -------
    Stress : ndarray, shape (n_points, 6)
        Free surface correction to stress tensor [Sxx, Syy, Szz, Sxy, Sxz, Syz]
    Strain : ndarray, shape (n_points, 6)
        Free surface correction to strain tensor [Exx, Eyy, Ezz, Exy, Exz, Eyz]
    """
    # Convert to arrays
    X = np.atleast_1d(X).flatten()
    Y = np.atleast_1d(Y).flatten()
    Z = np.atleast_1d(Z).flatten()

    PA = np.atleast_1d(PA).flatten()
    PB = np.atleast_1d(PB).flatten()

    n_points = len(X)

    # Calculate Poisson's ratio
    nu = 1.0 / (1.0 + lam/mu) / 2.0

    # Calculate TD side vector and the angle of the angular dislocation pair
    SideVec = PB - PA
    eZ = np.array([0.0, 0.0, 1.0])

    # Calculate beta angle
    norm_SideVec = np.linalg.norm(SideVec)
    if norm_SideVec < 1e-12:
        # Degenerate case: PA and PB are the same point
        Stress = np.zeros((n_points, 6))
        Strain = np.zeros((n_points, 6))
        return Stress, Strain

    beta = np.arccos(-np.dot(SideVec, eZ) / norm_SideVec)

    # Check for vertical sides
    if np.abs(beta) < np.finfo(float).eps or np.abs(np.pi - beta) < np.finfo(float).eps:
        # Vertical side - no contribution
        Stress = np.zeros((n_points, 6))
        Strain = np.zeros((n_points, 6))
        return Stress, Strain

    # Construct transformation matrix A
    # ey1 is horizontal projection of side vector
    ey1 = np.array([SideVec[0], SideVec[1], 0.0])
    norm_ey1 = np.linalg.norm(ey1)
    if norm_ey1 < 1e-12:
        # Side vector is vertical (already handled above, but safety check)
        Stress = np.zeros((n_points, 6))
        Strain = np.zeros((n_points, 6))
        return Stress, Strain
    ey1 = ey1 / norm_ey1

    ey3 = -eZ
    ey2 = np.cross(ey3, ey1)

    # Transformation matrix: columns are ey1, ey2, ey3
    A = np.column_stack([ey1, ey2, ey3])

    # Transform coordinates from EFCS to the first ADCS (point A)
    y1A, y2A, y3A = coord_trans(X - PA[0], Y - PA[1], Z - PA[2], A.T)

    # Transform side vector to ADCS
    y1AB, y2AB, y3AB = coord_trans(SideVec[0], SideVec[1], SideVec[2], A.T)

    # Coordinates in second ADCS (point B)
    y1B = y1A - y1AB
    y2B = y2A - y2AB
    y3B = y3A - y3AB

    # Transform slip vector components from EFCS to ADCS
    b1, b2, b3 = coord_trans(bX, bY, bZ, A.T)

    # Determine the best artifact-free configuration
    # Configuration based on (beta*y1A) >= 0
    I = (beta * y1A) >= 0

    # Initialize strain components
    v11A = np.zeros(n_points)
    v22A = np.zeros(n_points)
    v33A = np.zeros(n_points)
    v12A = np.zeros(n_points)
    v13A = np.zeros(n_points)
    v23A = np.zeros(n_points)

    v11B = np.zeros(n_points)
    v22B = np.zeros(n_points)
    v33B = np.zeros(n_points)
    v12B = np.zeros(n_points)
    v13B = np.zeros(n_points)
    v23B = np.zeros(n_points)

    # Configuration I: (beta*y1A) >= 0
    if np.any(I):
        v11A_I, v22A_I, v33A_I, v12A_I, v13A_I, v23A_I = ang_dis_strain_fsc(
            -y1A[I], -y2A[I], y3A[I],
            np.pi - beta, -b1, -b2, b3, nu, -PA[2]
        )
        v11A[I] = v11A_I
        v22A[I] = v22A_I
        v33A[I] = v33A_I
        v12A[I] = v12A_I
        v13A[I] = -v13A_I  # Note sign flip
        v23A[I] = -v23A_I  # Note sign flip

        v11B_I, v22B_I, v33B_I, v12B_I, v13B_I, v23B_I = ang_dis_strain_fsc(
            -y1B[I], -y2B[I], y3B[I],
            np.pi - beta, -b1, -b2, b3, nu, -PB[2]
        )
        v11B[I] = v11B_I
        v22B[I] = v22B_I
        v33B[I] = v33B_I
        v12B[I] = v12B_I
        v13B[I] = -v13B_I  # Note sign flip
        v23B[I] = -v23B_I  # Note sign flip

    # Configuration II: (beta*y1A) < 0
    if np.any(~I):
        v11A_II, v22A_II, v33A_II, v12A_II, v13A_II, v23A_II = ang_dis_strain_fsc(
            y1A[~I], y2A[~I], y3A[~I],
            beta, b1, b2, b3, nu, -PA[2]
        )
        v11A[~I] = v11A_II
        v22A[~I] = v22A_II
        v33A[~I] = v33A_II
        v12A[~I] = v12A_II
        v13A[~I] = v13A_II
        v23A[~I] = v23A_II

        v11B_II, v22B_II, v33B_II, v12B_II, v13B_II, v23B_II = ang_dis_strain_fsc(
            y1B[~I], y2B[~I], y3B[~I],
            beta, b1, b2, b3, nu, -PB[2]
        )
        v11B[~I] = v11B_II
        v22B[~I] = v22B_II
        v33B[~I] = v33B_II
        v12B[~I] = v12B_II
        v13B[~I] = v13B_II
        v23B[~I] = v23B_II

    # Calculate total Free Surface Correction to strains in ADCS (difference B - A)
    v11 = v11B - v11A
    v22 = v22B - v22A
    v33 = v33B - v33A
    v12 = v12B - v12A
    v13 = v13B - v13A
    v23 = v23B - v23A

    # Transform total Free Surface Correction to strains from ADCS to EFCS
    # Note: A' (transpose) is used for tensor transformation
    Exx, Eyy, Ezz, Exy, Exz, Eyz = tens_trans(v11, v22, v33, v12, v13, v23, A)

    # Calculate total Free Surface Correction to stresses in EFCS
    trace = Exx + Eyy + Ezz
    Sxx = 2*mu*Exx + lam*trace
    Syy = 2*mu*Eyy + lam*trace
    Szz = 2*mu*Ezz + lam*trace
    Sxy = 2*mu*Exy
    Sxz = 2*mu*Exz
    Syz = 2*mu*Eyz

    # Stack results
    Strain = np.column_stack([Exx, Eyy, Ezz, Exy, Exz, Eyz])
    Stress = np.column_stack([Sxx, Syy, Szz, Sxy, Sxz, Syz])

    return Stress, Strain
