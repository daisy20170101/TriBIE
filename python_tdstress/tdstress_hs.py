"""
Triangular Dislocation Stress calculation in elastic half-space.

Translation from MATLAB code by Nikkhoo & Walter (2015).

Reference:
Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
artefact-free solution. Geophysical Journal International.
"""

import numpy as np
from .tdstress_fs import tdstress_fs


def tdstress_hs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam):
    """
    Calculate stresses and strains for a triangular dislocation in elastic half-space.

    Parameters
    ----------
    X, Y, Z : array_like
        Coordinates of calculation points in EFCS (East, North, Up).
        X, Y and Z must have the same size. All Z coordinates must be negative!
    P1, P2, P3 : array_like, shape (3,)
        Coordinates of TD vertices in EFCS. All Z coordinates must be negative!
    Ss, Ds, Ts : float
        TD slip vector components (Strike-slip, Dip-slip, Tensile-slip)
    mu, lam : float
        Lame constants

    Returns
    -------
    Stress : ndarray, shape (n, 6)
        Stress tensor components [Sxx, Syy, Szz, Sxy, Sxz, Syz]
        in the same unit as Lame constants
    Strain : ndarray, shape (n, 6)
        Strain tensor components [Exx, Eyy, Ezz, Exy, Exz, Eyz]
        (dimensionless)

    Raises
    ------
    ValueError
        If any Z coordinates are positive (violates half-space condition)
    """
    # Convert to arrays
    X = np.atleast_1d(X).flatten()
    Y = np.atleast_1d(Y).flatten()
    Z = np.atleast_1d(Z).flatten()

    P1 = np.atleast_1d(P1).flatten()
    P2 = np.atleast_1d(P2).flatten()
    P3 = np.atleast_1d(P3).flatten()

    # Check half-space constraint
    if np.any(Z > 0) or P1[2] > 0 or P2[2] > 0 or P3[2] > 0:
        raise ValueError('Half-space solution: Z coordinates must be negative!')

    # Calculate main dislocation contribution to strains and stresses
    StsMS, StrMS = tdstress_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)

    # Calculate harmonic function contribution to strains and stresses
    StsFSC, StrFSC = tdstress_harfunc(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)

    # Calculate image dislocation contribution to strains and stresses
    P1_img = P1.copy()
    P2_img = P2.copy()
    P3_img = P3.copy()
    P1_img[2] = -P1_img[2]
    P2_img[2] = -P2_img[2]
    P3_img[2] = -P3_img[2]

    StsIS, StrIS = tdstress_fs(X, Y, Z, P1_img, P2_img, P3_img, Ss, Ds, Ts, mu, lam)

    # Special handling for horizontal triangles
    if P1[2] == 0 and P2[2] == 0 and P3[2] == 0:
        StsIS[:, 4] = -StsIS[:, 4]  # Sxz
        StsIS[:, 5] = -StsIS[:, 5]  # Syz
        StrIS[:, 4] = -StrIS[:, 4]  # Exz
        StrIS[:, 5] = -StrIS[:, 5]  # Eyz

    # Calculate complete stress and strain tensor components in EFCS
    Stress = StsMS + StsIS + StsFSC
    Strain = StrMS + StrIS + StrFSC

    return Stress, Strain


def tdstress_harfunc(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam):
    """
    Calculate harmonic function contribution to correct for free surface.

    This function calculates the harmonic function contribution to the
    strains and stresses associated with the main and image dislocations.

    Parameters
    ----------
    X, Y, Z : array_like
        Coordinates of calculation points
    P1, P2, P3 : array_like, shape (3,)
        Triangle vertices
    Ss, Ds, Ts : float
        Slip components
    mu, lam : float
        Lame constants

    Returns
    -------
    Stress : ndarray, shape (n, 6)
        Harmonic function stress contribution
    Strain : ndarray, shape (n, 6)
        Harmonic function strain contribution
    """
    # For now, this is a simplified implementation
    # The full implementation requires AngSetupFSC_S which is complex
    # This would need the angular dislocation pair calculations

    # Placeholder - returns zeros for now
    # Full implementation would calculate free surface corrections
    n_points = len(np.atleast_1d(X).flatten())

    Stress = np.zeros((n_points, 6))
    Strain = np.zeros((n_points, 6))

    # TODO: Implement full AngSetupFSC_S functionality
    # This requires:
    # 1. Angular dislocation pairs on each TD side
    # 2. Free surface correction calculations
    # 3. Transformation and summation of contributions

    return Stress, Strain
