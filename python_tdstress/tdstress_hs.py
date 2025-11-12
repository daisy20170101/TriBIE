"""
Triangular Dislocation Stress calculation in elastic half-space.

Translation from MATLAB code by Nikkhoo & Walter (2015).

Reference:
Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
artefact-free solution. Geophysical Journal International.
"""

import numpy as np
from .tdstress_fs import tdstress_fs
from .td_utils import coord_trans
from .ang_setup_fsc import ang_setup_fsc_s


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
    strains and stresses associated with a triangular dislocation in a
    half-space. The function cancels the surface normal tractions induced
    by the main and image dislocations.

    Parameters
    ----------
    X, Y, Z : array_like
        Coordinates of calculation points
    P1, P2, P3 : array_like, shape (3,)
        Triangle vertices
    Ss, Ds, Ts : float
        Slip components (Strike-slip, Dip-slip, Tensile-slip)
    mu, lam : float
        Lame constants

    Returns
    -------
    Stress : ndarray, shape (n, 6)
        Harmonic function stress contribution
    Strain : ndarray, shape (n, 6)
        Harmonic function strain contribution
    """
    # Map slip components from TDCS to local variables
    bx = Ts  # Tensile-slip
    by = Ss  # Strike-slip
    bz = Ds  # Dip-slip

    # Calculate unit strike, dip, and normal vectors for the TD
    # For a horizontal TD as an exception, if the normal vector points upward,
    # the strike and dip vectors point Northward and Westward, whereas if the
    # normal vector points downward, the strike and dip vectors point Southward
    # and Westward, respectively.

    # Calculate normal vector
    Vnorm = np.cross(P2 - P1, P3 - P1)
    norm_Vnorm = np.linalg.norm(Vnorm)
    if norm_Vnorm < 1e-12:
        # Degenerate triangle
        n_points = len(np.atleast_1d(X).flatten())
        Stress = np.zeros((n_points, 6))
        Strain = np.zeros((n_points, 6))
        return Stress, Strain
    Vnorm = Vnorm / norm_Vnorm

    # Calculate strike vector
    eY = np.array([0.0, 1.0, 0.0])
    eZ = np.array([0.0, 0.0, 1.0])
    Vstrike = np.cross(eZ, Vnorm)

    # Special case for horizontal triangles
    norm_Vstrike = np.linalg.norm(Vstrike)
    if norm_Vstrike < 1e-12:
        # Horizontal triangle: strike points North or South depending on normal direction
        Vstrike = eY * Vnorm[2]
        norm_Vstrike = np.linalg.norm(Vstrike)

    if norm_Vstrike > 1e-12:
        Vstrike = Vstrike / norm_Vstrike
    else:
        # Extremely degenerate case
        Vstrike = eY

    # Calculate dip vector
    Vdip = np.cross(Vnorm, Vstrike)

    # Transform slip vector components from TDCS into EFCS
    # A matrix has Vnorm, Vstrike, Vdip as columns
    A = np.column_stack([Vnorm, Vstrike, Vdip])
    bX, bY, bZ = coord_trans(bx, by, bz, A)

    # Calculate contribution of angular dislocation pair on each TD side
    # Side 1: P1-P2
    Stress1, Strain1 = ang_setup_fsc_s(X, Y, Z, bX, bY, bZ, P1, P2, mu, lam)

    # Side 2: P2-P3
    Stress2, Strain2 = ang_setup_fsc_s(X, Y, Z, bX, bY, bZ, P2, P3, mu, lam)

    # Side 3: P3-P1
    Stress3, Strain3 = ang_setup_fsc_s(X, Y, Z, bX, bY, bZ, P3, P1, mu, lam)

    # Calculate total harmonic function contribution to strains and stresses
    Stress = Stress1 + Stress2 + Stress3
    Strain = Strain1 + Strain2 + Strain3

    return Stress, Strain
