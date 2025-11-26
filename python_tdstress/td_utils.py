"""
Utility functions for Triangular Dislocation calculations.

Translation from MATLAB code by Nikkhoo & Walter (2015).

Reference:
Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
artefact-free solution. Geophysical Journal International.
"""

import numpy as np


def coord_trans(x1, x2, x3, A):
    """
    Transform coordinates from one system to another.

    CoordTrans transforms the coordinates of vectors from x1x2x3 coordinate
    system to X1X2X3 coordinate system. "A" is the transformation matrix,
    whose columns e1,e2 and e3 are the unit base vectors of the x1x2x3.
    The coordinates of e1,e2 and e3 in A must be given in X1X2X3.
    The transpose of A (i.e., A.T) will transform the coordinates from
    X1X2X3 into x1x2x3.

    Parameters
    ----------
    x1, x2, x3 : array_like
        Input coordinates
    A : ndarray, shape (3, 3)
        Transformation matrix

    Returns
    -------
    X1, X2, X3 : ndarray
        Transformed coordinates
    """
    x1 = np.atleast_1d(x1).flatten()
    x2 = np.atleast_1d(x2).flatten()
    x3 = np.atleast_1d(x3).flatten()

    # Stack coordinates and apply transformation
    r = A @ np.vstack([x1, x2, x3])

    X1 = r[0, :]
    X2 = r[1, :]
    X3 = r[2, :]

    return X1, X2, X3


def tens_trans(Txx1, Tyy1, Tzz1, Txy1, Txz1, Tyz1, A):
    """
    Transform tensor components from one coordinate system to another.

    TensTrans transforms the coordinates of tensors from x1y1z1 coordinate
    system to x2y2z2 coordinate system. "A" is the transformation matrix,
    whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1.
    The coordinates of e1,e2 and e3 in A must be given in x2y2z2.
    The transpose of A (i.e., A.T) does the transformation from x2y2z2 into x1y1z1.

    Parameters
    ----------
    Txx1, Tyy1, Tzz1, Txy1, Txz1, Tyz1 : float or ndarray
        Tensor components in coordinate system 1
    A : ndarray, shape (3, 3)
        Transformation matrix (column-major: columns are unit vectors)

    Returns
    -------
    Txx2, Tyy2, Tzz2, Txy2, Txz2, Tyz2 : ndarray
        Tensor components in coordinate system 2
    """
    # Flatten A to column-major order (like MATLAB)
    A_flat = A.T.flatten()  # A.T.flatten() gives column-major order

    Txx2 = (A_flat[0]**2 * Txx1 + 2*A_flat[0]*A_flat[3] * Txy1 +
            2*A_flat[0]*A_flat[6] * Txz1 + 2*A_flat[3]*A_flat[6] * Tyz1 +
            A_flat[3]**2 * Tyy1 + A_flat[6]**2 * Tzz1)

    Tyy2 = (A_flat[1]**2 * Txx1 + 2*A_flat[1]*A_flat[4] * Txy1 +
            2*A_flat[1]*A_flat[7] * Txz1 + 2*A_flat[4]*A_flat[7] * Tyz1 +
            A_flat[4]**2 * Tyy1 + A_flat[7]**2 * Tzz1)

    Tzz2 = (A_flat[2]**2 * Txx1 + 2*A_flat[2]*A_flat[5] * Txy1 +
            2*A_flat[2]*A_flat[8] * Txz1 + 2*A_flat[5]*A_flat[8] * Tyz1 +
            A_flat[5]**2 * Tyy1 + A_flat[8]**2 * Tzz1)

    Txy2 = (A_flat[0]*A_flat[1] * Txx1 +
            (A_flat[0]*A_flat[4] + A_flat[1]*A_flat[3]) * Txy1 +
            (A_flat[0]*A_flat[7] + A_flat[1]*A_flat[6]) * Txz1 +
            (A_flat[3]*A_flat[7] + A_flat[4]*A_flat[6]) * Tyz1 +
            A_flat[3]*A_flat[4] * Tyy1 + A_flat[6]*A_flat[7] * Tzz1)

    Txz2 = (A_flat[0]*A_flat[2] * Txx1 +
            (A_flat[0]*A_flat[5] + A_flat[2]*A_flat[3]) * Txy1 +
            (A_flat[0]*A_flat[8] + A_flat[2]*A_flat[6]) * Txz1 +
            (A_flat[3]*A_flat[8] + A_flat[5]*A_flat[6]) * Tyz1 +
            A_flat[3]*A_flat[5] * Tyy1 + A_flat[6]*A_flat[8] * Tzz1)

    Tyz2 = (A_flat[1]*A_flat[2] * Txx1 +
            (A_flat[1]*A_flat[5] + A_flat[2]*A_flat[4]) * Txy1 +
            (A_flat[1]*A_flat[8] + A_flat[2]*A_flat[7]) * Txz1 +
            (A_flat[4]*A_flat[8] + A_flat[5]*A_flat[7]) * Tyz1 +
            A_flat[4]*A_flat[5] * Tyy1 + A_flat[7]*A_flat[8] * Tzz1)

    return Txx2, Tyy2, Tzz2, Txy2, Txz2, Tyz2


def trimodefinder(x, y, z, p1, p2, p3):
    """
    Calculate normalized barycentric coordinates and determine configuration.

    trimodefinder calculates the normalized barycentric coordinates of
    the points with respect to the TD vertices and specifies the appropriate
    artefact-free configuration of the angular dislocations for the
    calculations. The input arrays x, y and z share the same size and
    correspond to the y, z and x coordinates in the TDCS, respectively. p1,
    p2 and p3 are two-component arrays representing the y and z coordinates
    of the TD vertices in the TDCS, respectively.

    The components of the output (trimode) corresponding to each calculation
    point are:
        1 for the first configuration
       -1 for the second configuration
        0 for calculation points that lie on the TD sides

    Parameters
    ----------
    x, y, z : array_like
        Coordinates in TDCS (note: x corresponds to y_td, y to z_td, z to x_td)
    p1, p2, p3 : array_like, shape (2,) or (3,)
        Triangle vertex coordinates (only first 2 components used: y and z in TDCS)

    Returns
    -------
    trimode : int or ndarray
        Configuration mode for each point
    """
    x = np.atleast_1d(x).flatten()
    y = np.atleast_1d(y).flatten()
    z = np.atleast_1d(z).flatten()

    # Extract 2D coordinates (y and z components)
    # Note: p1 = [x_TDCS, y_TDCS, z_TDCS], we need [y_TDCS, z_TDCS]
    p1_2d = np.array(p1[1:3])
    p2_2d = np.array(p2[1:3])
    p3_2d = np.array(p3[1:3])

    # Calculate barycentric coordinates
    denominator = ((p2_2d[1] - p3_2d[1]) * (p1_2d[0] - p3_2d[0]) +
                   (p3_2d[0] - p2_2d[0]) * (p1_2d[1] - p3_2d[1]))

    a = ((p2_2d[1] - p3_2d[1]) * (x - p3_2d[0]) +
         (p3_2d[0] - p2_2d[0]) * (y - p3_2d[1])) / denominator

    b = ((p3_2d[1] - p1_2d[1]) * (x - p3_2d[0]) +
         (p1_2d[0] - p3_2d[0]) * (y - p3_2d[1])) / denominator

    c = 1 - a - b

    # Initialize trimode to first configuration
    trimode = np.ones_like(x, dtype=int)

    # Check for second configuration (-1)
    trimode[(a <= 0) & (b > c) & (c > a)] = -1
    trimode[(b <= 0) & (c > a) & (a > b)] = -1
    trimode[(c <= 0) & (a > b) & (b > c)] = -1

    # Check for points on triangle sides (0)
    trimode[(a == 0) & (b >= 0) & (c >= 0)] = 0
    trimode[(a >= 0) & (b == 0) & (c >= 0)] = 0
    trimode[(a >= 0) & (b >= 0) & (c == 0)] = 0

    # Override for points not in the plane (z != 0)
    trimode[(trimode == 0) & (z != 0)] = 1

    # Return scalar if input was scalar
    if trimode.size == 1:
        return int(trimode[0])
    return trimode
