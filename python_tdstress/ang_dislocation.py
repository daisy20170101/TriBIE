"""
Angular dislocation strain calculations.

Translation from MATLAB code by Nikkhoo & Walter (2015).
"""

import numpy as np


def ang_dis_strain(x, y, z, alpha, bx, by, bz, nu):
    """
    Calculate strains associated with an angular dislocation in an elastic full-space.

    Parameters
    ----------
    x, y, z : array_like
        Coordinates of calculation points
    alpha : float
        Angular dislocation angle (radians)
    bx, by, bz : float
        Burgers vector components
    nu : float
        Poisson's ratio

    Returns
    -------
    Exx, Eyy, Ezz, Exy, Exz, Eyz : ndarray
        Strain tensor components
    """
    x = np.atleast_1d(x).flatten()
    y = np.atleast_1d(y).flatten()
    z = np.atleast_1d(z).flatten()

    sinA = np.sin(alpha)
    cosA = np.cos(alpha)
    eta = y * cosA - z * sinA
    zeta = y * sinA + z * cosA

    x2 = x * x
    y2 = y * y
    z2 = z * z

    r2 = x2 + y2 + z2
    r = np.sqrt(r2)
    r3 = r * r2

    rz = r * z
    r2z2 = r2 * z2
    r3z = r3 * z

    W = zeta - r
    W2 = W * W
    Wr = W * r
    W2r = W2 * r
    Wr3 = W * r3
    W2r2 = W2 * r2

    C = (r * cosA - z) / Wr
    S = (r * sinA - y) / Wr

    # Partial derivatives of Burgers' function
    rFi_rx = (eta / r / (r - zeta) - y / r / (r - z)) / (4 * np.pi)
    rFi_ry = (x / r / (r - z) - cosA * x / r / (r - zeta)) / (4 * np.pi)
    rFi_rz = (sinA * x / r / (r - zeta)) / (4 * np.pi)

    # Strain components
    Exx = (bx * rFi_rx +
           bx / (8 * np.pi * (1 - nu)) * (eta / Wr + eta * x2 / W2r2 -
                                           eta * x2 / Wr3 + y / rz - x2 * y / r2z2 - x2 * y / r3z) -
           by * x / (8 * np.pi * (1 - nu)) * (((2 * nu + 1) / Wr + x2 / W2r2 - x2 / Wr3) * cosA +
                                                (2 * nu + 1) / rz - x2 / r2z2 - x2 / r3z) +
           bz * x * sinA / (8 * np.pi * (1 - nu)) * ((2 * nu + 1) / Wr + x2 / W2r2 - x2 / Wr3))

    Eyy = (by * rFi_ry +
           bx / (8 * np.pi * (1 - nu)) * ((1 / Wr + S**2 - y2 / Wr3) * eta +
                                           (2 * nu + 1) * y / rz - y**3 / r2z2 - y**3 / r3z -
                                           2 * nu * cosA * S) -
           by * x / (8 * np.pi * (1 - nu)) * (1 / rz - y2 / r2z2 - y2 / r3z +
                                               (1 / Wr + S**2 - y2 / Wr3) * cosA) +
           bz * x * sinA / (8 * np.pi * (1 - nu)) * (1 / Wr + S**2 - y2 / Wr3))

    Ezz = (bz * rFi_rz +
           bx / (8 * np.pi * (1 - nu)) * (eta / W / r + eta * C**2 - eta * z2 / Wr3 +
                                           y * z / r3 + 2 * nu * sinA * C) -
           by * x / (8 * np.pi * (1 - nu)) * ((1 / Wr + C**2 - z2 / Wr3) * cosA + z / r3) +
           bz * x * sinA / (8 * np.pi * (1 - nu)) * (1 / Wr + C**2 - z2 / Wr3))

    Exy = (bx * rFi_ry / 2 + by * rFi_rx / 2 -
           bx / (8 * np.pi * (1 - nu)) * (x * y2 / r2z2 - nu * x / rz + x * y2 / r3z -
                                           nu * x * cosA / Wr + eta * x * S / Wr + eta * x * y / Wr3) +
           by / (8 * np.pi * (1 - nu)) * (x2 * y / r2z2 - nu * y / rz + x2 * y / r3z +
                                           nu * cosA * S + x2 * y * cosA / Wr3 + x2 * cosA * S / Wr) -
           bz * sinA / (8 * np.pi * (1 - nu)) * (nu * S + x2 * S / Wr + x2 * y / Wr3))

    Exz = (bx * rFi_rz / 2 + bz * rFi_rx / 2 -
           bx / (8 * np.pi * (1 - nu)) * (-x * y / r3 + nu * x * sinA / Wr +
                                           eta * x * C / Wr + eta * x * z / Wr3) +
           by / (8 * np.pi * (1 - nu)) * (-x2 / r3 + nu / r + nu * cosA * C +
                                           x2 * z * cosA / Wr3 + x2 * cosA * C / Wr) -
           bz * sinA / (8 * np.pi * (1 - nu)) * (nu * C + x2 * C / Wr + x2 * z / Wr3))

    Eyz = (by * rFi_rz / 2 + bz * rFi_ry / 2 +
           bx / (8 * np.pi * (1 - nu)) * (y2 / r3 - nu / r - nu * cosA * C + nu * sinA * S +
                                           eta * sinA * cosA / W2 - eta * (y * cosA + z * sinA) / W2r +
                                           eta * y * z / W2r2 - eta * y * z / Wr3) -
           by * x / (8 * np.pi * (1 - nu)) * (y / r3 + sinA * cosA**2 / W2 -
                                               cosA * (y * cosA + z * sinA) / W2r +
                                               y * z * cosA / W2r2 - y * z * cosA / Wr3) -
           bz * x * sinA / (8 * np.pi * (1 - nu)) * (y * z / Wr3 - sinA * cosA / W2 +
                                                       (y * cosA + z * sinA) / W2r - y * z / W2r2))

    return Exx, Eyy, Ezz, Exy, Exz, Eyz


def td_setup_s(x, y, z, alpha, bx, by, bz, nu, tri_vertex, side_vec):
    """
    Transform coordinates and slip components, then calculate strain.

    TDSetupS transforms coordinates of the calculation points as well as
    slip vector components from TDCS into ADCS. It then calculates the
    strains in ADCS and transforms them into TDCS.

    Parameters
    ----------
    x, y, z : array_like
        Coordinates in TDCS
    alpha : float
        Angular dislocation angle
    bx, by, bz : float
        Burgers vector components in TDCS
    nu : float
        Poisson's ratio
    tri_vertex : array_like, shape (3,)
        Triangle vertex coordinates in TDCS
    side_vec : array_like, shape (3,)
        Side vector in TDCS

    Returns
    -------
    exx, eyy, ezz, exy, exz, eyz : ndarray
        Strain components in TDCS
    """
    from .td_utils import tens_trans

    # Transformation matrix
    # A = [[SideVec(3);-SideVec(2)] SideVec(2:3)]'
    A = np.array([[side_vec[2], -side_vec[1]],
                  [side_vec[1], side_vec[2]]])

    # Transform coordinates from TDCS into ADCS
    y1 = A[0, 0] * (y - tri_vertex[1]) + A[0, 1] * (z - tri_vertex[2])
    z1 = A[1, 0] * (y - tri_vertex[1]) + A[1, 1] * (z - tri_vertex[2])

    # Transform slip vector components from TDCS into ADCS
    by1 = A[0, 0] * by + A[0, 1] * bz
    bz1 = A[1, 0] * by + A[1, 1] * bz

    # Calculate strains in ADCS
    exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs = ang_dis_strain(
        x, y1, z1, -np.pi + alpha, bx, by1, bz1, nu)

    # Transform strains from ADCS into TDCS
    # B = [[1 0 0];[zeros(2,1),A']]
    B = np.array([[1, 0, 0],
                  [0, A[0, 0], A[1, 0]],
                  [0, A[0, 1], A[1, 1]]])

    exx, eyy, ezz, exy, exz, eyz = tens_trans(
        exx_adcs, eyy_adcs, ezz_adcs, exy_adcs, exz_adcs, eyz_adcs, B)

    return exx, eyy, ezz, exy, exz, eyz
