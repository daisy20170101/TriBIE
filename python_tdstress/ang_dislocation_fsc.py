"""
Angular dislocation free surface correction functions.

These functions calculate the harmonic function contribution to strains
for angular dislocations in an elastic half-space.

Translation from MATLAB code by Nikkhoo & Walter (2015).
"""

import numpy as np


def ang_dis_strain_fsc(y1, y2, y3, beta, b1, b2, b3, nu, a):
    """
    Calculate harmonic function contribution to strains for angular dislocation.

    This function calculates the harmonic function contribution to the
    strains associated with an angular dislocation in an elastic half-space.

    Parameters
    ----------
    y1, y2, y3 : array_like
        Coordinates of calculation points in ADCS
        y1: along dip, y2: along strike (x in ADCS), y3: normal (z in ADCS)
    beta : float
        Dip angle in radians
    b1, b2, b3 : float
        Burgers vector components (b1=slip, b2=strike-slip, b3=tensile)
    nu : float
        Poisson's ratio
    a : float
        Vertical distance from free surface to angular dislocation origin

    Returns
    -------
    v11, v22, v33, v12, v13, v23 : ndarray
        Strain tensor components (harmonic contribution)
    """
    # Convert to arrays
    y1 = np.atleast_1d(y1).flatten()
    y2 = np.atleast_1d(y2).flatten()
    y3 = np.atleast_1d(y3).flatten()

    # Trig functions
    sinB = np.sin(beta)
    cosB = np.cos(beta)
    cotB = 1.0 / np.tan(beta) if np.abs(np.sin(beta)) > 1e-10 else 0.0

    # Modified coordinates for free surface
    y3b = y3 + 2*a
    z1b = y1*cosB + y3b*sinB
    z3b = -y1*sinB + y3b*cosB
    rb2 = y1**2 + y2**2 + y3b**2
    rb = np.sqrt(rb2)

    # Intermediate variables (W functions)
    W1 = rb*cosB + y3b
    W2 = cosB + a/rb
    W3 = cosB + y3b/rb
    W4 = nu + a/rb
    W5 = 2*nu + a/rb
    W6 = rb + y3b
    W7 = rb + z3b
    W8 = y3 + a
    W9 = 1 + a/rb/cosB

    N1 = 1 - 2*nu

    # Partial derivatives of Burgers' function
    rFib_ry2 = z1b/rb/(rb+z3b) - y1/rb/(rb+y3b)  # y2 = x in ADCS
    rFib_ry1 = y2/rb/(rb+y3b) - cosB*y2/rb/(rb+z3b)  # y1 = y in ADCS
    rFib_ry3 = -sinB*y2/rb/(rb+z3b)  # y3 = z in ADCS

    # Strain component v11 (Exx in ADCS)
    v11 = (
        b1*(1/4*(
            (-2+2*nu)*N1*rFib_ry1*cotB**2 -
            N1*y2/W6**2*((1-W5)*cotB - y1/W6*W4)/rb*y1 +
            N1*y2/W6*(a/rb**3*y1*cotB - 1/W6*W4 + y1**2/W6**2*W4/rb + y1**2/W6*a/rb**3) -
            N1*y2*cosB*cotB/W7**2*W2*(y1/rb - sinB) -
            N1*y2*cosB*cotB/W7*a/rb**3*y1 -
            3*a*y2*W8*cotB/rb**5*y1 -
            y2*W8/rb**3/W6*(-N1*cotB + y1/W6*W5 + a*y1/rb2)*y1 -
            y2*W8/rb2/W6**2*(-N1*cotB + y1/W6*W5 + a*y1/rb2)*y1 +
            y2*W8/rb/W6*(1/W6*W5 - y1**2/W6**2*W5/rb - y1**2/W6*a/rb**3 + a/rb2 - 2*a*y1**2/rb2**2) -
            y2*W8/rb**3/W7*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2)*y1 -
            y2*W8/rb/W7**2*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2)*(y1/rb - sinB) +
            y2*W8/rb/W7*(-cosB/W7**2*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB)*(y1/rb - sinB) +
                         cosB/W7*(1/rb*cosB*y1*(N1*cosB - a/rb)*cotB + W1*a/rb**3*y1*cotB + (2-2*nu)*(1/rb*sinB*y1 - 1)*cosB) +
                         2*a*y3b*cosB*cotB/rb2**2*y1)
        )/np.pi/(1-nu)) +
        b2*(1/4*(
            N1*(((2-2*nu)*cotB**2 + nu)/rb*y1/W6 - ((2-2*nu)*cotB**2 + 1)*cosB*(y1/rb - sinB)/W7) -
            N1/W6**2*(-N1*y1*cotB + nu*y3b - a + a*y1*cotB/rb + y1**2/W6*W4)/rb*y1 +
            N1/W6*(-N1*cotB + a*cotB/rb - a*y1**2*cotB/rb**3 + 2*y1/W6*W4 - y1**3/W6**2*W4/rb - y1**3/W6*a/rb**3) +
            N1*cotB/W7**2*(z1b*cosB - a*(rb*sinB - y1)/rb/cosB)*(y1/rb - sinB) -
            N1*cotB/W7*(cosB**2 - a*(1/rb*sinB*y1 - 1)/rb/cosB + a*(rb*sinB - y1)/rb**3/cosB*y1) -
            a*W8*cotB/rb**3 + 3*a*y1**2*W8*cotB/rb**5 -
            W8/W6**2*(2*nu + 1/rb*(N1*y1*cotB + a) - y1**2/rb/W6*W5 - a*y1**2/rb**3)/rb*y1 +
            W8/W6*(-1/rb**3*(N1*y1*cotB + a)*y1 + 1/rb*N1*cotB - 2*y1/rb/W6*W5 + y1**3/rb**3/W6*W5 + y1**3/rb2/W6**2*W5 + y1**3/rb2**2/W6*a - 2*a/rb**3*y1 + 3*a*y1**3/rb**5) -
            W8*cotB/W7**2*(-cosB*sinB + a*y1*y3b/rb**3/cosB + (rb*sinB - y1)/rb*((2-2*nu)*cosB - W1/W7*W9))*(y1/rb - sinB) +
            W8*cotB/W7*(a*y3b/rb**3/cosB - 3*a*y1**2*y3b/rb**5/cosB +
                        (1/rb*sinB*y1 - 1)/rb*((2-2*nu)*cosB - W1/W7*W9) -
                        (rb*sinB - y1)/rb**3*((2-2*nu)*cosB - W1/W7*W9)*y1 +
                        (rb*sinB - y1)/rb*(-1/rb*cosB*y1/W7*W9 + W1/W7**2*W9*(y1/rb - sinB) + W1/W7*a/rb**3/cosB*y1))
        )/np.pi/(1-nu)) +
        b3*(1/4*(
            N1*(-y2/W6**2*(1 + a/rb)/rb*y1 - y2/W6*a/rb**3*y1 + y2*cosB/W7**2*W2*(y1/rb - sinB) + y2*cosB/W7*a/rb**3*y1) +
            y2*W8/rb**3*(a/rb2 + 1/W6)*y1 -
            y2*W8/rb*(-2*a/rb2**2*y1 - 1/W6**2/rb*y1) -
            y2*W8*cosB/rb**3/W7*(W1/W7*W2 + a*y3b/rb2)*y1 -
            y2*W8*cosB/rb/W7**2*(W1/W7*W2 + a*y3b/rb2)*(y1/rb - sinB) +
            y2*W8*cosB/rb/W7*(1/rb*cosB*y1/W7*W2 - W1/W7**2*W2*(y1/rb - sinB) - W1/W7*a/rb**3*y1 - 2*a*y3b/rb2**2*y1)
        )/np.pi/(1-nu))
    )

    # Strain component v22 (Eyy in ADCS) - y2 derivatives
    v22 = (
        b1*(1/4*(
            N1*(((2-2*nu)*cotB**2 - nu)/rb*y2/W6 - ((2-2*nu)*cotB**2 + 1 - 2*nu)*cosB/rb*y2/W7) +
            N1/W6**2*(y1*cotB*(1 - W5) + nu*y3b - a + y2**2/W6*W4)/rb*y2 -
            N1/W6*(a*y1*cotB/rb**3*y2 + 2*y2/W6*W4 - y2**3/W6**2*W4/rb - y2**3/W6*a/rb**3) +
            N1*z1b*cotB/W7**2*W2/rb*y2 +
            N1*z1b*cotB/W7*a/rb**3*y2 +
            3*a*y2*W8*cotB/rb**5*y1 -
            W8/W6**2*(-2*nu + 1/rb*(N1*y1*cotB - a) + y2**2/rb/W6*W5 + a*y2**2/rb**3)/rb*y2 +
            W8/W6*(-1/rb**3*(N1*y1*cotB - a)*y2 + 2*y2/rb/W6*W5 - y2**3/rb**3/W6*W5 - y2**3/rb2/W6**2*W5 - y2**3/rb2**2/W6*a + 2*a/rb**3*y2 - 3*a*y2**3/rb**5) -
            W8/W7**2*(cosB**2 - 1/rb*(N1*z1b*cotB + a*cosB) + a*y3b*z1b*cotB/rb**3 - 1/rb/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1))/rb*y2 +
            W8/W7*(1/rb**3*(N1*z1b*cotB + a*cosB)*y2 - 3*a*y3b*z1b*cotB/rb**5*y2 +
                   1/rb**3/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)*y2 +
                   1/rb2/W7**2*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)*y2 -
                   1/rb/W7*(2*y2*cosB**2 + a*z1b*cotB/rb**3*W1*y2 - a*z1b*cotB/rb2*cosB*y2))
        )/np.pi/(1-nu)) +
        b2*(1/4*(
            (2-2*nu)*N1*rFib_ry2*cotB**2 +
            N1/W6*((W5 - 1)*cotB + y1/W6*W4) -
            N1*y2**2/W6**2*((W5 - 1)*cotB + y1/W6*W4)/rb +
            N1*y2/W6*(-a/rb**3*y2*cotB - y1/W6**2*W4/rb*y2 - y2/W6*a/rb**3*y1) -
            N1*cotB/W7*W9 +
            N1*y2**2*cotB/W7**2*W9/rb +
            N1*y2**2*cotB/W7*a/rb**3/cosB -
            a*W8*cotB/rb**3 + 3*a*y2**2*W8*cotB/rb**5 +
            W8/rb/W6*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6)) -
            y2**2*W8/rb**3/W6*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6)) -
            y2**2*W8/rb2/W6**2*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6)) +
            y2*W8/rb/W6*(2*nu*y1/W6**2/rb*y2 + a*y1/rb**3*(1/rb + 1/W6)*y2 - a*y1/rb*(-1/rb**3*y2 - 1/W6**2/rb*y2)) +
            W8*cotB/rb/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) -
            y2**2*W8*cotB/rb**3/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) -
            y2**2*W8*cotB/rb2/W7**2*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) +
            y2*W8*cotB/rb/W7*(1/rb*cosB*y2/W7*W9 - W1/W7**2*W9/rb*y2 - W1/W7*a/rb**3/cosB*y2 - 2*a*y3b/rb2**2/cosB*y2)
        )/np.pi/(1-nu)) +
        b3*(1/4*(
            N1*(-sinB/rb*y2/W7 + y2/W6**2*(1 + a/rb)/rb*y1 + y2/W6*a/rb**3*y1 - z1b/W7**2*W2/rb*y2 - z1b/W7*a/rb**3*y2) -
            y2*W8/rb**3*(a/rb2 + 1/W6)*y1 +
            y1*W8/rb*(-2*a/rb2**2*y2 - 1/W6**2/rb*y2) +
            W8/W7**2*(sinB*(cosB - a/rb) + z1b/rb*(1 + a*y3b/rb2) - 1/rb/W7*(y2**2*cosB*sinB - a*z1b/rb*W1))/rb*y2 -
            W8/W7*(sinB*a/rb**3*y2 - z1b/rb**3*(1 + a*y3b/rb2)*y2 - 2*z1b/rb**5*a*y3b*y2 +
                   1/rb**3/W7*(y2**2*cosB*sinB - a*z1b/rb*W1)*y2 +
                   1/rb2/W7**2*(y2**2*cosB*sinB - a*z1b/rb*W1)*y2 -
                   1/rb/W7*(2*y2*cosB*sinB + a*z1b/rb**3*W1*y2 - a*z1b/rb2*cosB*y2))
        )/np.pi/(1-nu))
    )

    # Strain component v33 (Ezz in ADCS) - y3 derivatives
    v33 = (
        b1*(1/4*(
            (2-2*nu)*(N1*rFib_ry3*cotB - y2/W6**2*W5*(y3b/rb + 1) -
                      1/2*y2/W6*a/rb**3*2*y3b +
                      y2*cosB/W7**2*W2*W3 +
                      1/2*y2*cosB/W7*a/rb**3*2*y3b) +
            y2/rb*(2*nu/W6 + a/rb2) -
            1/2*y2*W8/rb**3*(2*nu/W6 + a/rb2)*2*y3b +
            y2*W8/rb*(-2*nu/W6**2*(y3b/rb + 1) - a/rb2**2*2*y3b) +
            y2*cosB/rb/W7*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2) -
            1/2*y2*W8*cosB/rb**3/W7*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2)*2*y3b -
            y2*W8*cosB/rb/W7**2*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2)*W3 +
            y2*W8*cosB/rb/W7*(-(cosB*y3b/rb + 1)/W7*W2 + W1/W7**2*W2*W3 +
                              1/2*W1/W7*a/rb**3*2*y3b - a/rb2 + a*y3b/rb2**2*2*y3b)
        )/np.pi/(1-nu)) +
        b2*(1/4*(
            (-2 + 2*nu)*N1*cotB*((y3b/rb + 1)/W6 - cosB*W3/W7) +
            (2-2*nu)*y1/W6**2*W5*(y3b/rb + 1) +
            1/2*(2-2*nu)*y1/W6*a/rb**3*2*y3b +
            (2-2*nu)*sinB/W7*W2 -
            (2-2*nu)*z1b/W7**2*W2*W3 -
            1/2*(2-2*nu)*z1b/W7*a/rb**3*2*y3b +
            1/rb*(N1*cotB - 2*nu*y1/W6 - a*y1/rb2) -
            1/2*W8/rb**3*(N1*cotB - 2*nu*y1/W6 - a*y1/rb2)*2*y3b +
            W8/rb*(2*nu*y1/W6**2*(y3b/rb + 1) + a*y1/rb2**2*2*y3b) -
            1/W7*(cosB*sinB + W1*cotB/rb*((2-2*nu)*cosB - W1/W7) +
                  a/rb*(sinB - y3b*z1b/rb2 - z1b*W1/rb/W7)) +
            W8/W7**2*(cosB*sinB + W1*cotB/rb*((2-2*nu)*cosB - W1/W7) +
                      a/rb*(sinB - y3b*z1b/rb2 - z1b*W1/rb/W7))*W3 -
            W8/W7*((cosB*y3b/rb + 1)*cotB/rb*((2-2*nu)*cosB - W1/W7) -
                   1/2*W1*cotB/rb**3*((2-2*nu)*cosB - W1/W7)*2*y3b +
                   W1*cotB/rb*(-(cosB*y3b/rb + 1)/W7 + W1/W7**2*W3) -
                   1/2*a/rb**3*(sinB - y3b*z1b/rb2 - z1b*W1/rb/W7)*2*y3b +
                   a/rb*(-z1b/rb2 - y3b*sinB/rb2 + y3b*z1b/rb2**2*2*y3b -
                         sinB*W1/rb/W7 - z1b*(cosB*y3b/rb + 1)/rb/W7 +
                         1/2*z1b*W1/rb**3/W7*2*y3b + z1b*W1/rb/W7**2*W3))
        )/np.pi/(1-nu)) +
        b3*(1/4*(
            (2-2*nu)*rFib_ry3 -
            (2-2*nu)*y2*sinB/W7**2*W2*W3 -
            1/2*(2-2*nu)*y2*sinB/W7*a/rb**3*2*y3b +
            y2*sinB/rb/W7*(1 + W1/W7*W2 + a*y3b/rb2) -
            1/2*y2*W8*sinB/rb**3/W7*(1 + W1/W7*W2 + a*y3b/rb2)*2*y3b -
            y2*W8*sinB/rb/W7**2*(1 + W1/W7*W2 + a*y3b/rb2)*W3 +
            y2*W8*sinB/rb/W7*((cosB*y3b/rb + 1)/W7*W2 - W1/W7**2*W2*W3 -
                              1/2*W1/W7*a/rb**3*2*y3b + a/rb2 - a*y3b/rb2**2*2*y3b)
        )/np.pi/(1-nu))
    )

    # Shear strain v12 (Exy in ADCS)
    v12_term1 = b1/2*(1/4*(
        (-2 + 2*nu)*N1*rFib_ry2*cotB**2 +
        N1/W6*((1 - W5)*cotB - y1/W6*W4) -
        N1*y2**2/W6**2*((1 - W5)*cotB - y1/W6*W4)/rb +
        N1*y2/W6*(a/rb**3*y2*cotB + y1/W6**2*W4/rb*y2 + y2/W6*a/rb**3*y1) +
        N1*cosB*cotB/W7*W2 -
        N1*y2**2*cosB*cotB/W7**2*W2/rb -
        N1*y2**2*cosB*cotB/W7*a/rb**3 +
        a*W8*cotB/rb**3 - 3*a*y2**2*W8*cotB/rb**5 +
        W8/rb/W6*(-N1*cotB + y1/W6*W5 + a*y1/rb2) -
        y2**2*W8/rb**3/W6*(-N1*cotB + y1/W6*W5 + a*y1/rb2) -
        y2**2*W8/rb2/W6**2*(-N1*cotB + y1/W6*W5 + a*y1/rb2) +
        y2*W8/rb/W6*(-y1/W6**2*W5/rb*y2 - y2/W6*a/rb**3*y1 - 2*a*y1/rb2**2*y2) +
        W8/rb/W7*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2) -
        y2**2*W8/rb**3/W7*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2) -
        y2**2*W8/rb2/W7**2*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2) +
        y2*W8/rb/W7*(-cosB/W7**2*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB)/rb*y2 +
                     cosB/W7*(1/rb*cosB*y2*(N1*cosB - a/rb)*cotB + W1*a/rb**3*y2*cotB + (2-2*nu)/rb*sinB*y2*cosB) +
                     2*a*y3b*cosB*cotB/rb2**2*y2)
    )/np.pi/(1-nu))

    v12_term2 = b2/2*(1/4*(
        N1*(((2-2*nu)*cotB**2 + nu)/rb*y2/W6 - ((2-2*nu)*cotB**2 + 1)*cosB/rb*y2/W7) -
        N1/W6**2*(-N1*y1*cotB + nu*y3b - a + a*y1*cotB/rb + y1**2/W6*W4)/rb*y2 +
        N1/W6*(-a*y1*cotB/rb**3*y2 - y1**2/W6**2*W4/rb*y2 - y1**2/W6*a/rb**3*y2) +
        N1*cotB/W7**2*(z1b*cosB - a*(rb*sinB - y1)/rb/cosB)/rb*y2 -
        N1*cotB/W7*(-a/rb2*sinB*y2/cosB + a*(rb*sinB - y1)/rb**3/cosB*y2) +
        3*a*y2*W8*cotB/rb**5*y1 -
        W8/W6**2*(2*nu + 1/rb*(N1*y1*cotB + a) - y1**2/rb/W6*W5 - a*y1**2/rb**3)/rb*y2 +
        W8/W6*(-1/rb**3*(N1*y1*cotB + a)*y2 + y1**2/rb**3/W6*W5*y2 + y1**2/rb2/W6**2*W5*y2 + y1**2/rb2**2/W6*a*y2 + 3*a*y1**2/rb**5*y2) -
        W8*cotB/W7**2*(-cosB*sinB + a*y1*y3b/rb**3/cosB + (rb*sinB - y1)/rb*((2-2*nu)*cosB - W1/W7*W9))/rb*y2 +
        W8*cotB/W7*(-3*a*y1*y3b/rb**5/cosB*y2 + 1/rb2*sinB*y2*((2-2*nu)*cosB - W1/W7*W9) -
                    (rb*sinB - y1)/rb**3*((2-2*nu)*cosB - W1/W7*W9)*y2 +
                    (rb*sinB - y1)/rb*(-1/rb*cosB*y2/W7*W9 + W1/W7**2*W9/rb*y2 + W1/W7*a/rb**3/cosB*y2))
    )/np.pi/(1-nu))

    v12_term3 = b3/2*(1/4*(
        N1*(1/W6*(1 + a/rb) - y2**2/W6**2*(1 + a/rb)/rb - y2**2/W6*a/rb**3 - cosB/W7*W2 + y2**2*cosB/W7**2*W2/rb + y2**2*cosB/W7*a/rb**3) -
        W8/rb*(a/rb2 + 1/W6) +
        y2**2*W8/rb**3*(a/rb2 + 1/W6) -
        y2*W8/rb*(-2*a/rb2**2*y2 - 1/W6**2/rb*y2) +
        W8*cosB/rb/W7*(W1/W7*W2 + a*y3b/rb2) -
        y2**2*W8*cosB/rb**3/W7*(W1/W7*W2 + a*y3b/rb2) -
        y2**2*W8*cosB/rb2/W7**2*(W1/W7*W2 + a*y3b/rb2) +
        y2*W8*cosB/rb/W7*(1/rb*cosB*y2/W7*W2 - W1/W7**2*W2/rb*y2 - W1/W7*a/rb**3*y2 - 2*a*y3b/rb2**2*y2)
    )/np.pi/(1-nu))

    v12_term4 = b1/2*(1/4*(
        N1*(((2-2*nu)*cotB**2 - nu)/rb*y1/W6 - ((2-2*nu)*cotB**2 + 1 - 2*nu)*cosB*(y1/rb - sinB)/W7) +
        N1/W6**2*(y1*cotB*(1 - W5) + nu*y3b - a + y2**2/W6*W4)/rb*y1 -
        N1/W6*((1 - W5)*cotB + a*y1**2*cotB/rb**3 - y2**2/W6**2*W4/rb*y1 - y2**2/W6*a/rb**3*y1) -
        N1*cosB*cotB/W7*W2 +
        N1*z1b*cotB/W7**2*W2*(y1/rb - sinB) +
        N1*z1b*cotB/W7*a/rb**3*y1 -
        a*W8*cotB/rb**3 + 3*a*y1**2*W8*cotB/rb**5 -
        W8/W6**2*(-2*nu + 1/rb*(N1*y1*cotB - a) + y2**2/rb/W6*W5 + a*y2**2/rb**3)/rb*y1 +
        W8/W6*(-1/rb**3*(N1*y1*cotB - a)*y1 + 1/rb*N1*cotB - y2**2/rb**3/W6*W5*y1 - y2**2/rb2/W6**2*W5*y1 - y2**2/rb2**2/W6*a*y1 - 3*a*y2**2/rb**5*y1) -
        W8/W7**2*(cosB**2 - 1/rb*(N1*z1b*cotB + a*cosB) + a*y3b*z1b*cotB/rb**3 - 1/rb/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1))*(y1/rb - sinB) +
        W8/W7*(1/rb**3*(N1*z1b*cotB + a*cosB)*y1 - 1/rb*N1*cosB*cotB + a*y3b*cosB*cotB/rb**3 - 3*a*y3b*z1b*cotB/rb**5*y1 +
               1/rb**3/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)*y1 +
               1/rb/W7**2*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)*(y1/rb - sinB) -
               1/rb/W7*(-a*cosB*cotB/rb*W1 + a*z1b*cotB/rb**3*W1*y1 - a*z1b*cotB/rb2*cosB*y1))
    )/np.pi/(1-nu))

    v12_term5 = b2/2*(1/4*(
        (2-2*nu)*N1*rFib_ry1*cotB**2 -
        N1*y2/W6**2*((W5 - 1)*cotB + y1/W6*W4)/rb*y1 +
        N1*y2/W6*(-a/rb**3*y1*cotB + 1/W6*W4 - y1**2/W6**2*W4/rb - y1**2/W6*a/rb**3) +
        N1*y2*cotB/W7**2*W9*(y1/rb - sinB) +
        N1*y2*cotB/W7*a/rb**3/cosB*y1 +
        3*a*y2*W8*cotB/rb**5*y1 -
        y2*W8/rb**3/W6*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6))*y1 -
        y2*W8/rb2/W6**2*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6))*y1 +
        y2*W8/rb/W6*(-2*nu/W6 + 2*nu*y1**2/W6**2/rb - a/rb*(1/rb + 1/W6) + a*y1**2/rb**3*(1/rb + 1/W6) - a*y1/rb*(-1/rb**3*y1 - 1/W6**2/rb*y1)) -
        y2*W8*cotB/rb**3/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB)*y1 -
        y2*W8*cotB/rb/W7**2*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB)*(y1/rb - sinB) +
        y2*W8*cotB/rb/W7*(1/rb*cosB*y1/W7*W9 - W1/W7**2*W9*(y1/rb - sinB) - W1/W7*a/rb**3/cosB*y1 - 2*a*y3b/rb2**2/cosB*y1)
    )/np.pi/(1-nu))

    v12_term6 = b3/2*(1/4*(
        N1*(-sinB*(y1/rb - sinB)/W7 - 1/W6*(1 + a/rb) + y1**2/W6**2*(1 + a/rb)/rb + y1**2/W6*a/rb**3 + cosB/W7*W2 - z1b/W7**2*W2*(y1/rb - sinB) - z1b/W7*a/rb**3*y1) +
        W8/rb*(a/rb2 + 1/W6) -
        y1**2*W8/rb**3*(a/rb2 + 1/W6) +
        y1*W8/rb*(-2*a/rb2**2*y1 - 1/W6**2/rb*y1) +
        W8/W7**2*(sinB*(cosB - a/rb) + z1b/rb*(1 + a*y3b/rb2) - 1/rb/W7*(y2**2*cosB*sinB - a*z1b/rb*W1))*(y1/rb - sinB) -
        W8/W7*(sinB*a/rb**3*y1 + cosB/rb*(1 + a*y3b/rb2) - z1b/rb**3*(1 + a*y3b/rb2)*y1 - 2*z1b/rb**5*a*y3b*y1 +
               1/rb**3/W7*(y2**2*cosB*sinB - a*z1b/rb*W1)*y1 +
               1/rb/W7**2*(y2**2*cosB*sinB - a*z1b/rb*W1)*(y1/rb - sinB) -
               1/rb/W7*(-a*cosB/rb*W1 + a*z1b/rb**3*W1*y1 - a*z1b/rb2*cosB*y1))
    )/np.pi/(1-nu))

    v12 = v12_term1 + v12_term2 + v12_term3 + v12_term4 + v12_term5 + v12_term6

    # Shear strain v13 (Exz in ADCS)
    # Note: The formulas use factors like "1/2*...2*y3b" which simplifies to "...y3b"
    # but we keep the MATLAB pattern for accuracy
    v13_b1_half1 = b1/2*(1/4*(
        (-2 + 2*nu)*N1*rFib_ry3*cotB**2 -
        N1*y2/W6**2*((1 - W5)*cotB - y1/W6*W4)*(y3b/rb + 1) +
        N1*y2/W6*(1/2*a/rb**3*2*y3b*cotB + y1/W6**2*W4*(y3b/rb + 1) + 1/2*y1/W6*a/rb**3*2*y3b) -
        N1*y2*cosB*cotB/W7**2*W2*W3 -
        1/2*N1*y2*cosB*cotB/W7*a/rb**3*2*y3b +
        a/rb**3*y2*cotB -
        3/2*a*y2*W8*cotB/rb**5*2*y3b +
        y2/rb/W6*(-N1*cotB + y1/W6*W5 + a*y1/rb2) -
        1/2*y2*W8/rb**3/W6*(-N1*cotB + y1/W6*W5 + a*y1/rb2)*2*y3b -
        y2*W8/rb/W6**2*(-N1*cotB + y1/W6*W5 + a*y1/rb2)*(y3b/rb + 1) +
        y2*W8/rb/W6*(-y1/W6**2*W5*(y3b/rb + 1) - 1/2*y1/W6*a/rb**3*2*y3b - a*y1/rb2**2*2*y3b) +
        y2/rb/W7*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2) -
        1/2*y2*W8/rb**3/W7*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2)*2*y3b -
        y2*W8/rb/W7**2*(cosB/W7*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB) - a*y3b*cosB*cotB/rb2)*W3 +
        y2*W8/rb/W7*(-cosB/W7**2*(W1*(N1*cosB - a/rb)*cotB + (2-2*nu)*(rb*sinB - y1)*cosB)*W3 +
                     cosB/W7*((cosB*y3b/rb + 1)*(N1*cosB - a/rb)*cotB + 1/2*W1*a/rb**3*2*y3b*cotB + 1/2*(2-2*nu)/rb*sinB*2*y3b*cosB) -
                     a*cosB*cotB/rb2 + a*y3b*cosB*cotB/rb2**2*2*y3b)
    )/np.pi/(1-nu))

    v13_b2_half1 = b2/2*(1/4*(
        N1*(((2-2*nu)*cotB**2 + nu)*(y3b/rb + 1)/W6 - ((2-2*nu)*cotB**2 + 1)*cosB*W3/W7) -
        N1/W6**2*(-N1*y1*cotB + nu*y3b - a + a*y1*cotB/rb + y1**2/W6*W4)*(y3b/rb + 1) +
        N1/W6*(nu - 1/2*a*y1*cotB/rb**3*2*y3b - y1**2/W6**2*W4*(y3b/rb + 1) - 1/2*y1**2/W6*a/rb**3*2*y3b) +
        N1*cotB/W7**2*(z1b*cosB - a*(rb*sinB - y1)/rb/cosB)*W3 -
        N1*cotB/W7*(cosB*sinB - 1/2*a/rb2*sinB*2*y3b/cosB + 1/2*a*(rb*sinB - y1)/rb**3/cosB*2*y3b) -
        a/rb**3*y1*cotB +
        3/2*a*y1*W8*cotB/rb**5*2*y3b +
        1/W6*(2*nu + 1/rb*(N1*y1*cotB + a) - y1**2/rb/W6*W5 - a*y1**2/rb**3) -
        W8/W6**2*(2*nu + 1/rb*(N1*y1*cotB + a) - y1**2/rb/W6*W5 - a*y1**2/rb**3)*(y3b/rb + 1) +
        W8/W6*(-1/2/rb**3*(N1*y1*cotB + a)*2*y3b + 1/2*y1**2/rb**3/W6*W5*2*y3b + y1**2/rb/W6**2*W5*(y3b/rb + 1) + 1/2*y1**2/rb2**2/W6*a*2*y3b + 3/2*a*y1**2/rb**5*2*y3b) +
        cotB/W7*(-cosB*sinB + a*y1*y3b/rb**3/cosB + (rb*sinB - y1)/rb*((2-2*nu)*cosB - W1/W7*W9)) -
        W8*cotB/W7**2*(-cosB*sinB + a*y1*y3b/rb**3/cosB + (rb*sinB - y1)/rb*((2-2*nu)*cosB - W1/W7*W9))*W3 +
        W8*cotB/W7*(a/rb**3/cosB*y1 - 3/2*a*y1*y3b/rb**5/cosB*2*y3b +
                    1/2/rb2*sinB*2*y3b*((2-2*nu)*cosB - W1/W7*W9) -
                    1/2*(rb*sinB - y1)/rb**3*((2-2*nu)*cosB - W1/W7*W9)*2*y3b +
                    (rb*sinB - y1)/rb*(-(cosB*y3b/rb + 1)/W7*W9 + W1/W7**2*W9*W3 + 1/2*W1/W7*a/rb**3/cosB*2*y3b))
    )/np.pi/(1-nu))

    v13_b3_half1 = b3/2*(1/4*(
        N1*(-y2/W6**2*(1 + a/rb)*(y3b/rb + 1) - 1/2*y2/W6*a/rb**3*2*y3b + y2*cosB/W7**2*W2*W3 + 1/2*y2*cosB/W7*a/rb**3*2*y3b) -
        y2/rb*(a/rb2 + 1/W6) +
        1/2*y2*W8/rb**3*(a/rb2 + 1/W6)*2*y3b -
        y2*W8/rb*(-a/rb2**2*2*y3b - 1/W6**2*(y3b/rb + 1)) +
        y2*cosB/rb/W7*(W1/W7*W2 + a*y3b/rb2) -
        1/2*y2*W8*cosB/rb**3/W7*(W1/W7*W2 + a*y3b/rb2)*2*y3b -
        y2*W8*cosB/rb/W7**2*(W1/W7*W2 + a*y3b/rb2)*W3 +
        y2*W8*cosB/rb/W7*((cosB*y3b/rb + 1)/W7*W2 - W1/W7**2*W2*W3 - 1/2*W1/W7*a/rb**3*2*y3b + a/rb2 - a*y3b/rb2**2*2*y3b)
    )/np.pi/(1-nu))

    v13_b1_half2 = b1/2*(1/4*(
        (2-2*nu)*(N1*rFib_ry1*cotB - y1/W6**2*W5/rb*y2 - y2/W6*a/rb**3*y1 + y2*cosB/W7**2*W2*(y1/rb - sinB) + y2*cosB/W7*a/rb**3*y1) -
        y2*W8/rb**3*(2*nu/W6 + a/rb2)*y1 +
        y2*W8/rb*(-2*nu/W6**2/rb*y1 - 2*a/rb2**2*y1) -
        y2*W8*cosB/rb**3/W7*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2)*y1 -
        y2*W8*cosB/rb/W7**2*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2)*(y1/rb - sinB) +
        y2*W8*cosB/rb/W7*(-1/rb*cosB*y1/W7*W2 + W1/W7**2*W2*(y1/rb - sinB) + W1/W7*a/rb**3*y1 + 2*a*y3b/rb2**2*y1)
    )/np.pi/(1-nu))

    v13_b2_half2 = b2/2*(1/4*(
        (-2 + 2*nu)*N1*cotB*(1/rb*y1/W6 - cosB*(y1/rb - sinB)/W7) -
        (2-2*nu)/W6*W5 +
        (2-2*nu)*y1**2/W6**2*W5/rb +
        (2-2*nu)*y1**2/W6*a/rb**3 +
        (2-2*nu)*cosB/W7*W2 -
        (2-2*nu)*z1b/W7**2*W2*(y1/rb - sinB) -
        (2-2*nu)*z1b/W7*a/rb**3*y1 -
        W8/rb**3*(N1*cotB - 2*nu*y1/W6 - a*y1/rb2)*y1 +
        W8/rb*(-2*nu/W6 + 2*nu*y1**2/W6**2/rb - a/rb2 + 2*a*y1**2/rb2**2) +
        W8/W7**2*(cosB*sinB + W1*cotB/rb*((2-2*nu)*cosB - W1/W7) + a/rb*(sinB - y3b*z1b/rb2 - z1b*W1/rb/W7))*(y1/rb - sinB) -
        W8/W7*(1/rb2*cosB*y1*cotB*((2-2*nu)*cosB - W1/W7) - W1*cotB/rb**3*((2-2*nu)*cosB - W1/W7)*y1 +
               W1*cotB/rb*(-1/rb*cosB*y1/W7 + W1/W7**2*(y1/rb - sinB)) -
               a/rb**3*(sinB - y3b*z1b/rb2 - z1b*W1/rb/W7)*y1 +
               a/rb*(-y3b*cosB/rb2 + 2*y3b*z1b/rb2**2*y1 - cosB*W1/rb/W7 - z1b/rb2*cosB*y1/W7 + z1b*W1/rb**3/W7*y1 + z1b*W1/rb/W7**2*(y1/rb - sinB)))
    )/np.pi/(1-nu))

    v13_b3_half2 = b3/2*(1/4*(
        (2-2*nu)*rFib_ry1 -
        (2-2*nu)*y2*sinB/W7**2*W2*(y1/rb - sinB) -
        (2-2*nu)*y2*sinB/W7*a/rb**3*y1 -
        y2*W8*sinB/rb**3/W7*(1 + W1/W7*W2 + a*y3b/rb2)*y1 -
        y2*W8*sinB/rb/W7**2*(1 + W1/W7*W2 + a*y3b/rb2)*(y1/rb - sinB) +
        y2*W8*sinB/rb/W7*(1/rb*cosB*y1/W7*W2 - W1/W7**2*W2*(y1/rb - sinB) - W1/W7*a/rb**3*y1 - 2*a*y3b/rb2**2*y1)
    )/np.pi/(1-nu))

    v13 = v13_b1_half1 + v13_b2_half1 + v13_b3_half1 + v13_b1_half2 + v13_b2_half2 + v13_b3_half2

    # Shear strain v23 (Eyz in ADCS)
    # Read the MATLAB line 1019-1023 onward for the remaining part
    # This is split similar to v13
    v23_b1_half1 = b1/2*(1/4*(
        N1*(((2-2*nu)*cotB**2 - nu)*(y3b/rb + 1)/W6 - ((2-2*nu)*cotB**2 + 1 - 2*nu)*cosB*W3/W7) +
        N1/W6**2*(y1*cotB*(1 - W5) + nu*y3b - a + y2**2/W6*W4)*(y3b/rb + 1) -
        N1/W6*(1/2*a*y1*cotB/rb**3*2*y3b + nu - y2**2/W6**2*W4*(y3b/rb + 1) - 1/2*y2**2/W6*a/rb**3*2*y3b) -
        N1*sinB*cotB/W7*W2 +
        N1*z1b*cotB/W7**2*W2*W3 +
        1/2*N1*z1b*cotB/W7*a/rb**3*2*y3b -
        a/rb**3*y1*cotB +
        3/2*a*y1*W8*cotB/rb**5*2*y3b +
        1/W6*(-2*nu + 1/rb*(N1*y1*cotB - a) + y2**2/rb/W6*W5 + a*y2**2/rb**3) -
        W8/W6**2*(-2*nu + 1/rb*(N1*y1*cotB - a) + y2**2/rb/W6*W5 + a*y2**2/rb**3)*(y3b/rb + 1) +
        W8/W6*(-1/2/rb**3*(N1*y1*cotB - a)*2*y3b - 1/2*y2**2/rb**3/W6*W5*2*y3b - y2**2/rb/W6**2*W5*(y3b/rb + 1) - 1/2*y2**2/rb2**2/W6*a*2*y3b - 3/2*a*y2**2/rb**5*2*y3b) +
        1/W7*(cosB**2 - 1/rb*(N1*z1b*cotB + a*cosB) + a*y3b*z1b*cotB/rb**3 - 1/rb/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)) -
        W8/W7**2*(cosB**2 - 1/rb*(N1*z1b*cotB + a*cosB) + a*y3b*z1b*cotB/rb**3 - 1/rb/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1))*W3 +
        W8/W7*(1/2/rb**3*(N1*z1b*cotB + a*cosB)*2*y3b - 1/rb*N1*sinB*cotB + a*z1b*cotB/rb**3 + a*y3b*sinB*cotB/rb**3 - 3/2*a*y3b*z1b*cotB/rb**5*2*y3b +
               1/2/rb**3/W7*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)*2*y3b +
               1/rb/W7**2*(y2**2*cosB**2 - a*z1b*cotB/rb*W1)*W3 -
               1/rb/W7*(-a*sinB*cotB/rb*W1 + 1/2*a*z1b*cotB/rb**3*W1*2*y3b - a*z1b*cotB/rb*(cosB*y3b/rb + 1)))
    )/np.pi/(1-nu))

    v23_b2_half1 = b2/2*(1/4*(
        (2-2*nu)*N1*rFib_ry3*cotB**2 -
        N1*y2/W6**2*((W5 - 1)*cotB + y1/W6*W4)*(y3b/rb + 1) +
        N1*y2/W6*(-1/2*a/rb**3*2*y3b*cotB - y1/W6**2*W4*(y3b/rb + 1) - 1/2*y1/W6*a/rb**3*2*y3b) +
        N1*y2*cotB/W7**2*W9*W3 +
        1/2*N1*y2*cotB/W7*a/rb**3/cosB*2*y3b -
        a/rb**3*y2*cotB +
        3/2*a*y2*W8*cotB/rb**5*2*y3b +
        y2/rb/W6*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6)) -
        1/2*y2*W8/rb**3/W6*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6))*2*y3b -
        y2*W8/rb/W6**2*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6))*(y3b/rb + 1) +
        y2*W8/rb/W6*(2*nu*y1/W6**2*(y3b/rb + 1) + 1/2*a*y1/rb**3*(1/rb + 1/W6)*2*y3b - a*y1/rb*(-1/2/rb**3*2*y3b - 1/W6**2*(y3b/rb + 1))) +
        y2*cotB/rb/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) -
        1/2*y2*W8*cotB/rb**3/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB)*2*y3b -
        y2*W8*cotB/rb/W7**2*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB)*W3 +
        y2*W8*cotB/rb/W7*((cosB*y3b/rb + 1)/W7*W9 - W1/W7**2*W9*W3 - 1/2*W1/W7*a/rb**3/cosB*2*y3b + a/rb2/cosB - a*y3b/rb2**2/cosB*2*y3b)
    )/np.pi/(1-nu))

    v23_b3_half1 = b3/2*(1/4*(
        N1*(-sinB*W3/W7 + y1/W6**2*(1 + a/rb)*(y3b/rb + 1) + 1/2*y1/W6*a/rb**3*2*y3b + sinB/W7*W2 - z1b/W7**2*W2*W3 - 1/2*z1b/W7*a/rb**3*2*y3b) +
        y1/rb*(a/rb2 + 1/W6) -
        1/2*y1*W8/rb**3*(a/rb2 + 1/W6)*2*y3b +
        y1*W8/rb*(-a/rb2**2*2*y3b - 1/W6**2*(y3b/rb + 1)) -
        1/W7*(sinB*(cosB - a/rb) + z1b/rb*(1 + a*y3b/rb2) - 1/rb/W7*(y2**2*cosB*sinB - a*z1b/rb*W1)) +
        W8/W7**2*(sinB*(cosB - a/rb) + z1b/rb*(1 + a*y3b/rb2) - 1/rb/W7*(y2**2*cosB*sinB - a*z1b/rb*W1))*W3 -
        W8/W7*(1/2*sinB*a/rb**3*2*y3b + sinB/rb*(1 + a*y3b/rb2) - 1/2*z1b/rb**3*(1 + a*y3b/rb2)*2*y3b + z1b/rb*(a/rb2 - a*y3b/rb2**2*2*y3b) +
               1/2/rb**3/W7*(y2**2*cosB*sinB - a*z1b/rb*W1)*2*y3b +
               1/rb/W7**2*(y2**2*cosB*sinB - a*z1b/rb*W1)*W3 -
               1/rb/W7*(-a*sinB/rb*W1 + 1/2*a*z1b/rb**3*W1*2*y3b - a*z1b/rb*(cosB*y3b/rb + 1)))
    )/np.pi/(1-nu))

    # Second half of v23 - mixed derivatives
    v23_b1_half2 = b1/2*(1/4*(
        (2-2*nu)*(N1*rFib_ry2*cotB + 1/W6*W5 - y2**2/W6**2*W5/rb - y2**2/W6*a/rb**3 - cosB/W7*W2 + y2**2*cosB/W7**2*W2/rb + y2**2*cosB/W7*a/rb**3) +
        W8/rb*(2*nu/W6 + a/rb2) -
        y2**2*W8/rb**3*(2*nu/W6 + a/rb2) +
        y2*W8/rb*(-2*nu/W6**2/rb*y2 - 2*a/rb2**2*y2) +
        W8*cosB/rb/W7*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2) -
        y2**2*W8*cosB/rb**3/W7*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2) -
        y2**2*W8*cosB/rb2/W7**2*(1 - 2*nu - W1/W7*W2 - a*y3b/rb2) +
        y2*W8*cosB/rb/W7*(-1/rb*cosB*y2/W7*W2 + W1/W7**2*W2/rb*y2 + W1/W7*a/rb**3*y2 + 2*a*y3b/rb2**2*y2)
    )/np.pi/(1-nu))

    v23_b2_half2 = b2/2*(1/4*(
        (-2 + 2*nu)*N1*cotB*(1/W6 - y2**2/W6**2/rb - y2**2/W6*a/rb**3/rb - cosB/W7 + y2**2*cosB/W7**2/rb + y2**2*cosB/W7*a/rb**3) +
        1/rb*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6)) -
        y2**2/rb**3*(N1*cotB - 2*nu*y1/W6 - a*y1/rb*(1/rb + 1/W6)) +
        y2/rb*(2*nu*y1/W6**2/rb*y2 + a*y1/rb**3*(1/rb + 1/W6)*y2 - a*y1/rb*(-1/rb**3*y2 - 1/W6**2/rb*y2)) +
        cotB/rb/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) -
        y2**2*cotB/rb**3/W7*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) -
        y2**2*cotB/rb2/W7**2*((-2 + 2*nu)*cosB + W1/W7*W9 + a*y3b/rb2/cosB) +
        y2*cotB/rb/W7*(1/rb*cosB*y2/W7*W9 - W1/W7**2*W9/rb*y2 - W1/W7*a/rb**3/cosB*y2 - 2*a*y3b/rb2**2/cosB*y2)
    )/np.pi/(1-nu))

    v23_b3_half2 = b3/2*(1/4*(
        N1*(-sinB/W7 + y2**2*sinB/W7**2/rb + y2**2*sinB/W7*a/rb**3 + 1/W6 - y2**2/W6**2/rb - y2**2/W6*a/rb**3) +
        1/rb*(a/rb2 + 1/W6) -
        y2**2/rb**3*(a/rb2 + 1/W6) +
        y2/rb*(-2*a/rb2**2*y2 - 1/W6**2/rb*y2) +
        sinB/rb/W7*(1 + W1/W7*W2 + a*y3b/rb2) -
        y2**2*sinB/rb**3/W7*(1 + W1/W7*W2 + a*y3b/rb2) -
        y2**2*sinB/rb2/W7**2*(1 + W1/W7*W2 + a*y3b/rb2) +
        y2*sinB/rb/W7*(1/rb*cosB*y2/W7*W2 - W1/W7**2*W2/rb*y2 - W1/W7*a/rb**3*y2 - 2*a*y3b/rb2**2*y2)
    )/np.pi/(1-nu))

    v23 = v23_b1_half1 + v23_b2_half1 + v23_b3_half1 + v23_b1_half2 + v23_b2_half2 + v23_b3_half2

    return v11, v22, v33, v12, v13, v23
