"""
Test script for triangular dislocation stress calculations.

This script tests the Python implementation against known reference values
from the MATLAB implementation.
"""

import numpy as np
import sys
import os

# Add parent directory to path so we can import the package
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

from python_tdstress.tdstress_fs import tdstress_fs
from python_tdstress.tdstress_hs import tdstress_hs


def test_triangle_setup():
    """Test with the standard triangle configuration."""
    # Triangle vertices
    P1 = np.array([-1.0, -1.0, -5.0])
    P2 = np.array([1.0, -1.0, -5.0])
    P3 = np.array([-1.0, 1.0, -4.0])

    # Slip components
    Ss = 1.0   # Strike-slip
    Ds = -1.0  # Dip-slip
    Ts = 2.0   # Tensile-slip

    # Elastic parameters
    mu = 3.0e10
    lam = 3.0e10

    return P1, P2, P3, Ss, Ds, Ts, mu, lam


def test_reference_points():
    """Test against known reference values."""
    P1, P2, P3, Ss, Ds, Ts, mu, lam = test_triangle_setup()

    print("=" * 70)
    print("Testing Triangular Dislocation Implementation")
    print("=" * 70)
    print(f"Triangle vertices:")
    print(f"  P1 = {P1}")
    print(f"  P2 = {P2}")
    print(f"  P3 = {P3}")
    print(f"\nSlip components: Ss={Ss}, Ds={Ds}, Ts={Ts}")
    print(f"Elastic parameters: mu={mu:.2e}, lambda={lam:.2e}")
    print()

    # Test points with known reference values (all 15 points from Fortran test)
    test_cases = [
        {
            'name': 'Point 1 (center)',
            'coords': (-1.0/3.0, -1.0/3.0, -14.0/3.0),
            'expected_exx': 0.0481047005255181
        },
        {
            'name': 'Point 2',
            'coords': (0.0, 0.0, 0.0),
            'expected_exx': None  # Surface point - may be singular
        },
        {
            'name': 'Point 3',
            'coords': (0.0, 3.0, 0.0),
            'expected_exx': None  # Reference value not provided
        },
        {
            'name': 'Point 4',
            'coords': (7.0, -1.0, -5.0),
            'expected_exx': 0.000829157341339727
        },
        {
            'name': 'Point 5',
            'coords': (-7.0, -1.0, -5.0),
            'expected_exx': 0.00114439668841158
        },
        {
            'name': 'Point 6',
            'coords': (-1.0, 7.0, -5.0),
            'expected_exx': None
        },
        {
            'name': 'Point 7',
            'coords': (-1.0, -7.0, -5.0),
            'expected_exx': None
        },
        {
            'name': 'Point 8',
            'coords': (-1.0, -1.0, 7.0),
            'expected_exx': None
        },
        {
            'name': 'Point 9',
            'coords': (-1.0, -1.0, -12.0),
            'expected_exx': None
        },
        {
            'name': 'Point 10',
            'coords': (0.0, 0.0, -5.0),
            'expected_exx': None
        },
        {
            'name': 'Point 11',
            'coords': (0.0, -1.0, -5.0),
            'expected_exx': None
        },
        {
            'name': 'Point 12',
            'coords': (1.0, -1.0, -1.0),
            'expected_exx': 0.00441202690885827
        },
        {
            'name': 'Point 13',
            'coords': (-1.0, 1.0, -1.0),
            'expected_exx': None
        },
        {
            'name': 'Point 14',
            'coords': (-1.0, -1.0, -1.0),
            'expected_exx': None
        },
        {
            'name': 'Point 15',
            'coords': (1.0, -1.0, -8.0),
            'expected_exx': -0.000914111766849476
        },
    ]

    print("-" * 70)
    print("Testing TDstressFS (Full-Space)")
    print("-" * 70)

    for test_case in test_cases:
        X, Y, Z = test_case['coords']
        expected = test_case['expected_exx']

        try:
            stress, strain = tdstress_fs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)
            exx = strain[0, 0]

            print(f"\n{test_case['name']}: ({X:.3f}, {Y:.3f}, {Z:.3f})")
            print(f"  Got Exx:      {exx:.15e}")

            if expected is not None:
                error = abs(exx - expected)
                rel_error = abs(error / expected) * 100 if expected != 0 else np.nan

                print(f"  Expected Exx: {expected:.15e}")
                print(f"  Error:        {error:.15e}")
                if not np.isnan(rel_error):
                    print(f"  Rel. Error:   {rel_error:.6f}%")

                if np.isnan(exx):
                    print(f"  *** WARNING: Got NaN ***")
                elif rel_error < 0.1:
                    print(f"  *** PASS ***")
                elif rel_error < 1.0:
                    print(f"  *** CLOSE ***")
                else:
                    print(f"  *** FAIL ***")
            else:
                if np.isnan(exx):
                    print(f"  *** NaN (may be expected for singular points) ***")
                else:
                    print(f"  *** OK (no reference value) ***")

        except Exception as e:
            print(f"\n{test_case['name']}: ({X:.3f}, {Y:.3f}, {Z:.3f})")
            print(f"  *** ERROR: {str(e)} ***")

    print("\n" + "=" * 70)


def test_half_space():
    """Test half-space calculation with complete harmonic function."""
    P1, P2, P3, Ss, Ds, Ts, mu, lam = test_triangle_setup()

    print("\nTesting TDstressHS (Half-Space)")
    print("-" * 70)
    print("NOTE: Complete implementation with harmonic function")
    print()

    # Test with the same 5 points that have reference values
    test_cases = [
        {
            'name': 'Point 1 (center)',
            'coords': (-1.0/3.0, -1.0/3.0, -14.0/3.0),
            'expected_exx': 0.0481047005255181
        },
        {
            'name': 'Point 4',
            'coords': (7.0, -1.0, -5.0),
            'expected_exx': 0.000829157341339727
        },
        {
            'name': 'Point 5',
            'coords': (-7.0, -1.0, -5.0),
            'expected_exx': 0.00114439668841158
        },
        {
            'name': 'Point 12',
            'coords': (1.0, -1.0, -1.0),
            'expected_exx': 0.00441202690885827
        },
        {
            'name': 'Point 15',
            'coords': (1.0, -1.0, -8.0),
            'expected_exx': -0.000914111766849476
        },
    ]

    for test_case in test_cases:
        X, Y, Z = test_case['coords']
        expected = test_case['expected_exx']

        try:
            stress, strain = tdstress_hs(X, Y, Z, P1, P2, P3, Ss, Ds, Ts, mu, lam)
            exx = strain[0, 0]

            error = abs(exx - expected)
            rel_error = abs(error / expected) * 100 if expected != 0 else np.nan

            print(f"\n{test_case['name']}: ({X:.3f}, {Y:.3f}, {Z:.3f})")
            print(f"  Expected Exx: {expected:.15e}")
            print(f"  Got Exx:      {exx:.15e}")
            print(f"  Error:        {error:.15e}")
            if not np.isnan(rel_error):
                print(f"  Rel. Error:   {rel_error:.6f}%")

            if np.isnan(exx):
                print(f"  *** WARNING: Got NaN ***")
            elif rel_error < 0.1:
                print(f"  *** PASS ***")
            elif rel_error < 1.0:
                print(f"  *** CLOSE ***")
            else:
                print(f"  *** FAIL ***")

        except Exception as e:
            print(f"\n{test_case['name']}: ({X:.3f}, {Y:.3f}, {Z:.3f})")
            print(f"  *** ERROR: {str(e)} ***")

    print()


if __name__ == "__main__":
    test_reference_points()
    test_half_space()
    print("=" * 70)
    print("Testing complete!")
    print("=" * 70)
