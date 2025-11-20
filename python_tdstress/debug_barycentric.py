"""
Debug script to check barycentric coordinate calculation
"""
import numpy as np
import sys
sys.path.insert(0, '/home/user/TriBIE/python_tdstress')

from td_utils import coord_trans, trimodefinder

# Triangle vertices in EFCS (same as test)
P1 = np.array([-1., -1., -5.])
P2 = np.array([1., -1., -5.])
P3 = np.array([-1., 1., -4.])

# Center point in EFCS
X_center = -1.0/3.0
Y_center = -1.0/3.0
Z_center = -14.0/3.0

# Calculate unit normal, strike, and dip vectors (from tdstress_fs.py)
eY = np.array([0, 1, 0])
eZ = np.array([0, 0, 1])

Vnorm = np.cross(P2 - P1, P3 - P1)
Vnorm = Vnorm / np.linalg.norm(Vnorm)

Vstrike = np.cross(eZ, Vnorm)
if np.linalg.norm(Vstrike) == 0:
    Vstrike = eY * Vnorm[2]
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

x, y, z = coord_trans(X_center - P2[0], Y_center - P2[1], Z_center - P2[2], A)
p1[0], p1[1], p1[2] = coord_trans(P1[0] - P2[0], P1[1] - P2[1], P1[2] - P2[2], A)
p3[0], p3[1], p3[2] = coord_trans(P3[0] - P2[0], P3[1] - P2[1], P3[2] - P2[2], A)

print("=" * 70)
print("Barycentric Coordinate Debug")
print("=" * 70)
print(f"\nTriangle vertices in EFCS:")
print(f"  P1 = {P1}")
print(f"  P2 = {P2}")
print(f"  P3 = {P3}")

print(f"\nCenter point in EFCS:")
print(f"  ({X_center:.6f}, {Y_center:.6f}, {Z_center:.6f})")

print(f"\nTransformation vectors:")
print(f"  Vnorm = {Vnorm}")
print(f"  Vstrike = {Vstrike}")
print(f"  Vdip = {Vdip}")

print(f"\nTriangle vertices in TDCS:")
print(f"  p1 = {p1}")
print(f"  p2 = {p2}")
print(f"  p3 = {p3}")

print(f"\nCenter point in TDCS:")
print(f"  (x, y, z) = ({x[0]:.6f}, {y[0]:.6f}, {z[0]:.6f})")
print(f"  Note: MATLAB convention x=y_td, y=z_td, z=x_td")

# Extract 2D coordinates for barycentric calculation
p1_2d = np.array(p1[1:3])
p2_2d = np.array(p2[1:3])
p3_2d = np.array(p3[1:3])

print(f"\nExtracted 2D coordinates (should be [y_TDCS, z_TDCS]):")
print(f"  p1_2d = {p1_2d}")
print(f"  p2_2d = {p2_2d}")
print(f"  p3_2d = {p3_2d}")

# Calculate barycentric coordinates manually
denominator = ((p2_2d[1] - p3_2d[1]) * (p1_2d[0] - p3_2d[0]) +
               (p3_2d[0] - p2_2d[0]) * (p1_2d[1] - p3_2d[1]))

a = ((p2_2d[1] - p3_2d[1]) * (x - p3_2d[0]) +
     (p3_2d[0] - p2_2d[0]) * (y - p3_2d[1])) / denominator

b = ((p3_2d[1] - p1_2d[1]) * (x - p3_2d[0]) +
     (p1_2d[0] - p3_2d[0]) * (y - p3_2d[1])) / denominator

c = 1 - a - b

print(f"\nBarycentric coordinates (manual calculation):")
print(f"  a = {a[0]:.6f}")
print(f"  b = {b[0]:.6f}")
print(f"  c = {c[0]:.6f}")
print(f"  Sum = {(a+b+c)[0]:.6f} (should be 1.0)")

# Call trimodefinder
trimode = trimodefinder(y, z, x, p1, p2, p3)
print(f"\nTrimode from function: {trimode}")

print(f"\nExpected for center point:")
print(f"  a = b = c = 0.333333 (equilateral triangle)")
print(f"  trimode = 1 (inside triangle)")

print("\n" + "=" * 70)
