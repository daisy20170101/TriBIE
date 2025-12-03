"""
Trace the exact coordinate values through the barycentric calculation
"""
import numpy as np

# Triangle vertices in 2D (from debug output)
p1_2d = np.array([-2.0, 0.0])
p2_2d = np.array([0.0, 0.0])
p3_2d = np.array([-2.0, 2.23606798])

# Observation point (from debug output - these are y_TDCS, z_TDCS)
# From the debug output:
#   (x, y, z) = (-0.000000, -1.333333, 0.745356)
#   Note: MATLAB convention x=y_td, y=z_td, z=x_td
#
# But trimodefinder is called as: trimodefinder(y, z, x, p1, p2, p3)
# where y, z, x are arrays from TDCS

# From tdstress_fs.py line 98:
# Trimode = trimodefinder(y, z, x, p1, p2, p3)
#
# where y, z, x came from:
# x, y, z = coord_trans(X - P2[0], Y - P2[1], Z - P2[2], A)
#
# So inside trimodefinder:
# - first parameter (named 'x' in function) gets array 'y' = z_TDCS coords
# - second parameter (named 'y' in function) gets array 'z' = x_TDCS coords
# - third parameter (named 'z' in function) gets array 'x' = y_TDCS coords

x_param = -1.333333  # This is z_TDCS (the 'y' array from outside)
y_param = 0.745356   # This is x_TDCS (the 'z' array from outside)

print("=" * 70)
print("Barycentric Formula Trace")
print("=" * 70)

print(f"\nTriangle vertices (2D projection [y_TDCS, z_TDCS]):")
print(f"  p1_2d = {p1_2d}")
print(f"  p2_2d = {p2_2d}")
print(f"  p3_2d = {p3_2d}")

print(f"\nObservation point coordinates passed to formula:")
print(f"  x (parameter) = {x_param:.6f} (actually z_TDCS from coord_trans)")
print(f"  y (parameter) = {y_param:.6f} (actually x_TDCS from coord_trans)")

print(f"\nDenominator calculation:")
term1 = (p2_2d[1] - p3_2d[1]) * (p1_2d[0] - p3_2d[0])
term2 = (p3_2d[0] - p2_2d[0]) * (p1_2d[1] - p3_2d[1])
denominator = term1 + term2
print(f"  (p2_2d[1] - p3_2d[1]) * (p1_2d[0] - p3_2d[0])")
print(f"  = ({p2_2d[1]:.3f} - {p3_2d[1]:.3f}) * ({p1_2d[0]:.3f} - {p3_2d[0]:.3f})")
print(f"  = {p2_2d[1] - p3_2d[1]:.3f} * {p1_2d[0] - p3_2d[0]:.3f}")
print(f"  = {term1:.3f}")
print(f"  +")
print(f"  (p3_2d[0] - p2_2d[0]) * (p1_2d[1] - p3_2d[1])")
print(f"  = ({p3_2d[0]:.3f} - {p2_2d[0]:.3f}) * ({p1_2d[1]:.3f} - {p3_2d[1]:.3f})")
print(f"  = {p3_2d[0] - p2_2d[0]:.3f} * {p1_2d[1] - p3_2d[1]:.3f}")
print(f"  = {term2:.3f}")
print(f"  denominator = {denominator:.3f}")

print(f"\nBarycentric coordinate 'a' calculation:")
term1_a = (p2_2d[1] - p3_2d[1]) * (x_param - p3_2d[0])
term2_a = (p3_2d[0] - p2_2d[0]) * (y_param - p3_2d[1])
a = (term1_a + term2_a) / denominator
print(f"  (p2_2d[1] - p3_2d[1]) * (x - p3_2d[0])")
print(f"  = ({p2_2d[1]:.3f} - {p3_2d[1]:.3f}) * ({x_param:.3f} - {p3_2d[0]:.3f})")
print(f"  = {p2_2d[1] - p3_2d[1]:.3f} * {x_param - p3_2d[0]:.3f}")
print(f"  = {term1_a:.3f}")
print(f"  +")
print(f"  (p3_2d[0] - p2_2d[0]) * (y - p3_2d[1])")
print(f"  = ({p3_2d[0]:.3f} - {p2_2d[0]:.3f}) * ({y_param:.3f} - {p3_2d[1]:.3f})")
print(f"  = {p3_2d[0] - p2_2d[0]:.3f} * {y_param - p3_2d[1]:.3f}")
print(f"  = {term2_a:.3f}")
print(f"  a = ({term1_a:.3f} + {term2_a:.3f}) / {denominator:.3f} = {a:.6f}")

print(f"\nBarycentric coordinate 'b' calculation:")
term1_b = (p3_2d[1] - p1_2d[1]) * (x_param - p3_2d[0])
term2_b = (p1_2d[0] - p3_2d[0]) * (y_param - p3_2d[1])
b = (term1_b + term2_b) / denominator
print(f"  (p3_2d[1] - p1_2d[1]) * (x - p3_2d[0])")
print(f"  = ({p3_2d[1]:.3f} - {p1_2d[1]:.3f}) * ({x_param:.3f} - {p3_2d[0]:.3f})")
print(f"  = {p3_2d[1] - p1_2d[1]:.3f} * {x_param - p3_2d[0]:.3f}")
print(f"  = {term1_b:.3f}")
print(f"  +")
print(f"  (p1_2d[0] - p3_2d[0]) * (y - p3_2d[1])")
print(f"  = ({p1_2d[0]:.3f} - {p3_2d[0]:.3f}) * ({y_param:.3f} - {p3_2d[1]:.3f})")
print(f"  = {p1_2d[0] - p3_2d[0]:.3f} * {y_param - p3_2d[1]:.3f}")
print(f"  = {term2_b:.3f}")
print(f"  b = ({term1_b:.3f} + {term2_b:.3f}) / {denominator:.3f} = {b:.6f}")

c = 1 - a - b
print(f"\nBarycentric coordinate 'c':")
print(f"  c = 1 - a - b = 1 - {a:.6f} - {b:.6f} = {c:.6f}")

print(f"\nSum = {a + b + c:.6f} (should be 1.0)")

print(f"\nExpected for triangle center:")
print(f"  a = b = c = 0.333333")
print(f"\nActual values are WRONG!")
print("=" * 70)
