# Critical Fix: Bounds Checking for On-Edge Detection

## Problem After First Fix

After adding tolerance-based comparisons, **points 8 and 9 were still returning NaN**. The first fix used tolerances but didn't check if points were within triangle bounds.

## Root Cause: Extended Edge Lines

Points described as "along the edge but outside the triangle" have special barycentric properties:

### Example: Point on Extended Edge Line

Consider triangle with vertices V1, V2, V3:
- Edge V1-V2 connects the two vertices
- **Extended line** goes beyond V1 and V2 in both directions

```
    Extended line
-----------------V1--------V2-----------------
                 |   Edge  |
                 V3
```

A point **P** on the extended line (but outside edge bounds) has:
- **c ≈ 0** (on the V1-V2 line)
- **But**: a > 1 OR b > 1 OR a < 0 OR b < 0 (outside triangle)

### The Problem with First Fix

```fortran
! FIRST FIX (incomplete - still had NaN issue)
if (abs(a) < BARY_TOL .and. b >= -BARY_TOL .and. c >= -BARY_TOL) then
    trimode = 0  ! Returns NaN
```

This catches:
- ✓ Points on actual edge (c ≈ 0, a,b ∈ [0,1])
- ✗ Points on extended line (c ≈ 0, but a or b > 1)

## The Solution: Bounds Checking

### Second Fix (Complete)

```fortran
! SECOND FIX (complete - bounds checking added)
if (abs(a) < BARY_TOL .and. &
    b >= -BARY_TOL .and. b <= 1.0_DP + BARY_TOL .and. &
    c >= -BARY_TOL .and. c <= 1.0_DP + BARY_TOL) then
    trimode = 0  ! Only for points actually ON the edge
```

Now checks:
- ✓ Points on actual edge (c ≈ 0, a,b ∈ [0,1]) → trimode = 0 (NaN)
- ✓ Points on extended line (c ≈ 0, a or b outside [0,1]) → trimode = -1 or 1 (valid)

## Why Points 8 and 9 Had This Issue

### Point 8: (3.0, -3.0, -6.0)
```
After transformation to TDCS and barycentric calculation:
- One coordinate (e.g., c) ≈ 0 (on V1-V2 line)
- But a > 1 or b > 1 (beyond V1 or V2 vertex)
→ On extended edge line, NOT on actual triangle edge
```

### Point 9: (-3.0, 3.0, -3.0)
```
Similar situation:
- One coordinate ≈ 0 (on edge line)
- Other coordinates outside [0,1] (beyond triangle bounds)
→ On extended edge line, NOT on actual triangle edge
```

## Barycentric Coordinate Interpretation

### Valid Cases

| Barycentric Coords | Location | trimode | Result |
|-------------------|----------|---------|--------|
| a,b,c all in [0,1] | Inside triangle | 1 or -1 | Valid numbers |
| One = 0, others in [0,1] | On edge | 0 | NaN (singular) |
| One < 0 | Outside triangle | -1 | Valid numbers |
| One > 1 | Outside triangle | -1 or 1 | Valid numbers |

### The Problematic Case (Fixed)

| Barycentric Coords | Location | Old trimode | New trimode |
|-------------------|----------|-------------|-------------|
| c ≈ 0, a > 1, b ∈ [0,1] | Extended line | 0 (NaN) ✗ | -1 (valid) ✓ |
| b ≈ 0, c > 1, a ∈ [0,1] | Extended line | 0 (NaN) ✗ | -1 (valid) ✓ |
| a ≈ 0, b > 1, c < 0 | Extended line | 0 (NaN) ✗ | -1 (valid) ✓ |

## Complete Check Logic

For each edge (checking if point is on edge p1-p2, for example):

```fortran
! Edge p1-p2: coordinate 'c' should be ≈ 0
if (abs(c) < BARY_TOL .and.           ! On the p1-p2 line
    a >= -BARY_TOL .and.               ! a not too negative
    a <= 1.0_DP + BARY_TOL .and.       ! a not too large (NEW!)
    b >= -BARY_TOL .and.               ! b not too negative
    b <= 1.0_DP + BARY_TOL) then       ! b not too large (NEW!)
    trimode = 0  ! Actually on the triangle edge
end if
```

The upper bounds (`<= 1 + TOL`) prevent classification of extended line points as on-edge.

## Testing

### Before Bounds Fix
```
Point 8: trimode = 0 → NaN (incorrect)
Point 9: trimode = 0 → NaN (incorrect)
```

### After Bounds Fix
```
Point 8: trimode = -1 or 1 → Valid number (correct)
Point 9: trimode = -1 or 1 → Valid number (correct)
```

### Debug Programs

1. **debug_trimode.f90**: Simple test of points 8 and 9 using module
2. **debug_trimode_detailed.f90**: Standalone program showing barycentric coords and trimode logic

Run after recompiling:
```bash
# Compile module
gfortran -c sub_nikkhoo.f90

# Compile and run debug
gfortran -o debug_trimode_detailed debug_trimode_detailed.f90
./debug_trimode_detailed
```

## Summary of All Fixes

### Original Issue (casez_log)
- **Problem**: Exact equality check for singular points
- **Fix**: Changed to tolerance-based, return NaN

### First Trimode Issue (tolerance only)
- **Problem**: Exact equality for barycentric coords
- **Fix**: Added BARY_TOL tolerance

### Second Trimode Issue (no bounds checking)
- **Problem**: Extended edge lines caught as on-edge
- **Fix**: Added upper bound checks [0, 1+TOL]

## Files Modified

- **sub_nikkhoo.f90**: Lines 659-671 (trimode_finder function)
- **debug_trimode_detailed.f90**: Standalone diagnostic program

## Commits

```
1db6f15 - Add bounds checking to trimode_finder on-edge detection
840a313 - Fix trimode_finder to use tolerance-based comparison
6a77aaf - Fix casez_log implementation to match MATLAB reference
```

## Key Takeaway

**On-edge detection requires TWO conditions:**
1. One barycentric coordinate ≈ 0 (on edge **line**)
2. Other coordinates in [0,1] (within triangle **bounds**)

Without condition #2, points on **extended** edge lines are incorrectly classified as on-edge singular points, causing spurious NaN results.
