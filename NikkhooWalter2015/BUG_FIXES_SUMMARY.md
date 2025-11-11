# Summary of Critical Bug Fixes

## Bug 1: Uninitialized Variables (FIXED - Commit 5c49655)

### Problem
Variables declared but never calculated in `angdis_strain`:
- `W2`, `Wr`, `W2r`, `Wr3`, `W2r2`, `rz`, `r2z2`, `r3z`

### Symptom
NaN for ALL points due to using uninitialized (garbage) values in formulas

### Fix
Added missing calculations matching MATLAB (TDstressHS.m:614-618):
```fortran
W2 = W * W
Wr = W * r
W2r = W2 * r
Wr3 = W * r3
W2r2 = W2 * r2
rz = r * z
r2z2 = r2 * z2
r3z = r3 * z
```

---

## Bug 2: Wrong Array Indices in Barycentric Coordinates (FIXED - Commit 35c73f3)

### Problem
Incorrect array indexing when translating MATLAB's 2D barycentric formula to Fortran's 3D arrays.

**Array mapping:**
- MATLAB: `p = [p(1), p(2)]` = `[y_coord, z_coord]` (2D)
- Fortran: `p = [p(1), p(2), p(3)]` = `[x_coord, y_coord, z_coord]` (3D)
- Mapping: MATLAB's `p(1)` → Fortran's `p(2)`, MATLAB's `p(2)` → Fortran's `p(3)`

**Wrong formula** (using p(2) where it should be p(3)):
```fortran
denominator = (p2(2)-p3(2))*(p1(2)-p3(2)) + (p3(2)-p2(2))*(p1(3)-p3(3))
a = ((p2(2)-p3(2))*(x-p3(2)) + (p3(2)-p2(2))*(y-p3(3))) / denominator
b = ((p3(2)-p1(2))*(x-p3(2)) + (p1(2)-p3(2))*(y-p3(3))) / denominator
```

**Correct formula:**
```fortran
denominator = (p2(3)-p3(3))*(p1(2)-p3(2)) + (p3(2)-p2(2))*(p1(3)-p3(3))
a = ((p2(3)-p3(3))*(x-p3(2)) + (p3(2)-p2(2))*(y-p3(3))) / denominator
b = ((p3(3)-p1(3))*(x-p3(2)) + (p1(2)-p3(2))*(y-p3(3))) / denominator
```

### Symptom
- Triangle center (-0.333, -0.333, -4.667) classified as `trimode=0` (on edge)
- Barycentric coordinates: **(0.965, 0.0, 0.035)** instead of **(0.333, 0.333, 0.333)**
- `casez_log=TRUE` → returns NaN for that contribution
- Total result becomes NaN

### Root Cause
When `b=0.0` (due to wrong formula), the code thinks the point is on the edge P1-P3 (opposite to vertex P2), so it sets `trimode=0`, which triggers the singular case handling.

### Fix
Corrected array indices to match MATLAB's TDstressHS.m:457-460

---

## Testing After Fixes

### Test 1: Center Point
```bash
cd NikkhooWalter2015
./test_center.sh
```

**Expected:** Finite Exx value (NOT NaN)
**Reason:** Center with barycentric (0.333, 0.333, 0.333) should be `trimode=1` (inside)

### Test 2: Multiple Points
```bash
cd NikkhooWalter2015
gfortran -O0 -g -o test_multiple sub_nikkhoo.f90 test_multiple_points.f90
./test_multiple
```

**Expected:** All 7 points return finite values

### Test 3: Points 8 & 9 (Original Issue)
```bash
cd NikkhooWalter2015
./test_regularization.sh
```

**Expected results:**
- Point 8 (3.0, -3.0, -6.0): Exx = 7.064e-4
- Point 9 (-3.0, 3.0, -3.0): Exx = 2.113e-4

---

## Impact Analysis

### Before Both Fixes
- **Bug 1 alone**: NaN for ALL points (uninitialized variables)
- **Bug 2 alone**: Wrong trimode classification → wrong configuration → potential NaN

### After Both Fixes
- Barycentric coordinates calculated correctly
- Points classified correctly (inside/outside/on-edge)
- Intermediate variables properly initialized
- Code structure matches MATLAB implementation

---

## How These Bugs Occurred

### Bug 1: Accidental Deletion
When removing singularity handling (commit 6a9c7f1), the essential calculations were deleted along with the checks.

### Bug 2: Translation Error
When converting MATLAB (2D arrays) to Fortran (3D arrays), the index mapping was done incorrectly:
- Should have been systematic: MATLAB `p(i)` → Fortran `p(i+1)`
- Instead: Mixed up indices in the formula

---

## Prevention for Future

### Code Review Checklist
✓ Compare with reference implementation (MATLAB) line-by-line
✓ Verify all declared variables are initialized before use
✓ Check array index mappings when converting between languages
✓ Test with simple cases (e.g., triangle center) before edge cases
✓ Use compiler warnings (`-Wall -Wuninitialized`)

### Red Flags to Watch For
- Variables declared but no assignment statement before first use
- Different array dimensions between source and target languages
- "Simple" coordinate transformations with index arithmetic

---

## Verification Against MATLAB

Both fixes were verified by comparing with MATLAB source code:
1. **Bug 1**: TDstressHS.m lines 614-618 (W2, Wr, etc. calculations)
2. **Bug 2**: TDstressHS.m lines 457-460 (barycentric coordinate formula)

The Fortran code should now match MATLAB's behavior for trimode classification and angular dislocation calculations.
