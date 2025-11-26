# Array Extraction Bug Found: Bug #2 Identified

## Summary

**CRITICAL BUG FOUND:** Python `trimodefinder()` extracts the **wrong elements** from vertex arrays, causing incorrect barycentric coordinate calculations. This is **Bug #2** from the Fortran analysis.

## The Bug

### MATLAB Implementation (CORRECT)

**Call to trimodefinder:**
```matlab
Trimode = trimodefinder(y,z,x,p1(2:3),p2(2:3),p3(2:3));
```

Where:
- `p1 = [x_TDCS, y_TDCS, z_TDCS]` (3-element array)
- `p1(2:3)` extracts elements 2 and 3 = `[y_TDCS, z_TDCS]`
- Same for p2 and p3

**Inside trimodefinder:** Uses 2D projection onto y-z plane of TDCS
- `p1(1)` = y_TDCS component of vertex 1
- `p1(2)` = z_TDCS component of vertex 1

### Python Implementation (WRONG)

**Call to trimodefinder:**
```python
Trimode = trimodefinder(y, z, x, p1, p2, p3)
```

Where:
- `p1 = [x_TDCS, y_TDCS, z_TDCS]` (3-element array)
- Passes **full array** (not sliced)

**Inside trimodefinder (td_utils.py:144-146):**
```python
# Extract 2D coordinates (y and z components)
p1_2d = np.array(p1[:2])  # ❌ WRONG! Gets [x_TDCS, y_TDCS]
p2_2d = np.array(p2[:2])  # ❌ WRONG! Gets [x_TDCS, y_TDCS]
p3_2d = np.array(p3[:2])  # ❌ WRONG! Gets [x_TDCS, y_TDCS]
```

**Problem:**
- Python extracts `p1[:2]` = first 2 elements = `[p1[0], p1[1]]` = `[x_TDCS, y_TDCS]`
- Should extract `p1[1:3]` = last 2 elements = `[p1[1], p1[2]]` = `[y_TDCS, z_TDCS]`

## Impact

The barycentric coordinate formula (lines 152-153):
```python
a = ((p2_2d[1] - p3_2d[1]) * (x - p3_2d[0]) +
     (p3_2d[0] - p2_2d[0]) * (y - p3_2d[1])) / denominator
```

Uses:
- `p2_2d[0]` = x_TDCS of vertex 2 (WRONG, should be y_TDCS)
- `p2_2d[1]` = y_TDCS of vertex 2 (WRONG, should be z_TDCS)
- Same for p1_2d and p3_2d

This causes the barycentric coordinates to be calculated in the **wrong 2D plane** (x-y plane instead of y-z plane), leading to:
- Wrong trimode classification
- Wrong configuration selection
- Completely incorrect stress/strain calculations

### Evidence from Test Results

**Point 2 (triangle center):**
- Expected: -0.244189
- Got: -0.029210 (Full-space) and -0.032433 (Half-space)
- Error: 88.0% and 86.7%

The center point should have barycentric coordinates (1/3, 1/3, 1/3), but due to wrong plane projection, it gets incorrect coordinates, leading to massive errors.

## The Fix

### Option 1: Fix extraction in trimodefinder (Recommended)

**File:** `python_tdstress/td_utils.py` lines 144-146

**Current (WRONG):**
```python
# Extract 2D coordinates (y and z components)
p1_2d = np.array(p1[:2])  # Gets [x_TDCS, y_TDCS]
p2_2d = np.array(p2[:2])
p3_2d = np.array(p3[:2])
```

**Fixed (CORRECT):**
```python
# Extract 2D coordinates (y and z components)
p1_2d = np.array(p1[1:3])  # Gets [y_TDCS, z_TDCS]
p2_2d = np.array(p2[1:3])  # Gets [y_TDCS, z_TDCS]
p3_2d = np.array(p3[1:3])  # Gets [y_TDCS, z_TDCS]
```

### Option 2: Fix call sites to pass sliced arrays

**File:** `python_tdstress/tdstress_fs.py` line 98

**Current:**
```python
Trimode = trimodefinder(y, z, x, p1, p2, p3)
```

**Alternative fix:**
```python
Trimode = trimodefinder(y, z, x, p1[1:3], p2[1:3], p3[1:3])
```

**Recommendation:** Use Option 1 to fix it in one place and make it robust.

## Verification

After fixing, the center point (Point 2) should give correct barycentric coordinates:
- Should be (1/3, 1/3, 1/3)
- Should yield correct Exx ≈ -0.244189
- Error should drop from 88% to near 0%

## Related Files

- Bug location: `python_tdstress/td_utils.py:144-146`
- Call site: `python_tdstress/tdstress_fs.py:98`
- Test results: `python_tdstress/TEST_RESULTS.md`
- MATLAB reference: `/home/user/TriBIE/NikkhooWalter2015/TDstressHS.m:208`

## Priority

**CRITICAL** - This bug affects ALL calculations by using wrong barycentric coordinates, causing 88% error even at triangle center. Must be fixed immediately.

## Root Cause

The Python translator incorrectly assumed that passing full 3-element arrays and extracting `[:2]` would give the correct 2D projection. However, MATLAB explicitly passes `p1(2:3)` (last 2 elements), not `p1(1:2)` (first 2 elements).

This is a classic 0-based vs 1-based indexing translation error compounded by not recognizing that MATLAB was selecting elements 2-3, not 1-2.
