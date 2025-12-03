# Python TDstress Test Results

## Test Configuration
- Triangle vertices: P1=(-1,-1,-5), P2=(1,-1,-5), P3=(-1,1,-4)
- Slip: Ss=1.0, Ds=-1.0, Ts=2.0
- Material: mu=3e10, lambda=3e10

## Full-Space Results (TDstressFS)

All 15 test points from corrected Fortran reference:

| Point | Coordinates | Expected Exx | Got Exx | Rel. Error | Status |
|-------|-------------|--------------|---------|------------|--------|
| 1 | (-0.333, -0.333, -3.0) | 0.04810470 | 0.06317170 | 31.3% | ❌ FAIL |
| 2 (center) | (-0.333, -0.333, -4.667) | -0.24418898 | -0.02920989 | 88.0% | ❌ FAIL |
| 3 | (-0.333, -0.333, -6.0) | 0.05468314 | 2.22069477 | 3961% | ❌ FAIL |
| 4 | (7.0, -1.0, -5.0) | 0.00082916 | NaN | - | ❌ NaN |
| 5 | (-7.0, -1.0, -5.0) | 0.00114440 | NaN | - | ❌ NaN |
| 6 | (-1.0, -3.0, -6.0) | -0.00386292 | NaN | - | ❌ NaN |
| 7 | (-1.0, 3.0, -3.0) | -0.00243789 | 0.05612467 | 2402% | ❌ FAIL |
| 8 | (3.0, -3.0, -6.0) | 0.00070640 | 0.00548525 | 677% | ❌ FAIL |
| 9 | (-3.0, 3.0, -3.0) | 0.00021125 | NaN | - | ❌ NaN |
| 10 | (-1.0, -1.0, -1.0) | 0.00650801 | 0.00322825 | 50.4% | ❌ FAIL |
| 11 | (-1.0, 1.0, -1.0) | 0.00092245 | -0.00419984 | 555% | ❌ FAIL |
| 12 | (1.0, -1.0, -1.0) | 0.00441203 | NaN | - | ❌ NaN |
| 13 | (-1.0, -1.0, -8.0) | 0.00330232 | 0.02069253 | 527% | ❌ FAIL |
| 14 | (-1.0, 1.0, -8.0) | 0.00876399 | 0.16414011 | 1773% | ❌ FAIL |
| 15 | (1.0, -1.0, -8.0) | -0.00091411 | NaN | - | ❌ NaN |

**Summary:** 0 PASS, 9 FAIL, 6 NaN

## Half-Space Results (TDstressHS)

| Point | Coordinates | Expected Exx | Got Exx | Rel. Error | Status |
|-------|-------------|--------------|---------|------------|--------|
| 1 | (-0.333, -0.333, -3.0) | 0.04810470 | 0.06359737 | 32.2% | ❌ FAIL |
| 2 (center) | (-0.333, -0.333, -4.667) | -0.24418898 | -0.03162771 | 87.0% | ❌ FAIL |
| 3 | (-0.333, -0.333, -6.0) | 0.05468314 | 2.21736029 | 3955% | ❌ FAIL |
| 4 | (7.0, -1.0, -5.0) | 0.00082916 | NaN | - | ❌ NaN |
| 5 | (-7.0, -1.0, -5.0) | 0.00114440 | NaN | - | ❌ NaN |
| 6 | (-1.0, -3.0, -6.0) | -0.00386292 | NaN | - | ❌ NaN |
| 7 | (-1.0, 3.0, -3.0) | -0.00243789 | 0.04884027 | 2103% | ❌ FAIL |
| 8 | (3.0, -3.0, -6.0) | 0.00070640 | -0.01092830 | 1647% | ❌ FAIL |
| 9 | (-3.0, 3.0, -3.0) | 0.00021125 | NaN | - | ❌ NaN |
| 10 | (-1.0, -1.0, -1.0) | 0.00650801 | -0.01161042 | 278% | ❌ FAIL |
| 11 | (-1.0, 1.0, -1.0) | 0.00092245 | -0.01129609 | 1325% | ❌ FAIL |
| 12 | (1.0, -1.0, -1.0) | 0.00441203 | NaN | - | ❌ NaN |
| 13 | (-1.0, -1.0, -8.0) | 0.00330232 | 0.01341779 | 306% | ❌ FAIL |
| 14 | (-1.0, 1.0, -8.0) | 0.00876399 | 0.15831633 | 1706% | ❌ FAIL |
| 15 | (1.0, -1.0, -8.0) | -0.00091411 | NaN | - | ❌ NaN |

**Summary:** 0 PASS, 9 FAIL, 6 NaN

## Critical Issues Identified

### 1. Spurious NaNs (6 points)
**Points 4, 5, 6, 9, 12, 15** return NaN when they should be finite:
- Point 4: (7, -1, -5) - Far right of triangle
- Point 5: (-7, -1, -5) - Far left of triangle
- Point 6: (-1, -3, -6) - Below triangle
- Point 9: (-3, 3, -3) - Outside triangle
- Point 12: (1, -1, -1) - Above triangle
- Point 15: (1, -1, -8) - Below triangle

All these points are well away from triangle edges/vertices and should NOT be singular.

**Root causes:**
- Division by zero in angular dislocation calculations
- Incorrect trimode classification
- Overly strict edge detection (Bug #3 from Fortran)

### 2. Large Errors on All Valid Points
**All 9 non-NaN points have significant errors:**
- Point 3: 3961% error (worst case)
- Points 7, 8, 11, 13, 14: 500-2400% error
- Point 2 (center): 88% error (matches Bug #2 from Fortran)
- Points 1, 10: 31-50% error (smallest but still significant)

**Root causes:**
- Wrong barycentric coordinate formula (Bug #2 from Fortran)
- Incorrect matrix orientations (Bug #4 from Fortran)
- Sign errors in coordinate transformations
- Formula transcription errors from MATLAB

### 3. Triangle Center Point Wrong
**Point 2** (triangle center at barycentric coordinates 1/3, 1/3, 1/3):
- Full-space: -0.0292 vs expected -0.244 (88% error)
- Half-space: -0.0316 vs expected -0.244 (87% error)

This matches **Bug #2** in Fortran: incorrect barycentric coordinate formula in `trimodefinder()`.

## Comparison to Fortran Bugs

The Python implementation exhibits the **SAME BUGS** as the original unpatched Fortran code:

1. ✅ **Bug #1** (uninitialized variables) - Not applicable in Python
2. ❌ **Bug #2** (wrong barycentric indices) - **PRESENT** (center point wrong)
3. ❌ **Bug #3** (overly strict edge detection) - **PRESENT** (6 spurious NaNs)
4. ❌ **Bug #4** (matrix orientation) - **LIKELY PRESENT** (large errors everywhere)

## Required Fixes

Apply the same fixes as in corrected Fortran code:

### Fix #1: Barycentric Coordinate Formula
**File:** `python_tdstress/td_utils.py` function `trimodefinder()`

The barycentric coordinate calculation likely has wrong indices when translating from MATLAB 2D arrays to Python 3D arrays. Check the formula:
```python
a = ((p2_2d[1] - p3_2d[1]) * (x - p3_2d[0]) +
     (p3_2d[0] - p2_2d[0]) * (y - p3_2d[1])) / denom
```

### Fix #2: Edge Detection Logic
**File:** `python_tdstress/td_utils.py` function `trimodefinder()`

Simplify edge detection to match MATLAB exactly:
```python
# Current (overly strict):
if a < bary_tol and b >= 0 and c >= 0 and ...:

# Should be (MATLAB-like):
if a < bary_tol and b >= 0 and c >= 0:
```

### Fix #3: Numerical Stability
**File:** `python_tdstress/ang_dislocation.py`

Add epsilon checks before divisions to prevent spurious NaNs:
```python
if abs(Wr) < eps:
    # Handle near-zero case
```

### Fix #4: Matrix Orientations
**Files:** `ang_setup_fsc.py`, `tdstress_hs.py`, `td_utils.py`

Verify all coordinate transformations:
- Check if transformation matrices should be transposed
- Verify vector orientations (row vs column)
- Match MATLAB's matrix conventions exactly

## Warnings During Execution

Multiple runtime warnings indicate numerical instabilities:
- `RuntimeWarning: divide by zero encountered in divide`
- `RuntimeWarning: invalid value encountered in divide`
- `RuntimeWarning: invalid value encountered in subtract/multiply/add`

These warnings appear in:
- `td_utils.py:152, 155` - Barycentric coordinate calculation
- `ang_dislocation.py:58-110` - Multiple locations in strain calculations

## Next Steps

1. ✅ Run comprehensive test with all 15 reference points
2. ⏭️ Apply Bug #2 fix (barycentric coordinates) from Fortran to Python
3. ⏭️ Apply Bug #3 fix (edge detection) from Fortran to Python
4. ⏭️ Apply Bug #4 fix (matrix orientation) from Fortran to Python
5. ⏭️ Add numerical stability checks (epsilon guards)
6. ⏭️ Re-run tests to verify corrections
7. ⏭️ Add unit tests for barycentric coordinates
8. ⏭️ Document remaining discrepancies if any

## Test Environment

- Python version: (from environment)
- NumPy version: (from environment)
- Test date: 2025-11-12
- Reference: Corrected Fortran implementation with all 4 bugs fixed
