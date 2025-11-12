# Python TDstress Test Results

## Test Configuration
- Triangle vertices: P1=(-1,-1,-5), P2=(1,-1,-5), P3=(-1,1,-4)
- Slip: Ss=1.0, Ds=-1.0, Ts=2.0
- Material: mu=3e10, lambda=3e10

## Full-Space Results (TDstressFS)

**Coordinates now match Fortran test exactly**

| Point | Coordinates | Expected Exx | Got Exx | Status |
|-------|-------------|--------------|---------|--------|
| 1 | (-0.333, -0.333, -3.0) | - | 0.06317170 | ✓ OK |
| 2 (center) | (-0.333, -0.333, -4.667) | 0.04810470 | -0.02920989 | ❌ FAIL |
| 3 | (-0.333, -0.333, -6.0) | - | 2.22069477 | ✓ OK |
| 4 | (7.0, -1.0, -5.0) | 0.00082916 | NaN | ❌ NaN |
| 5 | (-7.0, -1.0, -5.0) | 0.00114440 | NaN | ❌ NaN |
| 6 | (-1.0, -3.0, -6.0) | - | NaN | ⚠️ NaN |
| 7 | (-1.0, 3.0, -3.0) | - | 0.05612467 | ✓ OK |
| 8 | (3.0, -3.0, -6.0) | - | 0.00548525 | ✓ OK |
| 9 | (-3.0, 3.0, -3.0) | - | NaN | ⚠️ NaN |
| 10 | (-1.0, -1.0, -1.0) | - | 0.00322825 | ✓ OK |
| 11 | (-1.0, 1.0, -1.0) | - | -0.00419984 | ✓ OK |
| 12 | (1.0, -1.0, -1.0) | 0.00441203 | NaN | ❌ NaN |
| 13 | (-1.0, -1.0, -8.0) | - | 0.02069253 | ✓ OK |
| 14 | (-1.0, 1.0, -8.0) | - | 0.16414011 | ✓ OK |
| 15 | (1.0, -1.0, -8.0) | -0.00091411 | NaN | ❌ NaN |

## Half-Space Results (TDstressHS)

| Point | Coordinates | Expected Exx | Got Exx | Status |
|-------|-------------|--------------|---------|--------|
| 2 (center) | (-0.333, -0.333, -4.667) | 0.04810470 | -0.03162771 | ❌ FAIL |
| 4 | (7.0, -1.0, -5.0) | 0.00082916 | NaN | ❌ NaN |
| 5 | (-7.0, -1.0, -5.0) | 0.00114440 | NaN | ❌ NaN |
| 12 | (1.0, -1.0, -1.0) | 0.00441203 | NaN | ❌ NaN |
| 15 | (1.0, -1.0, -8.0) | -0.00091411 | NaN | ❌ NaN |

## Issues Identified

### 1. Center Point Wrong Value
**Point 1** (triangle center) returns wrong value, not NaN:
- Full-space: -0.0292 (expected 0.0481) - 161% error
- Half-space: -0.0316 (expected 0.0481) - 166% error

This matches **Bug #2** in Fortran: incorrect barycentric coordinate formula.

### 2. Spurious NaNs at Valid Points
**Points 4, 5, 6, 9, 12, 15** return NaN when they should be finite:
- Point 4: (7, -1, -5) - Far right
- Point 5: (-7, -1, -5) - Far left
- Point 6: (-1, -3, -6) - Below triangle
- Point 9: (-3, 3, -3) - Outside triangle
- Point 12: (1, -1, -1) - Above triangle
- Point 15: (1, -1, -8) - Below triangle

All these points are well away from triangle edges/vertices and should NOT be singular.

Likely causes:
- Division by near-zero in angular dislocation calculations
- Incorrect trimode classification
- Overflow in intermediate calculations

This matches **Bug #3** in Fortran: overly strict edge detection.

## Comparison to Fortran Bugs

The Python implementation exhibits the **SAME BUGS** as the original unpatched Fortran code:

1. ✅ **Bug #1** (uninitialized variables) - Not applicable in Python (no uninitialized variables)
2. ❌ **Bug #2** (wrong barycentric indices) - **PRESENT** in Python (center point wrong)
3. ❌ **Bug #3** (overly strict edge detection) - **PRESENT** in Python (spurious NaNs)
4. ❓ **Bug #4** (matrix orientation) - May be present but hard to isolate

## Recommended Fixes for Python

Apply the same fixes as in Fortran:

1. **Fix barycentric coordinate formula** in `td_utils.py` `trimodefinder()`
   - Check indices when translating from 2D MATLAB to 3D Python arrays

2. **Simplify edge detection** in `td_utils.py` `trimodefinder()`
   - Remove overly strict bounds checking
   - Match MATLAB logic exactly: `a < BARY_TOL and b >= 0 and c >= 0`

3. **Add epsilon checks** in angular dislocation calculations
   - Prevent division by near-zero
   - Add range checks for distant points

4. **Verify matrix orientations** in coordinate transformations
   - Check if vectors should be rows vs columns
   - Ensure transpose operations match MATLAB

## Next Steps

1. Apply Fortran bug fixes to Python code
2. Re-run tests to verify corrections
3. Add unit tests for barycentric coordinates
4. Add numerical stability checks
