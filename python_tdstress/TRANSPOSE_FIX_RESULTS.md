# Results After Matrix Transpose Fixes

## Summary

Fixed 4 matrix transpose errors (removed `.T` from `coord_trans()` calls), but **large errors persist** indicating additional bugs beyond matrix transposes.

## Comparison: Before vs After Transpose Fixes

### Half-Space Results (TDstressHS)

| Point | Coordinates | Expected | Before Fix | After Fix | Error Before | Error After | Change |
|-------|-------------|----------|------------|-----------|--------------|-------------|---------|
| 1 | (-0.33, -0.33, -3.0) | 0.0481 | 0.0636 | 0.0625 | 32.2% | **29.9%** | ✓ Improved |
| 2 (center) | (-0.33, -0.33, -4.67) | -0.2442 | -0.0316 | -0.0324 | 87.0% | **86.7%** | ✓ Slight |
| 3 | (-0.33, -0.33, -6.0) | 0.0547 | 2.2174 | 2.2168 | 3955% | **3954%** | ≈ Same |
| 4 | (7.0, -1.0, -5.0) | 0.0008 | NaN | NaN | - | - | No change |
| 5 | (-7.0, -1.0, -5.0) | 0.0011 | NaN | NaN | - | - | No change |
| 6 | (-1.0, -3.0, -6.0) | -0.0039 | NaN | NaN | - | - | No change |
| 7 | (-1.0, 3.0, -3.0) | -0.0024 | 0.0488 | 0.0495 | 2103% | **2131%** | ✗ Worse |
| 8 | (3.0, -3.0, -6.0) | 0.0007 | -0.0109 | -0.0117 | 1647% | **1763%** | ✗ Worse |
| 9 | (-3.0, 3.0, -3.0) | 0.0002 | NaN | NaN | - | - | No change |
| 10 | (-1.0, -1.0, -1.0) | 0.0065 | -0.0116 | -0.0116 | 278% | **278%** | ≈ Same |
| 11 | (-1.0, 1.0, -1.0) | 0.0009 | -0.0113 | -0.0104 | 1325% | **1230%** | ✓ Improved |
| 12 | (1.0, -1.0, -1.0) | 0.0044 | NaN | NaN | - | - | No change |
| 13 | (-1.0, -1.0, -8.0) | 0.0033 | 0.0134 | 0.0131 | 306% | **296%** | ✓ Improved |
| 14 | (-1.0, 1.0, -8.0) | 0.0088 | 0.1583 | 0.1582 | 1706% | **1705%** | ≈ Same |
| 15 | (1.0, -1.0, -8.0) | -0.0009 | NaN | NaN | - | - | No change |

### Summary of Changes

- **Improvements:** Points 1, 2, 11, 13 (4 points improved)
- **Degradations:** Points 7, 8 (2 points worse)
- **No change:** Points 3, 10, 14 and all NaN points (9 points unchanged)

## Key Observations

### 1. Transpose Fixes Alone Are Insufficient

The matrix transpose fixes provided only **marginal improvements**:
- Best improvement: Point 1 (32.2% → 29.9%, only 2.3% reduction)
- Point 2 (center): Still 86.7% error despite being inside triangle
- Point 3: Still **3954% error** (essentially unchanged)
- Points 7, 8: Actually got **worse** after fixes

### 2. NaN Issues Persist

All 6 points returning NaN still return NaN after transpose fixes:
- Points 4, 5, 6, 9, 12, 15 unchanged
- This confirms NaNs are caused by **Bug #2** (barycentric) and **Bug #3** (edge detection)
- Not related to matrix transposes

### 3. Full-Space Unchanged

Full-space (TDstressFS) results are identical before and after transpose fixes:
- Expected: transpose errors were only in half-space harmonic function
- Confirms the fix targeted the correct code sections

### 4. Errors Still Massive

After transpose fixes, errors remain **unacceptably large**:
- 8 points: >100% error
- 6 points: >1000% error
- 1 point: >3000% error
- Only points 1 and 2 have <100% error (but still 30-87% wrong)

## Remaining Bugs to Fix

### Priority 1: Bug #2 - Barycentric Coordinate Formula

**Evidence:** Point 2 (triangle center) still has 86.7% error
- Center point should have barycentric coords (1/3, 1/3, 1/3)
- Getting -0.0324 instead of expected -0.2442
- Same error pattern as original Fortran Bug #2

**Location:** `td_utils.py:152-161` in `trimodefinder()`

**Fix needed:** Correct the barycentric coordinate calculation formula

### Priority 2: Bug #3 - Edge Detection

**Evidence:** 6 spurious NaNs at valid points far from triangle
- Points 4, 5, 6, 9, 12, 15 all return NaN
- All are >1 unit away from triangle edges
- Should NOT be classified as singular

**Location:** `td_utils.py:162-176` in `trimodefinder()`

**Fix needed:** Simplify edge detection logic (remove overly strict conditions)

### Priority 3: Additional Formula Errors

**Evidence:** Point 3 has 3954% error even after transpose fixes
- Cannot be explained by transposes or barycentric coords alone
- May indicate:
  - Sign errors in formulas
  - Missing terms in strain calculations
  - Incorrect tensor transformations
  - Wrong boundary condition handling

**Requires:** Detailed line-by-line comparison with MATLAB code

## Conclusion

The matrix transpose fixes were **necessary but not sufficient**:
1. ✅ Correctly identified and fixed
2. ✅ Provided marginal improvements (2-10% in some cases)
3. ❌ Did not resolve the fundamental calculation errors
4. ❌ Large errors (30-3954%) persist

**Next steps:**
1. Fix Bug #2 (barycentric coordinates) - should eliminate center point error
2. Fix Bug #3 (edge detection) - should eliminate 6 spurious NaNs
3. Deep comparison with MATLAB for remaining formula errors
4. Re-test after each fix to isolate impact

## Test Command

```bash
cd /home/user/TriBIE/python_tdstress
python test_tdstress.py
```

## Related Files

- Analysis: `MATRIX_TRANSPOSE_ISSUES.md`
- Test script: `test_tdstress.py`
- Fixed files: `tdstress_hs.py`, `ang_setup_fsc.py`
- Bug location: `td_utils.py` (barycentric and edge detection)
