# Matrix Transpose Inconsistencies Between MATLAB and Python

## Summary

**CRITICAL BUG FOUND:** The Python implementation incorrectly uses matrix transposes in multiple locations, contradicting the MATLAB reference code. This is likely **Bug #4** causing large errors in test results.

## Background: CoordTrans Function

Both MATLAB and Python implement the same coordinate transformation:

**MATLAB** (`TDstressHS.m:436`):
```matlab
r = A*[x1';x2';x3'];
```

**Python** (`td_utils.py:41`):
```python
r = A @ np.vstack([x1, x2, x3])
```

Both are equivalent: matrix A multiplies stacked coordinate vectors.

**Key principle from MATLAB documentation** (line 428-431):
> "A" is the transformation matrix, whose **columns e1,e2 and e3** are the unit base vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given in X1X2X3. **The transpose of A (i.e., A') will transform the coordinates from X1X2X3 into x1x2x3.**

## Issue #1: TDstress_HarFunc - WRONG TRANSPOSE

### MATLAB Code (TDstressHS.m:380-381)
```matlab
% Transform slip vector components from TDCS into EFCS
A = [Vnorm Vstrike Vdip];  % Columns are Vnorm, Vstrike, Vdip
[bX,bY,bZ] = CoordTrans(bx,by,bz,A);  % Use A directly (NO transpose)
```

### Python Code (tdstress_hs.py:162-163)
```python
# Transform slip vector components from TDCS into EFCS
# A matrix has Vnorm, Vstrike, Vdip as columns
A = np.column_stack([Vnorm, Vstrike, Vdip])  # Columns are Vnorm, Vstrike, Vdip
bX, bY, bZ = coord_trans(bx, by, bz, A.T)  # ❌ WRONG! Uses A.T (transpose)
```

**Problem:** Python uses `A.T` when MATLAB uses `A` directly.

**Impact:** This transforms in the OPPOSITE direction! Since we're converting slip from TDCS to EFCS, this bug completely reverses the transformation.

**Fix Required:**
```python
# CORRECT:
bX, bY, bZ = coord_trans(bx, by, bz, A)  # Remove .T
```

## Issue #2: AngSetupFSC_S - WRONG TRANSPOSE (3 locations)

### MATLAB Code (TDstressHS.m:517-528)
```matlab
A = [ey1,ey2,ey3]; % Transformation matrix, columns are ey1, ey2, ey3

% Transform coordinates from EFCS to the first ADCS
[y1A,y2A,y3A] = CoordTrans(X-PA(1),Y-PA(2),Z-PA(3),A);  % Use A directly

% Transform coordinates from EFCS to the second ADCS
[y1AB,y2AB,y3AB] = CoordTrans(SideVec(1),SideVec(2),SideVec(3),A);  % Use A directly

% Transform slip vector components from EFCS to ADCS
[b1,b2,b3] = CoordTrans(bX,bY,bZ,A);  % Use A directly
```

### Python Code (ang_setup_fsc.py:89-103)
```python
# Transformation matrix: columns are ey1, ey2, ey3
A = np.column_stack([ey1, ey2, ey3])

# Transform coordinates from EFCS to the first ADCS (point A)
y1A, y2A, y3A = coord_trans(X - PA[0], Y - PA[1], Z - PA[2], A.T)  # ❌ WRONG!

# Transform side vector to ADCS
y1AB, y2AB, y3AB = coord_trans(SideVec[0], SideVec[1], SideVec[2], A.T)  # ❌ WRONG!

# Transform slip vector components from EFCS to ADCS
b1, b2, b3 = coord_trans(bX, bY, bZ, A.T)  # ❌ WRONG!
```

**Problem:** Python uses `A.T` in all three locations when MATLAB uses `A` directly.

**Impact:** All coordinate transformations in the free surface correction are reversed! This affects:
- Observation point positions in ADCS
- Side vector in ADCS
- Slip vector in ADCS

This completely invalidates the harmonic function calculations.

**Fix Required:**
```python
# CORRECT (remove all .T):
y1A, y2A, y3A = coord_trans(X - PA[0], Y - PA[1], Z - PA[2], A)
y1AB, y2AB, y3AB = coord_trans(SideVec[0], SideVec[1], SideVec[2], A)
b1, b2, b3 = coord_trans(bX, bY, bZ, A)
```

## Correct Usage: TDstressFS

For comparison, here's where the Python code is CORRECT:

### MATLAB Code (TDstressHS.m:183-184)
```matlab
A = [Vnorm Vstrike Vdip]';  % Note the transpose! Rows become Vnorm, Vstrike, Vdip
[x,y,z] = CoordTrans(X'-P2(1),Y'-P2(2),Z'-P2(3),A);  % Use A directly
```

### Python Code (tdstress_fs.py:76-83)
```python
# Transformation matrix (rows are unit vectors)
A = np.array([Vnorm, Vstrike, Vdip])  # Rows are Vnorm, Vstrike, Vdip
# Transform coordinates from EFCS into TDCS
x, y, z = coord_trans(X - P2[0], Y - P2[1], Z - P2[2], A)  # ✅ CORRECT! Uses A directly
```

**Why this is correct:**
- MATLAB: Stacks as columns then transposes → rows are [Vnorm, Vstrike, Vdip]
- Python: Directly creates rows as [Vnorm, Vstrike, Vdip]
- Both use A directly (no transpose) in CoordTrans/coord_trans

## Root Cause Analysis

The confusion likely arose from Python's common convention of using `A.T` for inverse transformations. However, the MATLAB code is **inconsistent** in how it constructs A:

1. **TDstressFS**: Uses `A = [...]'` (transpose after construction)
2. **TDstress_HarFunc**: Uses `A = [...]` (no transpose)
3. **AngSetupFSC_S**: Uses `A = [...]` (no transpose)

The Python translator apparently assumed "if A has columns as unit vectors, then use A.T", but this is WRONG. The MATLAB code uses A directly when the columns are the unit vectors.

## Impact on Test Results

These matrix transpose errors explain the large errors observed in testing:

From `TEST_RESULTS.md`:
- Point 3: **3961% error** (worst case)
- Point 14: **1773% error**
- Point 7: **2402% error**
- Point 8: **677% error**

When coordinate transformations are reversed:
- Slip vectors point in wrong directions
- Observation points are in wrong positions relative to dislocation
- Stress/strain calculations use completely wrong geometric configuration

## Files Requiring Fixes

1. **python_tdstress/tdstress_hs.py:163** - Remove `.T` from slip vector transformation
2. **python_tdstress/ang_setup_fsc.py:92** - Remove `.T` from point A transformation
3. **python_tdstress/ang_setup_fsc.py:95** - Remove `.T` from side vector transformation
4. **python_tdstress/ang_setup_fsc.py:103** - Remove `.T` from slip vector transformation

## Testing Strategy

After fixing these issues:
1. Re-run `test_tdstress.py` with all 15 points
2. Compare results to expected values from corrected Fortran
3. Verify that errors decrease dramatically
4. Check that spurious NaNs remain (those are from Bug #2 and Bug #3, separate issues)

## References

- MATLAB source: `/home/user/TriBIE/NikkhooWalter2015/TDstressHS.m`
- Python source: `/home/user/TriBIE/python_tdstress/`
- Test results: `python_tdstress/TEST_RESULTS.md`

## Priority

**CRITICAL** - This bug affects all half-space calculations and causes errors up to 3961%. Must be fixed before any other debugging.
