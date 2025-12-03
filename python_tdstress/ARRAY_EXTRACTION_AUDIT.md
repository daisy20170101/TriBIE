# Array Extraction Audit: MATLAB vs Python

## Summary

Comprehensive audit of matrix and array element extractions between MATLAB and Python implementations. **One critical bug found** (Bug #2).

## Audit Results

### ✅ CORRECT Implementations

#### 1. SideVec extraction in AngSetupFSC_S

**MATLAB** (line 513):
```matlab
ey1 = [SideVec(1:2);0];
```
Extracts first 2 elements: `[SideVec(1), SideVec(2), 0]` = `[x, y, 0]`

**Python** (ang_setup_fsc.py:76):
```python
ey1 = np.array([SideVec[0], SideVec[1], 0.0])
```
Extracts first 2 elements: `[SideVec[0], SideVec[1], 0]` = `[x, y, 0]`

**Status:** ✅ **CORRECT** - Matches MATLAB

---

#### 2. Transformation matrix A in TDSetupS

**MATLAB** (line 479):
```matlab
A = [[SideVec(3);-SideVec(2)] SideVec(2:3)]';
```

Breakdown:
- `[SideVec(3);-SideVec(2)]` = column `[z; -y]`
- `SideVec(2:3)` = column `[y; z]`
- Horizontal concatenation: `[z, y; -y, z]`
- Transpose: `[z, -y; y, z]`

Result:
```
A = [ SideVec(3)  -SideVec(2) ]
    [ SideVec(2)   SideVec(3) ]
```

**Python** (ang_dislocation.py:147-148):
```python
A = np.array([[side_vec[2], -side_vec[1]],
              [side_vec[1], side_vec[2]]])
```

With 0-based indexing:
```
A = [ side_vec[2]  -side_vec[1] ]  # [z, -y]
    [ side_vec[1]   side_vec[2] ]  # [y,  z]
```

**Status:** ✅ **CORRECT** - Matches MATLAB

---

#### 3. Transformation matrix B in TDSetupS

**MATLAB** (line 495):
```matlab
B = [[1 0 0];[zeros(2,1),A']];
```

Result:
```
B = [ 1              0             0           ]
    [ 0   SideVec(3)   SideVec(2)  ]
    [ 0  -SideVec(2)   SideVec(3)  ]
```

**Python** (ang_dislocation.py:164-166):
```python
B = np.array([[1, 0, 0],
              [0, A[0, 0], A[1, 0]],
              [0, A[0, 1], A[1, 1]]])
```

Where `A[0,0]=side_vec[2]`, `A[0,1]=-side_vec[1]`, `A[1,0]=side_vec[1]`, `A[1,1]=side_vec[2]`

Result:
```
B = [ 1            0               0            ]
    [ 0   side_vec[2]      side_vec[1]  ]
    [ 0  -side_vec[1]      side_vec[2]  ]
```

**Status:** ✅ **CORRECT** - Matches MATLAB

---

#### 4. Barycentric coordinate formula structure

**MATLAB** (lines 457-458):
```matlab
a = ((p2(2)-p3(2)).*(x-p3(1))+(p3(1)-p2(1)).*(y-p3(2)))./...
    ((p2(2)-p3(2)).*(p1(1)-p3(1))+(p3(1)-p2(1)).*(p1(2)-p3(2)));
```

**Python** (td_utils.py:152-153):
```python
a = ((p2_2d[1] - p3_2d[1]) * (x - p3_2d[0]) +
     (p3_2d[0] - p2_2d[0]) * (y - p3_2d[1])) / denominator
```

**Status:** ✅ **CORRECT** - Formula structure is correct (but see Bug #2 below for wrong input)

---

### ❌ INCORRECT Implementation: Bug #2

#### Vertex array extraction in trimodefinder

**MATLAB** (line 208):
```matlab
Trimode = trimodefinder(y,z,x,p1(2:3),p2(2:3),p3(2:3));
```

Where `p1 = [x_TDCS, y_TDCS, z_TDCS]` (3 elements)
- `p1(2:3)` extracts elements 2 and 3 = `[y_TDCS, z_TDCS]`
- Projects triangle onto **y-z plane** of TDCS

**Python** (tdstress_fs.py:98):
```python
Trimode = trimodefinder(y, z, x, p1, p2, p3)
```

Where `p1 = [x_TDCS, y_TDCS, z_TDCS]` (3 elements)

**Then inside trimodefinder** (td_utils.py:144-146):
```python
# Extract 2D coordinates (y and z components)
p1_2d = np.array(p1[:2])  # ❌ WRONG!
p2_2d = np.array(p2[:2])  # ❌ WRONG!
p3_2d = np.array(p3[:2])  # ❌ WRONG!
```

- `p1[:2]` extracts first 2 elements = `[x_TDCS, y_TDCS]`
- Projects triangle onto **x-y plane** of TDCS (WRONG!)

**Status:** ❌ **INCORRECT** - **BUG #2 IDENTIFIED**

**Impact:**
- Barycentric coordinates calculated in wrong plane
- All trimode classifications incorrect
- Center point: 88% error
- All test points affected

**Fix:**
```python
# CORRECT extraction:
p1_2d = np.array(p1[1:3])  # Gets [y_TDCS, z_TDCS]
p2_2d = np.array(p2[1:3])  # Gets [y_TDCS, z_TDCS]
p3_2d = np.array(p3[1:3])  # Gets [y_TDCS, z_TDCS]
```

---

## Summary Table

| Component | MATLAB | Python | Status |
|-----------|--------|--------|--------|
| ey1 in AngSetupFSC_S | `SideVec(1:2)` | `SideVec[0:2]` | ✅ Correct |
| Matrix A in TDSetupS | Complex transpose | Direct construction | ✅ Correct |
| Matrix B in TDSetupS | With A' embedded | Direct construction | ✅ Correct |
| Barycentric formula | Standard formula | Standard formula | ✅ Correct |
| **Vertex extraction** | **`p1(2:3)`** | **`p1[:2]`** | ❌ **BUG #2** |

## Conclusion

Out of 5 major array/matrix extraction patterns audited:
- ✅ **4 are correct** (80%)
- ❌ **1 is incorrect** (20%) - **Bug #2**

The incorrect extraction is **critical** because it affects the fundamental barycentric coordinate calculation used by all stress/strain computations.

## Priority

**CRITICAL** - Bug #2 must be fixed before any other improvements. This single bug causes 88% error at the triangle center point.

## Next Steps

1. Fix Bug #2 (td_utils.py:144-146) by changing `p1[:2]` to `p1[1:3]`
2. Re-run tests to verify error reduction
3. Continue with Bug #3 audit (edge detection logic)

## Related Documentation

- Bug #2 detailed analysis: `BUG2_ARRAY_EXTRACTION.md`
- Matrix transpose fixes: `MATRIX_TRANSPOSE_ISSUES.md`
- Test results: `TEST_RESULTS.md`
- Post-transpose results: `TRANSPOSE_FIX_RESULTS.md`
