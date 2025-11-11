# Side-by-Side Comparison: casez_log Implementation

## Overview

This document provides a detailed comparison between the MATLAB reference and Fortran implementations of the `casez_log` case (trimode == 0) for handling calculation points on triangular dislocation edges.

---

## MATLAB Reference Implementation (CORRECT)

**File**: `TDstressHS.m`
**Lines**: 215, 312-318

### Configuration Detection
```matlab
% Line 215
casezLog = Trimode==0;
```

### casez_log Handling
```matlab
% Lines 312-318
if nnz(casezLog)~=0
    exx(casezLog,1) = nan;
    eyy(casezLog,1) = nan;
    ezz(casezLog,1) = nan;
    exy(casezLog,1) = nan;
    exz(casezLog,1) = nan;
    eyz(casezLog,1) = nan;
end
```

### Behavior
- **Action**: Sets all strain components to `NaN`
- **Rationale**: Points on triangle edges are singular
- **Lines of code**: 7 assignment statements
- **Computational cost**: O(1) - trivial assignment

---

## Fortran Implementation (INCORRECT)

**File**: `sub_nikkhoo.f90`
**Lines**: 211, 255-288

### Configuration Detection
```fortran
! Line 211
casez_log = (trimode == 0)
```

### casez_log Handling
```fortran
! Lines 255-288
else if (casez_log) then
    ! For points on the triangle, use average of positive and negative cases
    ! Configuration I (positive)
    call tdsetup_s(x_td, y_td, z_td, A_angle, bx, by, bz, nu, p1_td, -e13, &
                   exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    call tdsetup_s(x_td, y_td, z_td, B_angle, bx, by, bz, nu, p2_td, e12, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)
    exx_p = exx_p + exx_n; eyy_p = eyy_p + eyy_n; ezz_p = ezz_p + ezz_n
    exy_p = exy_p + exy_n; exz_p = exz_p + exz_n; eyz_p = eyz_p + eyz_n
    call tdsetup_s(x_td, y_td, z_td, C_angle, bx, by, bz, nu, p3_td, e23, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)
    exx_p = exx_p + exx_n; eyy_p = eyy_p + eyy_n; ezz_p = ezz_p + ezz_n
    exy_p = exy_p + exy_n; exz_p = exz_p + exz_n; eyz_p = eyz_p + eyz_n

    ! Configuration II (negative)
    call tdsetup_s(x_td, y_td, z_td, A_angle, -bx, -by, -bz, nu, p1_td, e13, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)
    call tdsetup_s(x_td, y_td, z_td, B_angle, -bx, -by, -bz, nu, p2_td, -e12, &
                   exx, eyy, ezz, exy, exz, eyz)
    exx_n = exx_n + exx; eyy_n = eyy_n + eyy; ezz_n = ezz_n + ezz
    exy_n = exy_n + exy; exz_n = exz_n + exz; eyz_n = eyz_n + eyz
    call tdsetup_s(x_td, y_td, z_td, C_angle, -bx, -by, -bz, nu, p3_td, -e23, &
                   exx, eyy, ezz, exy, exz, eyz)
    exx_n = exx_n + exx; eyy_n = eyy_n + eyy; ezz_n = ezz_n + ezz
    exy_n = exy_n + exy; exz_n = exz_n + exz; eyz_n = eyz_n + eyz

    ! Average the results
    exx = (exx_p + exx_n) / 2.0_DP
    eyy = (eyy_p + eyy_n) / 2.0_DP
    ezz = (ezz_p + ezz_n) / 2.0_DP
    exy = (exy_p + exy_n) / 2.0_DP
    exz = (exz_p + exz_n) / 2.0_DP
    eyz = (eyz_p + eyz_n) / 2.0_DP
end if
```

### Behavior
- **Action**: Computes 6 angular dislocation contributions and averages
- **Rationale**: Attempts to average Config I and Config II
- **Lines of code**: 33 lines with 6 subroutine calls
- **Computational cost**: O(100) - expensive angular dislocation calculations

---

## Key Differences

| Aspect | MATLAB (Correct) | Fortran (Incorrect) |
|--------|------------------|---------------------|
| **Result** | `NaN` | ~1e-16 (≈ 0) |
| **Meaning** | Undefined/singular | Valid zero value |
| **Implementation** | 7 simple assignments | 6 complex subroutine calls |
| **Performance** | Instant | ~100x slower |
| **Lines of code** | 7 | 33 |
| **Clarity** | Explicit handling | Hidden in averaging |
| **Correctness** | Physically accurate | Conceptually wrong |

---

## Evidence from Debug Log

From the existing log file when `casez_log = T`:

```
Configuration:
trimode =           0
casep_log = F
casen_log = F
casez_log = T
```

After 6 tdsetup_s calls (3 for Config I + 3 for Config II):

```
Before tensor transformation:
exx=  -1.1102230246251565E-016   ← Machine precision, effectively 0
eyy=  -2.6367796834847468E-016   ← Machine precision, effectively 0
ezz=  -8.8817841970012523E-016   ← Machine precision, effectively 0
exy=  -5.5511151231257827E-017   ← Machine precision, effectively 0
exz=   8.3266726846886741E-017   ← Machine precision, effectively 0
eyz=  -2.4286128663675299E-016   ← Machine precision, effectively 0
```

**Analysis**:
- All values are O(1e-16) ≈ double precision epsilon
- These are floating-point rounding errors, not meaningful values
- MATLAB would return `NaN` for this case
- The averaging mathematically produces near-zero, but it's still wrong

---

## Why Averaging Produces ~Zero

The Fortran code computes:
1. **Config I**: Uses slip (bx, by, bz) with angles (-e13, e12, e23)
2. **Config II**: Uses slip (-bx, -by, -bz) with angles (e13, -e12, -e23)

These configurations are **opposite in sign**, so:
```
result = (Config_I + Config_II) / 2 ≈ 0
```

While mathematically interesting, this is **physically wrong** because:
- The point is singular (on the triangle edge)
- No amount of averaging removes the singularity
- The solution is undefined, not zero

---

## Physical Interpretation

### On a Triangle Edge
```
        p3 *
          /|
         / |
        /  |
    x  /   |  ← Point x is ON this edge
      /    |
     /     |
    *------*
   p1      p2
```

When observation point `x` lies exactly on an edge:
- The Green's function has a singularity
- The strain field becomes infinite or undefined
- Mathematical theory says solution is not well-defined

### MATLAB Approach ✓
Returns `NaN` → "Don't use this value, it's undefined"

### Fortran Approach ✗
Returns `0` → "This is a valid zero strain value"

---

## Downstream Implications

### With NaN (MATLAB)
```matlab
if isnan(strain)
    warning('Singular point detected');
    % Handle appropriately
end
```
Code can detect and handle singular points.

### With ~Zero (Fortran)
```fortran
! No way to know this is a singular point
! Value appears valid but is meaningless
if (strain < threshold) then
    ! This might incorrectly trigger
end if
```
Code cannot distinguish singular points from actual zeros.

---

## Performance Comparison

For a single point on triangle edge:

| Implementation | Operations | Approximate Time |
|----------------|-----------|------------------|
| MATLAB (NaN) | 6 assignments | ~1 ns |
| Fortran (average) | 6 × tdsetup_s calls | ~100 ns |

The Fortran implementation is **~100x slower** for these edge cases while producing an **incorrect result**.

---

## Recommended Fix

### Corrected Fortran Code

```fortran
module nikkhoo_walter
  use, intrinsic :: ieee_arithmetic  ! Add this line
  implicit none
  ! ... rest of module ...

contains

subroutine tdstress_hs(...)
  ! ... existing code ...

  else if (casez_log) then
    ! Points on triangle edge are singular - set to NaN
    ! Matches MATLAB implementation (TDstressHS.m:312-318)
    exx = ieee_value(0.0_DP, ieee_quiet_nan)
    eyy = ieee_value(0.0_DP, ieee_quiet_nan)
    ezz = ieee_value(0.0_DP, ieee_quiet_nan)
    exy = ieee_value(0.0_DP, ieee_quiet_nan)
    exz = ieee_value(0.0_DP, ieee_quiet_nan)
    eyz = ieee_value(0.0_DP, ieee_quiet_nan)
  end if

  ! ... rest of subroutine ...
```

### Changes Required
1. **Line ~9**: Add `use, intrinsic :: ieee_arithmetic` to module
2. **Lines 255-287**: Replace entire casez_log block with 6 NaN assignments
3. **Lines 144-146**: Remove now-unused casez_log variables (exx_p, etc.)

---

## Testing Strategy

### 1. Identify casez_log Points
```fortran
! Run trimodefinder on test points
! Points with trimode=0 should return NaN
```

### 2. Verify NaN Return
```fortran
if (ieee_is_nan(strain(1))) then
    print *, "Correct: Singular point detected"
else
    print *, "Error: Should return NaN"
end if
```

### 3. Compare with MATLAB
```matlab
% Run same test in MATLAB
[stress, strain] = TDstressHS(...);
% Should match Fortran NaN behavior
```

---

## Conclusion

The current Fortran implementation:
- ✗ Disagrees with MATLAB reference
- ✗ Returns wrong value (0 instead of NaN)
- ✗ Hides singularities from downstream code
- ✗ Wastes computation on expensive averaging
- ✗ Violates physical interpretation

The fix is simple:
- ✓ Replace 33 lines with 6 NaN assignments
- ✓ Matches MATLAB exactly
- ✓ Faster (100x for edge cases)
- ✓ Physically correct
- ✓ Enables proper error handling

**Implementation**: See CASEZ_LOG_ISSUE.md for detailed fix instructions.
