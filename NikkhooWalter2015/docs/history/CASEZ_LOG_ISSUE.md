# casez_log Implementation Issue

## Summary

The Fortran implementation of `casez_log` (trimode == 0) differs from the MATLAB reference implementation in TDstressHS.m. This document details the discrepancy and provides the correct fix.

## Background

When a calculation point lies exactly on the edge of a triangular dislocation, it represents a **singularity** where the analytical solution is undefined. The trimode classifier returns `0` for such points.

## Current Implementations

### MATLAB Reference (TDstressHS.m:312-318) ✓ CORRECT

```matlab
if nnz(casezLog)~=0
    exx(casezLog,1) = nan;
    eyy(casezLog,1) = nan;
    ezz(casezLog,1) = nan;
    exy(casezLog,1) = nan;
    exz(casezLog,1) = nan;
    eyz(casezLog,1) = nan;
end
```

**Behavior**: Returns `NaN` (Not-a-Number) for all strain components when point is on triangle edge.

**Rationale**: The solution is undefined/singular at these locations.

### Fortran Implementation (sub_nikkhoo.f90:255-287) ✗ INCORRECT

```fortran
else if (casez_log) then
    ! For points on the triangle, use average of positive and negative cases
    ! Configuration I (positive)
    call tdsetup_s(x_td, y_td, z_td, A_angle, bx, by, bz, nu, p1_td, -e13, &
                   exx_p, eyy_p, ezz_p, exy_p, exz_p, eyz_p)
    call tdsetup_s(x_td, y_td, z_td, B_angle, bx, by, bz, nu, p2_td, e12, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)
    ! ... [sums contributions]

    ! Configuration II (negative)
    call tdsetup_s(x_td, y_td, z_td, A_angle, -bx, -by, -bz, nu, p1_td, e13, &
                   exx_n, eyy_n, ezz_n, exy_n, exz_n, eyz_n)
    ! ... [sums contributions]

    ! Average the results
    exx = (exx_p + exx_n) / 2.0_DP
    eyy = (eyy_p + eyy_n) / 2.0_DP
    ezz = (ezz_p + ezz_n) / 2.0_DP
    exy = (exy_p + exy_n) / 2.0_DP
    exz = (exz_p + exz_n) / 2.0_DP
    eyz = (eyz_p + eyz_n) / 2.0_DP
end if
```

**Behavior**: Computes average of Configuration I and Configuration II.

**Result**: Returns values near machine precision (~1e-16), effectively zero.

## Evidence from Debug Log

From `/home/user/TriBIE/NikkhooWalter2015/log` lines 303-421:

```
Configuration:
trimode =           0
casep_log = F
casen_log = F
casez_log = T
```

After computing 6 angular dislocation contributions and averaging:

```
Before tensor transformation:
exx=  -1.1102230246251565E-016
eyy=  -2.6367796834847468E-016
ezz=  -8.8817841970012523E-016
exy=  -5.5511151231257827E-017
exz=   8.3266726846886741E-017
eyz=  -2.4286128663675299E-016
```

All values are O(1e-16) ≈ machine precision, effectively zero.

## Why This Matters

### Semantic Difference

- **NaN**: "This point is singular/undefined - don't use this value"
- **Zero**: "This is a valid result with zero strain"

### Implications

1. **Numerical**: The Fortran code gives ~0 instead of NaN
2. **Physical**: Point on edge is singular, solution is undefined
3. **Downstream**: Code using these results might not know the value is invalid
4. **Scientific**: Violates the theoretical foundation of the method

### When Does This Occur?

Points exactly on triangle edges are rare in practice but can occur when:
- Using structured grids aligned with triangle boundaries
- Placing observation points at triangle vertices
- Using coarse meshes where elements share edges
- Testing/validation with synthetic geometries

## Test Case

Using the triangle:
- p1 = [-1.0, -1.0, -5.0]
- p2 = [1.0, -1.0, -5.0]
- p3 = [-1.0, 1.0, -4.0]

Test points that may trigger casez_log:
```fortran
! Points 10-12 are likely on edges (z = -1.0 is above triangle)
x(10) = -1.0, y(10) = -1.0, z(10) = -1.0  ! Possibly on edge p1-p3 projection
x(11) = -1.0, y(11) =  1.0, z(11) = -1.0  ! Near vertex p3
x(12) =  1.0, y(12) = -1.0, z(12) = -1.0  ! Near vertex p2
```

## Recommended Fix

### Option 1: Set to NaN (Matches MATLAB)

```fortran
else if (casez_log) then
    ! Points on triangle edge are singular - set to NaN
    exx = ieee_value(0.0_DP, ieee_quiet_nan)
    eyy = ieee_value(0.0_DP, ieee_quiet_nan)
    ezz = ieee_value(0.0_DP, ieee_quiet_nan)
    exy = ieee_value(0.0_DP, ieee_quiet_nan)
    exz = ieee_value(0.0_DP, ieee_quiet_nan)
    eyz = ieee_value(0.0_DP, ieee_quiet_nan)
end if
```

**Note**: Requires `use, intrinsic :: ieee_arithmetic` at module level.

### Option 2: Set to Large Value (If NaN causes issues)

```fortran
else if (casez_log) then
    ! Points on triangle edge are singular - set to large value
    real(DP), parameter :: SINGULAR_VALUE = huge(1.0_DP)
    exx = SINGULAR_VALUE
    eyy = SINGULAR_VALUE
    ezz = SINGULAR_VALUE
    exy = SINGULAR_VALUE
    exz = SINGULAR_VALUE
    eyz = SINGULAR_VALUE
end if
```

### Option 3: Set to Zero with Warning (Least recommended)

```fortran
else if (casez_log) then
    ! Points on triangle edge are singular - set to zero
    ! WARNING: This is a singular point, value is technically undefined
    exx = 0.0_DP
    eyy = 0.0_DP
    ezz = 0.0_DP
    exy = 0.0_DP
    exz = 0.0_DP
    eyz = 0.0_DP
end if
```

## Recommendation

**Use Option 1 (NaN)** to match the MATLAB reference implementation exactly. This is the most scientifically correct approach and clearly identifies singular points.

## Files to Modify

1. **Primary**: `/home/user/TriBIE/NikkhooWalter2015/sub_nikkhoo.f90`
   - Lines 255-287 (casez_log block in tdstress_hs subroutine)
   - Add `use, intrinsic :: ieee_arithmetic` to module

2. **Documentation**: Update README_NIKKHOO.md to note singular point handling

## Testing

After fix, verify:
1. Run debug_nikkhoo with test points
2. Compare with MATLAB TDstressHS.m results
3. Confirm NaN propagation in downstream calculations
4. Check that NaN handling doesn't cause crashes

## References

- Nikkhoo & Walter (2015), GJI: Original paper
- TDstressHS.m: MATLAB reference implementation (lines 312-318)
- sub_nikkhoo.f90: Fortran implementation (lines 255-287)
- Debug log: Evidence of current behavior
