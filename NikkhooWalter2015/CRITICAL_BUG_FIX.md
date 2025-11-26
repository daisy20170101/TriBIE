# Critical Bug Fix: Uninitialized Variables in angdis_strain

## The Problem

**Symptom**: NaN values for ALL points, including the triangle center at (-0.333, -0.333, -4.667)

**Root Cause**: Intermediate variables were declared but never calculated

## What Was Wrong

In `sub_nikkhoo.f90`, the `angdis_strain` subroutine had variables declared on line 766:

```fortran
real(DP) :: W, W2, Wr, W2r, Wr3, W2r2
real(DP) :: x2, y2, z2, r2, r, r3, rz, r2z2, r3z
```

But only `W` was calculated (line 785):

```fortran
W = zeta - r
```

The other variables (`W2`, `Wr`, `W2r`, `Wr3`, `W2r2`, `rz`, `r2z2`, `r3z`) were **never assigned values**.

These uninitialized variables were then used in formulas like:

```fortran
C = (r * cosA - z) / Wr    ! Wr was uninitialized!
S = (r * sinA - y) / Wr    ! Wr was uninitialized!
```

And throughout the strain calculations:

```fortran
exx = bx * rFi_rx +
      bx / (8.0_DP * PI * (1.0_DP - nu)) * (eta / Wr + eta * x2 / W2r2 - ...
      ! Using uninitialized Wr, W2r2, etc.
```

**Result**: Using uninitialized variables (which contain garbage values) in calculations produced NaN for every point.

## The Fix

Added the missing calculations (matching MATLAB code, TDstressHS.m lines 614-618):

```fortran
! W calculations
W = zeta - r
W2 = W * W
Wr = W * r
W2r = W2 * r
Wr3 = W * r3
W2r2 = W2 * r2

! Additional calculations
rz = r * z
r2z2 = r2 * z2
r3z = r3 * z
```

## How This Bug Occurred

The bug was introduced when singularity handling code was removed in commit 6a9c7f1:
- **Goal**: Remove ALL singularity handling to match MATLAB exactly
- **What happened**: The essential variable calculations were accidentally deleted along with the singularity checks
- **Why it wasn't caught**: The removal seemed logical at the time ("remove all special handling")

## Testing the Fix

### Test 1: Triangle Center (Should NEVER be singular)

```bash
cd NikkhooWalter2015
chmod +x test_center.sh
./test_center.sh
```

**Expected result**: Finite Exx value for center point (-0.333, -0.333, -4.667)

### Test 2: Multiple Points

```bash
cd NikkhooWalter2015
gfortran -O0 -g -o test_multiple test_multiple_points.f90 sub_nikkhoo.f90
./test_multiple
```

**Expected result**: All 7 test points should return finite values:
1. Center of triangle
2. Inside triangle (offset from center)
3. Far below triangle
4. Far to the side
5. Above triangle
6. Near vertex P1
7. Near edge P1-P2

### Test 3: Original Problematic Points 8 & 9

```bash
cd NikkhooWalter2015
./test_regularization.sh
```

**Expected results**:
- Point 8 (3.0, -3.0, -6.0): Exx = 7.064e-4 (or close)
- Point 9 (-3.0, 3.0, -3.0): Exx = 2.113e-4 (or close)

These may still have singularity issues (they're on extended edge lines), but at least the general case should work now.

## Verification

The fix was verified by:
1. **Comparing with MATLAB code**: The calculations now match TDstressHS.m exactly
2. **Mathematical correctness**: All intermediate variables are now properly defined before use
3. **Physical reasoning**: The center of a triangle cannot be singular

## Impact

**Before fix**:
- ALL points returned NaN (even well-behaved points like the center)
- Implementation was completely broken
- Issue appeared to be related to singularities but was actually basic initialization

**After fix**:
- General points should return finite values
- True singularities (if any) will be isolated to specific edge cases
- Code matches MATLAB implementation structure

## Next Steps

After verifying this fix works:

1. **Test the center**: Confirm finite value for triangle center
2. **Test multiple points**: Verify general case works
3. **Re-test Points 8 & 9**: See if they still have NaN or if it was just this bug
4. **If Points 8 & 9 still have NaN**: Then investigate true singularity handling
5. **Compare with MATLAB**: Run the comparison tests created earlier

## Lessons Learned

When removing code:
- ✗ Don't remove entire blocks without careful analysis
- ✓ Identify which parts are "checks/handling" vs "essential calculations"
- ✓ Verify all declared variables are properly initialized
- ✓ Compare structure with reference implementation (MATLAB)

**Red flag missed**: If a variable is declared but has no assignment statement before use, that's a critical bug.
