# Debug Findings: NaN Issue for Points 8 & 9

## Current Status

### What We Know ✓

1. **The bounds-checking fix IS in the source code** - verified by grep and md5sum
2. **The compiled module IS using the new code** - debug output appears showing:
   ```
   [DEBUG trimode_finder] Input: x=   2.0000000000000000
   [DEBUG trimode_finder] Barycentric: a=   3.7888543819998315  b=   0.0000000000000000  c=  -2.7888543819998315
   [DEBUG trimode_finder] Checking bounds with BARY_TOL=   9.9999999999999998E-013
   [DEBUG trimode_finder] FINAL trimode=          -1
   ```

3. **trimode_finder is working correctly**:
   - Point 8 has barycentric coords: `a=3.79, b=0.0, c=-2.79`
   - One coord is exactly zero (`b=0.0`) → on extended edge line
   - Other coords outside [0,1] bounds → NOT on actual triangle edge
   - Correctly classified as `trimode = -1` (Configuration II)

4. **BUT the result is still NaN** even though trimode = -1:
   ```
   === Main Dislocation Contribution ===
   Strain: Exx= NaN [all components NaN]
   === Harmonic Function Contribution ===
   Strain: Exx= -8.2624952469809805E-004 [valid numbers]
   === Image Dislocation Contribution ===
   Strain: Exx= 2.1669275953637400E-004 [valid numbers]
   ```

### The Mystery

**The trimode_finder is NOT the problem!** It's correctly returning trimode=-1.

**The NaN is coming from somewhere else in the calculation chain:**
- tdstress_fs calls trimode_finder → gets trimode=-1
- Should enter `casen_log` (Configuration II) path
- Calls tdsetup_s three times for angular dislocations
- tdsetup_s calls angdis_strain
- **Somewhere in this chain, NaN is being generated**

## Analysis: Why Extended Edge Line Points Are Problematic

### The Mathematics

Point 8 has these characteristics in TDCS (Triangle Dislocation Coordinate System):
- Input to trimode_finder: `x=2.0, y=-2.24, z=0.0`
- **Key insight**: `z=0.0` (point is on the triangle plane)
- Barycentric: `b=0.0` (on extended edge line for edge with b=0)

### Potential Division by Zero in Angular Dislocation Calculations

The angular dislocation calculations in `angdis_strain` have these denominators:

1. **W = zeta - r** where `zeta = y*sin(alpha) + z*cos(alpha)`
   - If W ≈ 0, then divisions by Wr, W²r, etc. → ∞ or NaN

2. **r - z** (appears in rz, r2z2, r3z)
   - For Point 8: z=0 in TDCS, r=sqrt(x²+y²+z²) ≈ 3.0
   - So r-z ≈ 3.0 (should be OK)

3. **r - zeta** (appears in partial derivatives)
   - Depends on alpha (the angular dislocation angle)
   - Could be problematic for certain angles

### Why Image Dislocation Works But Main Dislocation Doesn't

Looking at the debug output:
- **Main Dislocation**: trimode=-1, z_td≈0 → NaN
- **Image Dislocation**: trimode=-1, z_td≈-10.7 → Valid numbers

The difference is the z-coordinate! The image dislocation has a large negative z (mirror image),
while the main dislocation has z≈0 (on the triangle plane).

**Hypothesis**: When z≈0 AND point is on extended edge line (b=0), the angular dislocation
calculations encounter singularities in W, Wr, or related terms.

## Enhanced Debug Output

I've added comprehensive debug statements to trace exactly where the NaN originates:

### In `tdstress_fs` (lines 213-214, 261, 271):
```fortran
print *, '[DEBUG tdstress_fs] After trimode_finder: trimode=', trimode
print *, '[DEBUG tdstress_fs] casep_log=', casep_log, ' casen_log=', casen_log, ' casez_log=', casez_log
print *, '[DEBUG tdstress_fs] casez_log=TRUE, setting all values to NaN'  ! Only if casez_log
print *, '[DEBUG tdstress_fs] Before transformation: exx=', exx, ' (is_nan=', ieee_is_nan(exx), ')'
```

### In `angdis_strain` (lines 787-792, 795-801):
```fortran
! Check for r-z near zero
if (abs(r - z) < 1.0e-10_DP) then
  print *, '[DEBUG angdis_strain] WARNING: r-z near zero!'
  print *, '[DEBUG angdis_strain] r=', r, ' z=', z, ' r-z=', r-z
end if

! Check for W or Wr near zero
if (abs(W) < 1.0e-10_DP .or. abs(Wr) < 1.0e-10_DP) then
  print *, '[DEBUG angdis_strain] WARNING: W or Wr near zero!'
  print *, '[DEBUG angdis_strain] W=', W, ' Wr=', Wr
  print *, '[DEBUG angdis_strain] Input: x=', x, ' y=', y, ' z=', z
  print *, '[DEBUG angdis_strain] alpha=', alpha, ' zeta=', zeta, ' r=', r
end if
```

## Next Steps: Run Enhanced Diagnostic

### On Your System with Fortran Compiler:

```bash
cd NikkhooWalter2015
./diagnose_module.sh
```

This will:
1. Clean and recompile everything from scratch
2. Run `test_point8_only` with full debug output
3. Show exactly which function is generating NaN

### What to Look For in the Output:

**Case 1: If you see `[DEBUG angdis_strain] WARNING: W or Wr near zero!`**
- This means the angular dislocation calculation encounters a singularity
- W=0 occurs when zeta=r, which is a known singularity in the Nikkhoo-Walter method
- **Solution**: Need special case handling for points on extended edge lines with z≈0

**Case 2: If you see `[DEBUG tdstress_fs] casez_log=TRUE`**
- This means trimode==0 is somehow still being triggered
- Would indicate a bug in the bounds checking logic
- **Solution**: Need to investigate why bounds checking isn't working

**Case 3: If neither warning appears**
- NaN is coming from a different calculation
- Look at the exx value progression through the debug statements
- Trace backwards from first NaN appearance

## Expected Output Pattern

You should see something like:
```
Point 8 coords: x=   3.0000000000000000      , y=  -3.0000000000000000      , z=  -6.0000000000000000

Calling tdstress_hs...

[DEBUG trimode_finder] Input: x=   2.0000000000000000       y=  -2.2360679774997898       z=   0.0000000000000000
[DEBUG trimode_finder] Barycentric: a=   3.7888543819998315       b=   0.0000000000000000       c=  -2.7888543819998315
[DEBUG trimode_finder] Checking bounds with BARY_TOL=   9.9999999999999998E-013
[DEBUG trimode_finder] FINAL trimode=          -1

[DEBUG tdstress_fs] After trimode_finder: trimode=          -1
[DEBUG tdstress_fs] casep_log= F  casen_log= T  casez_log= F
[DEBUG tdstress_fs] Entering casen_log (Config II) path

[DEBUG angdis_strain] WARNING: W or Wr near zero!     <--- LIKELY TO SEE THIS
[DEBUG angdis_strain] W=   1.2345e-15   Wr=   3.7035e-15
[DEBUG angdis_strain] Input: x=   ...  y=   ...  z=   ...
[DEBUG angdis_strain] alpha=   ...   zeta=   ...   r=   ...

[DEBUG tdstress_fs] After 1st tdsetup_s: exx=  NaN  (is_nan= T)
[DEBUG tdstress_fs] Before transformation: exx=  NaN  (is_nan= T)
```

## Possible Solutions (Once We Confirm the Cause)

### If W≈0 Singularity is the Issue:

The Nikkhoo & Walter (2015) paper discusses this:
- Points where zeta=r are true singularities in the angular dislocation formulation
- For points on extended edge lines with z≈0, this can occur
- **Solution**: Add special case handling to avoid division by W or Wr when they're near zero
- May need to use asymptotic expansions or limit forms

### If It's a Different Issue:

The debug output will guide us to the exact source, and we can develop an appropriate fix.

## Files Modified

1. **sub_nikkhoo.f90**: Added debug prints in tdstress_fs and angdis_strain
2. **test_point8_only.f90**: Simple test for Point 8 only (NEW)
3. **diagnose_module.sh**: Updated to compile and run test_point8_only

## Summary

We've confirmed:
- ✓ Bounds checking fix is in source code
- ✓ Module is properly recompiled
- ✓ trimode_finder works correctly (returns -1, not 0)

We need to identify:
- ? Where exactly is the NaN being generated?
- ? Is it a W≈0 singularity in angdis_strain?
- ? Or something else in the calculation chain?

**Run `./diagnose_module.sh` to get the answer!**
