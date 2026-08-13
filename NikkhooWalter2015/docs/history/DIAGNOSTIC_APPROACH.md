# Diagnostic Approach for Persistent NaN Issue

## Problem Summary

Despite the bounds-checking fix being verified in `sub_nikkhoo.f90` source code:
- **Standalone program** (`debug_trimode_detailed`) works correctly → returns valid numbers for points 8 & 9
- **Module-based program** (`test_trimode_module`) still returns NaN → compiled module not reflecting source code

## Root Cause Analysis

This pattern indicates one of these issues:

1. **Stale compiled module files** - Old `.mod` or `.o` files are being used
2. **Multiple module files** - Compiler finding old module in different directory
3. **Compiler caching** - Fortran compiler not actually recompiling the module
4. **Module search path** - `nikkhoo_walter.mod` from unexpected location being used

## Diagnostic Solution

I've added **debug print statements** directly inside the `trimode_finder` function in `sub_nikkhoo.f90`.

### What Was Added

```fortran
! Lines 647-649: Print barycentric coordinates
print *, '[DEBUG trimode_finder] Input: x=', x, ' y=', y, ' z=', z
print *, '[DEBUG trimode_finder] Barycentric: a=', a, ' b=', b, ' c=', c

! Line 666: Print tolerance value
print *, '[DEBUG trimode_finder] Checking bounds with BARY_TOL=', BARY_TOL

! Lines 669, 673, 677: Print which edge case triggers
print *, '[DEBUG trimode_finder] Edge case A/B/C: ...'

! Line 684: Print z!=0 override
print *, '[DEBUG trimode_finder] z!=0 override: trimode 0->1'

! Line 688: Print final result
print *, '[DEBUG trimode_finder] FINAL trimode=', trimode
```

### Why This Works

When you compile the module and run the test:

**If you see debug output:**
- ✓ The NEW code with bounds checking IS being executed
- ✓ The module was successfully recompiled
- We can analyze the barycentric coordinates to see why NaN persists

**If you see NO debug output:**
- ✗ OLD code without debug statements is being executed
- ✗ The compiler is using cached/stale module files
- Need to investigate module search paths and compiler behavior

## Running the Diagnostic

### Option 1: Use the Automated Script

```bash
cd /home/user/TriBIE/NikkhooWalter2015
./diagnose_module.sh
```

This script will:
1. Verify the fix and debug statements are in source code
2. Delete ALL old compiled files (*.o, *.mod, executables)
3. Compile `sub_nikkhoo.f90` module from scratch
4. Verify module files were created with timestamps and checksums
5. Compile `test_trimode_module.f90` linked against the new module
6. Run the test and capture all debug output

### Option 2: Manual Steps

```bash
cd NikkhooWalter2015

# Clean everything
rm -f *.o *.mod test_* debug_*

# Compile module
gfortran -c -O2 sub_nikkhoo.f90

# Verify files created
ls -lh nikkhoo_walter.mod sub_nikkhoo.o

# Compile test
gfortran -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o

# Run test
./test_trimode_module
```

## Expected Output

### If Module Is Properly Recompiled

You should see extensive debug output like:

```
============================================================
Testing module version (using nikkhoo_walter module)
============================================================

This test uses the ACTUAL compiled module, not standalone code.
If this shows NaN, the module was not properly recompiled.

========== POINT 8 ==========
Coords: x=   3.0000000000000000      , y=  -3.0000000000000000      , z=  -6.0000000000000000
 [DEBUG trimode_finder] Input: x=  -4.5825756949558398       y=   5.0000000000000000       z=   3.5355339059327378
 [DEBUG trimode_finder] Barycentric: a=  -1.8284271247461898       b=   2.9142135623730949       c= -8.5786437934211137E-002
 [DEBUG trimode_finder] Checking bounds with BARY_TOL=   1.0000000000000000E-012
 [DEBUG trimode_finder] FINAL trimode=          -1

Strain(1) Exx =   -some valid number (not NaN)
  ✓ PASS: Got valid number (module has the fix)
```

### If Old Module Is Still Being Used

You'll see:

```
========== POINT 8 ==========
Coords: x=   3.0000000000000000      , y=  -3.0000000000000000      , z=  -6.0000000000000000
Strain(1) Exx =                        NaN
  ✗ FAIL: Got NaN (module still has old code without bounds checking)
```

**No debug output at all** - this means the old module without debug statements is being executed.

## Next Steps Based on Results

### Case 1: Debug Output Shows, But Still NaN

If you see the debug output with `BARY_TOL=1e-12` and bounds checking, but still get NaN:
- The trimode classification logic itself may need adjustment
- We can analyze the actual barycentric coordinates from the debug output
- May need to investigate the casez_log handling or other parts of the code

### Case 2: No Debug Output (Old Module Being Used)

If you see NO debug output at all:
1. Check for multiple module files: `find .. -name "nikkhoo_walter.mod"`
2. Check Fortran module search path: `gfortran -v`
3. Try absolute path compilation:
   ```bash
   cd NikkhooWalter2015
   gfortran -o test_trimode_module test_trimode_module.f90 $(pwd)/sub_nikkhoo.o
   ```
4. Check if old executables are in your PATH

## Files Modified

1. **sub_nikkhoo.f90**: Added debug print statements (lines 647-689)
2. **diagnose_module.sh**: Automated diagnostic script (NEW)
3. **DIAGNOSTIC_APPROACH.md**: This documentation (NEW)

## Rollback Instructions

To remove debug output after diagnosis:

```bash
cd NikkhooWalter2015
git checkout sub_nikkhoo.f90  # Restore version without debug prints
```

Or manually remove the `print *, '[DEBUG trimode_finder]'` lines.

## Summary

This diagnostic approach will **definitively show** whether:
- The new bounds-checking code is being executed
- The Fortran compiler is properly recompiling the module
- Where the problem actually lies (compilation vs. logic)

Run `./diagnose_module.sh` and examine the output carefully.
