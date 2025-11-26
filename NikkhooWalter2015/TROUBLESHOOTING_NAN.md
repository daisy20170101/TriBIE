# TROUBLESHOOTING: Points 8 and 9 Still Return NaN

## The Problem

You've confirmed the fix is in the source code (`verify_fix.sh` passed), you've run `compile_tests.sh`, but test programs still show NaN for points 8 and 9.

## Why This Happens

The most common cause is **stale compiled files** being used from somewhere else, or **compilation happening in the wrong directory**.

## Step-by-Step Solution

### Step 1: Verify You're in the Right Directory

```bash
cd /path/to/TriBIE/NikkhooWalter2015
pwd  # Should show .../TriBIE/NikkhooWalter2015
ls sub_nikkhoo.f90  # Should exist
```

### Step 2: Verify the Fix is in the Source Code

```bash
grep -c "b <= 1.0_DP + BARY_TOL" sub_nikkhoo.f90
```

**Expected output**: `2` or more

If you see `0`, you need to:
```bash
git pull origin claude/open-tribi-011CV1WxGz1Q8daU8f4udfNA
```

### Step 3: Complete Clean and Recompile

```bash
# Nuclear option - remove ALL compiled files
rm -f *.o *.mod test_* debug_* *.out

# Verify they're gone
ls *.o *.mod test_* debug_* 2>/dev/null  # Should show "No such file"

# Recompile
./compile_tests.sh
```

### Step 4: Run the Simple Module Test

```bash
./test_trimode_module
```

**Expected output**:
```
========== POINT 8 ==========
Strain(1) Exx =  <some number, not NaN>
  ✓ PASS: Got valid number (module has the fix)

========== POINT 9 ==========
Strain(1) Exx =  <some number, not NaN>
  ✓ PASS: Got valid number (module has the fix)
```

**If you still see NaN**, continue to Step 5.

### Step 5: Manual Compilation with Verbose Output

```bash
# Compile module with verbose output
gfortran -c -v sub_nikkhoo.f90 2>&1 | tee compile.log

# Check if module was created
ls -l nikkhoo_walter.mod sub_nikkhoo.o

# Compile simple test
gfortran -v -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o 2>&1 | tee link.log

# Run it
./test_trimode_module
```

Look at `compile.log` and `link.log` for any warnings or errors.

### Step 6: Check for Multiple Versions

Sometimes old compiled files exist in parent directories or system paths:

```bash
# Search for old module files
find .. -name "nikkhoo_walter.mod" 2>/dev/null
find .. -name "sub_nikkhoo.o" 2>/dev/null

# Check system module path (if any)
echo $GFORTRAN_MODULE_PATH
```

If you find old versions, delete them:
```bash
find .. -name "nikkhoo_walter.mod" -delete
find .. -name "sub_nikkhoo.o" -delete
```

### Step 7: Force Module Recompilation with Different Flag

Try compiling without optimization:

```bash
rm -f *.o *.mod
gfortran -c sub_nikkhoo.f90  # No -O2 flag
gfortran -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o
./test_trimode_module
```

### Step 8: Check GFortran Version

Some very old gfortran versions might have issues:

```bash
gfortran --version
```

**Minimum recommended**: gfortran 4.8 or later

### Step 9: Verify Module Contents

Check if the compiled module actually has the fix:

```bash
# Dump module interface (not all gfortran versions support this)
gfortran -fdump-fortran-original sub_nikkhoo.f90 2>&1 | grep -A 5 "trimode_finder"

# Or check object file symbols
nm sub_nikkhoo.o | grep trimode
```

## Common Issues and Solutions

### Issue 1: "Permission Denied" or "Text File Busy"

```bash
# Kill any running processes using the files
fuser -k test_casep test_casez

# Remove files
rm -f test_casep test_casez
```

### Issue 2: Working Directory Problems

Make sure you're running the executables from the SAME directory where you compiled:

```bash
# BAD - might run old version from PATH
test_casep

# GOOD - runs local version
./test_casep
```

### Issue 3: Module File (.mod) Incompatibility

If you previously compiled with a different compiler version:

```bash
# Remove all .mod files
rm -f *.mod

# Recompile EVERYTHING from scratch
./compile_tests.sh
```

### Issue 4: Cache Issues with Make/CMake

If you're using a build system:

```bash
make clean
rm -rf build/
./compile_tests.sh  # Use our script instead
```

## Ultimate Test: Inline Compilation

If nothing else works, compile everything in one command:

```bash
gfortran -o test_inline sub_nikkhoo.f90 test_trimode_module.f90
./test_inline
```

This bypasses any caching or module path issues.

## What Each Test Should Show

### test_trimode_module (simplest)
- Point 8: Valid number (not NaN)
- Point 9: Valid number (not NaN)

### debug_trimode_detailed (standalone, always works)
- Point 8: trimode = -1, valid numbers
- Point 9: trimode = 1, valid numbers

### test_casez and test_casep (should work after fix)
- Point 8: Valid strain value
- Point 9: Valid strain value

## If Still Getting NaN

If you've tried everything above and STILL get NaN, send me:

1. Output of `./verify_fix.sh`
2. Output of `./compile_tests.sh`
3. Output of `./test_trimode_module`
4. Output of `ls -lt *.f90 *.o *.mod test_* debug_*`
5. Output of `gfortran --version`

This will help diagnose any unusual environment issues.

## Quick Reference Commands

```bash
# Complete reset and recompile
cd NikkhooWalter2015
rm -f *.o *.mod test_* debug_*
./compile_tests.sh
./test_trimode_module

# If that passes, run full tests
./test_casep
```

## Success Criteria

✓ `test_trimode_module` shows valid numbers (not NaN) for points 8 and 9
✓ `test_casep` shows complete stress/strain tensors without NaN
✓ `debug_trimode_detailed` confirms trimode = -1 or 1 (not 0)

The fix IS in the code - it's a matter of ensuring your compiled binaries are using it!
