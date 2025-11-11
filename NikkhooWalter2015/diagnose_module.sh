#!/bin/bash
# Diagnostic script to verify module compilation and debug NaN issue

echo "========================================================"
echo "DIAGNOSTIC: Module Compilation and NaN Issue"
echo "========================================================"
echo ""

# Step 1: Verify the fix is in source code
echo "Step 1: Checking source code for bounds checking fix..."
BOUNDS_CHECK=$(grep -c "b <= 1.0_DP + BARY_TOL" sub_nikkhoo.f90)
if [ $BOUNDS_CHECK -ge 2 ]; then
    echo "  ✓ Bounds checking found in source ($BOUNDS_CHECK occurrences)"
else
    echo "  ✗ ERROR: Bounds checking NOT in source!"
    exit 1
fi

DEBUG_CHECK=$(grep -c "DEBUG trimode_finder" sub_nikkhoo.f90)
if [ $DEBUG_CHECK -ge 5 ]; then
    echo "  ✓ Debug statements found in source ($DEBUG_CHECK occurrences)"
else
    echo "  ✗ ERROR: Debug statements NOT in source!"
    exit 1
fi

# Step 2: Clean ALL old compiled files
echo ""
echo "Step 2: Removing ALL old compiled files..."
echo "  Searching for .mod files in current directory..."
find . -maxdepth 1 -name "*.mod" -exec ls -lh {} \;
echo "  Searching for .o files in current directory..."
find . -maxdepth 1 -name "*.o" -exec ls -lh {} \;

echo "  Deleting all compiled files..."
rm -fv *.o *.mod test_casez test_casep debug_trimode debug_trimode_detailed test_trimode_module

# Step 3: Compile module with verbose output
echo ""
echo "Step 3: Compiling sub_nikkhoo.f90 module..."
echo "Command: gfortran -c -O2 sub_nikkhoo.f90"
gfortran -c -O2 sub_nikkhoo.f90
COMPILE_STATUS=$?

if [ $COMPILE_STATUS -ne 0 ]; then
    echo "  ✗ ERROR: Module compilation failed!"
    exit 1
fi

# Step 4: Verify module files were created
echo ""
echo "Step 4: Verifying module files..."
if [ -f "nikkhoo_walter.mod" ]; then
    echo "  ✓ nikkhoo_walter.mod created ($(stat -c%s nikkhoo_walter.mod) bytes, $(date -r nikkhoo_walter.mod '+%Y-%m-%d %H:%M:%S'))"
    md5sum nikkhoo_walter.mod
else
    echo "  ✗ ERROR: nikkhoo_walter.mod NOT created!"
    exit 1
fi

if [ -f "sub_nikkhoo.o" ]; then
    echo "  ✓ sub_nikkhoo.o created ($(stat -c%s sub_nikkhoo.o) bytes, $(date -r sub_nikkhoo.o '+%Y-%m-%d %H:%M:%S'))"
    md5sum sub_nikkhoo.o
else
    echo "  ✗ ERROR: sub_nikkhoo.o NOT created!"
    exit 1
fi

# Step 5: Compile test program
echo ""
echo "Step 5: Compiling test_trimode_module.f90..."
echo "Command: gfortran -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o"
gfortran -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o
COMPILE_STATUS=$?

if [ $COMPILE_STATUS -ne 0 ]; then
    echo "  ✗ ERROR: Test program compilation failed!"
    exit 1
fi

if [ -f "test_trimode_module" ]; then
    echo "  ✓ test_trimode_module created ($(stat -c%s test_trimode_module) bytes, $(date -r test_trimode_module '+%Y-%m-%d %H:%M:%S'))"
else
    echo "  ✗ ERROR: test_trimode_module NOT created!"
    exit 1
fi

# Step 6: Run test and capture output
echo ""
echo "========================================================"
echo "Step 6: Running test_trimode_module with DEBUG output"
echo "========================================================"
echo ""
echo "If the debug output shows bounds checking (BARY_TOL), then"
echo "the NEW code is being used. If not, old code is cached."
echo ""
echo "Looking for these key indicators:"
echo "  - '[DEBUG trimode_finder]' messages"
echo "  - 'BARY_TOL=' showing tolerance value"
echo "  - Barycentric coordinates for points 8 and 9"
echo ""
echo "-------- OUTPUT BEGINS --------"
./test_trimode_module
echo "-------- OUTPUT ENDS --------"
echo ""

echo "========================================================"
echo "DIAGNOSTIC COMPLETE"
echo "========================================================"
echo ""
echo "If you see debug output above, the new code IS running."
echo "If you see NO debug output, then:"
echo "  1. Old module file from different location is being used"
echo "  2. Compiler is caching old version"
echo "  3. Module search path issue"
echo ""
