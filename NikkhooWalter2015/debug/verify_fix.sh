#!/bin/bash
# Verification script to check if the bounds checking fix is in place

echo "========================================================"
echo "Checking for bounds checking fix in sub_nikkhoo.f90..."
echo "========================================================"
echo ""

FILE="sub_nikkhoo.f90"

if [ ! -f "$FILE" ]; then
    echo "ERROR: $FILE not found!"
    echo "Make sure you're in the NikkhooWalter2015 directory"
    exit 1
fi

# Check for the bounds checking fix
BOUNDS_CHECK=$(grep -c "b <= 1.0_DP + BARY_TOL" $FILE)

if [ $BOUNDS_CHECK -ge 2 ]; then
    echo "✓ GOOD: Bounds checking fix found in source code"
    echo "  Found $BOUNDS_CHECK occurrences of bounds check"
else
    echo "✗ BAD: Bounds checking fix NOT found in source code"
    echo "  You need to pull the latest code from git"
    exit 1
fi

echo ""
echo "Source code is correct. Now checking for old compiled files..."
echo ""

# Check for old compiled files
OLD_FILES=""
if [ -f "nikkhoo_walter.mod" ]; then
    OLD_FILES="$OLD_FILES nikkhoo_walter.mod"
fi
if [ -f "sub_nikkhoo.o" ]; then
    OLD_FILES="$OLD_FILES sub_nikkhoo.o"
fi
if [ -f "test_casez" ]; then
    OLD_FILES="$OLD_FILES test_casez"
fi
if [ -f "test_casep" ]; then
    OLD_FILES="$OLD_FILES test_casep"
fi

if [ -n "$OLD_FILES" ]; then
    echo "⚠ WARNING: Old compiled files found:"
    for f in $OLD_FILES; do
        echo "    - $f ($(date -r $f '+%Y-%m-%d %H:%M:%S'))"
    done
    echo ""
    echo "These files need to be deleted and recompiled!"
    echo "Run: rm -f *.o *.mod test_casez test_casep"
else
    echo "✓ GOOD: No old compiled files found"
fi

echo ""
echo "========================================================"
echo "Summary:"
echo "========================================================"
echo ""
if [ -n "$OLD_FILES" ]; then
    echo "ACTION REQUIRED:"
    echo "1. Delete old files: rm -f *.o *.mod test_casez test_casep"
    echo "2. Recompile: ./compile_tests.sh"
    echo ""
    echo "The fix is in the source code, but you're running old"
    echo "compiled binaries that don't have the fix!"
else
    echo "Source code has the fix."
    echo "If test programs don't exist, compile with:"
    echo "  ./compile_tests.sh"
    echo ""
    echo "If test programs still show NaN for points 8 & 9,"
    echo "then you may be running old binaries from a different"
    echo "directory. Make sure you're running the newly compiled"
    echo "versions in THIS directory!"
fi
echo ""
