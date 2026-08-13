#!/bin/bash
# Test script to verify the angular dislocation singularity fix

echo "========================================================"
echo "Testing Angular Dislocation Singularity Fix"
echo "========================================================"
echo ""

# Clean and compile
echo "Step 1: Cleaning old files..."
rm -f *.o *.mod test_point8_only

echo ""
echo "Step 2: Compiling sub_nikkhoo.f90..."
gfortran -c -O2 sub_nikkhoo.f90
if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi
echo "SUCCESS: Module compiled"

echo ""
echo "Step 3: Compiling test_point8_only.f90..."
gfortran -o test_point8_only test_point8_only.f90 sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Test compilation failed!"
    exit 1
fi
echo "SUCCESS: Test program compiled"

echo ""
echo "Step 4: Running test..."
echo "========================================================"
./test_point8_only 2>&1 | grep -A2 "RESULTS for Point 8"
echo "========================================================"
echo ""

echo "Expected output:"
echo "  ✓ PASS: Got valid number"
echo ""
echo "If you see NaN, the fix didn't work."
echo "If you see a valid number, the fix is successful!"
echo ""
