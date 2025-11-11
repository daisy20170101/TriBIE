#!/bin/bash
# Test regularization approach for angular dislocation singularities

echo "========================================================"
echo "Testing Regularization Approach"
echo "========================================================"
echo ""

# Clean
echo "Cleaning old files..."
rm -f *.o *.mod quick_test

# Compile
echo "Compiling sub_nikkhoo.f90..."
gfortran -c -O2 sub_nikkhoo.f90
if [ $? -ne 0 ]; then
    echo "ERROR: Module compilation failed!"
    exit 1
fi

echo "Compiling quick_test.f90..."
gfortran -o quick_test quick_test.f90 sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Test compilation failed!"
    exit 1
fi

echo ""
echo "========================================================"
echo "Running test..."
echo "========================================================"
./quick_test

echo ""
echo "========================================================"
echo "Expected values:"
echo "  Point 8: Exx = 7.064e-4"
echo "  Point 9: Exx = 2.113e-4"
echo "========================================================"
