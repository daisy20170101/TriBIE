#!/bin/bash

echo "========================================================="
echo "Testing CENTER of Triangle"
echo "========================================================="
echo ""

# Clean up
rm -f test_center *.mod

# Compile
echo "Compiling..."
gfortran -O0 -g -fcheck=all -Wall -o test_center sub_nikkhoo.f90 test_center.f90

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi

echo "Compilation successful!"
echo ""

# Run
echo "Running test..."
echo ""
./test_center

echo ""
echo "========================================================="
echo "Test complete"
echo "========================================================="
