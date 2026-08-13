#!/bin/bash

echo "========================================================="
echo "Testing Problematic Points 4, 5, 12, 15"
echo "========================================================="
echo ""

# Clean up
rm -f test_problem_points *.mod

# Compile
echo "Compiling..."
gfortran -O0 -g -fcheck=all -Wall -o test_problem_points sub_nikkhoo.f90 test_problem_points.f90

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi

echo "Compilation successful!"
echo ""

# Run
echo "Running test..."
echo ""
./test_problem_points

echo ""
echo "========================================================="
echo ""
echo "Also running full test_casep to see all 15 points..."
echo ""
gfortran -O0 -g -o test_casep sub_nikkhoo.f90 test_casep.f90 2>/dev/null
if [ $? -eq 0 ]; then
    ./test_casep 2>&1 | grep -A1 "Point\|Exx"
fi

echo ""
echo "========================================================="
echo "Test complete"
echo "========================================================="
