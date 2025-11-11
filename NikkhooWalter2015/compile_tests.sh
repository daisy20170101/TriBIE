#!/bin/bash
# Script to recompile sub_nikkhoo module and test programs

echo "==========================================="
echo "Cleaning old compiled files..."
echo "==========================================="
rm -f *.o *.mod test_casez test_casep debug_trimode debug_trimode_detailed

echo ""
echo "==========================================="
echo "Compiling sub_nikkhoo.f90 module..."
echo "==========================================="
gfortran -c -O2 sub_nikkhoo.f90
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile sub_nikkhoo.f90"
    exit 1
fi
echo "SUCCESS: sub_nikkhoo.f90 compiled"

echo ""
echo "==========================================="
echo "Compiling test programs..."
echo "==========================================="

# Compile test_casez
echo "Compiling test_casez..."
gfortran -o test_casez test_casez.f90 sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile test_casez"
    exit 1
fi
echo "SUCCESS: test_casez compiled"

# Compile test_casep
echo "Compiling test_casep..."
gfortran -o test_casep test_casep.f90 sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile test_casep"
    exit 1
fi
echo "SUCCESS: test_casep compiled"

# Compile debug_trimode
echo "Compiling debug_trimode..."
gfortran -o debug_trimode debug_trimode.f90 sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile debug_trimode"
    exit 1
fi
echo "SUCCESS: debug_trimode compiled"

echo ""
echo "==========================================="
echo "All programs compiled successfully!"
echo "==========================================="
echo ""
echo "Run tests with:"
echo "  ./test_casez"
echo "  ./test_casep"
echo "  ./debug_trimode"
echo ""
