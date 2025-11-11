#!/bin/bash
# Script to recompile sub_nikkhoo module and test programs

echo "==========================================="
echo "Cleaning old compiled files..."
echo "==========================================="
rm -f *.o *.mod test_casez test_casep debug_trimode debug_trimode_detailed test_trimode_module

echo ""
echo "==========================================="
echo "Compiling sub_nikkhoo.f90 module..."
echo "==========================================="
gfortran -c -O2 sub_nikkhoo.f90
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile sub_nikkhoo.f90"
    exit 1
fi

# Verify module file was created
if [ ! -f "nikkhoo_walter.mod" ]; then
    echo "ERROR: nikkhoo_walter.mod was not created!"
    exit 1
fi
if [ ! -f "sub_nikkhoo.o" ]; then
    echo "ERROR: sub_nikkhoo.o was not created!"
    exit 1
fi

echo "SUCCESS: sub_nikkhoo.f90 compiled"
echo "  Created: nikkhoo_walter.mod ($(stat -c%s nikkhoo_walter.mod) bytes)"
echo "  Created: sub_nikkhoo.o ($(stat -c%s sub_nikkhoo.o) bytes)"

echo ""
echo "==========================================="
echo "Compiling test programs..."
echo "==========================================="

# Compile test_trimode_module (simple test using module)
echo "Compiling test_trimode_module..."
gfortran -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile test_trimode_module"
    exit 1
fi
echo "SUCCESS: test_trimode_module compiled"

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
echo "IMPORTANT: Run the simple module test first:"
echo "  ./test_trimode_module"
echo ""
echo "If that shows valid numbers (not NaN), then run:"
echo "  ./test_casez"
echo "  ./test_casep"
echo "  ./debug_trimode"
echo ""
