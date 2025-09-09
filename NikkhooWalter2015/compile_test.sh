#!/bin/bash
# Comprehensive compilation test for the Nikkhoo & Walter Fortran module

echo "=========================================="
echo "Testing Nikkhoo & Walter Fortran Module"
echo "=========================================="

# Check if gfortran is available
if ! command -v gfortran &> /dev/null; then
    echo "❌ gfortran not found. Please install gfortran first."
    echo "   On macOS: brew install gcc"
    echo "   Or download from: https://gcc.gnu.org/wiki/GFortranBinaries#MacOS"
    exit 1
fi

echo "✅ gfortran found: $(gfortran --version | head -n1)"
echo ""

# Clean up any existing object files
rm -f *.o test_nikkhoo

echo "Step 1: Compiling sub_nikkhoo.f90 module..."
gfortran -c sub_nikkhoo.f90 -o sub_nikkhoo.o

if [ $? -eq 0 ]; then
    echo "✅ Module compilation successful!"
    echo ""
    
    echo "Step 2: Compiling test_nikkhoo.f90..."
    gfortran -c test_nikkhoo.f90 -o test_nikkhoo.o
    
    if [ $? -eq 0 ]; then
        echo "✅ Test program compilation successful!"
        echo ""
        
        echo "Step 3: Linking executable..."
        gfortran sub_nikkhoo.o test_nikkhoo.o -o test_nikkhoo
        
        if [ $? -eq 0 ]; then
            echo "✅ Linking successful!"
            echo ""
            echo "🎉 All compilation tests passed!"
            echo ""
            echo "You can now run the test program:"
            echo "  ./test_nikkhoo"
            echo ""
            echo "Or run it directly:"
            echo "  make -f Makefile.nikkhoo run"
        else
            echo "❌ Linking failed"
            echo "Check the error messages above"
            exit 1
        fi
    else
        echo "❌ Test program compilation failed"
        echo "Check the error messages above"
        exit 1
    fi
else
    echo "❌ Module compilation failed"
    echo "Check the error messages above"
    exit 1
fi

# Clean up object files
rm -f *.o

echo ""
echo "Compilation test completed successfully! 🚀"
