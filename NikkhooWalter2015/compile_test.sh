#!/bin/bash
# Simple compilation test for the Nikkhoo & Walter Fortran module

echo "Testing compilation of sub_nikkhoo.f90..."

# Try to compile the module
gfortran -c sub_nikkhoo.f90 -o sub_nikkhoo.o

if [ $? -eq 0 ]; then
    echo "✅ Module compilation successful!"
    
    # Try to compile the test program
    echo "Testing compilation of test_nikkhoo.f90..."
    gfortran -c test_nikkhoo.f90 -o test_nikkhoo.o
    
    if [ $? -eq 0 ]; then
        echo "✅ Test program compilation successful!"
        
        # Try to link
        echo "Testing linking..."
        gfortran sub_nikkhoo.o test_nikkhoo.o -o test_nikkhoo
        
        if [ $? -eq 0 ]; then
            echo "✅ Linking successful!"
            echo "✅ All compilation tests passed!"
            echo ""
            echo "You can now run: ./test_nikkhoo"
        else
            echo "❌ Linking failed"
        fi
    else
        echo "❌ Test program compilation failed"
    fi
else
    echo "❌ Module compilation failed"
    echo "Please check the error messages above"
fi

# Clean up object files
rm -f *.o
