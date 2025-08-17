#!/bin/bash

# Compilation script for SIMD-optimized 3dtri_BP5
# This version includes advanced loop tiling and cache optimization

echo "Compiling SIMD-optimized 3dtri_BP5..."

# Compiler and flags for SIMD optimization
MPIF90="mpif90"
FFLAGS="-O3 -march=native -mtune=native -funroll-loops -ftree-vectorize -ffast-math"
OPENMP_FLAGS="-fopenmp"
SIMD_FLAGS="-faggressive-loop-optimizations -floop-nest-optimize -fprefetch-loop-arrays"

# Compile the SIMD-optimized version
echo "Using compiler: $MPIF90"
echo "Optimization flags: $FFLAGS"
echo "OpenMP flags: $OPENMP_FLAGS"
echo "SIMD flags: $SIMD_FLAGS"

# First compile the required module
echo "Compiling phy3d_module_non..."
$MPIF90 $FFLAGS $OPENMP_FLAGS $SIMD_FLAGS -c phy3d_module_non.f90

if [ $? -ne 0 ]; then
    echo "Failed to compile phy3d_module_non.f90"
    exit 1
fi

# Then compile the main program and link with the module
echo "Compiling 3dtri_BP5_simd..."
$MPIF90 $FFLAGS $OPENMP_FLAGS $SIMD_FLAGS -o 3dtri_BP5_simd 3dtri_BP5_simd.f90 phy3d_module_non.o

if [ $? -eq 0 ]; then
    echo "Successfully compiled 3dtri_BP5_simd"
    echo "Executable: 3dtri_BP5_simd"
    
    # Clean up object files
    echo "Cleaning up object files..."
    rm -f phy3d_module_non.o
    
    echo ""
    echo "To run:"
    echo "  mpirun -np <number_of_processes> ./3dtri_BP5_simd"
    echo ""
    echo "Performance features:"
    echo "  - Multi-level cache blocking (L1: 32, L2: 256, L3: 1024)"
    echo "  - SIMD vectorization with OpenMP"
    echo "  - Advanced loop tiling for cache optimization"
    echo "  - Non-blocking MPI communications"
    echo "  - Optimized stiffness matrix reading"
    echo "  - SIMD-optimized physics calculations"
else
    echo "Compilation failed!"
    # Clean up object files even on failure
    rm -f phy3d_module_non.o
    exit 1
fi
