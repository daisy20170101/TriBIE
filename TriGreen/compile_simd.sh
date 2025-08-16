#!/bin/bash

# Compilation script for SIMD-optimized calc_trigreen
# This version includes advanced loop tiling and cache optimization

echo "Compiling SIMD-optimized calc_trigreen..."

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

$MPIF90 $FFLAGS $OPENMP_FLAGS $SIMD_FLAGS -o calc_trigreen_simd calc_trigreen_simd.f90

if [ $? -eq 0 ]; then
    echo "Successfully compiled calc_trigreen_simd"
    echo "Executable: calc_trigreen_simd"
    echo ""
    echo "To run:"
    echo "  mpirun -np <number_of_processes> ./calc_trigreen_simd"
    echo ""
    echo "Performance features:"
    echo "  - Multi-level cache blocking (L1: 32, L2: 256, L3: 1024)"
    echo "  - SIMD vectorization with OpenMP"
    echo "  - Advanced loop tiling for cache optimization"
    echo "  - Reduced error checking overhead in inner loops"
else
    echo "Compilation failed!"
    exit 1
fi
