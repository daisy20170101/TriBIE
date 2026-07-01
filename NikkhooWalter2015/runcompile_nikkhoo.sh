#!/bin/bash

#===============================================================================
# Nikkhoo-Walter Stiffness Calculation - Compilation Script
# Compiles calc_nikkhoo.f90 and its dependencies with MPI/OpenMP support
#
# Usage:
#   ./runcompile_nikkhoo.sh          # Use mpif90 (default, MPI+OpenMP)
#   ./runcompile_nikkhoo.sh gfortran # Use gfortran (OpenMP only, no MPI)
#   ./runcompile_nikkhoo.sh mpi      # Use mpif90 explicitly
#===============================================================================

echo "=========================================="
echo "Nikkhoo-Walter Stiffness Compilation"
echo "=========================================="

# Parse command line argument for compiler selection
COMPILER_OPTION="${1:-mpi}"

# Set compiler based on option
case "$COMPILER_OPTION" in
    gfortran|gnu|serial)
        COMPILER="gfortran"
        echo "Mode: Serial/OpenMP (gfortran)"
        ;;
    mpi|mpif90|*)
        COMPILER="mpif90"
        echo "Mode: MPI/OpenMP (mpif90)"
        ;;
esac

# Set compiler flags
OPTIMIZATION_FLAGS="-O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize"
DEBUG_FLAGS="-g -fbacktrace -Wall"
OPENMP_FLAGS="-fopenmp"
FLOATING_POINT_FLAGS="-fno-trapping-math"
LINE_LENGTH_FLAGS="-ffree-line-length-none -ffree-form"

# Combine all flags
ALL_FLAGS="$OPTIMIZATION_FLAGS $DEBUG_FLAGS $OPENMP_FLAGS $FLOATING_POINT_FLAGS $LINE_LENGTH_FLAGS"

echo "Compiler: $COMPILER"
echo "Flags: $ALL_FLAGS"
echo ""

# Check if source files exist
echo "Checking source files..."
if [ ! -f "m_nikkhoo_green.f90" ]; then
    echo "ERROR: m_nikkhoo_green.f90 not found!"
    exit 1
fi

if [ ! -f "sub_nikkhoo.f90" ]; then
    echo "ERROR: sub_nikkhoo.f90 not found!"
    exit 1
fi

if [ ! -f "calc_nikkhoo.f90" ]; then
    echo "ERROR: calc_nikkhoo.f90 not found!"
    exit 1
fi

echo "All source files found."
echo ""

# Clean previous object files
echo "Cleaning previous object files..."
rm -f *.o *.mod calc_nikkhoo
echo "Cleanup complete."
echo ""

# Compilation steps
echo "Starting compilation..."

# Step 1: Compile m_nikkhoo_green.f90 (parameter module)
echo "Step 1: Compiling m_nikkhoo_green.f90..."
$COMPILER $ALL_FLAGS -c m_nikkhoo_green.f90 -o m_nikkhoo_green.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile m_nikkhoo_green.f90"
    exit 1
fi
echo "  m_nikkhoo_green.f90 compiled successfully"

# Step 2: Compile sub_nikkhoo.f90 (nikkhoo_walter module)
echo "Step 2: Compiling sub_nikkhoo.f90..."
$COMPILER $ALL_FLAGS -c sub_nikkhoo.f90 -o sub_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile sub_nikkhoo.f90"
    exit 1
fi
echo "  sub_nikkhoo.f90 compiled successfully"

# Step 3: Compile calc_nikkhoo.f90 (main program)
echo "Step 3: Compiling calc_nikkhoo.f90..."
$COMPILER $ALL_FLAGS -c calc_nikkhoo.f90 -o calc_nikkhoo.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile calc_nikkhoo.f90"
    exit 1
fi
echo "  calc_nikkhoo.f90 compiled successfully"

# Step 4: Link all object files
echo "Step 4: Linking object files..."
$COMPILER $ALL_FLAGS -o calc_nikkhoo calc_nikkhoo.o sub_nikkhoo.o m_nikkhoo_green.o

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to link object files"
    exit 1
fi
echo "  Linking completed successfully"

# Final cleanup
echo ""
echo "Cleaning up object files..."
rm -f *.o *.mod

# Check if executable was created
if [ -f "calc_nikkhoo" ]; then
    echo ""
    echo "=========================================="
    echo "COMPILATION SUCCESSFUL!"
    echo "=========================================="
    echo "Executable: calc_nikkhoo"
    echo "Size: $(ls -lh calc_nikkhoo | awk '{print $5}')"
    echo ""
    echo "Usage:"
    if [ "$COMPILER" = "gfortran" ]; then
        echo "  Single process:  ./calc_nikkhoo"
        echo "  OpenMP parallel: OMP_NUM_THREADS=<threads> ./calc_nikkhoo"
    else
        echo "  Single process:  ./calc_nikkhoo"
        echo "  MPI parallel:    mpirun -np <nprocs> ./calc_nikkhoo"
        echo "  Hybrid MPI+OMP:  OMP_NUM_THREADS=<threads> mpirun -np <nprocs> ./calc_nikkhoo"
    fi
    echo ""
    echo "Input:  triangular_mesh.gts (GTS format mesh file)"
    echo "Output: trigreen_<rank>.bin, position.bin"
    echo ""
else
    echo ""
    echo "=========================================="
    echo "COMPILATION FAILED!"
    echo "=========================================="
    echo "Executable was not created"
    exit 1
fi

echo "Compilation script completed."
