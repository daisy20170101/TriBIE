#!/bin/bash

# Compilation script for improved TriBIE with hybrid MPI+OpenMP parallelization
# Enhanced version with dynamic load balancing and performance monitoring

echo "Compiling improved TriBIE with hybrid MPI+OpenMP parallelization..."

# Compiler options for hybrid parallelization
MPIF90="mpif90"
FFLAGS="-O3 -fopenmp -march=native -mtune=native -ffast-math -funroll-loops"
DEBUG_FLAGS="-g -fcheck=all -fbacktrace -Wall -Wextra"
OPTIMIZATION_FLAGS="-O3 -fopenmp -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize"

# Choose compilation mode
if [ "$1" = "debug" ]; then
    echo "Compiling in DEBUG mode..."
    FFLAGS="$DEBUG_FLAGS"
elif [ "$1" = "optimize" ]; then
    echo "Compiling in OPTIMIZED mode..."
    FFLAGS="$OPTIMIZATION_FLAGS"
else
    echo "Compiling in STANDARD mode..."
    FFLAGS="-O2 -fopenmp"
fi

# Compile the module first
echo "Compiling phy3d_module_non.f90..."
$MPIF90 $FFLAGS -c phy3d_module_non.f90 -o phy3d_module_non.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile phy3d_module_non.f90"
    exit 1
fi

# Compile the main program
echo "Compiling 3dtriBIE_v1.f90..."
$MPIF90 $FFLAGS -c 3dtriBIE_v1.f90 -o 3dtriBIE_v1.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile 3dtriBIE_v1.f90"
    exit 1
fi

# Link the executable
echo "Linking executable..."
$MPIF90 $FFLAGS -o 3dtriBIE_improved phy3d_module_non.o 3dtriBIE_v1.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to link executable"
    exit 1
fi

echo "Compilation successful!"
echo "Executable: 3dtriBIE_improved"
echo ""
echo "To run with hybrid parallelization:"
echo "  mpirun -np <MPI_processes> -x OMP_NUM_THREADS=<threads_per_process> ./3dtriBIE_improved"
echo ""
echo "Example for 2 nodes with 20 MPI processes and 8 OpenMP threads each:"
echo "  mpirun -np 40 -x OMP_NUM_THREADS=8 ./3dtriBIE_improved"
echo ""
echo "For Slurm submission, use the provided submit_job.slurm script"
