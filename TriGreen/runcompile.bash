#!/bin/bash

# Compilation script for TriGreen with improved error handling
# Enhanced version with floating-point exception handling

echo "Compiling TriGreen with improved error handling..."

# Compiler options for hybrid parallelization and error handling
MPIF90="mpif90"
FFLAGS="-O3 -fopenmp -march=native -mtune=native -ffast-math -funroll-loops"
DEBUG_FLAGS="-g -fcheck=all -fbacktrace -Wall -Wextra"
OPTIMIZATION_FLAGS="-O3 -fopenmp -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize"

# Floating-point exception handling flags
FP_FLAGS="-fno-trapping-math -fno-signaling-nans -fno-unsafe-math-optimizations"

# Line length and formatting flags
FORMAT_FLAGS="-ffree-line-length-none -ffree-form"

# Choose compilation mode
if [ "$1" = "debug" ]; then
    echo "Compiling in DEBUG mode with error checking..."
    FFLAGS="$DEBUG_FLAGS $FP_FLAGS $FORMAT_FLAGS"
elif [ "$1" = "optimize" ]; then
    echo "Compiling in OPTIMIZED mode..."
    FFLAGS="$OPTIMIZATION_FLAGS $FP_FLAGS $FORMAT_FLAGS"
else
    echo "Compiling in STANDARD mode with error handling..."
    FFLAGS="-O2 -fopenmp $FP_FLAGS $FORMAT_FLAGS"
fi

# Compile the module first
echo "Compiling m_calc_green.f90..."
$MPIF90 $FFLAGS -c m_calc_green.f90 -o m_calc_green.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile m_calc_green.f90"
    exit 1
fi

# Compile the main program
echo "Compiling calc_trigreen.f90..."
$MPIF90 $FFLAGS -c calc_trigreen.f90 -o calc_trigreen.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile calc_trigreen.f90"
    exit 1
fi

# Link the executable
echo "Linking executable..."
$MPIF90 $FFLAGS -o calc_trigreen m_calc_green.o calc_trigreen.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to link executable"
    exit 1
fi

echo "Compilation successful!"
echo "Executable: calc_trigreen"
echo ""
echo "To run with hybrid parallelization:"
echo "  mpirun -np <MPI_processes> -x OMP_NUM_THREADS=<threads_per_process> ./calc_trigreen"
echo ""
echo "Example for 2 nodes with 2 MPI processes and 4 OpenMP threads each:"
echo "  mpirun -np 2 -x OMP_NUM_THREADS=4 ./calc_trigreen"
echo ""
echo "For Slurm submission, use the provided submit_job.slurm script"
echo ""
echo "Note: This version includes improved error handling and floating-point exception handling"
echo "to prevent crashes due to numerical issues."
echo ""
echo "Line length issues have been resolved with -ffree-line-length-none flag."
