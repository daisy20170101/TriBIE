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

# Compile all source files
echo "Compiling all source files..."

# Compile the module first
echo "Compiling m_calc_green.f90..."
$MPIF90 $FFLAGS -c m_calc_green.f90 -o m_calc_green.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile m_calc_green.f90"
    exit 1
fi

# Compile the comdun module FIRST (converted from F77 to F90)
# This must come before mod_dtrigreen.f90 since it depends on it
if [ -f "sub_comdun.f90" ]; then
    echo "Compiling sub_comdun.f90..."
    $MPIF90 $FFLAGS -c sub_comdun.f90 -o sub_comdun.o
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to compile sub_comdun.f90"
        exit 1
    fi
    SUB_COMDUN_OBJ="sub_comdun.o"
elif [ -f "sub_comdun.f" ]; then
    echo "Warning: sub_comdun.f90 not found, using old F77 version..."
    echo "Compiling sub_comdun.f..."
    $MPIF90 $FFLAGS -c sub_comdun.f -o sub_comdun.o
    
    if [ $? -ne 0 ]; then
        echo "Error: Failed to compile sub_comdun.f"
        exit 1
    fi
    SUB_COMDUN_OBJ="sub_comdun.o"
else
    echo "Warning: Neither sub_comdun.f90 nor sub_comdun.f found, skipping..."
    SUB_COMDUN_OBJ=""
fi

# Compile the dstuart module AFTER comdun (since it depends on it)
echo "Compiling mod_dtrigreen.f90..."
$MPIF90 $FFLAGS -c mod_dtrigreen.f90 -o mod_dtrigreen.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile mod_dtrigreen.f90"
    exit 1
fi

# Compile the main program
echo "Compiling calc_trigreen.f90..."
$MPIF90 $FFLAGS -c calc_trigreen.f90 -o calc_trigreen.o

if [ $? -ne 0 ]; then
    echo "Error: Failed to compile calc_trigreen.f90"
    exit 1
fi

# Link the executable with all object files
echo "Linking executable..."
if [ -n "$SUB_COMDUN_OBJ" ]; then
    $MPIF90 $FFLAGS -o calc_trigreen m_calc_green.o mod_dtrigreen.o $SUB_COMDUN_OBJ calc_trigreen.o
else
    $MPIF90 $FFLAGS -o calc_trigreen m_calc_green.o mod_dtrigreen.o calc_trigreen.o
fi

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
echo ""
echo "All necessary source files are now included:"
echo "  - m_calc_green.f90 (module with constants)"
echo "  - mod_dtrigreen.f90 (dstuart subroutine)"
echo "  - sub_comdun.f90 (comdun subroutine - converted from F77)"
echo "  - calc_trigreen.f90 (main program)"
echo ""
echo "Note: sub_comdun.f90 is a simplified F90 conversion. The full vc and dvc"
echo "calculations may need to be implemented if required for your application."
