#!/bin/bash

#===============================================================================
# TriGreen Compilation Script
# Compiles calc_trigreen.f90 and its dependencies
#===============================================================================

echo "=========================================="
echo "TriGreen Compilation Script"
echo "=========================================="

# Set compiler and flags
COMPILER="mpif90"
OPTIMIZATION_FLAGS="-O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize"
DEBUG_FLAGS="-g -fcheck=all -fbacktrace -Wall -Wextra"
OPENMP_FLAGS="-fopenmp"
FLOATING_POINT_FLAGS="-fno-trapping-math -fno-signaling-nans -fno-unsafe-math-optimizations"
LINE_LENGTH_FLAGS="-ffree-line-length-none -ffree-form"

# Combine all flags
ALL_FLAGS="$OPTIMIZATION_FLAGS $DEBUG_FLAGS $OPENMP_FLAGS $FLOATING_POINT_FLAGS $LINE_LENGTH_FLAGS"

echo "Compiler: $COMPILER"
echo "Flags: $ALL_FLAGS"
echo ""

# Check if source files exist
echo "Checking source files..."
if [ ! -f "mod_dtrigreen.f90" ]; then
    echo "ERROR: mod_dtrigreen.f90 not found!"
    exit 1
fi

if [ ! -f "m_calc_green.f90" ]; then
    echo "ERROR: m_calc_green.f90 not found!"
    exit 1
fi

if [ ! -f "calc_trigreen.f90" ]; then
    echo "ERROR: calc_trigreen.f90 not found!"
    exit 1
fi

echo "All source files found."
echo ""

# Clean previous object files
echo "Cleaning previous object files..."
rm -f *.o *.mod calc_trigreen
echo "Cleanup complete."
echo ""

# Compilation steps
echo "Starting compilation..."

# Step 1: Compile sub_comdun.f90 FIRST (contains comdun subroutine)
if [ -f "sub_comdun.f90" ]; then
    echo "Step 1: Compiling sub_comdun.f90 (FIRST - contains comdun subroutine)..."
    $COMPILER $ALL_FLAGS -c sub_comdun.f90 -o sub_comdun.o
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to compile sub_comdun.f90"
        exit 1
    fi
    echo "✓ sub_comdun.f90 compiled successfully"
    SUB_COMDUN_OBJ="sub_comdun.o"
else
    echo "ERROR: sub_comdun.f90 not found! This file is required."
    echo "It contains the comdun subroutine used by other modules."
    exit 1
fi

# Step 2: Compile mod_dtrigreen.f90 (depends on sub_comdun)
echo "Step 2: Compiling mod_dtrigreen.f90 (depends on sub_comdun)..."
$COMPILER $ALL_FLAGS -c mod_dtrigreen.f90 -o mod_dtrigreen.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile mod_dtrigreen.f90"
    exit 1
fi
echo "✓ mod_dtrigreen.f90 compiled successfully"

# Step 3: Compile m_calc_green.f90
echo "Step 3: Compiling m_calc_green.f90..."
$COMPILER $ALL_FLAGS -c m_calc_green.f90 -o m_calc_green.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile m_calc_green.f90"
    exit 1
fi
echo "✓ m_calc_green.f90 compiled successfully"

# Step 4: Compile calc_trigreen.f90
echo "Step 4: Compiling calc_trigreen.f90..."
$COMPILER $ALL_FLAGS -c calc_trigreen.f90 -o calc_trigreen.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile calc_trigreen.f90"
    exit 1
fi
echo "✓ calc_trigreen.f90 compiled successfully"

# Step 5: Link all object files in dependency order
echo "Step 5: Linking object files in dependency order..."
$COMPILER $ALL_FLAGS -o calc_trigreen calc_trigreen.o m_calc_green.o mod_dtrigreen.o $SUB_COMDUN_OBJ

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to link object files"
    exit 1
fi
echo "✓ Linking completed successfully"

# Final cleanup
echo ""
echo "Cleaning up object files..."
rm -f *.o *.mod

# Check if executable was created
if [ -f "calc_trigreen" ]; then
    echo ""
    echo "=========================================="
    echo "✓ COMPILATION SUCCESSFUL!"
    echo "=========================================="
    echo "Executable: calc_trigreen"
    echo "Size: $(ls -lh calc_trigreen | awk '{print $5}')"
    echo ""
    echo "You can now run: ./calc_trigreen"
    echo ""
else
    echo ""
    echo "=========================================="
    echo "✗ COMPILATION FAILED!"
    echo "=========================================="
    echo "Executable was not created"
    exit 1
fi

echo "Compilation script completed."
