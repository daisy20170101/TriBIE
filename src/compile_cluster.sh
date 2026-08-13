#!/bin/bash

#===============================================================================
# HPC Cluster Compilation Script for TriBIE with EasyBuild HDF5
# Optimized for clusters using EBROOTHDF5 environment variable
#===============================================================================

echo "=========================================="
echo "HPC Cluster Compilation Script for TriBIE"
echo "=========================================="

# Check for EasyBuild HDF5
if [ -z "$EBROOTHDF5" ]; then
    echo "ERROR: EBROOTHDF5 environment variable not set!"
    echo "Please load the HDF5 module first:"
    echo "  module load HDF5"
    echo "  # or"
    echo "  module load hdf5"
    echo ""
    echo "Then check the module:"
    echo "  module list | grep hdf5"
    echo "  echo \$EBROOTHDF5"
    exit 1
fi

echo "EasyBuild HDF5 detected: $EBROOTHDF5"
echo "HDF5 version: $(basename $EBROOTHDF5)"

# Set compiler and flags optimized for HPC clusters
COMPILER="mpif90"
OPTIMIZATION_FLAGS="-O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize"
DEBUG_FLAGS="-g -fcheck=all -fbacktrace -Wall -Wextra"
OPENMP_FLAGS="-fopenmp"
FLOATING_POINT_FLAGS="-fno-trapping-math -fno-signaling-nans -fno-unsafe-math-optimizations"
LINE_LENGTH_FLAGS="-ffree-line-length-none -ffree-form"

# HDF5 paths from EasyBuild
HDF5_FLAGS="-I$EBROOTHDF5/include"
HDF5_LIBS="-L$EBROOTHDF5/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"

# Check if HDF5 files exist
if [ ! -d "$EBROOTHDF5/include" ]; then
    echo "ERROR: HDF5 include directory not found: $EBROOTHDF5/include"
    exit 1
fi

if [ ! -d "$EBROOTHDF5/lib" ]; then
    echo "ERROR: HDF5 library directory not found: $EBROOTHDF5/lib"
    exit 1
fi

# Check for required HDF5 files
if [ ! -f "$EBROOTHDF5/include/hdf5.mod" ] && [ ! -f "$EBROOTHDF5/include/hdf5.mod" ]; then
    echo "WARNING: hdf5.mod not found in $EBROOTHDF5/include"
    echo "This might cause compilation issues"
fi

# Combine all flags
ALL_FLAGS="$OPTIMIZATION_FLAGS $DEBUG_FLAGS $OPENMP_FLAGS $FLOATING_POINT_FLAGS $LINE_LENGTH_FLAGS"

echo "Compiler: $COMPILER"
echo "Flags: $ALL_FLAGS"
echo "HDF5 Include: $HDF5_FLAGS"
echo "HDF5 Libraries: $HDF5_LIBS"
echo ""

# Check if source files exist
echo "Checking source files..."
if [ ! -f "phy3d_module_non.f90" ]; then
    echo "ERROR: phy3d_module_non.f90 not found!"
    exit 1
fi

if [ ! -f "3dtri_BP5.f90" ]; then
    echo "ERROR: 3dtri_BP5.f90 not found!"
    exit 1
fi

echo "All source files found."
echo ""

# Clean previous object files
echo "Cleaning previous object files..."
rm -f *.o *.mod 3dtri_BP5
echo "Cleanup complete."
echo ""

# Compilation steps
echo "Starting compilation..."

# Step 1: Compile phy3d_module_non.f90 first
echo "Step 1: Compiling phy3d_module_non.f90..."
$COMPILER $ALL_FLAGS -c phy3d_module_non.f90 -o phy3d_module_non.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile phy3d_module_non.f90"
    exit 1
fi
echo "✓ phy3d_module_non.f90 compiled successfully"

# Step 2: Compile 3dtri_BP5.f90 with HDF5 support
echo "Step 2: Compiling 3dtri_BP5.f90 with HDF5 support..."
$COMPILER $ALL_FLAGS $HDF5_FLAGS -c 3dtri_BP5.f90 -o 3dtri_BP5.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile 3dtri_BP5.f90"
    echo "HDF5 compilation failed. Check:"
    echo "1. HDF5 module is loaded: module list | grep hdf5"
    echo "2. EBROOTHDF5 is set: echo \$EBROOTHDF5"
    echo "3. HDF5 files exist: ls -la \$EBROOTHDF5/include/hdf5*"
    echo "4. HDF5 libraries exist: ls -la \$EBROOTHDF5/lib/libhdf5*"
    exit 1
fi
echo "✓ 3dtri_BP5.f90 compiled successfully"

# Step 3: Link all object files with HDF5 libraries
echo "Step 3: Linking object files with HDF5 libraries..."
$COMPILER $ALL_FLAGS -o 3dtri_BP5 3dtri_BP5.o phy3d_module_non.o $HDF5_LIBS

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to link object files"
    echo "HDF5 linking failed. Check:"
    echo "1. HDF5 libraries are accessible: ls -la \$EBROOTHDF5/lib/libhdf5*"
    echo "2. Library paths are correct: echo \$EBROOTHDF5/lib"
    echo "3. Try manual linking: mpif90 -o 3dtri_BP5 3dtri_BP5.o phy3d_module_non.o $HDF5_LIBS"
    exit 1
fi
echo "✓ Linking completed successfully"

# Final cleanup
echo ""
echo "Cleaning up object files..."
rm -f *.o *.mod

# Check if executable was created
if [ -f "3dtri_BP5" ]; then
    echo ""
    echo "=========================================="
    echo "✓ COMPILATION SUCCESSFUL!"
    echo "=========================================="
    echo "Executable: 3dtri_BP5"
    echo "Size: $(ls -lh 3dtri_BP5 | awk '{print $5}')"
    echo ""
    echo "HDF5 Support: Enabled"
    echo "HDF5 Path: $EBROOTHDF5"
    echo "You can now run: ./3dtri_BP5"
    echo ""
    echo "Note: This executable now supports HDF5 output for time-series data"
    echo "      including cosine slip, SSE, and other monitoring variables"
    echo ""
    
    # Verify HDF5 linking
    echo "Verifying HDF5 linking..."
    if ldd 3dtri_BP5 2>/dev/null | grep -q hdf5; then
        echo "✓ HDF5 libraries successfully linked"
    else
        echo "⚠️  HDF5 libraries not found in ldd output (this may be normal on some systems)"
    fi
else
    echo ""
    echo "=========================================="
    echo "✗ COMPILATION FAILED!"
    echo "=========================================="
    echo "Executable was not created"
    exit 1
fi

echo "Compilation script completed."
