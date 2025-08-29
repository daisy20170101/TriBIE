#!/bin/bash

#===============================================================================
# TriBIE Compilation Script with HDF5 Support
# Compiles 3dtri_BP5.f90 and its dependencies
#===============================================================================

echo "=========================================="
echo "TriBIE Compilation Script with HDF5 Support"
echo "=========================================="

# Set compiler and flags
COMPILER="mpif90"
OPTIMIZATION_FLAGS="-O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize"
DEBUG_FLAGS="-g -fcheck=all -fbacktrace -Wall -Wextra"
OPENMP_FLAGS="-fopenmp"
FLOATING_POINT_FLAGS="-fno-trapping-math -fno-signaling-nans -fno-unsafe-math-optimizations"
LINE_LENGTH_FLAGS="-ffree-line-length-none -ffree-form"

# HDF5 specific flags
HDF5_FLAGS=""
HDF5_LIBS="-lhdf5_fortran -lhdf5 -lz -ldl -lm"

# Check for EasyBuild HDF5 environment variable (common on HPC clusters)
if [ ! -z "$EBROOTHDF5" ]; then
    echo "EasyBuild HDF5 detected: $EBROOTHDF5"
    HDF5_FLAGS="-I$EBROOTHDF5/include"
    HDF5_LIBS="-L$EBROOTHDF5/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
    echo "Using HDF5 from: $EBROOTHDF5"
elif [ ! -z "$HDF5_ROOT" ]; then
    echo "HDF5_ROOT detected: $HDF5_ROOT"
    HDF5_FLAGS="-I$HDF5_ROOT/include"
    HDF5_LIBS="-L$HDF5_ROOT/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
    echo "Using HDF5 from: $HDF5_ROOT"
# Alternative HDF5 paths (common on different systems)
elif [ -d "/usr/include/hdf5/openmpi" ]; then
    HDF5_FLAGS="-I/usr/include/hdf5/openmpi -I/usr/include/hdf5/serial"
    HDF5_LIBS="-lhdf5_fortran -lhdf5 -lz -ldl -lm"
elif [ -d "/opt/hdf5/include" ]; then
    HDF5_FLAGS="-I/opt/hdf5/include"
    HDF5_LIBS="-L/opt/hdf5/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
elif [ -d "/usr/local/hdf5/include" ]; then
    HDF5_FLAGS="-I/usr/local/hdf5/include"
    HDF5_LIBS="-L/usr/local/hdf5/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
elif [ -d "/sw/hdf5/include" ]; then
    HDF5_FLAGS="-I/sw/hdf5/include"
    HDF5_LIBS="-L/sw/hdf5/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
fi

# Check for HDF5 module (common on HPC clusters)
if command -v module >/dev/null 2>&1; then
    if module list 2>&1 | grep -q hdf5; then
        echo "HDF5 module detected, using system HDF5"
        HDF5_FLAGS=""
        HDF5_LIBS="-lhdf5_fortran -lhdf5 -lz -ldl -lm"
    fi
fi

# Combine all flags
ALL_FLAGS="$OPTIMIZATION_FLAGS $DEBUG_FLAGS $OPENMP_FLAGS $FLOATING_POINT_FLAGS $LINE_LENGTH_FLAGS"

echo "Compiler: $COMPILER"
echo "Flags: $ALL_FLAGS"
echo "HDF5 Flags: $HDF5_FLAGS"
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
    echo "This might be due to missing HDF5 libraries or incorrect paths"
    echo "Try installing HDF5 development libraries or loading HDF5 module"
    exit 1
fi
echo "✓ 3dtri_BP5.f90 compiled successfully"

# Step 3: Link all object files with HDF5 libraries
echo "Step 3: Linking object files with HDF5 libraries..."
$COMPILER $ALL_FLAGS -o 3dtri_BP5 3dtri_BP5.o phy3d_module_non.o $HDF5_LIBS

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to link object files"
    echo "This might be due to missing HDF5 libraries"
    echo "Check HDF5 installation and library paths"
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
    echo "You can now run: ./3dtri_BP5"
    echo ""
    echo "Note: This executable now supports HDF5 output for time-series data"
    echo "      including cosine slip, SSE, and other monitoring variables"
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
