#!/bin/bash

#===============================================================================
# Simple HDF5 Compilation Script for TriBIE
# Use this if the main compile.sh doesn't work
# Supports EBROOTHDF5 environment variable
#===============================================================================

echo "=========================================="
echo "Simple HDF5 Compilation Script for TriBIE"
echo "=========================================="

# Basic compilation flags
COMPILER="mpif90"
FLAGS="-O3 -fopenmp -ffree-line-length-none -ffree-form"

# HDF5 libraries - check environment variables first
if [ ! -z "$EBROOTHDF5" ]; then
    echo "Using EasyBuild HDF5: $EBROOTHDF5"
    HDF5_INCLUDE="-I$EBROOTHDF5/include"
    HDF5_LIBS="-L$EBROOTHDF5/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
elif [ ! -z "$HDF5_ROOT" ]; then
    echo "Using HDF5_ROOT: $HDF5_ROOT"
    HDF5_INCLUDE="-I$HDF5_ROOT/include"
    HDF5_LIBS="-L$HDF5_ROOT/lib -lhdf5_fortran -lhdf5 -lz -ldl -lm"
else
    # Default system paths
    HDF5_INCLUDE="-I/usr/include/hdf5/openmpi"
    HDF5_LIBS="-lhdf5_fortran -lhdf5 -lz -ldl -lm"
fi

echo "Compiler: $COMPILER"
echo "Flags: $FLAGS"
echo "HDF5 Include: $HDF5_INCLUDE"
echo "HDF5 Libraries: $HDF5_LIBS"
echo ""

# Clean up
rm -f *.o *.mod 3dtri_BP5

# Compile phy3d_module_non.f90
echo "Compiling phy3d_module_non.f90..."
$COMPILER $FLAGS -c phy3d_module_non.f90 -o phy3d_module_non.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile phy3d_module_non.f90"
    exit 1
fi

# Compile 3dtri_BP5.f90 with HDF5
echo "Compiling 3dtri_BP5.f90 with HDF5..."
$COMPILER $FLAGS $HDF5_INCLUDE -c 3dtri_BP5.f90 -o 3dtri_BP5.o
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to compile 3dtri_BP5.f90"
    echo "Check HDF5 installation and paths"
    echo "Current HDF5_INCLUDE: $HDF5_INCLUDE"
    exit 1
fi

# Link
echo "Linking with HDF5 libraries..."
$COMPILER $FLAGS -o 3dtri_BP5 3dtri_BP5.o phy3d_module_non.o $HDF5_LIBS
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to link"
    echo "Check HDF5 libraries and paths"
    echo "Current HDF5_LIBS: $HDF5_LIBS"
    exit 1
fi

# Clean up
rm -f *.o *.mod

echo "Compilation successful! Executable: 3dtri_BP5"
echo "HDF5 Support: Enabled"
