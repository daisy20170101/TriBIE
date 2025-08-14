#!/bin/bash

#===============================================================================
# COMPILATION SCRIPT FOR OPTIMIZED 3DTRI_BP5.F90
# 
# This script compiles the earthquake simulation code with MPI and OpenMP
# parallelization optimizations.
#
# Author: Optimized compilation for improved earthquake simulation code
# Last Modified: Current compilation setup
#===============================================================================

# Set error handling
set -e

# Configuration variables
PROGRAM_NAME="3dtri_BP5"
SOURCE_FILE="3dtri_BP5.f90"
MODULE_FILE="phy3d_module_non.f90"
OUTPUT_DIR="build"

# Compiler selection (modify as needed for your system)
FC="mpif90"
CC="mpicc"

# Compiler flags for optimization
FFLAGS="-O3 -march=native -mtune=native -funroll-loops -ftree-vectorize"
FFLAGS="$FFLAGS -fopenmp -fno-stack-protector -fno-range-check"
FFLAGS="$FFLAGS -std=f2008 -Wall -Wextra"

# MPI flags
MPI_FLAGS="-DMPI"

# OpenMP flags
OMP_FLAGS="-fopenmp"

# Linker flags
LDFLAGS="-fopenmp"

# Debug flags (uncomment for debugging)
# FFLAGS="$FFLAGS -g -fcheck=all -fbounds-check -fbacktrace"

# Create build directory
echo "Creating build directory..."
mkdir -p $OUTPUT_DIR

# Compilation steps
echo "Compiling module file: $MODULE_FILE"
$FC $FFLAGS $MPI_FLAGS -c $MODULE_FILE -o $OUTPUT_DIR/phy3d_module_non.o

echo "Compiling main program: $SOURCE_FILE"
$FC $FFLAGS $MPI_FLAGS -c $SOURCE_FILE -o $OUTPUT_DIR/3dtri_BP5.o

echo "Linking executable: $PROGRAM_NAME"
$FC $FFLAGS $MPI_FLAGS -o $OUTPUT_DIR/$PROGRAM_NAME \
    $OUTPUT_DIR/phy3d_module_non.o \
    $OUTPUT_DIR/3dtri_BP5.o \
    $LDFLAGS

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Compilation successful!"
    echo "Executable: $OUTPUT_DIR/$PROGRAM_NAME"
    echo "=========================================="
    
    # Display executable information
    echo "Executable details:"
    ls -la $OUTPUT_DIR/$PROGRAM_NAME
    
    # Optional: Run a quick test if parameter file exists
    if [ -f "../input/parameter1.txt" ]; then
        echo ""
        echo "Parameter file found. You can run the simulation with:"
        echo "cd $OUTPUT_DIR && mpirun -np <number_of_processes> ./$PROGRAM_NAME"
    else
        echo ""
        echo "Note: Parameter file not found. Make sure to create parameter1.txt before running."
    fi
else
    echo "=========================================="
    echo "Compilation failed!"
    echo "=========================================="
    exit 1
fi

# Optional: Create a run script
echo ""
echo "Creating run script..."
cat > $OUTPUT_DIR/run.sh << 'EOF'
#!/bin/bash
# Run script for 3dtri_BP5 simulation

# Default number of processes
NP=${1:-4}

# Check if executable exists
if [ ! -f "./3dtri_BP5" ]; then
    echo "Error: Executable not found. Please compile first."
    exit 1
fi

# Check if parameter file exists
if [ ! -f "../input/parameter1.txt" ]; then
    echo "Error: Parameter file not found in ../input/"
    echo "Please create parameter1.txt before running."
    exit 1
fi

echo "Starting simulation with $NP processes..."
echo "Using executable: $(pwd)/3dtri_BP5"
echo "Parameter file: ../input/parameter1.txt"
echo ""

# Run the simulation
mpirun -np $NP ./3dtri_BP5

echo ""
echo "Simulation completed."
EOF

chmod +x $OUTPUT_DIR/run.sh
echo "Run script created: $OUTPUT_DIR/run.sh"
echo "Usage: cd $OUTPUT_DIR && ./run.sh [number_of_processes]"

echo ""
echo "Compilation script completed successfully!"
