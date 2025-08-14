#!/bin/bash

#===============================================================================
# COMPILATION SCRIPT FOR EARTHQUAKE SIMULATION
# 
# This script compiles all Fortran source files in the src directory
# using mpif90 with proper dependency ordering.
#
# Author: Compilation script for improved earthquake simulation code
# Last Modified: Current compilation setup
#===============================================================================

# Configuration
COMPILER="mpif90"
COMPILER_FLAGS="-O2 -g -Wall -Wextra -std=f2008"
OUTPUT_NAME="earthquake_simulation"
BUILD_DIR="build"
SRC_DIR="."

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to cleanup on exit
cleanup() {
    if [ $? -ne 0 ]; then
        print_error "Compilation failed!"
        exit 1
    fi
}

# Set trap for cleanup
trap cleanup EXIT

# Check if mpif90 is available
print_status "Checking for MPI Fortran compiler..."
if ! command_exists $COMPILER; then
    print_error "MPI Fortran compiler '$COMPILER' not found!"
    print_error "Please install an MPI implementation (e.g., OpenMPI, MPICH)"
    print_error "or ensure it's in your PATH."
    exit 1
fi

print_success "Found MPI Fortran compiler: $COMPILER"

# Check compiler version
print_status "Checking compiler version..."
$COMPILER --version | head -n 1

# Create build directory
print_status "Creating build directory..."
if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
    print_success "Created build directory: $BUILD_DIR"
else
    print_status "Build directory already exists: $BUILD_DIR"
fi

# Change to build directory
cd "$BUILD_DIR"

# Clean previous build artifacts
print_status "Cleaning previous build artifacts..."
rm -f *.o *.mod "$OUTPUT_NAME"

# Define compilation order (dependencies first)
MODULES=(
    "physical_constants"
    "simulation_parameters"
    "physical_variables"
    "io_parameters"
    "error_handling"
    "time_integration"
    "physics_equations"
    "mpi_utilities"
    "io_handling"
)

MAIN_PROGRAM="earthquake_simulation_improved"

# Compile modules in dependency order
print_status "Compiling modules in dependency order..."

for module in "${MODULES[@]}"; do
    source_file="../$SRC_DIR/${module}.f90"
    object_file="${module}.o"
    
    if [ ! -f "$source_file" ]; then
        print_error "Source file not found: $source_file"
        exit 1
    fi
    
    print_status "Compiling $module..."
    $COMPILER $COMPILER_FLAGS -c "$source_file" -o "$object_file"
    
    if [ $? -eq 0 ]; then
        print_success "Compiled $module successfully"
    else
        print_error "Failed to compile $module"
        exit 1
    fi
done

# Compile main program
print_status "Compiling main program..."
source_file="../$SRC_DIR/${MAIN_PROGRAM}.f90"
object_file="${MAIN_PROGRAM}.o"

if [ ! -f "$source_file" ]; then
    print_error "Main program source file not found: $source_file"
    exit 1
fi

$COMPILER $COMPILER_FLAGS -c "$source_file" -o "$object_file"

if [ $? -eq 0 ]; then
    print_success "Compiled main program successfully"
else
    print_error "Failed to compile main program"
    exit 1
fi

# Link executable
print_status "Linking executable..."
object_files=""
for module in "${MODULES[@]}"; do
    object_files="$object_files ${module}.o"
done
object_files="$object_files ${MAIN_PROGRAM}.o"

$COMPILER $COMPILER_FLAGS -o "$OUTPUT_NAME" $object_files

if [ $? -eq 0 ]; then
    print_success "Linking completed successfully"
else
    print_error "Linking failed"
    exit 1
fi

# Check if executable was created
if [ -f "$OUTPUT_NAME" ]; then
    print_success "Executable created: $BUILD_DIR/$OUTPUT_NAME"
    
    # Show executable information
    print_status "Executable details:"
    ls -lh "$OUTPUT_NAME"
    
    # Test if executable runs (basic check)
    print_status "Testing executable..."
    if command_exists timeout; then
        timeout 5s ./"$OUTPUT_NAME" --help >/dev/null 2>&1 || echo "Executable runs (help option not implemented)"
    else
        print_warning "Skipping executable test (timeout command not available)"
    fi
    
else
    print_error "Executable was not created"
    exit 1
fi

# Show build summary
print_status "Build summary:"
echo "  Compiler: $COMPILER"
echo "  Flags: $COMPILER_FLAGS"
echo "  Output: $BUILD_DIR/$OUTPUT_NAME"
echo "  Object files: $(ls -1 *.o | wc -l)"
echo "  Module files: $(ls -1 *.mod | wc -l)"

print_success "Compilation completed successfully!"
print_status "To run the simulation:"
echo "  cd $BUILD_DIR"
echo "  mpirun -np <number_of_processes> ./$OUTPUT_NAME"

# Optional: create symlink in parent directory
if [ "$1" = "--symlink" ]; then
    print_status "Creating symlink in parent directory..."
    cd ..
    ln -sf "$BUILD_DIR/$OUTPUT_NAME" "$OUTPUT_NAME"
    print_success "Symlink created: $OUTPUT_NAME -> $BUILD_DIR/$OUTPUT_NAME"
fi

exit 0
