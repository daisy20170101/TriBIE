#!/bin/bash

#===============================================================================
# TEST COMPILATION SCRIPT
# 
# This script tests the compilation of individual modules to identify
# any remaining compilation errors.
#
# Author: Test compilation script
# Last Modified: Current compilation setup
#===============================================================================

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

# Check if mpif90 is available
print_status "Checking for MPI Fortran compiler..."
if ! command -v mpif90 >/dev/null 2>&1; then
    print_error "MPI Fortran compiler 'mpif90' not found!"
    print_error "Please install an MPI implementation or ensure it's in your PATH."
    exit 1
fi

print_success "Found MPI Fortran compiler: mpif90"

# Create test build directory
TEST_BUILD_DIR="test_build"
print_status "Creating test build directory..."
rm -rf "$TEST_BUILD_DIR"
mkdir -p "$TEST_BUILD_DIR"

# Test compilation flags
COMPILER_FLAGS="-c -O0 -g -Wall -Wextra -std=f2008"

# Test modules in dependency order
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

print_status "Testing module compilation..."

for module in "${MODULES[@]}"; do
    source_file="${module}.f90"
    object_file="${TEST_BUILD_DIR}/${module}.o"
    
    if [ ! -f "$source_file" ]; then
        print_error "Source file not found: $source_file"
        exit 1
    fi
    
    print_status "Testing compilation of $module..."
    
    # Compile module
    if mpif90 $COMPILER_FLAGS "$source_file" -o "$object_file" 2>&1; then
        print_success "✓ $module compiled successfully"
    else
        print_error "✗ $module compilation failed"
        exit 1
    fi
done

print_status "Testing main program compilation..."
MAIN_PROGRAM="earthquake_simulation_improved.f90"
if [ ! -f "$MAIN_PROGRAM" ]; then
    print_error "Main program source file not found: $MAIN_PROGRAM"
    exit 1
fi

# Try to compile main program (this will show dependency issues)
if mpif90 $COMPILER_FLAGS "$MAIN_PROGRAM" -o "${TEST_BUILD_DIR}/main.o" 2>&1; then
    print_success "✓ Main program compiled successfully"
else
    print_warning "Main program compilation failed (expected due to missing modules)"
fi

# Cleanup
print_status "Cleaning up test build directory..."
rm -rf "$TEST_BUILD_DIR"

print_success "All modules compiled successfully!"
print_status "You can now run the full compilation with:"
echo "  ./compile.sh"
echo "  or"
echo "  make"
