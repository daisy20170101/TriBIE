#===============================================================================
# Makefile for Nikkhoo-Walter Triangular Stiffness Calculation
# MPI/OpenMP hybrid parallel implementation
#
# Usage:
#   make                    # Build with mpif90 (default, MPI+OpenMP)
#   make FC=gfortran        # Build with gfortran (OpenMP only, no MPI)
#   make gfortran           # Same as above (convenience target)
#===============================================================================

# Compiler settings (can be overridden: make FC=gfortran)
FC ?= mpif90
FFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize
FFLAGS += -fopenmp
FFLAGS += -ffree-line-length-none -ffree-form
FFLAGS += -g -fbacktrace

# Debug build (uncomment for debugging)
# FFLAGS = -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -ffree-line-length-none -ffree-form -fopenmp

# Source files (in dependency order)
SOURCES = m_nikkhoo_green.f90 sub_nikkhoo.f90 calc_nikkhoo.f90
OBJECTS = $(SOURCES:.f90=.o)
TARGET = calc_nikkhoo

# Module dependencies
calc_nikkhoo.o: m_nikkhoo_green.o sub_nikkhoo.o

#===============================================================================
# Targets
#===============================================================================

# Default target
all: $(TARGET)

# Link executable
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# Compile Fortran source files
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET) *.mod

# Clean all (including output files)
cleanall: clean
	rm -f trigreen_*.bin position.bin

# Run single process
run: $(TARGET)
	./$(TARGET)

# Run with MPI (4 processes by default)
run-mpi: $(TARGET)
	mpirun -np 4 ./$(TARGET)

# Run with custom MPI processes
run-mpi-%: $(TARGET)
	mpirun -np $* ./$(TARGET)

# Build with gfortran (convenience target)
gfortran: clean
	$(MAKE) FC=gfortran all

# Help
help:
	@echo "Nikkhoo-Walter Stiffness Calculation Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  all       - Build the executable (default, uses mpif90)"
	@echo "  gfortran  - Build with gfortran (OpenMP only, no MPI)"
	@echo "  clean     - Remove object files, modules, and executable"
	@echo "  cleanall  - Remove all build and output files"
	@echo "  run       - Build and run with single process"
	@echo "  run-mpi   - Build and run with 4 MPI processes"
	@echo "  run-mpi-N - Build and run with N MPI processes (e.g., run-mpi-8)"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Compiler selection:"
	@echo "  make                 # Use mpif90 (default)"
	@echo "  make FC=gfortran    # Use gfortran"
	@echo "  make gfortran       # Same as above (convenience)"
	@echo ""
	@echo "Environment variables:"
	@echo "  OMP_NUM_THREADS - Number of OpenMP threads per MPI process"
	@echo ""
	@echo "Example usage:"
	@echo "  make                          # Build with MPI"
	@echo "  make gfortran                 # Build without MPI"
	@echo "  make run                      # Run single process"
	@echo "  make run-mpi                  # Run with 4 MPI processes"
	@echo "  OMP_NUM_THREADS=4 make run-mpi-8  # 8 MPI x 4 OpenMP threads"
	@echo ""
	@echo "Input:  triangular_mesh.gts (GTS format mesh file)"
	@echo "Output: trigreen_<rank>.bin, position.bin"

.PHONY: all clean cleanall run run-mpi gfortran help
