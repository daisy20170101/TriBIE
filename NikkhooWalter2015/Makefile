#===============================================================================
# Makefile for Nikkhoo-Walter Triangular Stiffness Calculation
# MPI/OpenMP hybrid parallel implementation
#===============================================================================

# Compiler settings
FC = mpif90
FFLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize
FFLAGS += -fopenmp
FFLAGS += -ffree-line-length-none -ffree-form
FFLAGS += -g -fbacktrace

# Debug build (uncomment for debugging)
# FFLAGS = -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -ffree-line-length-none -ffree-form -fopenmp

# Source files (in dependency order)
SOURCES = m_nikkhoo_green.f90 calc_nikkhoo.f90
OBJECTS = $(SOURCES:.f90=.o)
TARGET = calc_nikkhoo

# Module dependencies
calc_nikkhoo.o: m_nikkhoo_green.o

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
	rm -f nikkhoo_*.bin position.bin

# Run single process
run: $(TARGET)
	./$(TARGET)

# Run with MPI (4 processes by default)
run-mpi: $(TARGET)
	mpirun -np 4 ./$(TARGET)

# Run with custom MPI processes
run-mpi-%: $(TARGET)
	mpirun -np $* ./$(TARGET)

# Help
help:
	@echo "Nikkhoo-Walter Stiffness Calculation Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  all       - Build the executable (default)"
	@echo "  clean     - Remove object files, modules, and executable"
	@echo "  cleanall  - Remove all build and output files"
	@echo "  run       - Build and run with single process"
	@echo "  run-mpi   - Build and run with 4 MPI processes"
	@echo "  run-mpi-N - Build and run with N MPI processes (e.g., run-mpi-8)"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Environment variables:"
	@echo "  OMP_NUM_THREADS - Number of OpenMP threads per MPI process"
	@echo ""
	@echo "Example usage:"
	@echo "  make                          # Build"
	@echo "  make run                      # Run single process"
	@echo "  make run-mpi                  # Run with 4 MPI processes"
	@echo "  OMP_NUM_THREADS=4 make run-mpi-8  # 8 MPI x 4 OpenMP threads"
	@echo ""
	@echo "Input:  triangular_mesh.gts (GTS format mesh file)"
	@echo "Output: nikkhoo_<rank>.bin, position.bin"

.PHONY: all clean cleanall run run-mpi help
