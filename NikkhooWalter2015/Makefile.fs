#===============================================================================
# Makefile.fs -- builds calc_nikkhoo_fs (full-space stiffness driver, BP8)
#
# Usage:
#   make -f Makefile.fs                # Build with mpif90 (default)
#   make -f Makefile.fs run-mpi-N       # Build and run with N MPI processes
#===============================================================================

FC = mpif90
FFLAGS = -O3 -fopenmp -ffree-line-length-none -ffree-form -g -fbacktrace

SOURCES = sub_nikkhoo.f90 calc_nikkhoo_fs.f90
OBJECTS = $(SOURCES:.f90=.o)
TARGET = calc_nikkhoo_fs

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

# Ensures sub_nikkhoo.o is built before calc_nikkhoo_fs.o (module dependency).
# Placed after 'all' so it doesn't silently become the default goal (GNU
# Make picks the first target in the file as the default when none is
# given on the command line -- a prerequisite-only rule still counts).
calc_nikkhoo_fs.o: sub_nikkhoo.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) *.mod

cleanall: clean
	rm -f trigreen_*.bin position.bin

run-mpi-%: $(TARGET)
	mpirun -np $* ./$(TARGET)

help:
	@echo "calc_nikkhoo_fs: full-space (whole-space) stiffness driver for BP8"
	@echo "  make -f Makefile.fs               - build"
	@echo "  make -f Makefile.fs run-mpi-N      - build and run with N MPI ranks"
	@echo "Input:  triangular_mesh.gts (flat mesh, see example4/make_bp8_mesh.py)"
	@echo "Output: trigreen_{22,23,32,33}_<rank>.bin, position.bin"

.PHONY: all clean cleanall help
