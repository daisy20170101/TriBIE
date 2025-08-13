#!/bin/bash
# Enhanced compilation script with OpenMP and optimization flags
mpiifort -O3 -xHost -ipo -qopenmp -parallel \
         -D_OPENMP \
         -qopt-zmm-usage=high \
         -qopt-report=2 \
         sub_comdun.f mod_dtrigreen.f90 m_calc_green.f90 calc_trigreen.f90 \
         -o calc_trigreen_optimized

# Alternative compilation for systems without Intel compiler
# mpif90 -O3 -fopenmp -march=native -mtune=native \
#        -ffast-math -funroll-loops \
#        sub_comdun.f mod_dtrigreen.f90 m_calc_green.f90 calc_trigreen.f90 \
#        -o calc_trigreen_optimized

echo "Compilation completed. Executable: calc_trigreen_optimized"
echo "To run with OpenMP: export OMP_NUM_THREADS=<number_of_threads>"
echo "Example: mpirun -np 4 ./calc_trigreen_optimized"
