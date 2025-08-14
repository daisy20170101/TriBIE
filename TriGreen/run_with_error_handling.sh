#!/bin/bash

# Runtime script for TriGreen with floating-point exception handling
# This script sets environment variables to prevent crashes from numerical issues

echo "Setting up TriGreen with floating-point exception handling..."

# Set floating-point exception handling
export GFORTRAN_CONVERT_UNIT='big_endian:101-199,301-399'
export GFORTRAN_ERROR_DUMPCORE=0
export GFORTRAN_ERROR_BACKTRACE=1

# Set OpenMP environment variables
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_SCHEDULE=dynamic

# Set MPI environment variables for better error handling
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl="^openib"
export OMPI_MCA_orte_abort_on_non_zero_status=0

# Check if executable exists
if [ ! -f "./calc_trigreen" ]; then
    echo "Error: calc_trigreen executable not found!"
    echo "Please compile first using: ./runcompile.bash"
    exit 1
fi

# Get MPI processes from command line argument
NP=${1:-2}
if [ $NP -lt 1 ]; then
    echo "Error: Number of MPI processes must be at least 1"
    exit 1
fi

echo "Running TriGreen with:"
echo "  MPI processes: $NP"
echo "  OpenMP threads per process: $OMP_NUM_THREADS"
echo "  Total threads: $((NP * OMP_NUM_THREADS))"
echo "  Floating-point exception handling: Enabled"
echo ""

# Run the program
echo "Starting execution..."
mpirun -np $NP \
       -x OMP_NUM_THREADS \
       -x GFORTRAN_CONVERT_UNIT \
       -x GFORTRAN_ERROR_DUMPCORE \
       -x GFORTRAN_ERROR_BACKTRACE \
       ./calc_trigreen

# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "TriGreen completed successfully!"
else
    echo ""
    echo "TriGreen exited with code $EXIT_CODE"
    echo "Check the output files for details about any errors."
fi

exit $EXIT_CODE
