#!/bin/bash
#SBATCH --job-name=trigreen_improved
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=trigreen_%j.out
#SBATCH --error=trigreen_%j.err
#SBATCH --partition=compute
#SBATCH --exclusive

# Load necessary modules (adjust for your cluster)
module load openmpi
module load gcc

# Set OpenMP threads per MPI process
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Performance tuning
export OMP_SCHEDULE=dynamic
export OMP_DYNAMIC=TRUE
export OMP_NESTED=TRUE

# MPI tuning
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"
export OMPI_MCA_btl="^openib"

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Running on $SLURM_NNODES nodes"
echo "MPI processes per node: $SLURM_NTASKS_PER_NODE"
echo "OpenMP threads per MPI process: $OMP_NUM_THREADS"
echo "Total parallel threads: $((SLURM_NNODES * SLURM_NTASKS_PER_NODE * OMP_NUM_THREADS))"
echo "Node list: $SLURM_NODELIST"

# Create output directory
mkdir -p output_${SLURM_JOB_ID}

# Run the improved executable
echo "Starting TriBIE calculation with improved parallelization..."
mpirun --bind-to core:overload-allowed \
       --map-by ppr:${SLURM_NTASKS_PER_NODE}:node \
       -np $((SLURM_NNODES * SLURM_NTASKS_PER_NODE)) \
       ./3dtriBIE_improved

# Check exit status
if [ $? -eq 0 ]; then
    echo "Job completed successfully!"
    echo "Results saved in output_${SLURM_JOB_ID}/"
else
    echo "Job failed with exit code $?"
    exit 1
fi
