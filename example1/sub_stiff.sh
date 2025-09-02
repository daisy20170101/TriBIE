#!/bin/bash

#SBATCH -J calc_trigreen
#SBATCH --time=50:00
#SBATCH --mem-per-cpu=15G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --output=slog_%j.out
#SBATCH --error=slog_%j.err


# User specific aliases and functions


#srun --mpi=pmix /home/DuoL/seissol_bin/bin/SeisSol_Release_dskx_3_elastic parameters.par
export OMP_NUM_THREADS=12
mpirun -np 1 /scratch/duol/GitHub/TriBIE/TriGreen/calc_trigreen



