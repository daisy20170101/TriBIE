#!/bin/bash

#SBATCH -J 3dtri
#SBATCH --time=50:00
#SBATCH --mem-per-cpu=15G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --output=slog_%j.out
#SBATCH --error=slog_%j.err


# User specific aliases and functions


#srun --mpi=pmix /home/DuoL/seissol_bin/bin/SeisSol_Release_dskx_3_elastic parameters.par
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close
export OMP_PLACES=cores

export OMP_DISPLAY_ENV=TRUE
export OMP_DISPLAY_AFFINITY=True

mpirun -np 1 /scratch/duol/GitHub/TriBIE_imp5/src/3dtri_BP5



