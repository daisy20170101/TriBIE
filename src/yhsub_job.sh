#!/bin/bash
#SBATCH --output=job-%j.out
#SBATCH --error=error-%j.out
#SBATCH --nodes=5
#SBATCH --ntasks=120
#SBATCH --partition=TH_LONG
yhrun ./calc_oblique_ns
