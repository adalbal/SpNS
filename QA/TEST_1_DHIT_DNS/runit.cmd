#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --output=output.out
#SBATCH --error=output.err
#SBATCH --cpus-per-task=4
#SBATCH --job-name=DHIT
#SBATCH --time=00:10:00

if which srun >/dev/null; then
  srun ../../build/a.out params1_DHIT_DNS
else
  mpirun -np 4 ../../build/a.out params1_DHIT_DNS
fi
