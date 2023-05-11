#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --output=output.out
#SBATCH --error=output.err
#SBATCH --job-name=DHIT
#SBATCH --time=00:10:00

if which srun >/dev/null; then
  srun ../../build/a.out params8_DHIT_SMAGORINSKY_ROTATING_NONCUBIC
else
  mpirun -np 4 ../../build/a.out params8_DHIT_SMAGORINSKY_ROTATING_NONCUBIC
fi
