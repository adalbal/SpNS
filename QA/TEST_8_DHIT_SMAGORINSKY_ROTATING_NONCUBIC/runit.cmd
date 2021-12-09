#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --output=output.out
#SBATCH --error=output.err
#SBATCH --job-name=DHIT
#SBATCH --time=00:10:00

srun ../../build/a.out params8_DHIT_SMAGORINSKY_ROTATING_NONCUBIC