#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --output=output.out
#SBATCH --error=output.err
#SBATCH --cpus-per-task=4
#SBATCH --job-name=DHIT
#SBATCH --time=00:10:00

srun ../../build/a.out params1_DHIT_DNS
