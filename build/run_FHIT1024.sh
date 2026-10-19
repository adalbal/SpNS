#!/bin/bash

#SBATCH --job-name=FHIT-1024                    # Job name
#SBATCH --nodes=60                              # Number of nodes
#SBATCH --cpus-per-task=1                       # Number of cores per MPI rank 
#SBATCH --ntasks-per-node=112                   # How many tasks on each node
#SBATCH --ntasks-per-socket=56                  # How many tasks on each CPU or socket
#SBATCH --time=3-00:00:00                       # Time limit hrs:min:sec
#SBATCH --output=FHIT1024/log_%j.out            # Standard output log
#SBATCH --error=FHIT1024/log_%j.err             # Standard error log

pwd; hostname; date

module load fftw

export OMP_NUM_THREADS=1;
cat params_FHIT1024
srun ./a.out params_FHIT1024

date
