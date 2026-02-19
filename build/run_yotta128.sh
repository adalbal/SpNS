#!/bin/bash

#SBATCH --job-name=yotta128                     # Job name
#SBATCH --mail-type=ALL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adel.alsalti@upc.edu        # Where to send mail	
#SBATCH --account=upc61
#SBATCH --qos=gp_resa
#SBATCH --partition=gpp
#SBATCH --nodes=1                               # Number of nodes
##SBATCH --ntasks=56                            # Number of MPI ranks
#SBATCH --cpus-per-task=1                       # Number of cores per MPI rank 
#SBATCH --ntasks-per-node=112                   # How many tasks on each node
#SBATCH --ntasks-per-socket=56                  # How many tasks on each CPU or socket
#SBATCH --time=72:00:00                         # Time limit hrs:min:sec
#SBATCH --output=yotta128/log_%j.out            # Standard output log
#SBATCH --error=yotta128/log_%j.err             # Standard error log

pwd; hostname; date

module load fftw

export OMP_NUM_THREADS=1;
cat params_yotta128nu
srun ./a.out params_yotta128nu

date
