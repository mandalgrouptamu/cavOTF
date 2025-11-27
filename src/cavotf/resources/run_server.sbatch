#!/bin/bash
#SBATCH --job-name=server
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:15:00
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=server_output.txt


module load Anaconda3/2024.02-1  # if needed
echo "$(hostname)"
echo "$(hostname)" > server_hostname.txt
srun python3 server_DFTB.py $1
