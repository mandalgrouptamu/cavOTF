#!/bin/bash

#SBATCH --partition=cpu                                                                                                                            
#SBATCH --job-name=getdmu         #Set the job name to "JobExample1"                                                                                    
#SBATCH --time=0-0:05:00          #Set the wall clock limit to 1hr and 30min                                                                         
#SBATCH --ntasks=1               #Request 1 task                                                                                                      
#SBATCH --nodes=1                                                                                                                                     
#SBATCH --cpus-per-task=1                                                                                                                             
#SBATCH --mem=2000M                                                                                                                                  
#SBATCH --output=out.%j        

# Clean environment
ml purge
ml Anaconda3
source "$(conda info --base)"/etc/profile.d/conda.sh
conda activate sci-env
ml purge
ml Anaconda3
source "$(conda info --base)"/etc/profile.d/conda.sh
conda activate dftbplus
module purge


conda activate dftb


# Debug info
python -c "import sys, numpy, ase; print(sys.version); print(numpy.__version__); print(ase.__version__)"
which dftb+

python get_mu.py
