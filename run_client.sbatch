#!/bin/bash
#SBATCH --job-name=cmd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:20:00
#SBATCH --mem=4G  
#SBATCH --output=cmd_%j.txt
## array job    
module purge
# Activate the Python virtual environment
source /scratch/user/u.sw216206/sachithpy3/bin/activate

python -c "import sys; import numpy; print(sys.version); print(numpy.__version__); print(numpy.__file__)"


module load GCC/11.3.0 OpenMPI/4.1.4 DFTB+/24.1
# give the server a few seconds to start listening
sleep 10

/sw/eb/sw/Python/3.10.4-GCCcore-11.3.0/bin/python client_DFTB.py $1
