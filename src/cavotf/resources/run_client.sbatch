#!/bin/bash
#SBATCH --job-name=cmd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:20:00
#SBATCH --mem=4G  
#SBATCH --output=cmd_%j.txt
## array job    
# Optionally activate your environment, e.g.:
# [ -n "$CAVOTF_ACTIVATE" ] && source "$CAVOTF_ACTIVATE"

# give the server a few seconds to start listening
sleep 10

python client_DFTB.py $1
