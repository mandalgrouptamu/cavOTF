#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=cavityDFTB    # Set the job name to "md"
#SBATCH --time=2-00:00:00        # Set the wall clock limit to 10 minutes
#SBATCH --ntasks=1               # Request 1 task
#SBATCH --mem=1000M              # Request 1000MB of memory
#SBATCH --output=MDout.%j        # Output file (the %j will be replaced by the job ID)

# Activate the Python virtual environment
source /scratch/user/u.sw216206/sachithpy3/bin/activate

# Run the Python script
python dftbLM.py
