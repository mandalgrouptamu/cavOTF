#!/bin/bash

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=cavityDFTB    # Set the job name to "md"
#SBATCH --time=7-00:00:00        # Set the wall clock limit to 10 minutes
#SBATCH --ntasks=1               # Request 1 task
#SBATCH --mem=1000M              # Request 1000MB of memory
#SBATCH --output=MDout.%j        # Output file (the %j will be replaced by the job ID)

# Activate any environment you need here, for example:
# [ -n "$CAVOTF_ACTIVATE" ] && source "$CAVOTF_ACTIVATE"

# Run the Python script
python dftbLM.py
