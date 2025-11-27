#!/bin/bash
#SBATCH --job-name=muDFTB
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=500M
#SBATCH --output=mu.%j.out

set -euo pipefail

# This template simply writes a placeholder mu value. Replace the body with the
# real workflow your cluster uses to produce the dipole moment file.

# Activate any environment you need here, for example:
# [ -n "$CAVOTF_ACTIVATE" ] && source "$CAVOTF_ACTIVATE"

python get_mu.py
