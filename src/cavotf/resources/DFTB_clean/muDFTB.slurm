#!/bin/bash
#SBATCH --job-name=muDFTB
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=500M
#SBATCH --output=mu.%j.out

set -euo pipefail

# This template simply writes a placeholder mu value. Replace the body with the
# real workflow your cluster uses to produce the dipole moment file.
python - <<'PY'
import numpy as np
from pathlib import Path

outfile = Path("dmu.dat")
if not outfile.exists():
    np.savetxt(outfile, np.array([0.0]))
print(f"Wrote placeholder mu to {outfile.resolve()}")
PY
