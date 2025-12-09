#!/bin/bash
#SBATCH --job-name=10box
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=48:00:00
#SBATCH -A tam116
#SBATCH --output=workflow_%j.out
#SBATCH --error=workflow_%j.err

# Load environment
module purge
source /expanse/lustre/projects/tam116/aamini1/mambaforge/etc/profile.d/conda.sh
conda activate dftb

# Path to geometry folders
GEOM_DIR="/home/aamini1/geometries"

# List of 40 geometries
GEOMS=(
    MDLM_2_10 MDLM_2_11 MDLM_2_14 MDLM_2_15 MDLM_2_16
    MDLM_2_18 MDLM_2_2  MDLM_2_20 MDLM_2_3  MDLM_2_4
    MDLM_2_7  MDLM_2_8  MDLM_2_9  MDLM_3_10 MDLM_3_14
    MDLM_3_18 MDLM_4_1  MDLM_4_10 MDLM_4_11 MDLM_4_12
    MDLM_4_15 MDLM_4_18 MDLM_4_19 MDLM_4_2  MDLM_4_20
    MDLM_4_4  MDLM_4_5  MDLM_4_8  MDLM_5_1  MDLM_5_11
    MDLM_5_15 MDLM_5_5  MDLM_5_8  MDLM_6_12 MDLM_6_13
    MDLM_6_14 MDLM_6_17 MDLM_6_2  MDLM_6_3  MDLM_6_7
)

# Enable debug mode so every command is echoed
set -x

echo "==============================="
echo " Workflow started at $(date)"
echo " Running on host: $(hostname)"
echo " SLURM_JOBID: $SLURM_JOB_ID"
echo "==============================="

# Loop through geometries
for geom in "${GEOMS[@]}"; do
    echo "---------------------------------"
    echo "=== Starting geometry $geom at $(date) ==="
    echo "---------------------------------"

    # Step 1: run_first with geometry
    echo ">>> Running run_first.py on $geom"
    python -u run_first.py "$geom"
    echo ">>> Finished run_first for $geom at $(date)"

    # Step 2: run_second
    echo ">>> Running run_second.py"
    python -u run_second.py
    echo ">>> Finished run_second at $(date)"

    # Step 3: run
    echo ">>> Running run.py"
    python -u run.py
    echo ">>> Finished run.py at $(date)"

    echo "=== Completed geometry $geom at $(date) ==="
done

echo "==============================="
echo " Workflow completed at $(date)"
echo "==============================="
