# cavOTF packaged workflow

This repository now ships the cavOTF workflow as an installable Python package with a single command-line entry point. The pipeline previously split between `run_first.py`, `run_second.py`, and `run.py` is executed end-to-end from a single configuration file.

## Installation

```bash
pip install .
```

The package installs a console script named `cavotf` and also supports `python -m cavotf`.

## Configuration

All user-configurable options live in one configuration file (default: `input.txt`). The file is parsed as simple `key=value` pairs with optional blank lines and `#` comments. Relative paths are resolved with respect to the configuration file location. Use `DEFAULT` to fall back to built-in resources for `clean_template_dir`, `client_sbatch`, or `server_sbatch`.

Example configuration:

```
# Paths
working_dir=.
data_source_dir=/scratch/user/u.sw216206/vsc-fluxside/project1/cavityDFTB/with_dmu/
clean_template_dir=DEFAULT
client_sbatch=DEFAULT
server_sbatch=DEFAULT
mu_submission_script=muDFTB.slurm

# Runtime
run_prefix=run-
server_wait_seconds=30
mu_results_filename=dmu.dat
mu_poll_interval=60
mu_timeout_seconds=0
thermalization_multiplier=50
initial_positions_name=thermaliz_water__InTheBox.dat
initial_velocities_name=thermaliz_water_vel__InTheBox.dat
initial_position_target=initXYZ.dat
initial_velocity_target=initPxPyPz.dat
```

Notes:
- `clean_template_dir` defaults to the packaged `cavotf/resources/DFTB_clean` directory. The templates are copied into each `run-*` folder before submitting jobs.
- `client_sbatch` and `server_sbatch` default to the packaged submission scripts. They are copied into the relevant working directories automatically.
- `mu_submission_script` **must** point to a valid Slurm script that produces the `mu_results_filename` (defaults to `dmu.dat`) inside each `run-*` directory. The workflow will raise a clear error if the script cannot be found.

## Usage

Run the complete workflow with a single command:

```bash
cavotf --config input.txt
# or
python -m cavotf --config input.txt
```

The command performs the following steps:
1. Cleans previous `run-*` folders and `*.txt` files in `working_dir`.
2. Copies the clean template, submission scripts, and starting coordinates from `data_source_dir` into `run-*` folders (one per `param.nk` from `funcLM.py`).
3. Submits the `mu_submission_script` job in each folder and waits for the `mu_results_filename` to appear (poll interval and optional timeout are configurable).
4. Generates `initial.dat` files via the thermalization loop originally in `run_second.py`.
5. Submits the server Slurm job followed by one client Slurm job per `run-*` folder after the configured wait.

Be mindful of HPC environments and file paths: the workflow runs commands inside each `run-*` directory using `sbatch`, so ensure the submission scripts are compatible with your cluster.
