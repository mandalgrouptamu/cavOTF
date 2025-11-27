# cavOTF packaged workflow

This repository now ships the cavOTF workflow as an installable Python package with a single command-line entry point. The pipeline previously split between `run_first.py`, `run_second.py`, and `run.py` is executed end-to-end from a single configuration file.

## Installation

```bash
pip install .
```

The package installs a console script named `cavotf` and also supports `python -m cavotf`.

## Configuration

All user-configurable options live in one INI-style configuration file (default: `input.txt`). Relative paths are resolved with respect to the configuration file location. Use `DEFAULT` to fall back to built-in resources for `clean_template_dir`, `client_sbatch`, or `server_sbatch`.

Example configuration showing every configurable field:

```
[general]
n_trajectories = 4
geometry_path = /scratch/user/u.aa271894/VSC/geometries/
working_dir = .
clean_template_dir = DEFAULT
client_sbatch = DEFAULT
server_sbatch = DEFAULT
trajectory_prefix = trajectory_
init_xyz_file = thermaliz_water__InTheBox.dat
init_vel_file = thermaliz_water_vel__InTheBox.dat
initial_position_target = initXYZ.dat
initial_velocity_target = initPxPyPz.dat
mu_results_filename = dmu.dat
mu_poll_interval = 60
mu_timeout_seconds = 0
thermalization_multiplier = 50
use_thermostat = yes
thermostat_type = andersen
thermostat_steps = 250
collision_frequency = 0.001
thermostat_reassign_particles = 16

[geometry]
geom_root_head = /scratch/user/u.aa271894/VSC/geometries/
geom_root      = /scratch/user/u.aa271894/VSC/geometries/
use_head_root  = true

[physics]
use_param_nk   = true       ; if true, N = param.nk, else use n_trajectories
nk             = 25
beta           = 1052.8
lambda         = 0.001
omega_c        = 0.190/27.2114
eta_b          = 0.0003
use_thermostat = true
thermostat_type = andersen
thermostat_steps = 250
thermostat_reassign_particles = 16

[hpc]
cpus_per_job   = 8
partition      = compute
account        = myaccount
mu_submission_script = muDFTB.slurm
server_wait_seconds = 30
run_get_mu     = true
run_dynamics   = true
dftb_prefix = /scratch/user/u.sw216206/dftb_sk_files/mio-1-1
dftb_command = /scratch/user/u.aa271894/.conda/envs/dftbplus/bin/dftb+ > PREFIX.out
```

Notes:
- `clean_template_dir` defaults to the packaged `cavotf/resources/DFTB_clean` directory. The templates (including the placeholder `get_mu.py`) are copied into each `trajectory_*` folder before submitting jobs.
- `client_sbatch` and `server_sbatch` default to the packaged submission scripts. They are copied into the relevant working directories automatically.
- `mu_submission_script` **must** point to a valid Slurm script that produces the `mu_results_filename` (defaults to `dmu.dat`) inside each run directory. The workflow will raise a clear error if the script cannot be found. A minimal template `muDFTB.slurm` is shipped in `cavotf/resources/DFTB_clean` that simply calls `get_mu.py` to write a placeholder `dmu.dat`; replace its contents with your actual μ job.
- `run_get_mu` toggles whether the μ jobs are submitted during staging; `run_dynamics` toggles whether the final server/clients are launched.
- `cpus_per_job`, `partition`, and `account` are passed directly to `sbatch`; `DFTB_PREFIX` and `DFTB_COMMAND` are exported to the job environment if provided.
- The default Slurm templates no longer try to activate a hard-coded conda environment. If you need to activate something, set an environment variable like `CAVOTF_ACTIVATE=/path/to/activate_script.sh` before launching `cavotf` (each script checks for it before activation).

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
