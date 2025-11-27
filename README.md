# cavOTF packaged workflow

This repository now ships the cavOTF workflow as an installable Python package with a single command-line entry point. The pipeline previously split between `run_first.py`, `run_second.py`, and `run.py` is executed end-to-end from a single configuration file.

## Installation

```bash
pip install .
```

The package installs a console script named `cavotf` and also supports `python -m cavotf`.

## Configuration

All user-configurable options live in one configuration file (default: `input.txt`). The file is parsed as an INI with sections. Relative paths are resolved with respect to the configuration file location. Use `DEFAULT` to fall back to built-in resources for `clean_template_dir`, `client_sbatch`, or `server_sbatch`.

Example configuration (with the requested section layout):

```
[general]
n_trajectories = 4
working_dir = .
geometry_path = /scratch/user/u.aa271894/VSC/geometries/
init_xyz_file = thermaliz_water__InTheBox.dat
init_vel_file = thermaliz_water_vel__InTheBox.dat
trajectory_prefix = trajectory_
use_thermostat = yes
thermostat_type = andersen
thermostat_steps = 250
collision_frequency = 0.001
thermostat_reassign_particles = 16
clean_template_dir = DEFAULT
client_sbatch = DEFAULT
server_sbatch = DEFAULT
mu_submission_script = muDFTB.slurm
mu_results_filename = dmu.dat
server_wait_seconds = 30
mu_poll_interval = 60
mu_timeout_seconds = 0
thermalization_multiplier = 50
initial_position_target = initXYZ.dat
initial_velocity_target = initPxPyPz.dat

[geometry]
geom_root_head = /scratch/user/u.aa271894/VSC/geometries/
geom_root      = /scratch/user/u.aa271894/VSC/geometries/
use_head_root  = true       ; keep your location_heaed vs location logic

[physics]
# These mirror funcLM.param defaults; we can override later if needed
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
run_get_mu     = true
run_dynamics   = true
dftb_prefix = /scratch/user/u.sw216206/dftb_sk_files/mio-1-1
dftb_command = /scratch/user/u.aa271894/.conda/envs/dftbplus/bin/dftb+ > PREFIX.out
```

Notes:
- `clean_template_dir` defaults to the packaged `cavotf/resources/DFTB_clean` directory. The templates are copied into each `trajectory_*` folder before submitting jobs.
- `client_sbatch` and `server_sbatch` default to the packaged submission scripts. They are copied into the relevant working directories automatically and can be overridden per input file.
- `mu_submission_script` **must** point to a valid Slurm script that produces the `mu_results_filename` (defaults to `dmu.dat`) inside each `trajectory_*` directory. The workflow will raise a clear error if the script cannot be found.
- HPC options feed directly into `sbatch` invocations: `partition`, `account`, and `cpus_per_job` become CLI options, while `dftb_prefix` and `dftb_command` are exported as environment variables.

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
