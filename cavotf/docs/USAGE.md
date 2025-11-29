# cavOTF usage guide

This guide explains how cavOTF parses configuration files, prepares run directories, and interacts with Slurm. Use it alongside `input.txt.example` when customizing your own runs.

## Command line interface

The CLI defined in `cavotf/cli.py` exposes two subcommands:

- `cavotf --config input.txt validate` performs a dry run. It checks geometry folders, renders sbatch scripts, and logs what would be submitted without writing files or launching jobs.
- `cavotf --config input.txt run` executes the full workflow. It creates run directories, writes multiprog configuration files, initializes cavity coordinates, and submits sbatch scripts through `sbatch`.

`python -m cavotf` is equivalent to invoking the `cavotf` console script.

## Configuration parsing

`cavotf.config` reads INI-style configuration files into dataclasses. All paths are resolved relative to the config file directory unless already absolute. The following sections are required:

### [general]
- `geometry_path`: directory containing one subfolder per trajectory seed. Each subfolder must include the geometry and velocity files specified below.
- `init_xyz_file`: filename (inside each geometry subfolder) holding atomic positions for the starting frame.
- `init_vel_file`: filename holding velocities for the same frame.
- `trajectory_prefix`: prefix used when generating trajectory names.
- `clean_template_dir` (optional): override path for the bundled `cavotf/resources/DFTB_clean` templates; defaults to the packaged copy.
- Thermostat settings (`use_thermostat`, `thermostat_type`, `thermostat_steps`, `collision_frequency`, `thermostat_reassign_particles`) mirror the defaults used in the packaged templates.

### [physics]
- `nk`: number of trajectories/clients to spawn; must not exceed the number of geometry subfolders present under `geometry_path`.
- `beta`, `lambda`, `omega_c`, `eta_b`: cavity parameters forwarded into the bundled parameter object.
- Thermostat options: same meanings as in `[general]`; can be overridden per-physics section.
- `calculate_dipole_derivatives` and `dipole_derivative_interval`: control whether dipole derivatives are computed and how frequently.
- `thermal_steps`: number of warmup steps used when initializing cavity coordinates.

### [outputs]
Booleans controlling which artifacts are produced during dynamics:
- `write_logfile`, `write_results`: enable logging and final results files.
- `record_k_space`, `print_k_space`: control whether k-space coordinates are recorded or echoed.
- `write_xyz_trajectory`, `write_histogram`, `write_output_client`, `write_midpoint_snapshots`: toggle client-side outputs written by the dynamics routines.

### [dftb]
Provides overrides for ASE Dftb calculator keywords. Values are parsed intelligently:
- Numeric strings are converted to integers or floats when possible.
- `off`, `none`, or `false` (case-insensitive) remove an option from the defaults.
- String values are cleaned to avoid stray newlines in generated input.

### [hpc]
Slurm submission settings used by `cavotf.hpc`:
- `cpus_per_job`, `partition`, `account`: scheduler resource requests.
- `run_get_mu`, `run_dynamics`: toggle submission of the dipole-collection and dynamics phases.
- `dftb_prefix`: SK path exported as `DFTB_PREFIX` for ASE/DFTB+.
- `dftb_command`: command template run by the server/client processes (defaults to `dftb+ > PREFIX.out`).
- `sbatch_template`: optional path to a custom sbatch template; falls back to the built-in template if omitted.
- `walltime`, `memory`: resource limits injected into the rendered script.

## Workflow steps

1. **Geometry preparation (`cavotf.geometry`)**
   - Discovers all subfolders under `geometry_path` and ensures at least `nk` exist.
   - Randomly selects `nk` cases and copies `init_xyz_file` and `init_vel_file` into `run-<index>` directories.

2. **DFTB and multiprog setup (`cavotf.dftb`)**
   - Builds configuration files for the `get_mu` and dynamics phases using parameters from `[dftb]`, `[outputs]`, and `[hpc]`.
   - Writes multiprog definitions consumed by `srun --multi-prog` in the sbatch script.

3. **Cavity initialization (`cavotf.dynamics`)**
   - Imports the parameter module bundled in the clean template directory and overrides key fields from `[physics]`.
   - Reads dipole values (`dmu.dat`) produced by the `get_mu` jobs, thermalizes cavity coordinates, and writes `initial.dat` to each run directory.

4. **Slurm integration (`cavotf.hpc`)**
   - Renders sbatch scripts using either the default template or a custom one. Context variables include job name, walltime, CPUs, memory, partition/account, DFTB paths, and the multiprog config.
   - Writes scripts to disk and submits them with `sbatch` unless the workflow is running in validation mode.

## Preparing geometry seeds

Each entry under `geometry_path` must look like:

```
geometry_path/
├─ case-000/
│  ├─ <init_xyz_file>
│  └─ <init_vel_file>
├─ case-001/
│  ├─ <init_xyz_file>
│  └─ <init_vel_file>
└─ ...
```

The names for `<init_xyz_file>` and `<init_vel_file>` come directly from the `[general]` section. If fewer than `nk` subfolders are present, the workflow raises an error during validation before any jobs are launched.

## Customizing sbatch scripts

If you need site-specific scheduler directives, point `sbatch_template` at a template file. Any placeholder surrounded by `{{...}}` in the template will be replaced with values from the workflow context. The default template sets `DFTB_PREFIX`, `DFTB_COMMAND`, and runs `srun --multi-prog` against the generated configuration file.

## Troubleshooting tips

- Use the `validate` command before `run` to confirm geometry counts, template resolution, and sbatch rendering.
- Ensure `sbatch` and `srun` are available in your environment; otherwise `run` will raise an error when trying to submit jobs.
- When overriding `clean_template_dir`, confirm the replacement directory contains the expected scripts (e.g., `funcLM.py`) required by `cavotf.dynamics.load_param`.