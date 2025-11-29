# cavOTF

cavOTF orchestrates cavity molecular dynamics workflows that couple ASE/DFTB calculations with cavity initialization and Slurm job management. The package exposes both a Python API and a CLI entry point so users can validate inputs, render sbatch scripts, and launch production runs with consistent configuration handling.

## Installation

The project uses a standard setuptools build backend and ships all required Python modules and resource templates. Install from a local checkout, a wheel, or a Git URL:

```bash
pip install .
pip install cavotf                 # from PyPI or a built wheel
pip install git+https://<repo-url>
```

After installation, the following entry points are available:

- `cavotf --config input.txt validate` — dry-run validation without submitting jobs.
- `cavotf --config input.txt run` — full workflow execution, including sbatch submission.
- `python -m cavotf --config input.txt ...` — module invocation equivalent to the CLI.

## Quickstart

1. Copy `input.txt.example` and adjust the paths and parameters for your system.
2. Make sure your geometry directory contains one subfolder per trajectory seed with matching coordinate (`init_xyz_file`) and velocity (`init_vel_file`) files.
3. Run `cavotf --config input.txt validate` to confirm sbatch scripts and run directories will be generated correctly.
4. Run `cavotf --config input.txt run` to create run folders, initialize cavity coordinates, and submit get_mu and dynamics jobs through Slurm.

During execution, cavOTF will:

- Read configuration sections into strongly typed dataclasses (see `cavotf.config`).
- Prepare run directories and copy seed geometries (`cavotf.geometry`).
- Build DFTB input and multiprog definitions for dipole collection and dynamics (`cavotf.dftb`).
- Render sbatch scripts with any custom template provided (`cavotf.hpc`).
- Initialize cavity coordinates from generated dipoles before launching dynamics (`cavotf.dynamics`).

## Configuration overview

Configuration files use INI syntax. The example file documents all supported keys. Key sections include:

- **[general]** — geometry path, initial coordinate/velocity filenames, trajectory prefix, and optional override for bundled `DFTB_clean` templates.
- **[physics]** — number of trajectories (`nk`), cavity parameters (`beta`, `lambda`, `omega_c`, `eta_b`), thermostat controls, dipole derivative options, and thermalization steps for cavity initialization.
- **[outputs]** — controls for logfile/results generation, k-space tracking, and client-side artifacts such as `initial.dat` and trajectory histograms.
- **[dftb]** — overrides for ASE Dftb calculator keywords; set values to `off`/`none` to drop a default line entirely.
- **[hpc]** — Slurm integration: CPUs per task, partition/account, DFTB prefix/command, optional sbatch template path, walltime, and memory.

See `docs/USAGE.md` for a deeper walkthrough of how these settings are parsed and used across the workflow.

## Bundled resources

The package ships `cavotf/resources/DFTB_clean` and `cavotf/server_DFTB.py` inside the wheel so default templates and the server script are available even when running from an installed package. You can override the template directory via `clean_template_dir` in `[general]` if you need to supply custom DFTB inputs.

## Development notes

- The CLI banner and option parsing live in `cavotf/cli.py`.
- High-level workflow validation and execution are defined in `cavotf.workflow`.
- sbatch rendering and submission helpers are in `cavotf.hpc`.
- Geometry preparation and cavity initialization routines are in `cavotf.geometry` and `cavotf.dynamics`.

Contributions should avoid changing simulation algorithms unless intentionally extending scientific behavior; most user-facing tweaks can be made through configuration and resource templates.
