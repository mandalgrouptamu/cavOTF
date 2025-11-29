
# cavOTF

**cavOTF.py** is a highly customizable molecular dynamics package based on DFTB+ (Density Functional Tight Binding) for simulating cavity systems. Due to the computational intensity and long runtimes, this package is designed for use on High Performance Computing (HPC) systems.

## Installation

Before installing **cavOTF.py**, ensure that you have a clean DFTB+ installation in your environment. For instructions, refer to the [DFTB+ Installation Guide](https://dftbplus-recipes.readthedocs.io/en/latest/introduction.html). We recommend installing DFTB+ using the **Anaconda** package manager within a separate environment dedicated solely to DFTB+. Avoid using pre-compiled DFTB+ binaries as they may significantly reduce the performance of your calculations. **cavOTF** has been tested with the **no MPI** version of DFTB+, so it is advisable to use this version, as `numpy` and/or `scipy` may have compatibility issues with some MPI installations.

The package uses the **setuptools** build system and includes most required Python modules and resource templates. Ensure you have an up-to-date version of `setuptools` before installation.

To install **cavOTF**, clone the repository and run:

```bash
pip install .
```

Alternatively, you can install directly from the repository:

```bash
pip install git+https://<repo-url>
```

### Available Entry Points

Once installed, you can use the following commands:

- `cavotf --config input.txt validate` — Perform a dry-run validation without submitting jobs.
- `cavotf --config input.txt run` — Execute the full workflow, including sbatch job submissions.
- `python -m cavotf --config input.txt ...` — Invoke the module programmatically, equivalent to the CLI.

## Quickstart

1. Copy the `input.txt.example` file and modify the paths and parameters according to your system. Detailed descriptions of input options are provided within the file.
2. Ensure your geometry directory contains one subfolder per trajectory seed, with corresponding coordinate (`init_xyz_file`) and velocity (`init_vel_file`) files.
3. Run the following command to validate the sbatch scripts and ensure that the run directories will be created correctly:
   ```bash
   cavotf --config input.txt validate
   ```
4. Execute the full workflow with:
   ```bash
   cavotf --config input.txt run
   ```
   This will create the necessary run folders, initialize the cavity coordinates, and submit the required jobs through Slurm.
5. Optionally, modify the provided sbatch template to fit your cluster’s requirements. We recommend using the maximum allowable wall time to ensure the simulations run smoothly.

<!-- During execution, **cavOTF** will:

- Parse the configuration file into strongly-typed dataclasses (located in `cavotf.config`).
- Prepare the run directories and copy the seed geometries (handled by `cavotf.geometry`).
- Build the DFTB+ input files and multiprog definitions for dipole collection and dynamics (using `cavotf.dftb`).
- Render Slurm sbatch scripts, incorporating any custom template provided (`cavotf.hpc`).
- Initialize cavity coordinates based on generated dipoles before launching dynamics simulations (`cavotf.dynamics`). -->

## Configuration Overview

The `input.txt.example` file provides an example configuration, including all supported keys. Key sections include:

- **[general]** — Specifies the geometry path, initial coordinate/velocity filenames and optional overrides for bundled `DFTB_clean` templates. Modifications to the `DFTB_clean` folder should be done cautiously.
- **[physics]** — Defines the number of clients (`nk`), cavity parameters (`beta`, `lambda`, `omega_c`, `eta_b`), thermostat controls, dipole derivative options, and thermalization steps for cavity initialization.
- **[outputs]** — Controls for logfile and results generation, k-space tracking, and client-side outputs such as `qt.dat` and trajectory histograms.
- **[dftb]** — Override ASE DFTB+ calculator settings. Set values to `off` or `none` to omit default lines. The provided defaults are configured for a 99-atom water simulation, so modify these as needed.
- **[hpc]** — Slurm integration details, including CPUs per task, partition/account, DFTB+ prefix/command, sbatch template path, walltime, and memory requirements. Ensure sufficient memory allocation to prevent out-of-memory (OOM) issues.

<!-- Refer to `docs/USAGE.md` for a more detailed walkthrough of how these settings are parsed and utilized throughout the workflow. -->

<!-- ## Bundled Resources

The package includes the default `DFTB_clean` templates and `server_DFTB.py`, making these resources available even when running from an installed package. You can override the template directory via the `clean_template_dir` option in the `[general]` section if custom DFTB+ input files are required. -->

---

**Authors:**  
Sachith Wickramasinghe, Amir H. Amini, Arkajit Mandal  
November 2025, Chemistry Department, Texas A&M University
