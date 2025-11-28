"""Orchestration of the CAVOTF workflow."""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Optional

from .config import Config
from .dftb import prepare_get_mu, prepare_run
from .dynamics import initialize_cavity
from .geometry import prepare_run_directories
from .hpc import render_sbatch, submit_job, write_sbatch

LOGGER = logging.getLogger(__name__)


def validate_workflow(config: Config) -> None:
    """Perform a dry-run validation and log planned actions."""
    LOGGER.info("Validating workflow from %s", config.path)
    run_dirs = prepare_run_directories(config, dry_run=True)
    get_mu_conf = prepare_get_mu(config, run_dirs, dry_run=True)
    run_conf = prepare_run(config, run_dirs, dry_run=True)
    LOGGER.info("Would create %s run directories", len(run_dirs))
    LOGGER.info("get_mu multiprog: %s", get_mu_conf)
    LOGGER.info("run multiprog: %s", run_conf)

    if config.hpc.run_get_mu:
        script_text = render_sbatch(config, Path("get_mu.conf"), job_name="getMU", ntasks=len(run_dirs))
        LOGGER.info("Would submit get_mu sbatch script:\n%s", script_text)
    if config.hpc.run_dynamics:
        script_text = render_sbatch(config, Path("run.conf"), job_name="Run", ntasks=len(run_dirs) + 1)
        LOGGER.info("Would submit dynamics sbatch script:\n%s", script_text)


def run_workflow(config: Config) -> None:
    """Execute the workflow following the legacy scripts."""
    LOGGER.info("Running workflow from %s", config.path)
    run_dirs = prepare_run_directories(config, dry_run=False)

    get_mu_conf = prepare_get_mu(config, run_dirs, dry_run=False)
    if config.hpc.run_get_mu:
        script_text = render_sbatch(config, get_mu_conf, job_name="getMU", ntasks=len(run_dirs))
        script_path = write_sbatch(script_text, config.path.parent / "get_mu.sh")
        submit_job(script_path, dry_run=False)
        _wait_for_dipoles(run_dirs)

    # Initialize cavity coordinates using dmu.dat produced by get_mu
    initialize_cavity(config, run_dirs, dry_run=False)

    run_conf = prepare_run(config, run_dirs, dry_run=False)
    if config.hpc.run_dynamics:
        script_text = render_sbatch(config, run_conf, job_name="Run", ntasks=len(run_dirs) + 1)
        script_path = write_sbatch(script_text, config.path.parent / "run.sh")
        submit_job(script_path, dry_run=False)


def _wait_for_dipoles(run_dirs: list[Path], poll_interval: int = 30, timeout: Optional[int] = None) -> None:
    """Block until every ``dmu.dat`` file exists in the provided run directories."""

    missing = _missing_dipoles(run_dirs)
    if not missing:
        return

    LOGGER.info("Waiting for %s dmu.dat files to be produced", len(missing))
    start = time.monotonic()
    while missing:
        if timeout is not None and (time.monotonic() - start) > timeout:
            missing_list = ", ".join(str(path) for path in sorted(missing))
            raise TimeoutError(f"Timed out waiting for dipole files: {missing_list}")

        time.sleep(poll_interval)
        missing = _missing_dipoles(run_dirs)

    LOGGER.info("All dmu.dat files detected; proceeding to cavity initialization")


def _missing_dipoles(run_dirs: list[Path]) -> set[Path]:
    return {run_dir / "dmu.dat" for run_dir in run_dirs if not (run_dir / "dmu.dat").is_file()}

