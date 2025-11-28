# =============================================================================
#  Project:     cavOTF.py
#  File:        workflow.py
#  Author:      Amir H. Amini <amiramini@tamu.edu>
#  Last update: 11/28/2025
#
#  Description:
#      Management for sbatch script generation and job submission to HPC systems.
# =============================================================================

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict

from .config import Config

LOGGER = logging.getLogger(__name__)

DEFAULT_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={{job_name}}
#SBATCH --time={{walltime}}
#SBATCH --nodes=1
#SBATCH --ntasks={{ntasks}}
#SBATCH --cpus-per-task={{cpus}}
#SBATCH --mem={{memory}}
#SBATCH --partition={{partition}}
#SBATCH --output=out.%j
{{account_line}}

export DFTB_PREFIX={{dftb_prefix}}
export DFTB_COMMAND="{{dftb_command}}"

srun --multi-prog {{conf_path}}
"""


def _render_template(template: str, context: Dict[str, str]) -> str:
    rendered = template
    for key, value in context.items():
        rendered = rendered.replace(f"{{{{{key}}}}}", value)
    return rendered


def render_sbatch(config: Config, conf_path: Path, job_name: str, ntasks: int) -> str:
    # This routine generates an sbatch script text based on the provided 
    # configuration and template.
    template_path = config.hpc.sbatch_template
    if template_path:
        template_text = template_path.read_text()
    else:
        template_text = DEFAULT_TEMPLATE

    account_line = f"#SBATCH --account={config.hpc.account}" if config.hpc.account else ""
    context = {
        "job_name": job_name,
        "walltime": config.hpc.walltime,
        "ntasks": str(ntasks),
        "cpus": str(config.hpc.cpus_per_job),
        "memory": config.hpc.memory,
        "partition": config.hpc.partition or "compute",
        "account_line": account_line,
        "dftb_prefix": str(config.hpc.dftb_prefix or ""),
        "dftb_command": config.hpc.dftb_command or "dftb+ > PREFIX.out",
        "conf_path": str(conf_path),
    }
    return _render_template(template_text, context)


def write_sbatch(script_text: str, path: Path, dry_run: bool = False) -> Path:
    LOGGER.info("Writing sbatch script to %s", path)
    if not dry_run:
        path.write_text(script_text)
        os.chmod(path, 0o755)
    return path


def submit_job(script_path: Path, dry_run: bool = False) -> None:
    # This routine submits the sbatch script to the job scheduler.
    cmd = ["sbatch", str(script_path)]
    LOGGER.info("Submitting job: %s", " ".join(cmd))
    if dry_run:
        return
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:  # noqa: BLE001
        raise RuntimeError("sbatch command not found in PATH") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"sbatch submission failed: {exc}") from exc

