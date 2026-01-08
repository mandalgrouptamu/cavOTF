# =============================================================================
#  Project:     cavOTF.py
#  File:        dftb.py
#  Author:      Amir H. Amini <amiramini@tamu.edu>
#  Last update: 11/28/2025
#
#  Description:
#      DFTB+ integration routines.
#      Note: I will consider renaming this script later to avoid conflict.
# =============================================================================
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, List

from .config import Config
from .resources import server_script_path

LOGGER = logging.getLogger(__name__)


def build_multiprog_commands(commands: Iterable[str]) -> List[str]:
    # Build multiprog commands with indices.
    return [f"{idx}  {cmd}" for idx, cmd in enumerate(commands)]


def write_conf(commands: Iterable[str], path: Path, dry_run: bool = False) -> Path:
    # Write the multiprog configuration file.
    lines = list(commands)
    LOGGER.debug("Writing multiprog configuration to %s", path)
    if not dry_run:
        path.write_text("\n".join(lines) + "\n")
    return path


def prepare_get_mu(config: Config, run_dirs: List[Path], dry_run: bool = False) -> Path:
    # Generate the get_mu multiprog file.
    base = config.path.parent
    script = config.general.clean_template_dir / "get_mu.py"
    config_arg = f"--config {config.path}"
    commands = [f"python {script} --workdir {run_dir} {config_arg}" for run_dir in run_dirs]
    prefixed = build_multiprog_commands(commands)
    return write_conf(prefixed, base / "get_mu.conf", dry_run=dry_run)


def prepare_run(config: Config, run_dirs: List[Path], dry_run: bool = False) -> Path:
    # Generate the main run multiprog file (server + clients).
    base = config.path.parent
    config_arg = f"--config {config.path}"
    server_script = server_script_path()
    commands = [f"python {server_script} {len(run_dirs)} {config_arg}"]
    client_script = config.general.clean_template_dir / "client_DFTB.py"
    commands.extend(
        f"python {client_script} {idx} --workdir {run_dir} --base {base} {config_arg}" for idx, run_dir in enumerate(run_dirs)
    )
    prefixed = build_multiprog_commands(commands)
    return write_conf(prefixed, base / "run.conf", dry_run=dry_run)

