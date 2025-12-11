# =============================================================================
#  Project:     cavOTF.py
#  File:        geometry.py
#  Author:      Amir H. Amini <amiramini@tamu.edu>
#  Last update: 11/28/2025
#
#  Description:
#      Geometry preparation utilities for cavOTF.py simulations.
# =============================================================================
from __future__ import annotations

import logging
import random
import shutil
from pathlib import Path
from typing import Iterable, List

from .config import Config

LOGGER = logging.getLogger(__name__)


def _copy_tree(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Template directory not found: {src}")
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def collect_geometry_cases(geometry_path: Path) -> List[Path]:
    # Retrieve all subdirectories under the specified geometry path.
    if not geometry_path.is_dir():
        raise FileNotFoundError(f"Geometry path does not exist: {geometry_path}")
    entries = [p for p in geometry_path.iterdir() if p.is_dir()]
    if not entries:
        raise FileNotFoundError(f"No geometry subfolders found under {geometry_path}")
    return entries


def prepare_run_directories(config: Config, dry_run: bool = False) -> List[Path]:
    # Create run directories, each initialized with geometry and velocity files.
    logger = LOGGER.getChild("prepare")
    cases = collect_geometry_cases(config.general.geometry_path)
    if len(cases) < config.physics.nk:
        raise ValueError(
            f"Found {len(cases)} geometry folders but need {config.physics.nk} for nk trajectories"
        )

    random.shuffle(cases)
    selected = cases[: config.physics.nk]
    logger.info("Preparing %s run directories from %s geometries", config.physics.nk, len(cases))

    run_dirs: List[Path] = []
    base = config.path.parent
    for idx, case in enumerate(selected):
        run_dir = base / f"run-{idx}"
        run_dirs.append(run_dir)
        logger.debug("Setting up %s from source %s", run_dir, case)
        if dry_run:
            continue
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)
        xyz_src = case / config.general.init_xyz_file
        vel_src = case / config.general.init_vel_file
        if not xyz_src.is_file():
            raise FileNotFoundError(f"Initial geometry file not found: {xyz_src}")
        if not vel_src.is_file():
            raise FileNotFoundError(f"Initial velocity file not found: {vel_src}")
        shutil.copy2(xyz_src, run_dir / "initXYZ.dat")
        shutil.copy2(vel_src, run_dir / "initPxPyPz.dat")
    return run_dirs

