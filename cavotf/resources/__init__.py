"""Resource helpers for bundled simulation templates."""
from __future__ import annotations

from importlib import resources
from pathlib import Path


def default_template_dir() -> Path:
    """Return the bundled clean template directory bundled with the package."""
    return Path(resources.files(__package__).joinpath("DFTB_clean"))


def server_script_path() -> Path:
    """Return the path to the packaged server script."""
    return Path(resources.files("cavotf").joinpath("server_DFTB.py"))
