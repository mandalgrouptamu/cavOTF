from __future__ import annotations

from dataclasses import dataclass
import importlib.resources as resources
from pathlib import Path
from typing import Dict

RESOURCE_PACKAGE = "cavotf.resources"


def _resource_path(name: str) -> Path:
    return Path(resources.files(RESOURCE_PACKAGE) / name)


@dataclass
class CavOTFConfig:
    working_dir: Path
    data_source_dir: Path
    clean_template_dir: Path
    run_prefix: str = "run-"
    client_sbatch: Path = _resource_path("run_client.sbatch")
    server_sbatch: Path = _resource_path("run_server.sbatch")
    mu_submission_script: Path | str = "muDFTB.slurm"
    server_wait_seconds: int = 30
    mu_results_filename: str = "dmu.dat"
    mu_poll_interval: int = 60
    mu_timeout_seconds: int = 0
    thermalization_multiplier: int = 50
    initial_positions_name: str = "thermaliz_water__InTheBox.dat"
    initial_velocities_name: str = "thermaliz_water_vel__InTheBox.dat"
    initial_position_target: str = "initXYZ.dat"
    initial_velocity_target: str = "initPxPyPz.dat"

    @classmethod
    def from_file(cls, path: Path | str) -> "CavOTFConfig":
        path = Path(path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        base_dir = path.parent
        values = _parse_key_values(path)

        working_dir = _resolve_path(values.get("working_dir", base_dir), base_dir)
        data_source_dir = _resolve_path(
            values.get(
                "data_source_dir",
                "/scratch/user/u.sw216206/vsc-fluxside/project1/cavityDFTB/with_dmu/",
            ),
            base_dir,
        )

        clean_template_dir = _resolve_path(
            values.get("clean_template_dir", "DEFAULT"), base_dir, default=_resource_path("DFTB_clean")
        )

        client_sbatch = _resolve_path(
            values.get("client_sbatch", "DEFAULT"), base_dir, default=_resource_path("run_client.sbatch")
        )
        server_sbatch = _resolve_path(
            values.get("server_sbatch", "DEFAULT"), base_dir, default=_resource_path("run_server.sbatch")
        )

        mu_submission_script_raw = Path(values.get("mu_submission_script", "muDFTB.slurm"))
        if not mu_submission_script_raw.is_absolute():
            candidate = (base_dir / mu_submission_script_raw).resolve()
            mu_submission_script = candidate if candidate.exists() else mu_submission_script_raw
        else:
            mu_submission_script = mu_submission_script_raw

        return cls(
            working_dir=working_dir,
            data_source_dir=data_source_dir,
            clean_template_dir=clean_template_dir,
            run_prefix=values.get("run_prefix", "run-"),
            client_sbatch=client_sbatch,
            server_sbatch=server_sbatch,
            mu_submission_script=mu_submission_script,
            server_wait_seconds=int(values.get("server_wait_seconds", 30)),
            mu_results_filename=values.get("mu_results_filename", "dmu.dat"),
            mu_poll_interval=int(values.get("mu_poll_interval", 60)),
            mu_timeout_seconds=int(values.get("mu_timeout_seconds", 0)),
            thermalization_multiplier=int(values.get("thermalization_multiplier", 50)),
            initial_positions_name=values.get(
                "initial_positions_name", "thermaliz_water__InTheBox.dat"
            ),
            initial_velocities_name=values.get(
                "initial_velocities_name", "thermaliz_water_vel__InTheBox.dat"
            ),
            initial_position_target=values.get("initial_position_target", "initXYZ.dat"),
            initial_velocity_target=values.get("initial_velocity_target", "initPxPyPz.dat"),
        )


def _parse_key_values(path: Path) -> Dict[str, str]:
    pairs: Dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"Invalid line in config: '{raw_line}'")
        key, value = line.split("=", 1)
        pairs[key.strip()] = value.strip()
    return pairs


def _resolve_path(value: str | Path, base_dir: Path, default: Path | None = None) -> Path:
    if isinstance(value, str) and value.upper() == "DEFAULT":
        if default is None:
            raise ValueError("DEFAULT specified but no default provided")
        return default

    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()

    if default is not None and not candidate.exists():
        return default
    return candidate
