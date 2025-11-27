from __future__ import annotations

import configparser
from dataclasses import dataclass
import importlib.resources as resources
from pathlib import Path
from typing import Optional

RESOURCE_PACKAGE = "cavotf.resources"


def _resource_path(name: str) -> Path:
    return Path(resources.files(RESOURCE_PACKAGE) / name)


def _load_parser(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(path)
    return parser


def _get_bool(parser: configparser.ConfigParser, section: str, option: str, default: bool) -> bool:
    if parser.has_option(section, option):
        return parser.getboolean(section, option)
    return default


def _get_float(parser: configparser.ConfigParser, section: str, option: str, default: float) -> float:
    if parser.has_option(section, option):
        raw_value = parser.get(section, option)
        try:
            return float(raw_value)
        except ValueError:
            try:
                return float(eval(raw_value, {"__builtins__": {}}, {}))  # type: ignore[arg-type]
            except Exception as exc:  # noqa: BLE001
                raise ValueError(f"Could not parse numeric value '{raw_value}' for {section}.{option}") from exc
    return default


def _resolve_path(value: str | Path, base_dir: Path, default: Optional[Path] = None) -> Path:
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


@dataclass
class GeneralConfig:
    working_dir: Path
    data_source_dir: Path
    clean_template_dir: Path
    client_sbatch: Path
    server_sbatch: Path
    mu_submission_script: Path | str
    run_prefix: str = "run-"
    server_wait_seconds: int = 30
    mu_results_filename: str = "dmu.dat"
    mu_poll_interval: int = 60
    mu_timeout_seconds: int = 0
    thermalization_multiplier: int = 50
    initial_positions_name: str = "thermaliz_water__InTheBox.dat"
    initial_velocities_name: str = "thermaliz_water_vel__InTheBox.dat"
    initial_position_target: str = "initXYZ.dat"
    initial_velocity_target: str = "initPxPyPz.dat"
    n_trajectories: int = 25
    geometry_path: Optional[Path] = None
    init_xyz_file: str = "thermaliz_water__InTheBox.dat"
    init_vel_file: str = "thermaliz_water_vel__InTheBox.dat"
    trajectory_prefix: str = "run-"
    use_thermostat: bool = False
    thermostat_type: str = "andersen"
    thermostat_steps: int = 250
    collision_frequency: float = 0.0
    thermostat_reassign_particles: int = 0


@dataclass
class GeometryConfig:
    geom_root_head: Path
    geom_root: Path
    use_head_root: bool = False

    @property
    def active_root(self) -> Path:
        return self.geom_root_head if self.use_head_root else self.geom_root


@dataclass
class PhysicsConfig:
    use_param_nk: bool = True
    nk: int = 25
    beta: float = 1052.8
    lambda_: float = 0.001
    omega_c: float = 0.190 / 27.2114
    eta_b: float = 0.0003
    use_thermostat: bool = True
    thermostat_type: str = "andersen"
    thermostat_steps: int = 250
    thermostat_reassign_particles: int = 16


@dataclass
class HPCConfig:
    cpus_per_job: int = 1
    partition: str | None = None
    account: str | None = None
    run_get_mu: bool = True
    run_dynamics: bool = True
    dftb_prefix: Path | None = None
    dftb_command: str | None = None


class CavOTFConfig:
    def __init__(self, general: GeneralConfig, geometry: GeometryConfig, physics: PhysicsConfig, hpc: HPCConfig):
        self.general = general
        self.geometry = geometry
        self.physics = physics
        self.hpc = hpc

    @classmethod
    def from_file(cls, path: Path | str) -> "CavOTFConfig":
        path = Path(path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        base_dir = path.parent
        parser = _load_parser(path)

        general = _build_general(parser, base_dir)
        geometry = _build_geometry(parser, base_dir, general)
        physics = _build_physics(parser)
        hpc = _build_hpc(parser, base_dir)

        # If the user provided a geometry_path, prefer it for data sources
        if general.geometry_path is not None:
            general.data_source_dir = general.geometry_path

        # Honor geometry section selection
        general.data_source_dir = geometry.active_root

        return cls(general=general, geometry=geometry, physics=physics, hpc=hpc)

    # Properties used throughout the workflow for backwards compatibility
    @property
    def working_dir(self) -> Path:
        return self.general.working_dir

    @property
    def data_source_dir(self) -> Path:
        return self.general.data_source_dir

    @property
    def clean_template_dir(self) -> Path:
        return self.general.clean_template_dir

    @property
    def client_sbatch(self) -> Path:
        return self.general.client_sbatch

    @property
    def server_sbatch(self) -> Path:
        return self.general.server_sbatch

    @property
    def mu_submission_script(self) -> Path | str:
        return self.general.mu_submission_script

    @property
    def run_prefix(self) -> str:
        return self.general.run_prefix

    @property
    def server_wait_seconds(self) -> int:
        return self.general.server_wait_seconds

    @property
    def mu_results_filename(self) -> str:
        return self.general.mu_results_filename

    @property
    def mu_poll_interval(self) -> int:
        return self.general.mu_poll_interval

    @property
    def mu_timeout_seconds(self) -> int:
        return self.general.mu_timeout_seconds

    @property
    def thermalization_multiplier(self) -> int:
        return self.general.thermalization_multiplier

    @property
    def initial_positions_name(self) -> str:
        return self.general.initial_positions_name

    @property
    def initial_velocities_name(self) -> str:
        return self.general.initial_velocities_name

    @property
    def initial_position_target(self) -> str:
        return self.general.initial_position_target

    @property
    def initial_velocity_target(self) -> str:
        return self.general.initial_velocity_target


def _build_general(parser: configparser.ConfigParser, base_dir: Path) -> GeneralConfig:
    working_dir = _resolve_path(parser.get("general", "working_dir", fallback=base_dir), base_dir)

    clean_template_dir = _resolve_path(
        parser.get("general", "clean_template_dir", fallback="DEFAULT"),
        base_dir,
        default=_resource_path("DFTB_clean"),
    )
    client_sbatch = _resolve_path(
        parser.get("general", "client_sbatch", fallback="DEFAULT"),
        base_dir,
        default=_resource_path("run_client.sbatch"),
    )
    server_sbatch = _resolve_path(
        parser.get("general", "server_sbatch", fallback="DEFAULT"),
        base_dir,
        default=_resource_path("run_server.sbatch"),
    )

    mu_submission_script_raw = Path(parser.get("hpc", "mu_submission_script", fallback="muDFTB.slurm"))
    if not mu_submission_script_raw.is_absolute():
        candidate = (base_dir / mu_submission_script_raw).resolve()
        mu_submission_script = candidate if candidate.exists() else mu_submission_script_raw
    else:
        mu_submission_script = mu_submission_script_raw

    geometry_path = parser.get("general", "geometry_path", fallback=None)
    geometry_path_resolved: Optional[Path] = None
    if geometry_path:
        geometry_path_resolved = _resolve_path(geometry_path, base_dir)

    return GeneralConfig(
        working_dir=working_dir,
        data_source_dir=_resolve_path(
            parser.get("general", "data_source_dir", fallback=geometry_path_resolved or working_dir),
            base_dir,
        ),
        clean_template_dir=clean_template_dir,
        client_sbatch=client_sbatch,
        server_sbatch=server_sbatch,
        mu_submission_script=mu_submission_script,
        run_prefix=parser.get("general", "trajectory_prefix", fallback=parser.get("general", "run_prefix", fallback="run-")),
        server_wait_seconds=parser.getint("hpc", "server_wait_seconds", fallback=30),
        mu_results_filename=parser.get("general", "mu_results_filename", fallback="dmu.dat"),
        mu_poll_interval=parser.getint("general", "mu_poll_interval", fallback=60),
        mu_timeout_seconds=parser.getint("general", "mu_timeout_seconds", fallback=0),
        thermalization_multiplier=parser.getint("general", "thermalization_multiplier", fallback=50),
        initial_positions_name=parser.get("general", "init_xyz_file", fallback="thermaliz_water__InTheBox.dat"),
        initial_velocities_name=parser.get("general", "init_vel_file", fallback="thermaliz_water_vel__InTheBox.dat"),
        initial_position_target=parser.get("general", "initial_position_target", fallback="initXYZ.dat"),
        initial_velocity_target=parser.get("general", "initial_velocity_target", fallback="initPxPyPz.dat"),
        n_trajectories=parser.getint("general", "n_trajectories", fallback=25),
        geometry_path=geometry_path_resolved,
        init_xyz_file=parser.get("general", "init_xyz_file", fallback="thermaliz_water__InTheBox.dat"),
        init_vel_file=parser.get("general", "init_vel_file", fallback="thermaliz_water_vel__InTheBox.dat"),
        trajectory_prefix=parser.get("general", "trajectory_prefix", fallback="run-"),
        use_thermostat=_get_bool(parser, "general", "use_thermostat", False),
        thermostat_type=parser.get("general", "thermostat_type", fallback="andersen"),
        thermostat_steps=parser.getint("general", "thermostat_steps", fallback=250),
        collision_frequency=_get_float(parser, "general", "collision_frequency", 0.0),
        thermostat_reassign_particles=parser.getint("general", "thermostat_reassign_particles", fallback=0),
    )


def _build_geometry(parser: configparser.ConfigParser, base_dir: Path, general: GeneralConfig) -> GeometryConfig:
    geom_root_head = _resolve_path(
        parser.get("geometry", "geom_root_head", fallback=general.geometry_path or general.data_source_dir),
        base_dir,
    )
    geom_root = _resolve_path(
        parser.get("geometry", "geom_root", fallback=general.geometry_path or general.data_source_dir),
        base_dir,
    )

    return GeometryConfig(
        geom_root_head=geom_root_head,
        geom_root=geom_root,
        use_head_root=_get_bool(parser, "geometry", "use_head_root", False),
    )


def _build_physics(parser: configparser.ConfigParser) -> PhysicsConfig:
    return PhysicsConfig(
        use_param_nk=_get_bool(parser, "physics", "use_param_nk", True),
        nk=parser.getint("physics", "nk", fallback=25),
        beta=_get_float(parser, "physics", "beta", 1052.8),
        lambda_=_get_float(parser, "physics", "lambda", 0.001),
        omega_c=_get_float(parser, "physics", "omega_c", 0.190 / 27.2114),
        eta_b=_get_float(parser, "physics", "eta_b", 0.0003),
        use_thermostat=_get_bool(parser, "physics", "use_thermostat", True),
        thermostat_type=parser.get("physics", "thermostat_type", fallback="andersen"),
        thermostat_steps=parser.getint("physics", "thermostat_steps", fallback=250),
        thermostat_reassign_particles=parser.getint("physics", "thermostat_reassign_particles", fallback=16),
    )


def _build_hpc(parser: configparser.ConfigParser, base_dir: Path) -> HPCConfig:
    dftb_prefix = parser.get("hpc", "dftb_prefix", fallback=None)
    dftb_prefix_path = _resolve_path(dftb_prefix, base_dir) if dftb_prefix else None

    return HPCConfig(
        cpus_per_job=parser.getint("hpc", "cpus_per_job", fallback=1),
        partition=parser.get("hpc", "partition", fallback=None),
        account=parser.get("hpc", "account", fallback=None),
        run_get_mu=_get_bool(parser, "hpc", "run_get_mu", True),
        run_dynamics=_get_bool(parser, "hpc", "run_dynamics", True),
        dftb_prefix=dftb_prefix_path,
        dftb_command=parser.get("hpc", "dftb_command", fallback=None),
    )
