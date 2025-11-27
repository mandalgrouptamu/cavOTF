from __future__ import annotations

from configparser import ConfigParser
from dataclasses import dataclass
import importlib.resources as resources
from pathlib import Path

RESOURCE_PACKAGE = "cavotf.resources"


def _resource_path(name: str) -> Path:
    return Path(resources.files(RESOURCE_PACKAGE) / name)


@dataclass
class GeneralSettings:
    working_dir: Path
    n_trajectories: int
    geometry_path: Path | None
    init_xyz_file: str
    init_vel_file: str
    trajectory_prefix: str
    use_thermostat: bool
    thermostat_type: str
    thermostat_steps: int
    collision_frequency: float
    thermostat_reassign_particles: int
    clean_template_dir: Path
    client_sbatch: Path
    server_sbatch: Path
    mu_submission_script: Path | str
    mu_results_filename: str
    server_wait_seconds: int
    mu_poll_interval: int
    mu_timeout_seconds: int
    thermalization_multiplier: int
    initial_position_target: str
    initial_velocity_target: str


@dataclass
class GeometrySettings:
    geom_root_head: Path
    geom_root: Path
    use_head_root: bool


@dataclass
class PhysicsSettings:
    use_param_nk: bool
    nk: int
    beta: float
    lambda_: float
    omega_c: float
    eta_b: float
    use_thermostat: bool
    thermostat_type: str
    thermostat_steps: int
    thermostat_reassign_particles: int


@dataclass
class HPCSettings:
    cpus_per_job: int
    partition: str
    account: str
    run_get_mu: bool
    run_dynamics: bool
    dftb_prefix: str
    dftb_command: str


@dataclass
class CavOTFConfig:
    general: GeneralSettings
    geometry: GeometrySettings
    physics: PhysicsSettings
    hpc: HPCSettings
    base_dir: Path

    @property
    def working_dir(self) -> Path:
        return self.general.working_dir

    @property
    def data_source_dir(self) -> Path:
        if self.general.geometry_path:
            return self.general.geometry_path
        return self.geometry.geom_root_head if self.geometry.use_head_root else self.geometry.geom_root

    @property
    def clean_template_dir(self) -> Path:
        return self.general.clean_template_dir

    @property
    def run_prefix(self) -> str:
        return self.general.trajectory_prefix

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
    def mu_results_filename(self) -> str:
        return self.general.mu_results_filename

    @property
    def server_wait_seconds(self) -> int:
        return self.general.server_wait_seconds

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
        return self.general.init_xyz_file

    @property
    def initial_velocities_name(self) -> str:
        return self.general.init_vel_file

    @property
    def initial_position_target(self) -> str:
        return self.general.initial_position_target

    @property
    def initial_velocity_target(self) -> str:
        return self.general.initial_velocity_target

    @classmethod
    def from_file(cls, path: Path | str) -> "CavOTFConfig":
        path = Path(path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Configuration file not found: {path}")

        base_dir = path.parent
        parser = ConfigParser(inline_comment_prefixes=("#", ";"))
        parser.optionxform = str  # preserve case
        parser.read(path)

        general = _load_general(parser, base_dir)
        geometry = _load_geometry(parser, base_dir, general.geometry_path)
        physics = _load_physics(parser)
        hpc = _load_hpc(parser)

        return cls(general=general, geometry=geometry, physics=physics, hpc=hpc, base_dir=base_dir)


def _load_general(parser: ConfigParser, base_dir: Path) -> GeneralSettings:
    section = "general"
    if section not in parser:
        raise ValueError("[general] section is required in the configuration file")
    get = parser[section].get

    working_dir = _resolve_path(get("working_dir", base_dir), base_dir)
    geometry_raw = get("geometry_path", "").strip()
    geometry_path = _resolve_path(geometry_raw, base_dir) if geometry_raw else None

    clean_template_dir = _resolve_path(
        get("clean_template_dir", "DEFAULT"), base_dir, default=_resource_path("DFTB_clean")
    )
    client_sbatch = _resolve_path(
        get("client_sbatch", "DEFAULT"), base_dir, default=_resource_path("run_client.sbatch")
    )
    server_sbatch = _resolve_path(
        get("server_sbatch", "DEFAULT"), base_dir, default=_resource_path("run_server.sbatch")
    )

    mu_submission_script_raw = Path(get("mu_submission_script", "muDFTB.slurm"))
    if not mu_submission_script_raw.is_absolute():
        candidate = (base_dir / mu_submission_script_raw).resolve()
        mu_submission_script = candidate if candidate.exists() else mu_submission_script_raw
    else:
        mu_submission_script = mu_submission_script_raw

    return GeneralSettings(
        working_dir=working_dir,
        n_trajectories=int(_as_number(get("n_trajectories", "1"))),
        geometry_path=geometry_path,
        init_xyz_file=get("init_xyz_file", "thermaliz_water__InTheBox.dat"),
        init_vel_file=get("init_vel_file", "thermaliz_water_vel__InTheBox.dat"),
        trajectory_prefix=get("trajectory_prefix", "run-"),
        use_thermostat=_as_bool(get("use_thermostat", "yes")),
        thermostat_type=get("thermostat_type", "andersen"),
        thermostat_steps=int(_as_number(get("thermostat_steps", "250"))),
        collision_frequency=float(_as_number(get("collision_frequency", "0.001"))),
        thermostat_reassign_particles=int(_as_number(get("thermostat_reassign_particles", "16"))),
        clean_template_dir=clean_template_dir,
        client_sbatch=client_sbatch,
        server_sbatch=server_sbatch,
        mu_submission_script=mu_submission_script,
        mu_results_filename=get("mu_results_filename", "dmu.dat"),
        server_wait_seconds=int(_as_number(get("server_wait_seconds", "30"))),
        mu_poll_interval=int(_as_number(get("mu_poll_interval", "60"))),
        mu_timeout_seconds=int(_as_number(get("mu_timeout_seconds", "0"))),
        thermalization_multiplier=int(_as_number(get("thermalization_multiplier", "50"))),
        initial_position_target=get("initial_position_target", "initXYZ.dat"),
        initial_velocity_target=get("initial_velocity_target", "initPxPyPz.dat"),
    )


def _load_geometry(parser: ConfigParser, base_dir: Path, geometry_path: Path | None) -> GeometrySettings:
    section = "geometry"
    if section not in parser:
        fallback = geometry_path if geometry_path is not None else base_dir
        return GeometrySettings(geom_root_head=fallback, geom_root=fallback, use_head_root=True)

    get = parser[section].get
    default_path = geometry_path if geometry_path is not None else base_dir
    geom_root_head = _resolve_path(get("geom_root_head", default_path), base_dir)
    geom_root = _resolve_path(get("geom_root", default_path), base_dir)
    use_head_root = _as_bool(get("use_head_root", "true"))
    return GeometrySettings(geom_root_head=geom_root_head, geom_root=geom_root, use_head_root=use_head_root)


def _load_physics(parser: ConfigParser) -> PhysicsSettings:
    section = "physics"
    if section not in parser:
        raise ValueError("[physics] section is required in the configuration file")
    get = parser[section].get
    return PhysicsSettings(
        use_param_nk=_as_bool(get("use_param_nk", "true")),
        nk=int(_as_number(get("nk", "25"))),
        beta=float(_as_number(get("beta", "1052.8"))),
        lambda_=float(_as_number(get("lambda", "0.001"))),
        omega_c=float(_as_number(get("omega_c", "0.190/27.2114"))),
        eta_b=float(_as_number(get("eta_b", "0.0003"))),
        use_thermostat=_as_bool(get("use_thermostat", "true")),
        thermostat_type=get("thermostat_type", "andersen"),
        thermostat_steps=int(_as_number(get("thermostat_steps", "250"))),
        thermostat_reassign_particles=int(_as_number(get("thermostat_reassign_particles", "16"))),
    )


def _load_hpc(parser: ConfigParser) -> HPCSettings:
    section = "hpc"
    if section not in parser:
        raise ValueError("[hpc] section is required in the configuration file")
    get = parser[section].get
    return HPCSettings(
        cpus_per_job=int(_as_number(get("cpus_per_job", "1"))),
        partition=get("partition", ""),
        account=get("account", ""),
        run_get_mu=_as_bool(get("run_get_mu", "true")),
        run_dynamics=_as_bool(get("run_dynamics", "true")),
        dftb_prefix=get("dftb_prefix", ""),
        dftb_command=get("dftb_command", ""),
    )


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


def _as_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def _as_number(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return float(eval(value, {}, {}))
