# =============================================================================
#  Project:     cavOTF.py
#  File:        config.py
#  Author:      Amir H. Amini <amiramini@tamu.edu>
#  Last update: 11/28/2025
#
#  Description:
#      Configuration management for cavOTF.py simulations, including parameter
#      handling and validation. This file controls the input options for the
#      cavOTF.py workflow.
# =============================================================================

from __future__ import annotations

import configparser
import dataclasses
import logging
from pathlib import Path
from typing import Optional

from .resources import default_template_dir

LOGGER = logging.getLogger(__name__)


def _parse_bool(value: str) -> bool:
    # Make sure to parse boolean values robustly! Include every common variant.
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"Invalid boolean value: {value}")


def _parse_float(value: str) -> float:
    # Gracefully handle float parsing, including simple expressions! 
    # Maybe useful for legacy inputs, maybe not.
    # Might consider removing this later and enforcing strict float inputs.
    try:
        return float(value)
    except ValueError:
        try:
            # Allow simple expressions like 0.190/27.2114 used in legacy inputs
            return float(eval(value, {"__builtins__": {}}, {}))
        except Exception as exc:  # noqa: BLE001
            raise ValueError(f"Invalid float expression: {value}") from exc


@dataclasses.dataclass
class GeneralConfig:
    # General settings for the simulation.
    # I should consider permanently moving thermostat settings to PhysicsConfig later
    # it does feel more appropriate there!

    geometry_path: Path
    init_xyz_file: str
    init_vel_file: str
    trajectory_prefix: str = "trajectory_"
    clean_template_dir: Path = default_template_dir()
    use_thermostat: bool = True
    thermostat_type: str = "andersen"
    thermostat_steps: int = 250
    thermostat_BetaScale: float = 0.75
    collision_frequency: float = 0.001
    thermostat_reassign_particles: int = 16


@dataclasses.dataclass
class PhysicsConfig:
    # Cavity parameters and physical settings.

    nk: int
    beta: float
    beta_run0: float
    lambda_: float
    omega_c: float
    eta_b: float
    use_thermostat: bool = True
    thermostat_type: str = "andersen"
    thermostat_steps: int = 250
    thermostat_BetaScale: float = 0.75
    thermostat_reassign_particles: int = 16
    calculate_dipole_derivatives: bool = True
    dipole_derivative_interval: int = 10
    thermal_steps: int = 1000


@dataclasses.dataclass
class HPCConfig:
    # HPC-related settings for job submission.
    # I might consider adding support for different sbatch templates for each
    # job type later if needed.

    cpus_per_job: int = 1
    partition: Optional[str] = None
    account: Optional[str] = None
    run_get_mu: bool = True
    run_dynamics: bool = True
    dftb_prefix: Optional[Path] = None
    dftb_command: Optional[str] = None
    sbatch_template: Optional[Path] = None

    walltime: str = "2-00:00:00"
    memory: str = "80G"


@dataclasses.dataclass
class OutputConfig:
    # Output control settings.

    write_logfile: bool = True
    write_results: bool = True
    record_k_space: bool = True
    print_k_space: bool = False
    write_xyz_trajectory: bool = True
    write_histogram: bool = True
    write_output_client: bool = True
    write_midpoint_snapshots: bool = True


@dataclasses.dataclass
class DFTBConfig:
    # DFTB+ calculator options
    # Generates DFTB+ input.

    parameters: dict[str, object]

    @staticmethod
    def defaults() -> "DFTBConfig":
        # Default DFTB+ parameters for cavOTF.py simulations suited for O33H66 system ONLY!
        return DFTBConfig(
            parameters={
                "label": "O33H66",
                "Hamiltonian_SCC": "Yes",
                "Hamiltonian_SCCTolerance": 1e-5,
                "Hamiltonian_MaxSCCIterations": 400,
                "Hamiltonian_Mixer": "Anderson { MixingParameter = 0.026 }",
                "Hamiltonian_ConvergentSccOnly": "No",
                "Hamiltonian_MaxAngularMomentum_": "",
                "Hamiltonian_MaxAngularMomentum_O": "p",
                "Hamiltonian_MaxAngularMomentum_H": "s",
                "Hamiltonian_Charge": 0.0,
                "Hamiltonian_ReadInitialCharges": "Yes",
                "Hamiltonian_SpinConstants": "{O = { -0.035 -0.030 -0.030 -0.028 } H = { -0.072 }}",
                "Options_WriteDetailedXml": "No",
                "Options_WriteEigenvectors": "No",
                "Options_WriteResultsTag": "Yes",
                "kpts": (3, 3, 3),
            }
        )

    def merged(self, overrides: dict[str, object]) -> "DFTBConfig":
        # Return a new DFTBConfig with overridden parameters.
        # Gracefully handle None values to remove parameters.
        merged = dict(self.parameters)
        for key, value in overrides.items():
            if value is None:
                merged.pop(key, None)
            else:
                merged[key] = value
        return DFTBConfig(parameters=merged)


@dataclasses.dataclass
class Config:
    # Complete configuration for a cavOTF.py simulation.

    path: Path
    general: GeneralConfig
    physics: PhysicsConfig
    hpc: HPCConfig
    outputs: OutputConfig
    dftb: DFTBConfig


class ConfigError(RuntimeError):
    """Parsing fails."""
    
class MissingOptionError(ConfigError):
    """Required option is missing."""
    
class InvalidOptionError(ConfigError):
    """Option has invalid value."""


def _require(section: configparser.SectionProxy, option: str) -> str:
    if option not in section:
        raise MissingOptionError(f"Missing required option '{option}' in section [{section.name}]")
    return section[option]


def _resolve_path(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path)


def _load_outputs(config: configparser.ConfigParser) -> OutputConfig:
    section = config["outputs"] if "outputs" in config else {}

    def _get_bool(option: str, default: bool) -> bool:
        return _parse_bool(section.get(option, str(default)))

    if "write_output_client" in section:
        write_output_client = _get_bool("write_output_client", True)
    elif "write_phase_space" in section:
        write_output_client = _parse_bool(section["write_phase_space"])
    else:
        write_output_client = True

    return OutputConfig(
        write_logfile=_get_bool("write_logfile", True),
        write_results=_get_bool("write_results", True),
        record_k_space=_get_bool("record_k_space", True),
        print_k_space=_get_bool("print_k_space", False),
        write_xyz_trajectory=_get_bool("write_xyz_trajectory", True),
        write_histogram=_get_bool("write_histogram", True),
        write_output_client=write_output_client,
        write_midpoint_snapshots=_get_bool("write_midpoint_snapshots", True),
    )


def _parse_value(raw: str) -> object:
    # Parse a raw string value into an appropriate Python type.
    import ast

    normalized = raw.strip()
    if normalized.lower() in {"off", "none", "false", "no", ""}:
        return None
    for parser in (int, float):
        try:
            return parser(normalized)
        except ValueError:
            continue
    try:
        return ast.literal_eval(normalized)
    except Exception:  # noqa: BLE001
        return raw


def _clean_string_value(value: str) -> str:
    # Clean up string values by collapsing whitespace and removing newlines
    # Avoid unintended formatting issues in DFTB+ input.

    collapsed = value.replace("\\n", " ").replace("\\\n", " ")
    collapsed = " ".join(collapsed.split())
    return collapsed.strip()


def _load_dftb(config: configparser.ConfigParser) -> DFTBConfig:
    defaults = DFTBConfig.defaults()
    if "dftb" not in config:
        return defaults
    section = config["dftb"]
    overrides: dict[str, object] = {}
    for key, raw in section.items():
        value = _parse_value(raw)
        if isinstance(value, str):
            value = _clean_string_value(value)
        overrides[key] = value
    return defaults.merged(overrides)


def _load_general(config: configparser.ConfigParser, base: Path) -> GeneralConfig:
    if "general" not in config:
        raise MissingOptionError("Missing required [general] section")
    section = config["general"]

    geometry_path = _resolve_path(base, _require(section, "geometry_path"))
    init_xyz_file = _require(section, "init_xyz_file")
    init_vel_file = _require(section, "init_vel_file")
    trajectory_prefix = section.get("trajectory_prefix", "trajectory_")

    template_override = section.get("clean_template_dir")
    clean_template_dir = default_template_dir()
    if template_override:
        candidate = _resolve_path(base, template_override)
        if candidate.is_dir():
            clean_template_dir = candidate
        else:
            LOGGER.warning(
                "clean_template_dir '%s' not found; falling back to packaged templates",
                candidate,
            )
    use_thermostat = _parse_bool(section.get("use_thermostat", "true"))
    thermostat_type = section.get("thermostat_type", "andersen")
    thermostat_steps = int(section.get("thermostat_steps", 250))
    thermostat_BetaScale = _parse_float(section.get("thermostat_BetaScale", "0.75"))
    collision_frequency = _parse_float(section.get("collision_frequency", "0.001"))
    thermostat_reassign_particles = int(section.get("thermostat_reassign_particles", 16))

    return GeneralConfig(
        geometry_path=geometry_path,
        init_xyz_file=init_xyz_file,
        init_vel_file=init_vel_file,
        trajectory_prefix=trajectory_prefix,
        clean_template_dir=clean_template_dir,
        use_thermostat=use_thermostat,
        thermostat_type=thermostat_type,
        thermostat_steps=thermostat_steps,
        thermostat_BetaScale=thermostat_BetaScale,
        collision_frequency=collision_frequency,
        thermostat_reassign_particles=thermostat_reassign_particles,
    )


def _load_physics(config: configparser.ConfigParser) -> PhysicsConfig:
    if "physics" not in config:
        raise MissingOptionError("Missing required [physics] section")
    section = config["physics"]

    nk = int(_require(section, "nk"))
    beta = _parse_float(_require(section, "beta"))
    beta_run0_raw = section.get("beta_run0")
    beta_run0 = _parse_float(beta_run0_raw) if beta_run0_raw is not None else beta
    lambda_ = _parse_float(_require(section, "lambda"))
    omega_c = _parse_float(_require(section, "omega_c"))
    eta_b = _parse_float(_require(section, "eta_b"))
    use_thermostat = _parse_bool(section.get("use_thermostat", "true"))
    thermostat_type = section.get("thermostat_type", "andersen")
    thermostat_steps = int(section.get("thermostat_steps", 250))
    thermostat_BetaScale = _parse_float(section.get("thermostat_BetaScale", "0.75"))
    thermostat_reassign_particles = int(section.get("thermostat_reassign_particles", 16))
    calculate_dipole_derivatives = _parse_bool(section.get("calculate_dipole_derivatives", "true"))
    dipole_derivative_interval = int(section.get("dipole_derivative_interval", 10))
    thermal_steps = int(section.get("thermal_steps", 1000))

    return PhysicsConfig(
        nk=nk,
        beta=beta,
        beta_run0=beta_run0,
        lambda_=lambda_,
        omega_c=omega_c,
        eta_b=eta_b,
        use_thermostat=use_thermostat,
        thermostat_type=thermostat_type,
        thermostat_steps=thermostat_steps,
        thermostat_BetaScale=thermostat_BetaScale,
        thermostat_reassign_particles=thermostat_reassign_particles,
        calculate_dipole_derivatives=calculate_dipole_derivatives,
        dipole_derivative_interval=dipole_derivative_interval,
        thermal_steps=thermal_steps,
    )


def _load_hpc(config: configparser.ConfigParser, base: Path) -> HPCConfig:
    if "hpc" not in config:
        raise MissingOptionError("Missing required [hpc] section")
    section = config["hpc"]

    cpus_per_job = int(section.get("cpus_per_job", 1))
    partition = section.get("partition")
    account = section.get("account")
    run_get_mu = _parse_bool(section.get("run_get_mu", "true"))
    run_dynamics = _parse_bool(section.get("run_dynamics", "true"))
    dftb_prefix_raw = section.get("dftb_prefix")
    dftb_prefix = _resolve_path(base, dftb_prefix_raw) if dftb_prefix_raw else None
    dftb_command = section.get("dftb_command")
    sbatch_template_raw = section.get("sbatch_template")
    sbatch_template = _resolve_path(base, sbatch_template_raw) if sbatch_template_raw else None
    walltime = section.get("walltime", "2-00:00:00")
    memory = section.get("memory", "80G")

    return HPCConfig(
        cpus_per_job=cpus_per_job,
        partition=partition,
        account=account,
        run_get_mu=run_get_mu,
        run_dynamics=run_dynamics,
        dftb_prefix=dftb_prefix,
        dftb_command=dftb_command,
        sbatch_template=sbatch_template,
        walltime=walltime,
        memory=memory,
    )


def load_config(path: Path) -> Config:
    # Load and parse the configuration file at the given path.
    absolute_path = path.expanduser().resolve()
    if not absolute_path.is_file():
        raise ConfigError(f"Config file not found: {absolute_path}")

    parser = configparser.ConfigParser()
    parser.optionxform = str  # preserve case
    with absolute_path.open() as handle:
        parser.read_file(handle)

    base = absolute_path.parent
    try:
        general = _load_general(parser, base)
        physics = _load_physics(parser)
        hpc = _load_hpc(parser, base)
        outputs = _load_outputs(parser)
        dftb = _load_dftb(parser)
    except (MissingOptionError, ValueError) as exc:
        raise InvalidOptionError(str(exc)) from exc

    LOGGER.debug("Configuration loaded from %s", absolute_path)
    return Config(
        path=absolute_path,
        general=general,
        physics=physics,
        hpc=hpc,
        outputs=outputs,
        dftb=dftb,
    )

from .resources import default_template_dir