import configparser
from dataclasses import dataclass
from pathlib import Path

from .funcLM import param as ParamClass  # your original param


@dataclass
class GeneralConfig:
    job_name: str
    n_trajectories: int
    clean_folder: str
    fold_prefix: str


@dataclass
class GeometryConfig:
    geom_root_head: str
    geom_root: str
    use_head_root: bool


@dataclass
class PhysicsConfig:
    use_param_nk: bool
    nk: int  # only used if use_param_nk = False


@dataclass
class HPCConfig:
    cpus_per_job: int
    partition: str
    account: str
    run_get_mu: bool
    run_dynamics: bool


@dataclass
class Config:
    general: GeneralConfig
    geometry: GeometryConfig
    physics: PhysicsConfig
    hpc: HPCConfig
    param: object          # your original param() object
    root_dir: Path         # where we launched cavOTF from


def load_config(path: str = "input.txt") -> Config:
    cp = configparser.ConfigParser()
    read_files = cp.read(path)
    if not read_files:
        raise FileNotFoundError(f"Could not read input file: {path}")

    root_dir = Path(".").absolute()

    # ---------- general ----------
    g = cp["general"]
    general = GeneralConfig(
        job_name=g.get("job_name", "cavotf_job"),
        n_trajectories=g.getint("n_trajectories", 1),
        clean_folder=g.get("clean_folder", "DFTB_clean"),
        fold_prefix=g.get("fold_prefix", "run-"),
    )

    # ---------- geometry ----------
    geo = cp["geometry"]
    geometry = GeometryConfig(
        geom_root_head=geo.get("geom_root_head", geo.get("geom_root", ".")),
        geom_root=geo.get("geom_root", "."),
        use_head_root=geo.getboolean("use_head_root", True),
    )

    # ---------- physics ----------
    phys = cp["physics"]
    physics = PhysicsConfig(
        use_param_nk=phys.getboolean("use_param_nk", True),
        nk=phys.getint("nk", 25),
    )

    # ---------- hpc ----------
    hpc_sec = cp["hpc"]
    hpc = HPCConfig(
        cpus_per_job=hpc_sec.getint("cpus_per_job", 1),
        partition=hpc_sec.get("partition", "compute"),
        account=hpc_sec.get("account", ""),
        run_get_mu=hpc_sec.getboolean("run_get_mu", True),
        run_dynamics=hpc_sec.getboolean("run_dynamics", True),
    )

    # your original param object
    param_obj = ParamClass()

    # Option: if physics says use_param_nk = True, force N from param.nk
    if physics.use_param_nk:
        physics.nk = param_obj.nk

    return Config(
        general=general,
        geometry=geometry,
        physics=physics,
        hpc=hpc,
        param=param_obj,
        root_dir=root_dir,
    )
