from .config import load_config
from .stage1_get_mu import stage1_prepare_folders
from .stage2_thermal import stage2_init_cavity
from .stage3_run_md import stage3_run_dynamics


def run_all(input_path: str = "input.txt"):
    cfg = load_config(input_path)

    # Stage 1: build folders, copy geometries, submit get_mu
    stage1_prepare_folders(cfg)

    # Stage 2: after µj are computed (may be a separate batch job in practice)
    stage2_init_cavity(cfg)

    # Stage 3: run DFTB dynamics
    stage3_run_dynamics(cfg)
