import os

from .config import Config

RED = "\033[31m"
GREEN = "\033[32m"
RESET = "\033[0m"


def stage3_run_dynamics(cfg: Config):
    """
    1. Write run.conf with server + client commands
    2. Submit run.sh via sbatch
    """

    param = cfg.param
    N = cfg.physics.nk if cfg.physics.use_param_nk else cfg.general.n_trajectories
    foldprefix = cfg.general.fold_prefix

    cwd = cfg.root_dir
    f = open("run.conf", "w")

    print(f"{GREEN}Server started with {N} clients.{RESET}")

    cmd_server = f"python {cwd}/server_DFTB.py {N}"
    print(0, " ", cmd_server, file=f)

    for i in range(N):
        cmd_client = f"python {cwd}/{foldprefix}{i}/client_DFTB.py {i}"
        print(i + 1, " ", cmd_client, file=f)

    f.close()

    if cfg.hpc.run_dynamics:
        os.system("sbatch run.sh")
        print(f"{RED} {N} Clients started.{RESET}")
    else:
        print(f"{RED}run_dynamics = false → not submitting run.sh.{RESET}")
