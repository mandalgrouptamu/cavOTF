import os
import random
from pathlib import Path

from .config import Config


RED = "\033[31m"
GREEN = "\033[32m"
RESET = "\033[0m"


def stage1_prepare_folders(cfg: Config):
    """
    1. Create run-i folders (i = 0 .. N-1)
    2. Copy DFTB_clean content + run_client.sbatch into each
    3. Copy initial geometry + velocity into each as initXYZ.dat / initPxPyPz.dat
    4. Write get_mu.conf with python commands
    5. Submit get_mu.sh
    """

    param = cfg.param
    # number of "clients" / trajectories
    N = cfg.physics.nk if cfg.physics.use_param_nk else cfg.general.n_trajectories

    foldprefix = cfg.general.fold_prefix
    cleanFolder = cfg.general.clean_folder

    root = cfg.root_dir
    cwd = root  # base directory for this job

    # open get_mu.conf in root dir
    get_mu_conf_path = root / "get_mu.conf"
    f = get_mu_conf_path.open("w")

    # clean up previous run files
    os.system("rm -rf *.txt")
    os.system(f"rm -rf {foldprefix}*")

    # geometry folders
    location_head = cfg.geometry.geom_root_head
    location = cfg.geometry.geom_root

    # choose folders randomly
    folders_head = [
        d for d in os.listdir(location_head)
        if os.path.isdir(os.path.join(location_head, d))
    ]
    random.shuffle(folders_head)

    # ---------- First trajectory (index 0) ----------
    os.system(f"mkdir -p {foldprefix}0")
    os.system(f"cp -r {cleanFolder}/* {foldprefix}0/")
    os.system(f"cp run_client.sbatch {foldprefix}0/")

    os.system(
        f"cp {location_head}{folders_head[0]}/thermaliz_water__InTheBox.dat "
        f"{foldprefix}0/initXYZ.dat"
    )
    os.system(
        f"cp {location_head}{folders_head[0]}/thermaliz_water_vel__InTheBox.dat "
        f"{foldprefix}0/initPxPyPz.dat"
    )

    cmd0 = f"python {cwd}/{foldprefix}0/get_mu.py"
    print(0, " ", cmd0, file=f)

    # ---------- Remaining trajectories ----------
    folders = [
        d for d in os.listdir(location)
        if os.path.isdir(os.path.join(location, d))
    ]
    random.shuffle(folders)

    for i in range(N - 1):
        run_dir = f"{foldprefix}{i+1}"

        os.system(f"mkdir -p {run_dir}")
        os.system(f"cp -r {cleanFolder}/* {run_dir}/")
        os.system(f"cp run_client.sbatch {run_dir}/")

        os.system(
            f"cp {location}{folders[i+1]}/thermaliz_water__InTheBox.dat "
            f"{run_dir}/initXYZ.dat"
        )
        os.system(
            f"cp {location}{folders[i+1]}/thermaliz_water_vel__InTheBox.dat "
            f"{run_dir}/initPxPyPz.dat"
        )

        cmd = f"python {cwd}/{run_dir}/get_mu.py"
        print(i + 1, " ", cmd, file=f)

    f.close()

    # submit get_mu.sh via sbatch
    if cfg.hpc.run_get_mu:
        os.system("sbatch get_mu.sh")
        print(f"{GREEN}Submitted get_mu.sh for {N} runs.{RESET}")
    else:
        print(f"{RED}run_get_mu = false → not submitting get_mu.sh.{RESET}")
