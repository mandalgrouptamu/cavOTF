"""Legacy get_mu script wired for shared templates and config overrides."""
from __future__ import annotations

import argparse
import os
import pathlib
import time

import numpy as np
from ase import Atoms
from ase.io import write

from dftb import getForcesCharges, getdµ, set_calculator_options
from funcLM import param

try:
    from cavotf.config import load_config
    from cavotf.dynamics import _recompute_mode_grid
except Exception:  # noqa: BLE001
    load_config = None
    _recompute_mode_grid = None


ATM_SYMBOLS = "O33H66"  # Retained to preserve legacy behavior


def apply_config_overrides(params: param, config_path: pathlib.Path | None) -> None:
    """Apply configuration overrides from ``input.txt`` if available."""
    if not config_path or load_config is None or _recompute_mode_grid is None:
        return

    cfg = load_config(config_path)
    overrides = {
        "nk": cfg.physics.nk,
        "β": cfg.physics.beta,
        "λ": cfg.physics.lambda_,
        "ωc": cfg.physics.omega_c,
        "ηb": cfg.physics.eta_b,
        "thermal_steps": cfg.physics.thermal_steps,
    }
    for key, value in overrides.items():
        if hasattr(params, key):
            setattr(params, key, value)

    # Keep cavity grid consistent with overrides
    _recompute_mode_grid(params)

    if cfg.hpc.dftb_prefix:
        os.environ["DFTB_PREFIX"] = str(cfg.hpc.dftb_prefix)
    if cfg.hpc.dftb_command:
        os.environ["DFTB_COMMAND"] = cfg.hpc.dftb_command
    set_calculator_options(cfg.dftb.parameters)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", type=str, default=None, help="Directory containing initXYZ.dat and initPxPyPz.dat")
    parser.add_argument("--config", type=str, default=None, help="Path to cavotf input.txt")
    args = parser.parse_args()

    workdir = pathlib.Path(args.workdir) if args.workdir else pathlib.Path.cwd()
    os.chdir(workdir)

    params = param()
    config_path = pathlib.Path(args.config) if args.config else None
    apply_config_overrides(params, config_path)

    print("Start")
    t0 = time.time()
    bhr = 1.8897259886
    evdivA = 27.2114 * bhr
    AngdivPs2AU = bhr / 41341.3733365614

    natoms = params.natoms
    dt = params.dt
    thermal_steps = params.thermal_steps

    coordina_initial = np.loadtxt("initXYZ.dat", usecols=(2, 3, 4))
    velocity_initial = np.loadtxt("initPxPyPz.dat", usecols=(0, 1, 2))
    coordinates = np.array(coordina_initial) * bhr  # a.u.

    atoms = Atoms(ATM_SYMBOLS, positions=coordinates / bhr)
    mass = atoms.get_masses() * 1822.8884
    box = params.box
    atoms.set_cell([box, box, box])
    atoms.set_pbc(True)

    coordinates = atoms.get_positions(wrap=False) * bhr
    atoms.set_positions(coordinates / bhr)

    write("WaterMD_Cavity.xyz", atoms, format="xyz", append=True)

    rxj, ryj, rzj = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]
    rj = np.concatenate((rxj, ryj, rzj))
    Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)

    fj, charges = getForcesCharges(rj, natoms, ATM_SYMBOLS, box)
    µj = np.sum(charges * rxj)

    np.savetxt("dmu.dat", [µj], fmt="%.6f")
    print(f"Completed in {time.time() - t0:.2f} s")


if __name__ == "__main__":
    main()
