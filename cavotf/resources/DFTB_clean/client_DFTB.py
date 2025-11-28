"""Shared client script for CAVOTF dynamics.

Runs from a centralized location without copying into each run directory.
Configuration overrides from ``input.txt`` are applied to keep cavity
parameters consistent across server and clients.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import socket
import sys
import time
from types import SimpleNamespace

import numpy as np
from ase import Atoms
from ase.io import write
from dftb import getForcesCharges, getCharges, getdµ, set_calculator_options
from funcLM import *  # noqa: F403,F401

try:
    from cavotf.config import OutputConfig, load_config
    from cavotf.dynamics import _recompute_mode_grid
except Exception:  # noqa: BLE001
    load_config = None
    _recompute_mode_grid = None
    OutputConfig = None


def comm(msg, host, port):
    msg["execTime"] = time.time()
    payload = json.dumps(msg)
    with socket.socket() as cli:
        cli.connect((host, port))
        cli.sendall(payload.encode())
        reply = cli.recv(4096)
    return reply.decode().strip()


def apply_config_overrides(params, cfg):
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
    _recompute_mode_grid(params)

    if cfg.hpc.dftb_prefix:
        os.environ["DFTB_PREFIX"] = str(cfg.hpc.dftb_prefix)
    if cfg.hpc.dftb_command:
        os.environ["DFTB_COMMAND"] = cfg.hpc.dftb_command
    set_calculator_options(cfg.dftb.parameters)


def _default_output_config():
    if OutputConfig:
        return OutputConfig()
    return SimpleNamespace(
        write_logfile=True,
        write_results=True,
        record_k_space=True,
        print_k_space=False,
        write_xyz_trajectory=True,
        write_histogram=True,
        write_phase_space=True,
        write_midpoint_snapshots=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("idx", type=str, nargs="?", default="0")
    parser.add_argument("--workdir", type=str, default=None, help="Run directory containing input files")
    parser.add_argument("--base", type=str, default=None, help="Base directory containing server_hostname.txt")
    parser.add_argument("--config", type=str, default=None, help="Path to cavotf input.txt")
    args = parser.parse_args()

    workdir = pathlib.Path(args.workdir) if args.workdir else pathlib.Path.cwd()
    os.chdir(workdir)

    base_dir = pathlib.Path(args.base) if args.base else workdir.parent
    with open(base_dir / "server_hostname.txt", "r") as fob:
        host = fob.readline().strip()
        port = int(fob.readline().strip())

    idx = args.idx
    output_cfg = _default_output_config()
    params = param()

    derivative_interval = 5
    if args.config and load_config and _recompute_mode_grid:
        try:
            cfg = load_config(pathlib.Path(args.config))
            apply_config_overrides(params, cfg)
            derivative_interval = cfg.physics.dipole_derivative_interval
            output_cfg = cfg.outputs
        except Exception as exc:  # noqa: BLE001
            print(f"Warning: failed to apply config overrides: {exc}")

    time.sleep(30)

    print(r"""
  ▄▄             █             ▗▄▖ ▗▄▄▄▖▗▄▄▄▖
 █▀▀▌            ▀   ▐▌        █▀█ ▝▀█▀▘▐▛▀▀▘
▐▛    ▟██▖▐▙ ▟▌ ██  ▐███ ▝█ █▌▐▌ ▐▌  █  ▐▌        ▐▙█▙ ▝█ █▌
▐▌    ▘▄▟▌ █ █   █   ▐▌   █▖█ ▐▌ ▐▌  █  ▐███      ▐▛ ▜▌ █▖█
▐▙   ▗█▀▜▌ ▜▄▛   █   ▐▌   ▐█▛ ▐▌ ▐▌  █  ▐▌        ▐▌ ▐▌ ▐█▛
 █▄▄▌▐▙▄█▌ ▐█▌ ▗▄█▄▖ ▐▙▄   █▌  █▄█   █  ▐▌     █  ▐█▄█▘  █▌
  ▀▀  ▀▀▝▘  ▀  ▝▀▀▀▘  ▀▀   █   ▝▀▘   ▀  ▝▘     ▀  ▐▌▀▘   █
                          █▌                      ▐▌    █▌
       Mandal Group, TAMU
""")
    print("══════════════════════════════════════════════════════")
    print("   CALCULATIONS EXECUTED BY AMIR H. AMINI (Nov–Dec 2025)")
    print("   Texas A&M University | College Station, TX")
    print("══════════════════════════════════════════════════════")
    print("   IMPORTANT NOTICE")
    print("   These results should be used with extreme caution.")
    print("   This code is experimental and not a final, polished release.")
    print("   If reused for independent simulations, ensure it is thoroughly")
    print("   reviewed, debugged, and validated prior to application.")
    print("══════════════════════════════════════════════════════")
    print("   INITIATING COMPUTATIONAL ROUTINE...")

    t0 = time.time()
    bhr = 1.8897259886
    evdivA = 27.2114 * bhr
    AngdivPs2AU = bhr / 41341.3733365614
    natoms = params.natoms
    steps = params.steps
    dt = params.dt
    dt2 = dt / 2
    thermal_steps = params.thermal_steps

    atm = "O33H66"  # Legacy default
    coordina_initial = np.loadtxt("initXYZ.dat", usecols=(2, 3, 4))
    velocity_initial = np.loadtxt("initPxPyPz.dat", usecols=(0, 1, 2))
    coordinates = np.array(coordina_initial) * bhr

    velocity = np.array(velocity_initial) * AngdivPs2AU
    vx, vy, vz = velocity[:, 0], velocity[:, 1], velocity[:, 2]

    atoms = Atoms(atm, positions=coordinates / bhr)
    mass = atoms.get_masses() * 1822.8884
    box = params.box
    atoms.set_cell([box, box, box])
    atoms.set_pbc(True)
    px, py, pz = vx * mass, vy * mass, vz * mass

    coordinates = atoms.get_positions(wrap=False) * bhr
    atoms.set_positions(coordinates / bhr)

    if output_cfg.write_xyz_trajectory:
        write("WaterMD_Cavity.xyz", atoms, format="xyz", append=True)

    pj = np.concatenate((px, py, pz))
    masses = np.concatenate((mass, mass, mass))
    rxj, ryj, rzj = coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]
    rj = np.concatenate((rxj, ryj, rzj))

    fj, charges = getForcesCharges(rj, natoms, atm, box)
    μj = np.sum(charges * rxj)

    xk, pk = init(μj, params)  # noqa: F405
    x_save = np.zeros((int(thermal_steps * 0.1)))

    if output_cfg.write_histogram:
        data = np.histogram(x_save, bins=100)
        np.savetxt("hist_coupled.txt", np.c_[data[1][1:], data[0]])
    print("Thermalization done")

    xk, pk = np.loadtxt("initial.dat").T
    xk = np.array([xk])
    pk = np.array([pk])

    pk += dpk(xk, μj, params) * dt2  # noqa: F405
    pk = pk * np.cos(params.ωc * dt2) - params.ωc * xk * np.sin(params.ωc * dt2)

    f = None
    if output_cfg.write_phase_space:
        f = open("qt.out", "w")
    x0 = -(2 / params.ωc) * μj * params.ηb
    dµ = getdµ(natoms, rj, μj, atm, box, dr=0.01)

    output_format = "{0: >5d} {1: >#016.8f} {2: >#016.8f} {3: >#016.8f} {4: >#016.8f} {5: >#016.8f} {6: >#016.8f}"
    fjt = dpj(xk, fj[:natoms], dµ, μj, params)  # noqa: F405
    fxt = dpk(xk, μj, params)  # noqa: F405
    Tk = np.sum(pj**2 / (2 * masses))
    if output_cfg.write_phase_space:
        print(output_format.format(0, xk[0], pk[0], np.sum(μj), fxt, fjt[0], Tk), file=f)

    def andersen_thermostat(Px, Py, Pz, mass, β, timestep, N_atoms):
        N = len(Px)
        mass = np.asarray(mass)
        Px_new, Py_new, Pz_new = Px.copy(), Py.copy(), Pz.copy()
        reassign = np.zeros(N, dtype=bool)
        selected = np.random.choice(N, size=N_atoms, replace=False)
        reassign[selected] = True
        std_dev = np.sqrt(mass[reassign] / β)
        Px_new[reassign] = np.random.normal(0, std_dev)
        Py_new[reassign] = np.random.normal(0, std_dev)
        Pz_new[reassign] = np.random.normal(0, std_dev)
        return Px_new, Py_new, Pz_new

    def calculation(rj, pj, xk, pk, fj, μj, dµ, f, params, output_cfg, i):
        pj[:natoms] += dpj(xk, fj[:natoms], dµ, μj, params) * dt2  # noqa: F405
        pj[natoms:3 * natoms] += fj[natoms:3 * natoms] * dt2
        rj += pj * dt / masses
        pk += dpk(xk, μj, params) * dt2  # noqa: F405

        fj, charges = getForcesCharges(rj, natoms, atm, box)
        Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)
        μj = np.sum(charges * (rj[:natoms] - Rcom))

        if i % derivative_interval == 0:
            dµ = getdµ(natoms, rj, μj, atm, box, dr=0.01)

        fjt = dpj(xk, fj[:natoms], dµ, μj, params)  # noqa: F405
        fxt = dpk(xk, μj, params)  # noqa: F405

        pj[:natoms] += fjt * dt2
        pj[natoms:3 * natoms] += fj[natoms:3 * natoms] * dt2
        pk += fxt * dt2

        if output_cfg.write_phase_space and i % 2 == 0:
            if f:
                f.close()
            f = open("qt.out", "a")

        if output_cfg.write_midpoint_snapshots and i % 2 == 0:
            with open("midpoint.dat", "w") as f2:
                print(i + 1, xk, pk, rj, pj, file=f2)

        Tk = np.sum(pj**2 / (2 * masses))
        if output_cfg.write_phase_space:
            print(output_format.format((i + 1), xk[0], pk[0], np.sum(μj), fxt, fjt[0], Tk), file=f)

        coordinates = np.column_stack((rj[:natoms], rj[natoms:2 * natoms], rj[2 * natoms:3 * natoms]))
        atoms.set_positions(coordinates / bhr)
        coordinates = atoms.get_positions(wrap=False) * bhr
        rj = np.concatenate((coordinates[:, 0], coordinates[:, 1], coordinates[:, 2]))

        if output_cfg.write_xyz_trajectory and i % derivative_interval == 0:
            rxj = rj[:natoms]
            ryj = rj[natoms:2 * natoms]
            rzj = rj[2 * natoms:3 * natoms]
            coordinates = np.column_stack((rxj, ryj, rzj))
            atoms.set_positions(coordinates / bhr)
            atoms.set_positions(atoms.get_positions(wrap=False))
            write("WaterMD_Cavity.xyz", atoms, format="xyz", append=True)
            print("Step xyz ", i + 1)

        return rj, pj, xk, pk, fj, μj, dµ, f

    q = xk[0]
    p = pk[0]
    sleeptime = 1.5
    loglevel = 2

    for i in range(steps):
        dat = {"q": q, "p": p, "idx": idx, "killed": False, "step": i}
        reply = json.loads(comm(dat, host, port))
        globalStep = reply["step"]

        while reply["step"] < i:
            time.sleep(sleeptime)
            if loglevel >= 1:
                print(f" Waiting..... (Server {globalStep} | Client {i})")
            reply = json.loads(comm(dat, host, port))
            globalStep = reply["step"]

        q = reply["q"]
        p = reply["p"]
        xk[0] = q
        pk[0] = p

        if loglevel >= 1:
            print(q, f"Step {i + 1} of {steps}")

        rj, pj, xk, pk, fj, μj, dµ, f = calculation(
            rj,
            pj,
            xk,
            pk,
            fj,
            μj,
            dµ,
            f,
            params,
            output_cfg,
            i,
        )
        dat["step"] = i + 1
        dat["q"] = xk[0]
        dat["p"] = pk[0]
        reply = json.loads(comm(dat, host, port))
        globalStep = reply["step"]

    dat["killed"] = True
    dat["step"] = steps
    dat["q"] = xk[0]
    dat["p"] = pk[0]
    comm(dat, host, port)

    if output_cfg.write_phase_space and f:
        f.close()


if __name__ == "__main__":
    main()
