# =============================================================================
#  Project:     cavOTF.py
#  File:        dynamics.py
#  Author:      Sachith Wickramasinghe
#  Modified by: Amir H. Amini <amiramini@tamu.edu>
#  Last update: 11/28/2025
#
#  Description:
#     Functions to initialize and propagate cavity dynamics.
# =============================================================================
from __future__ import annotations

import importlib
import logging
import sys
from pathlib import Path
from typing import List

import numpy as np
from numpy.random import normal as gran
from scipy.fft import fft, ifft

from .config import Config

LOGGER = logging.getLogger(__name__)


def load_param(clean_template_dir: Path, config: Config):
    # Load parameters
    sys.path.insert(0, str(clean_template_dir))
    funcLM = importlib.import_module("funcLM")
    param_obj = funcLM.param()
    overrides = {
        "nk": config.physics.nk,
        "β": config.physics.beta,
        "λ": config.physics.lambda_,
        "ωc": config.physics.omega_c,
        "ηb": config.physics.eta_b,
        "thermal_steps": config.physics.thermal_steps,
    }
    for key, value in overrides.items():
        if hasattr(param_obj, key):
            setattr(param_obj, key, value)

    _recompute_mode_grid(param_obj)
    return param_obj


def _recompute_mode_grid(param_obj) -> None:

    Lx = 200000 * 4
    param_obj.ω0 = param_obj.ωc
    param_obj.dL = Lx / param_obj.nk
    ky = np.fft.fftfreq(param_obj.nk) * param_obj.nk * 2 * np.pi / Lx
    param_obj.ωk = np.sqrt(param_obj.ω0**2 + (param_obj.c * ky) ** 2)
    param_obj.ωk[param_obj.ωk > 5 * param_obj.ω0] = 5 * param_obj.ω0


def jtok(qj, pj, ω, ωk):
    an = np.sqrt(ω / 2) * (qj + 1j * pj / ω)
    ak = fft(an, norm="ortho")
    akd = np.conj(ak)
    qk = (ak + akd) / np.sqrt(2 * ωk)
    pk = -1j * (ak - akd) * np.sqrt(ωk / 2)
    return qk.real, pk.real


def ktoj(qk, pk, ω, ωk):
    ak = np.sqrt(ωk / 2) * (qk + 1j * pk / ωk)
    aj = ifft(ak, norm="ortho")
    ajd = np.conj(aj)
    qj = (aj + ajd) / np.sqrt(2 * ω)
    pj = -1j * (aj - ajd) * np.sqrt(ω / 2)
    return qj.real, pj.real


def pdot(q, p, muj, param):
    dp = q * 0.0
    dp[:] = -muj * param.ηb
    return dp


def vvl(q, p, μj, param): 
    ndof = len(q)
    β = param.β
    λ = param.λ
    σ = (2.0 * λ / (β * param.m)) ** 0.5
    ξ = gran(0, 1, ndof)
    θ = gran(0, 1, ndof)
    const = 0.28867513459

    dt, nk, ω0, ωk = param.dt, param.nk, param.ω0, param.ωk
    dt2 = dt / 2

    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0
    q[:nk], p[:nk] = ktoj(qk, pk, ω0, ωk)

    f1 = pdot(q * 1.0, p * 1, μj, param)
    A = (0.5 * dt**2) * (f1 / param.m - λ * p) + (σ * dt ** (3.0 / 2.0)) * (0.5 * ξ + const * θ)

    q += A
    f2 = pdot(q * 1.0, p * 1, μj, param)
    p += (0.5 * dt * (f1 + f2) / param.m - dt * λ * p + σ * (dt**0.5) * ξ - A * λ)

    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0
    q[:nk], p[:nk] = ktoj(qk, pk, ω0, ωk)
    return q, p


def init(μj, param):
    β = param.β
    ωc = param.ωc
    ωk = param.ωk
    nk = param.nk

    σp = (1 / β) ** 0.5
    σK = σp / ωc
    x0 = -(2 / ωc) * μj * param.ηb

    xk = np.random.normal(0, σK, nk)
    pk = np.random.normal(0, σp, nk)
    return xk + x0, pk


def initialize_cavity(config: Config, run_dirs: List[Path], dry_run: bool = False) -> None:
    # Initialize cavity dynamics and write initial conditions to each run directory
    # This routine replaces the legacy run_second.py script.
    
    logger = LOGGER.getChild("init")
    param = load_param(config.general.clean_template_dir, config)
    mu_values = []
    for run_dir in run_dirs:
        dmu_file = run_dir / "dmu.dat"
        if not dmu_file.is_file():
            raise FileNotFoundError(f"Missing dipole file: {dmu_file}")
        mu_values.append(float(np.loadtxt(dmu_file)))
    μj = np.array(mu_values)

    q, p = init(μj, param)
    for _ in range(param.thermal_steps):
        q, p = vvl(q, p, μj, param)

    for idx, run_dir in enumerate(run_dirs):
        logger.debug("Writing initial.dat for %s", run_dir)
        if not dry_run:
            np.savetxt(run_dir / "initial.dat", np.c_[q[idx], p[idx]])

