import os
import numpy as np
from numpy.random import normal as gran
from operator import attrgetter
from scipy.fft import fft, ifft

from .config import Config


def vvl(q, p, μj, param):  # all degrees of freedom (unchanged)
    ndof = len(q)
    β = param.β
    λ = param.λ
    σ = (2.0 * λ / (β * param.m)) ** 0.5
    ξ = gran(0, 1, ndof)
    θ = gran(0, 1, ndof)
    const = 0.28867513459

    dt, nk, ω0, ωk = attrgetter("dt", "nk", "ω0", "ωk")(param)
    dt2 = dt / 2

    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    q[:nk], p[:nk] = ktoj(qk1, pk1, ω0, ωk)

    f1 = pdot(q * 1.0, p * 1, μj, param)

    A = (
        0.5 * dt**2 * (f1 / param.m - λ * p)
        + σ * dt ** (3.0 / 2.0) * (0.5 * ξ + const * θ)
    )

    q += A
    f2 = pdot(q * 1.0, p * 1, μj, param)
    p += (
        0.5 * dt * (f1 + f2) / param.m
        - dt * λ * p
        + σ * dt**0.5 * ξ
        - A * λ
    )

    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    q[:nk], p[:nk] = ktoj(qk1, pk1, ω0, ωk)

    return q, p


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


def qdot(q, p, param):
    pr = p * 1.0
    pr[: param.nk] = 0.0
    return pr


def init(µj, param):  # initialize the xk, pk
    β = param.β
    ωc = param.ωc
    nk = param.nk

    σp = (1 / β) ** 0.5
    σK = σp / ωc
    x0 = -(2 / ωc) * µj * param.ηb

    xk = np.random.normal(0, σK, nk)
    pk = np.random.normal(0, σp, nk)

    return xk + x0, pk


def stage2_init_cavity(cfg: Config):
    """
    1. Read µj[i] from each run-i/dmu.dat
    2. Thermalize cavity degrees of freedom
    3. Write initial.dat into each run-i
    """

    param = cfg.param
    N = cfg.physics.nk if cfg.physics.use_param_nk else cfg.general.n_trajectories
    foldprefix = cfg.general.fold_prefix

    µj = np.zeros(N)
    for i in range(N):
        os.chdir(f"{foldprefix}{i}")
        µj[i] = np.loadtxt("dmu.dat")
        os.chdir("..")

    q, p = init(µj, param)

    for t in range(param.thermal_steps):
        q, p = vvl(q, p, µj, param)

    print("q", q)
    print("p", p)

    for i in range(N):
        os.chdir(f"{foldprefix}{i}")
        np.savetxt("initial.dat", np.c_[q[i], p[i]])
        os.chdir("..")
