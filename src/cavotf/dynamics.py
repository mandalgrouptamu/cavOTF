from __future__ import annotations

import numpy as np
from numpy.random import normal as gran
from operator import attrgetter
from scipy.fft import fft, ifft


def vvl(q: np.ndarray, p: np.ndarray, mu_j: np.ndarray, param) -> tuple[np.ndarray, np.ndarray]:
    ndof = len(q)
    β = param.β
    λ = param.λ
    σ = (2.0 * λ / (β * param.m)) ** 0.5
    ξ = gran(0, 1, ndof)
    θ = gran(0, 1, ndof)
    const = 0.28867513459

    dt, nk, ω0, ωk = attrgetter("dt", "nk", "ω0", "ωk")(param)
    dt = 1.0
    dt2 = dt / 2

    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0
    q[:nk], p[:nk] = ktoj(qk, pk, ω0, ωk)

    f1 = pdot(q * 1.0, p * 1, mu_j, param)

    a_vec = (0.5 * dt**2) * (f1 / param.m - λ * p) + (σ * dt ** (3.0 / 2.0)) * (0.5 * ξ + const * θ)
    q += a_vec

    f2 = pdot(q * 1.0, p * 1, mu_j, param)
    p += (0.5 * dt * (f1 + f2) / param.m - dt * λ * p + σ * (dt**0.5) * ξ - a_vec * λ)

    dt, nk, ω0, ωk = attrgetter("dt", "nk", "ω0", "ωk")(param)
    dt2 = dt / 2
    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0
    q[:nk], p[:nk] = ktoj(qk, pk, ω0, ωk)

    return q, p


def jtok(qj: np.ndarray, pj: np.ndarray, ω: float, ωk: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    an = np.sqrt(ω / 2) * (qj + 1j * pj / ω)
    ak = fft(an, norm="ortho")
    akd = np.conj(ak)
    qk = (ak + akd) / np.sqrt(2 * ωk)
    pk = -1j * (ak - akd) * np.sqrt(ωk / 2)
    return qk.real, pk.real


def ktoj(qk: np.ndarray, pk: np.ndarray, ω: float, ωk: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ak = np.sqrt(ωk / 2) * (qk + 1j * pk / ωk)
    aj = ifft(ak, norm="ortho")
    ajd = np.conj(aj)
    qj = (aj + ajd) / np.sqrt(2 * ω)
    pj = -1j * (aj - ajd) * np.sqrt(ω / 2)
    return qj.real, pj.real


def pdot(q: np.ndarray, p: np.ndarray, mu_j: np.ndarray, param) -> np.ndarray:
    dp = q * 0.0
    dp[:] = -mu_j * param.ηb
    return dp


def qdot(q: np.ndarray, p: np.ndarray, param) -> np.ndarray:
    pr = p * 1.0
    pr[: param.nk] = 0.0
    return pr


def init(mu_j: np.ndarray, param) -> tuple[np.ndarray, np.ndarray]:
    β = param.β
    ωc = param.ωc
    ωk = param.ωk
    nk = param.nk

    σp = (1 / β) ** 0.5
    σK = σp / ωc
    x0 = -(2 / ωc) * mu_j * param.ηb

    xk = np.random.normal(0, σK, nk)
    pk = np.random.normal(0, σp, nk)

    return xk + x0, pk
