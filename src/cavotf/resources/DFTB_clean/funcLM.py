# funcLM.py
# Functions, Classes and Parameters

import numpy as np
from numpy.random import normal as gran

def dpk(x, µ, par):
    ηb = par.ηb
    ωc = par.ωc
    return  - ηb * µ

def dpkT(x, µ, par):
    ηb = par.ηb
    ωc = par.ωc
    return -ωc ** 2 * x  - ηb * µ

def dpj(x, fj, dµ, μ, par):
    ηb = par.ηb
    ωc = par.ωc
    return fj - ηb * dµ * x - (ηb ** 2 * dµ * µ / ωc ** 2)

def vvl(x, p, µ, param, f1): 
    ndof = 1
    β  = param.β
    v = p/param.m
    dt = param.dt
    λ = param.λ
    σ = (2.0 * λ / (β * param.m )) ** 0.5
    ξ = gran(0, 1, ndof)
    θ = gran(0, 1, ndof)
    c = 0.28867513459
    A = (0.5 * dt ** 2) * (f1 / param.m - λ * v) + (σ * dt ** (3.0/2.0)) * (0.5 * ξ + c * θ) 
    # ===== position update ======
    x += (v * dt + A) 
    f2 = dpkT(x * 1.0, µ, param)
    # ===== velocity update ======
    v += ( 0.5 * dt * (f1 + f2) / param.m - dt * λ * v + σ * (dt ** 0.5) * ξ - A * λ ) 
    return x, v * param.m, f2

def init(μ, param):
    β = param.β
    ωc = param.ωc

    σp = (1 / β) ** 0.5
    σK = σp / ωc
    x0 = - (2 / ωc) * μ * param.ηb

    xk = np.random.normal(0, σK)
    pk = np.random.normal(0, σp)

    return xk + x0, pk

class param:
    def __init__(self, ωc = 0.190/27.2114):
        self.ωc = ωc
        self.ω0 = ωc
        self.β  = 1052.8
        self.λ = 0.001
        self.natoms = 99 # No. of atoms
        self.box = 10.0 # Size of the periodic box in Angstrom
        self.c = 137.0

        totaltime    = 1000   # Total simulation time in fs
        thermal_time = 10000   # Total thermalization time in fs
        dt           = 0.5    # Time step in fs

        self.steps = int(totaltime // dt)
        self.dt = dt * 41.3413733365614
        self.thermal_steps = int(thermal_time // dt)

        self.dL = 20 # Box length in Angstrom
        self.nk = 25 # No. of boxes

        self.ky = np.fft.fftfreq(self.nk) * self.nk * 2 * np.pi / (self.dL * self.nk)
        self.ωk = np.sqrt(self.ωc**2 + (self.c * self.ky)**2)
 
        self.m = 1.0
          
        self.ηb = 0.0003  # Light-matter coupling strength


