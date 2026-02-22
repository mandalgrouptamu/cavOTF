# =============================================================================
#  Project:     cavOTF.py
#  File:        dftb.py
#  Author:      Sachith Wickramasinghe
#  Last update: 02/02/2026
#
#  Description:
#  parameters and functions.
# =============================================================================

import numpy as np
from numpy.random import normal as gran
from types import SimpleNamespace

def dpk(x, µ, par, idx, t=0):
    ηb = par.ηb
    ωc = par.ωc
    idx = int(idx)
    
    return  - ηb * µ - par.gl[idx]*np.sin(par.ωl * t )

def dpkT(x, µ, par):
    ηb = par.ηb
    ωc = par.ωc
    return -ωc**2 * x  - ηb * µ



def dpj(x, fj, dµ, μ, par, idx, t=0):
    ηb = par.ηb
    ωc = par.ωc
    idx = int(idx)

    return fj - ηb * dµ * x - (ηb**2 * dµ * µ /ωc**2) - (dµ * par.glm[idx]*np.sin(par.ωlm * t))

def vvl(x, p, µ, param, f1): #only for 1 cavity
    ndof = 1
    β  = param.β
    v = p/param.m
    dt = param.dt
    λ = param.λ #/ param.m
    σ = (2.0 * λ/(β * param.m )) ** 0.5
    ξ = gran(0, 1, ndof)  #np.array([0.5 * gran() for i in range(len(x))])
    θ = gran(0, 1, ndof) #np.array([gran() * 0.28867513459  for i in range(len(x))])
    c = 0.28867513459
    A = (0.5 * dt**2) * (f1/param.m - λ * v) + (σ * dt**(3.0/2.0)) * (0.5 * ξ + c * θ) 
    #---- X update -----------
    
    x += (v * dt + A) 
    #-------------------------
    f2 = dpkT(x * 1.0, µ, param)
    #---- V update ----------- 
    v += ( 0.5 * dt * (f1+f2)/param.m - dt * λ * v +  σ * (dt**0.5) * ξ - A * λ ) 
    #-------------------------
    return x, v * param.m, f2


def init(μ, param): # initialize the xk, pk
    β = param.β
    ωc = param.ωc

    σp = (1/β)**0.5
    σK = σp/ωc
    x0 = - (1/ωc**2) * μ * param.ηb 
  

    #-------- define initial positions and momennam for cavity ----------
    xk = np.random.normal(0,σK)  
    pk = np.random.normal(0,σp)

    return xk + x0, pk 

def default_physics():
    return SimpleNamespace(
        omega_c=0.19,    
        beta=1052.8,
        lambda_=0.001,
        nk=81,
        omega_l=0.0,
        gl_val=0.0,
        gl_n_active=0,
        omega_lm=0.0,
        gl_valm=0.0,
        gl_n_activem=0,
        eta_b=0.0002,
    )


class param:
    def __init__(self, physics = None):
        if physics is None:
            physics = default_physics()
        self.ωc = physics.omega_c
        self.ω0 = self.ωc
        self.β  = physics.beta #* (300.0/200.0)
        self.λ = physics.lambda_
        self.natoms = 99 # number of atoms in the MD simulation
        self.box = 10.0 # size of the periodic box
        self.c = 137.0 

        totaltime    = 6000   # total simulation time in fs
        thermal_time = 50000   # thermalization time in fs
        dt           = 0.3    #time step in fs

        self.steps = int(totaltime // dt)
        self.dt = dt * 41.3413733365614
        self.thermal_steps = int(thermal_time // dt)

        Lx = 200000 * 4
        self.nk = physics.nk
        self.dL = Lx/self.nk
        
        self.ωl = physics.omega_l
        gl_val_Ha = physics.gl_val / 27.2114
        gl = np.zeros(self.nk)
        gl[:min(physics.gl_n_active, self.nk)] = gl_val_Ha
        self.gl = gl
        
        self.ωlm = physics.omega_lm
        gl_val_Ham = physics.gl_valm / 27.2114
        glm = np.zeros(self.nk)
        gl[:min(physics.gl_n_activem, self.nk)] = gl_val_Ham
        self.glm = glm
        
        
        ky = np.fft.fftfreq(self.nk) * self.nk * 2 * np.pi / (Lx)
        self.ωk = np.sqrt(self.ω0**2 + (self.c * ky)**2)
        self.ωk[self.ωk> 5 * self.ω0] = 5 * self.ω0
        self.m = 1.0
          
        self.ηb = physics.eta_b  #0.007 


