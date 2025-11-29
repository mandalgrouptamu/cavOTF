# =============================================================================
#  Project:     cavOTF.py
#  File:        dftb.py
#  Author:      Sachith Wickramasinghe
#  Last update: 11/28/2025
#
#  Description:
#  parameters and functions.
# =============================================================================

import numpy as np
from numpy.random import normal as gran

def dpk(x, µ, par):
    ηb = par.ηb
    ωc = par.ωc
    return  - ηb * µ

def dpkT(x, µ, par):
    ηb = par.ηb
    ωc = par.ωc
    return -ωc**2 * x  - ηb * µ



def dpj(x, fj, dµ, μ, par):
    ηb = par.ηb
    ωc = par.ωc

    return fj - ηb * dµ * x - (ηb**2 * dµ * µ /ωc**2) 

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



class param:
    def __init__(self, ωc = 0.190/27.2114):
        self.ωc = ωc
        self.ω0 = ωc
        self.β  = 1052.8 #* (300.0/200.0)
        self.λ = 0.001
        self.natoms = 99 # number of atoms in the MD simulation
        self.box = 10.0 # size of the periodic box
        self.c = 137.0 

        totaltime    = 3000   # total simulation time in fs
        thermal_time = 50000   # thermalization time in fs
        dt           = 0.3    #time step in fs

        self.steps = int(totaltime // dt)
        self.dt = dt * 41.3413733365614
        self.thermal_steps = int(thermal_time // dt)

        Lx = 200000 * 4
        self.nk = 81
        self.dL = Lx/self.nk
 

        ky = np.fft.fftfreq(self.nk) * self.nk * 2 * np.pi / (Lx)
        self.ωk = np.sqrt(self.ω0**2 + (self.c * ky)**2)
        self.ωk[self.ωk> 5 * self.ω0] = 5 * self.ω0
        self.m = 1.0
          
        self.ηb = 0.0002  #0.007 


