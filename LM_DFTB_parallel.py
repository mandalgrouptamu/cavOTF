#!/sw/eb/sw/Anaconda3/2024.02-1/bin/python

##NECESSARY JOB SPECIFICATIONS                                                                                                                        
#SBATCH --partition=shared                                                                                                                            
#SBATCH --job-name=MD         #Set the job name to "JobExample1"                                                                                    
#SBATCH --time=0-20:20:00          #Set the wall clock limit to 1hr and 30min                                                                         
#SBATCH --ntasks=1               #Request 1 task                                                                                                      
#SBATCH --nodes=1                                                                                                                                     
#SBATCH --cpus-per-task=1                                                                                                                             
#SBATCH --mem=2000M                                                                                                                                                                                                                                                                      
#SBATCH --output=output/out.%j    

import numpy as np
from operator import attrgetter
import time 
import os
import sys
from numpy.random import normal as gran
from scipy.fft import fft, dct, ifft
#from ftt import k_to_j, j_to_k
pwd = os.getcwd()
sys.path.append(pwd)

class param:
    def __init__(self, ηb = 0.001):

        self.c = 137.0                           # Speed of light in atomic units
        self.ω0 = 0.10 / 27.2114                 # Base cavity frequency
        self.nk = 50000                          # Number of wave vectors
        self.nkPrint = 100
        self.nkl = 100                           # Number of k modes subjest to langevin dynamics
        self.dL = 200                            # Spacing in j direction
        self.ηb = ηb                             # Dissipation coefficient
        self.dt = 1                             # Time step
        self.steps =  int(1000 /self.dt)       # Number of time steps
        self.tsteps = 500 #150000                # Number of thermalization time steps
        self.ntraj = 1                           # Number of trajectories
        
        self.T = 298.0/1000                      # Temperature in kelvin
        self.β = 315774/self.T

        self.m = 1.0                             # Mass
        assert self.m == 1.0

        # Frequency and wave vector calculations
        self.nj = self.nk                        # Number of boxes
        self.ky = np.fft.fftfreq(self.nk) * self.nk * 2 * np.pi / (self.dL * self.nk)
        self.ωk = np.sqrt(self.ω0**2 + (self.c * self.ky)**2)

        self.ωkmax = 5 * self.ω0
        self.ωk = np.minimum(self.ωk, self.ωkmax)

        self.Rj = np.arange(self.nj) * self.dL

    
        self.nm = len(self.Cn)                        # Number of matter in a box
        self.ndof = (self.nj * self.nm ) + self.nk    # Number of matter degrees of freedom
        self.nmdof = self.nj * self.nm                # Number of matter degrees of freedom
        

        self.λ = np.zeros(self.ndof) + 0.0005#2  # λ for the Langevin dynamics

        self.nskip = 1

def jtok(qj, pj, ω, ωk):
    an = np.sqrt(ω/2) * (qj + 1j*pj/ω)
    ak = fft(an, norm='ortho')
    akd =  np.conj(ak)
    qk = (ak + akd)/np.sqrt(2*ωk)
    pk = -1j*(ak - akd)*np.sqrt(ωk/2)
     
    return qk.real, pk.real


def ktoj(qk, pk, ω, ωk):
    
    ak = np.sqrt(ωk/2) * (qk + 1j*pk/ωk)

    
    aj = ifft(ak, norm='ortho')
    ajd =  np.conj(aj)
    qj = (aj + ajd)/np.sqrt(2*ω)
    pj = -1j*(aj - ajd)*np.sqrt(ω/2)
    return qj.real, pj.real

def init(param):
    """Initializes the cavity."""
    β, ωk, ndof, nk  = attrgetter("β", "ωk", "ndof", "nk")(param)

    σp = (1/β)**0.5
    σK = σp/ωk

    q, p = np.zeros(nk), np.zeros(nk)   # initializing q and p for the cavity modes
    
    q = np.random.normal(0, σK, nk)     # initializing q for the cavity modes
    p = np.random.normal(0, σp, ndof)   # initializing p for the cavity modes
    return q, p

def vvlKC(q, p, param): #only for 1 cavity
    dt, nkl, ω0, ωk , λ = attrgetter("dt", "nkl",  "ω0", "ωk", "λ")(param)
    
    ndof = nkl
    β  = param.β
    dt2 = dt/2
    λ = λ[:nkl] #/ param.m
    
    σ = (2.0 * λ/(β * param.m )) ** 0.5
    ξ = gran(0, 1, ndof)  #np.array([0.5 * gran() for i in range(len(x))])
    θ = gran(0, 1, ndof) #np.array([gran() * 0.28867513459  for i in range(len(x))])
    c = 0.28867513459

    qk, pk = jtok(q[:], p[:], ω0, ωk)

    f1 = -ωk[:nkl] * qk[:nkl]

    A = (0.5 * dt2**2) * (f1 - λ * p[:nkl]) + (σ * dt2**(3.0/2.0)) * (0.5 * ξ + c * θ) 
    #---- q update -----------
    qk[:nkl] += (pk[:nkl] * dt2 + A) 
    #-------------------------
    f2 = -ωk[:nkl] * qk[:nkl]
    #---- p update ----------- 
    pk[:nkl] += ( 0.5 * dt2 * (f1+f2) - dt2 * λ * pk[:nkl] +  σ * (dt**0.5) * ξ - A * λ ) 
    #-------------------------

    
    qk1 = qk[nkl:] * np.cos(ωk[nkl:] * dt2) + pk[nkl:] * np.sin(ωk[nkl:] * dt2)/ωk[nkl:]
    pk1 = pk[nkl:] * np.cos(ωk[nkl:] * dt2) - ωk[nkl:] * qk[nkl:] * np.sin(ωk[nkl:] * dt2)
    qk[nkl:], pk[nkl:] = qk1 * 1.0, pk1 * 1.0
    q[:], p[:] = ktoj(qk, pk, ω0, ωk)


    return q, p

def vv2(q, p, param):
    dt, nk = attrgetter("dt", "nk")(param)

    # bare photon------------------------- 
    q[:nk], p[:nk] = vvlKC(q[:nk], p[:nk], param)

    #---- vv -----------
    #client update this part

 
    # bare photon------------------------- 
    q[:nk], p[:nk] = vvlKC(q[:nk], p[:nk], param)

    return q, p