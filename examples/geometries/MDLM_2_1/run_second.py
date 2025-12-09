import os 
import sys
import numpy as np
import time
import random
from operator import attrgetter
from numpy.random import normal as gran
from scipy.fft import fft, ifft

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'
sys.path.append(os.popen("pwd").read().replace("\n","")+"/DFTB_clean")
from funcLM import param
param = param()  # Initialize parameters
N = param.nk
foldprefix = "run-"
cleanFolder = "DFTB_clean"

µj = np.zeros(N)


for i in range(N):

    os.chdir(foldprefix + str(i))  # Change to the client's directory
    µj[i] = np.loadtxt('dmu.dat')  # Load the computed µj value
    os.chdir("..")  # Change back to the parent directory



def vvl(q, p, μj, param): # all degrees of freedom
    ndof = len(q)
    β  = param.β
    λ = param.λ #/ param.m
    σ = (2.0 * λ/(β * param.m )) ** 0.5
    ξ = gran(0, 1, ndof)  #np.array([0.5 * gran() for i in range(len(x))])
    θ = gran(0, 1, ndof) #np.array([gran() * 0.28867513459  for i in range(len(x))])
    const = 0.28867513459

    # # bare photon------------------------- 
    dt, nk, ω0, ωk = attrgetter("dt", "nk",  "ω0", "ωk")(param)
    dt = 1.0
    dt2 = dt/2

    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2)/ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0
    q[:nk], p[:nk] = ktoj(qk, pk, ω0, ωk)
    #--------------------------------------


    f1 = pdot(q * 1.0, p * 1, μj,param)

    A = (0.5 * dt**2) * (f1/param.m - λ * p) + (σ * dt**(3.0/2.0)) * (0.5 * ξ + const * θ)
   
    #---- X update -----------
    q += A 
    #-------------------------
    f2 = pdot(q * 1.0, p * 1, μj, param)
    #---- V update ----------- 
    p += ( 0.5 * dt * (f1+f2)/param.m - dt * λ * p +  σ * (dt**0.5) * ξ - A * λ ) 
 
    #-------------------------

    # # bare photon------------------------- 
    dt, nk, ω0, ωk = attrgetter("dt", "nk",  "ω0", "ωk")(param)
    dt2 = dt/2
    qk, pk = jtok(q[:nk], p[:nk], ω0, ωk)
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2)/ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0
    q[:nk], p[:nk] = ktoj(qk, pk, ω0, ωk)
    #--------------------------------------

    return q, p


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

def pdot(q,p, muj,param):
    """Calculates the force from the Hamiltonian."""
    dp = q * 0.0
    # photon digrees of freedom
    dp[:] = - muj*param.ηb # light-matter coupling
    return dp


def qdot(q, p, param):
    """Calculates the force from the Hamiltonian."""
    pr = p * 1.0
    pr[:param.nk] = 0.0
    return pr

def init(µj, param): # initialize the xk, pk
    β = param.β
    ωc = param.ωc
    ωk = param.ωk
    nk = param.nk

    σp = (1/β)**0.5
    σK = σp/ωc
    x0 = - (2/ωc) * µj * param.ηb 
  

    #-------- define initial positions and momennam for cavity ----------
    xk = np.random.normal(0,σK,nk)  
    pk = np.random.normal(0,σp,nk)

    return xk + x0, pk

q, p = init(µj, param)

for t in range(param.thermal_steps*50):

    q, p = vvl(q, p, µj, param)
    
print("q", q)
print("p", p)

for i in range(N):
    os.chdir(foldprefix + str(i))  # Change to the client's directory
    np.savetxt('initial.dat', np.c_[q[i], p[i]])
    os.chdir("..")  # Change back to the parent directory
