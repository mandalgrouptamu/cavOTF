# client_DFTB.py
# Client code 

import os
import sys
import time
import json
import socket

import numpy as np
from ase.io import write
from ase import Atoms

from dftb import getForcesCharges, getdµ
from funcLM import *

# Read server hostname and port from file

fob = open("../server_hostname.txt", "r")
HOST = fob.readline().strip()
PORT = int(fob.readline().strip())

idx = sys.argv[1] if len(sys.argv) > 1 else "0"

def comm(msg):
    msg['execTime'] = time.time()
    msg = json.dumps(msg)
    with socket.socket() as cli:
        cli.connect((HOST, PORT))
        cli.sendall(msg.encode())
        reply = cli.recv(4096)
    return reply.decode().strip()

# === Parameters and initializations ===

t0 = time.time()
bhr = 1.8897259886
evdivA = 27.2114 * bhr
AngdivPs2AU = bhr / 41341.3733365614 
params = param()

natoms = params.natoms
steps = params.steps
dt = params.dt
dt2 = dt/2
thermal_steps = params.thermal_steps

# Path to DFTB+ SK files
# e.g., mio-1-1

os.system('export DFTB_PREFIX=path_to_your_DFTB_SK_files') 
os.environ['DFTB_PREFIX'] = 'path_to_your_DFTB_SK_files'  

atm = 'atomic_structure' # Define the atomic structure, e.g., 'O33H66' for 33 water molecules

# Load the atomic coordinates from the 'finalstep_InTheBox.xyz' file
coordina_initial = np.loadtxt('initXYZ.dat', usecols=(2,3,4))
velocity_initial = np.loadtxt('initPxPyPz.dat', usecols=(0,1,2))
coordinates = np.array(coordina_initial) * bhr #a.u.

# --- velocity
velocity = np.array(velocity_initial) * AngdivPs2AU
vx, vy, vz  = velocity[:,0], velocity[:,1], velocity[:,2]
#--------------

atoms = Atoms(atm,  positions =   coordinates /bhr) # atomic structure
mass = atoms.get_masses()  * 1822.8884
box = params.box
atoms.set_cell([box, box, box])
atoms.set_pbc(True)
px, py, pz = vx * mass, vy * mass, vz * mass

coordinates = atoms.get_positions(wrap=False) * bhr
atoms.set_positions(coordinates/bhr)

write('coordinates.xyz',atoms,format='xyz', append = True)

pj = np.concatenate((px, py, pz))
masses = np.concatenate((mass, mass, mass))
rxj, ryj, rzj = coordinates[:,0], coordinates[:,1], coordinates[:,2]
rj = np.concatenate((rxj, ryj, rzj))
Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)

fj, charges = getForcesCharges(rj, natoms, atm, box)
µj = np.sum(charges * rxj)

xk, pk = init(µj, params)
x_save = np.zeros((int(thermal_steps * 0.1)))    

        
data = np.histogram(x_save, bins=100)
np.savetxt('hist_coupled.txt', np.c_[data[1][1:], data[0]])
print('Thermalization done')

xk, pk = np.loadtxt('initial.dat').T
xk = np.array([xk]) 
pk = np.array([pk])

#-----for Leapfrog -----------------
pk += dpk(xk, μj, params) * dt2
pk = pk * np.cos(params.ωc * dt2) - params.ωc * xk * np.sin(params.ωc * dt2)

f = open('qt.out', 'w') 
x0 = - (2/params.ωc) * μj * params.ηb 
dµ = getdµ(natoms, rj, atm, box, dr=0.0001)


output_format = '{0: >5d} {1: >#016.8f} {2: >#016.8f} {3: >#016.8f} {4: >#016.8f} {5: >#016.8f} {6: >#016.8f}'
fjt = dpj(xk, fj[:natoms], dµ, μj, params)
fxt = dpk(xk, μj, params)
Tk =np.sum(pj**2 / (2 * masses))
print(output_format.format(0,xk[0],pk[0], np.sum(μj), fxt, fjt[0], Tk), file=f)


tstep0 = time.time()

#------------------------------------------------------
#------------------------------------------------------
def andersen_thermostat(Px, Py, Pz, mass, β, timestep, collision_freq):
    N = len(Px)
    mass = np.asarray(mass)
    Px_new, Py_new, Pz_new = Px.copy(), Py.copy(), Pz.copy()
    collision_prob = 1 - np.exp(-collision_freq * timestep)

    #reassign = np.random.rand(N) < collision_prob

    reassign = np.zeros(N, dtype=bool)
    selected = np.random.choice(N, size=16, replace=False)
    reassign[selected] = True


    std_dev = np.sqrt(mass[reassign] / β )

    Px_new[reassign] = np.random.normal(0, std_dev)
    Py_new[reassign] = np.random.normal(0, std_dev)
    Pz_new[reassign] = np.random.normal(0, std_dev)

    return Px_new, Py_new, Pz_new


#------------------------------------------------------
#------------------------------------------------------



def calculation(rj, pj, xk, pk, fj, μj, dµ, f, params,i):
    # Update the momenta vv1
    pj[:natoms] += dpj(xk, fj[:natoms], dµ, μj, params) * dt2
    pj[natoms:3*natoms] += fj[natoms:3*natoms] * dt2

    # Update the positions vv2
    rj += pj * dt / masses

    # Update the cavity momenta  vv3
    pk += dpk(xk, μj, params) * dt2

    # get forces and charges
    fj, charges = getForcesCharges(rj, natoms, atm, box)  
    Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)
    μj = np.sum(charges * (rj[:natoms] - Rcom))

    if i % 5 == 0:
        dµ = getdµ(natoms, rj, atm, box, dr=0.0001)

    # Update the momenta vv4
    fjt = dpj(xk, fj[:natoms], dµ, μj, params)
    fxt = dpk(xk, μj, params)
    
    pj[:natoms] += fjt * dt2
    pj[natoms:3*natoms] += fj[natoms:3*natoms] * dt2


    pk += fxt * dt2

    if i < 250:
        # Apply Andersen thermostat
        collision_freq = 0.001
        pj[:natoms], pj[natoms:2*natoms], pj[2*natoms:3*natoms] = andersen_thermostat(
            pj[:natoms], pj[natoms:2*natoms], pj[2*natoms:3*natoms], mass, params.β, dt, collision_freq)

    if i % 2 == 0:
        f.close()
        f = open('qt.out' , 'a') # Open the file in write mode

    if i % 2 == 0:
        f2 = open('midpoint.dat' , 'w') # Open the file in write mode
        print(i+1,xk,pk,rj,pj, file=f2)
        f2.close()

    Tk = np.sum(pj**2 / (2 * masses))
    print(output_format.format((i+1),xk[0], pk[0], np.sum(μj), fxt, fjt[0] , Tk), file=f)

    coordinates = np.column_stack((rj[:natoms], rj[natoms:2*natoms], rj[2*natoms:3*natoms]))
    atoms.set_positions(coordinates/bhr)
    coordinates = atoms.get_positions(wrap=False) * bhr
    rj = np.concatenate((coordinates[:,0], coordinates[:,1], coordinates[:,2]))

    if i % 5 == 0:
        rxj = rj[:natoms]
        ryj = rj[natoms:2*natoms]
        rzj = rj[2*natoms:3*natoms]
    
        coordinates = np.column_stack((rxj, ryj, rzj))
        atoms.set_positions(coordinates/bhr)
        atoms.set_positions(atoms.get_positions(wrap=False))
        write('WaterMD_Cavity.xyz', atoms, format='xyz', append=True)
        print("Step xyz ", i+1 )
    return rj, pj, xk, pk, fj, μj, dµ, f



# --- main loop ---
q = xk[0]
p = pk[0]
steps = params.steps
sleeptime = 1.5
param = {}
loglevel = 2

for i in range(steps):
    # 1) announce arrival at step i
    dat = {'q': q, 'p': p, 'idx': idx, 'killed': False, 'step': i}
    reply = json.loads(comm(dat))
    globalStep = reply['step']

    # 2) barrier: wait until server.step > i
    while reply['step'] < i:
        time.sleep(sleeptime)
        if loglevel >= 1:
            print(f" Waiting..... (Server {globalStep} | Client {i})")
        reply = json.loads(comm(dat))
        globalStep = reply['step']

    # 3) barrier passed: read updated q
    q = reply['q']
    p = reply['p']
    xk[0] = q  # update xk
    pk[0] = p  # update pk

    if loglevel >= 1:
        print(q, f"Step {i+1} of {steps}")

    # 4) compute next q and send it for step i+1
    # q = client(q, param)
    
    rj, pj, xk, pk, fj, μj, dµ, f = calculation(rj, pj, xk, pk, fj, μj, dµ, f, params,i)
    dat['step'] = i + 1
    dat['q'] = xk[0]  # update q
    dat['p'] = pk[0]  # update p
    reply = json.loads(comm(dat))
    globalStep = reply['step']



# 5) tell server we're done
dat['killed'] = True
dat['step'] = steps
dat['q'] = xk[0]
dat['p'] = pk[0]
comm(dat)
