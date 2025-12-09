import os
import numpy as np
#from ase.build import molecule
#from ase.calculators.dftb import Dftb
from ase.io import write
from ase import Atoms
import time
from dftb import getForcesCharges, getCharges, getdµ
from funcLM import *


# client.py
import socket, sys, time, json
import numpy as np





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


#------------ INITIALIZATION ----------------

print('Start')
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

os.system('export DFTB_PREFIX=/home/sachith/vsc-fluxside/dftb_sk_files/mio-1-1/')
os.environ['DFTB_PREFIX'] = '/home/sachith/vsc-fluxside/dftb_sk_files/mio-1-1/'
os.system('export A="/home/sachith/vsc-fluxside/dftb_package/dftbplus-24.1.x86_64-linux/bin/dftb+ "')

os.environ["DFTB_COMMAND"] = "/home/sachith/vsc-fluxside/dftb_package/dftbplus-24.1.x86_64-linux/bin/dftb+"

atm = 'O33H66' # Define the atomic structure
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
#scaled_positions = atoms.get_positions(wrap=True)
px, py, pz = vx * mass, vy * mass, vz * mass

coordinates = atoms.get_positions(wrap=False) * bhr
atoms.set_positions(coordinates/bhr)

write('WaterMD_Cavity.xyz',atoms,format='xyz', append = True)


pj = np.concatenate((px, py, pz))
masses = np.concatenate((mass, mass, mass))
rxj, ryj, rzj = coordinates[:,0], coordinates[:,1], coordinates[:,2]
rj = np.concatenate((rxj, ryj, rzj))
Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)

fj, charges = getForcesCharges(rj, natoms, atm, box)
µj = np.sum(charges * rxj)


# initialize the cavity variables
xk, pk = init(µj, params) # initialize the one-dimensional xk, pk
x_save = np.zeros((int(thermal_steps * 0.1)))    

#------------------------------------------------------
# # Thermalize the cavity
# fk = dpkT(xk, μj, params)
# for ti in range(thermal_steps):
#     xk, pk, fk = vvl(xk, pk, µj, params, fk)
#     if ti > thermal_steps  - len(x_save):
#         x_save[ti - thermal_steps  + len(x_save)] = xk
        
data = np.histogram(x_save, bins=100)
np.savetxt('hist_coupled.txt', np.c_[data[1][1:], data[0]])
print('Thermalization done')
# xk ,pk = -19.13384715389926, 0.0180175741057247781  #BUG: set the cavity position and momentum to a specific value after thermalization
# xk, pk = np.array([xk]), np.array([pk])  # Convert to numpy arrays for consistency
# print('BUG. Cavity position and momentum set to:', xk, pk)
#------------------------------------------------------
xk, pk = np.loadtxt('initial.dat').T
xk = np.array([xk])  # Ensure xk is a 1D array
pk = np.array([pk])  # Ensure pk is a 1D array

#-----for Leapfrog -----------------
pk += dpk(xk, μj, params) * dt2
pk = pk * np.cos(params.ωc * dt2) - params.ωc * xk * np.sin(params.ωc * dt2)

f = open('qt.out', 'w') # Open the file in write mode
x0 = - (2/params.ωc) * μj * params.ηb 
dµ = getdµ(natoms, rj, μj, atm, box, dr=0.0001)


output_format = '{0: >5d} {1: >#016.8f} {2: >#016.8f} {3: >#016.8f} {4: >#016.8f} {5: >#016.8f} {6: >#016.8f}'
fjt = dpj(xk, fj[:natoms], dµ, μj, params)
fxt = dpk(xk, μj, params)
#print(output_format.format(0,xk[0], np.sum(μj), fxt[0], fjt[0]), file=f)
Tk =np.sum(pj**2 / (2 * masses))
# Print the initial state to the output file
print(output_format.format(0,xk[0],pk[0], np.sum(μj), fxt, fjt[0], Tk), file=f)
# Main MD loop
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
    # i = params.i
    # Update the momenta vv1
    pj[:natoms] += dpj(xk, fj[:natoms], dµ, μj, params) * dt2
    pj[natoms:3*natoms] += fj[natoms:3*natoms] * dt2
    # pj[:] += fj[:] * dt2

    # Update the positions vv2
    rj += pj * dt / masses

    # Update the cavity momenta  vv3
    pk += dpk(xk, μj, params) * dt2
    # xk += pk * dt

    # get forces and charges
    fj, charges = getForcesCharges(rj, natoms, atm, box)  
    Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)
    μj = np.sum(charges * (rj[:natoms] - Rcom))

    if i % 5 == 0:
        dµ = getdµ(natoms, rj, μj, atm, box, dr=0.0001)

    # Update the momenta vv4
    fjt = dpj(xk, fj[:natoms], dµ, μj, params)
    fxt = dpk(xk, μj, params)
    
    pj[:natoms] += fjt * dt2
    pj[natoms:3*natoms] += fj[natoms:3*natoms] * dt2
    # pj[:] += fj[:] * dt2


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

    #print(output_format.format((i+1), xk[0], np.sum(μj), fxt[0], fjt[0]), file=f)
    Tk = np.sum(pj**2 / (2 * masses))
    print(output_format.format((i+1),xk[0], pk[0], np.sum(μj), fxt, fjt[0] , Tk), file=f)

    #---------- some stuff -----------------
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
    #---------- some stuff -----------------
    return rj, pj, xk, pk, fj, μj, dµ, f


# def client(q, params):
#     rj, pj, xk, pk, fj, μj, dµ, f = calculation(rj, pj, xk, pk, fj, μj, dµ, f, params)
#     return q

# --- main loop ---
q = xk[0]
p = pk[0]
# steps = 10
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
