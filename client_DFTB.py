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
    msg  = json.dumps(msg) 
    with socket.socket() as cli:
        cli.connect((HOST, PORT))
        cli.sendall(msg.encode())
        reply = cli.recv(4096)
        # print(f"[CLIENT] Reply: {reply.decode().strip()}")
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

os.system('which dftb+')

os.system('export DFTB_PREFIX=/scratch/user/u.sw216206/dftb_sk_files/mio-1-1/')
os.environ['DFTB_PREFIX'] = '/scratch/user/u.sw216206/dftb_sk_files/mio-1-1/'

atm = 'O33H66' # Define the atomic structure
# Load the atomic coordinates from the 'finalstep_InTheBox.xyz' file
coordina_initial = np.loadtxt('init/initXYZ.dat', usecols=(2,3,4))
velocity_initial = np.loadtxt('init/initPxPyPz.dat', usecols=(0,1,2))
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
# Thermalize the cavity
fk = dpk(xk, μj, params)
for ti in range(thermal_steps):
    xk, pk, fk = vvl(xk, pk, µj, params, fk)
    if ti > thermal_steps  - len(x_save):
        x_save[ti - thermal_steps  + len(x_save)] = xk
        
data = np.histogram(x_save, bins=100)
np.savetxt('hist_coupled.txt', np.c_[data[1][1:], data[0]])
print('Thermalization done')
#------------------------------------------------------

f = open('qt.out', 'w') # Open the file in write mode
x0 = - (2/params.ωc) * μj * params.ηb 
dµ = getdµ(natoms, rj, μj, atm, box, False)

output_format = '{0: >5d} {1: >#016.8f} {2: >#016.8f} {3: >#016.8f} {4: >#016.8f}'
fjt = dpj(xk, fj[:natoms], dµ, μj, params)
fxt = dpk(xk, μj, params)
print(output_format.format(0,xk[0], np.sum(μj), fxt[0], fjt[0]), file=f)
# Main MD loop
tstep0 = time.time()

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
    xk += pk * dt

    # get forces and charges
    fj, charges = getForcesCharges(rj, natoms, atm, box)  
    Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)
    μj = np.sum(charges * (rj[:natoms] - Rcom))

    if i % 5 == 0:
        dµ = getdµ(natoms, rj, μj, atm, box, False)

    # Update the momenta vv4
    fjt = dpj(xk, fj[:natoms], dµ, μj, params)
    fxt = dpk(xk, μj, params)
    
    pj[:natoms] += fjt * dt2
    pj[natoms:3*natoms] += fj[natoms:3*natoms] * dt2
    pk += fxt * dt2

    if i % 2 == 0:
        f.close()
        f = open(f'qt_{i}.dat' , 'w') # Open the file in write mode

    print(output_format.format((i+1), xk[0], np.sum(μj), fxt[0], fjt[0]), file=f)

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
    #---------- some stuff -----------------
    return rj, pj, xk, pk, fj, μj, dµ, f





# Actual code
q, p = xk, pk
steps = params.steps
sleeptime = 1.5
 

for i in range(steps): # time loop
    dat = {'q': q[0], 'p': p[0], 'idx': idx, 'killed': False}
    reply = comm(dat)
    reply = json.loads(reply)
    runs = reply['runs']
    
    pause = True
    while pause:
        time.sleep(sleeptime)
        if runs:
            pause = False
            # Client CALCULATION
            rj, pj, xk, pk, fj, μj, dµ, f = calculation(rj, pj, xk, pk, fj, μj, dµ, f, params,i)
            q, p = xk, pk
            break
        reply = comm(dat)
        reply = json.loads(reply)
        
        
        #---- update from server----
        runs = reply['runs']
        q = float(reply['q'])
        p = float(reply['p'])

dat['killed'] = True
comm(dat)


f.close()  # Close the
print('End')
print('Time taken:',time.time()-t0)