import os
import numpy as np
import re
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


def load_midpoints(path="midpoint.dat"):
    # Each record: i  [xk]  [pk]  [rj]  [pj]
    pat = re.compile(r"""
        \s*(\d+)\s*           # i
        (\[[^\]]*\])\s*       # xk
        (\[[^\]]*\])\s*       # pk
        (\[[^\]]*\])\s*       # rj
        (\[[^\]]*\])          # pj
    """, re.S | re.X)

    out = []
    with open(path, "r") as f:
        text = f.read()

    for m in pat.finditer(text):
        i  = int(m.group(1))
        xk = np.fromstring(m.group(2)[1:-1], sep=' ')
        pk = np.fromstring(m.group(3)[1:-1], sep=' ')
        rj = np.fromstring(m.group(4)[1:-1], sep=' ')
        pj = np.fromstring(m.group(5)[1:-1], sep=' ')
        out.append((i, xk, pk, rj, pj))
    return out


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

# 1) Read all records from the file
records = load_midpoints("midpoint.dat")   # returns list of (i, xk, pk, rj, pj)


# 2) Get the first record (example)
i, xk, pk, rj, pj = records[0]

# 3) Get the last record
# i_last, xk_last, pk_last, rj_last, pj_last = records[-1]

# 4) Get the record with a specific i (e.g., i == 18)
# rec_by_i = {i: (xk, pk, rj, pj) for i, xk, pk, rj, pj in records}
# xk18, pk18, rj18, pj18 = rec_by_i[18]   # raises KeyError if 18 not present

i_past  = np.array([t[0] for t in records])[0]                     # shape (n_records,)
xk = np.array([t[1].item() for t in records])              # if xk is always length-1
pk = np.array([t[2].item() for t in records])              # if pk is always length-1
rj = np.vstack([t[3] for t in records])[0]                    # shape (n_records, len_rj)
pj = np.vstack([t[4] for t in records])[0]


rxj = rj[:natoms]
ryj = rj[natoms:2*natoms]
rzj = rj[2*natoms:3*natoms]

coordinates = np.column_stack((rxj, ryj, rzj))


atoms = Atoms(atm,  positions =   coordinates /bhr) # atomic structure

mass = atoms.get_masses()  * 1822.8884
box = params.box
atoms.set_cell([box, box, box])
atoms.set_pbc(True)
#scaled_positions = atoms.get_positions(wrap=True)

atoms.set_positions(coordinates/bhr)

write('WaterMD_Cavity.xyz',atoms,format='xyz', append = True)


masses = np.concatenate((mass, mass, mass))


Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)

fj, charges = getForcesCharges(rj, natoms, atm, box)

µj = np.sum(charges * rj[:natoms])





#-----for Leapfrog -----------------
pk += dpk(xk, μj, params) * dt2
pk = pk * np.cos(params.ωc * dt2) - params.ωc * xk * np.sin(params.ωc * dt2)

f = open('qt.out', 'a') # Open the file in write mode
x0 = - (2/params.ωc) * μj * params.ηb 
dµ = getdµ(natoms, rj, μj, atm, box, dr=0.0001)

output_format = '{0: >5d} {1: >#016.8f} {2: >#016.8f} {3: >#016.8f} {4: >#016.8f} {5: >#016.8f} {6: >#016.8f}'
fjt = dpj(xk, fj[:natoms], dµ, μj, params)
fxt = dpk(xk, μj, params)
#print(output_format.format(0,xk[0], np.sum(μj), fxt[0], fjt[0]), file=f)
Tk =np.sum(pj**2 / (2 * masses))
# Print the initial state to the output file
print(output_format.format(i_past,xk[0],pk[0], np.sum(μj), fxt, fjt[0], Tk), file=f)
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
    
    rj, pj, xk, pk, fj, μj, dµ, f = calculation(rj, pj, xk, pk, fj, μj, dµ, f, params,(i+i_past))
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
