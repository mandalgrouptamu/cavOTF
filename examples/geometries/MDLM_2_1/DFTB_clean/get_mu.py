import os
import numpy as np
#from ase.build import molecule
#from ase.calculators.dftb import Dftb
from ase.io import write
from ase import Atoms
import time
from dftb import getForcesCharges, getCharges, getdµ
from funcLM import *

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

os.system('export DFTB_PREFIX=/home/aamini1/dftb_sk_files/mio-1-1/')
os.environ['DFTB_PREFIX'] = '/home/aamini1/dftb_sk_files/mio-1-1/'
os.system('export A="/home/aamini1/DFTB/dftbplus-24.1.x86_64-linux/bin/dftb+ "')

os.environ["DFTB_COMMAND"] = "/home/aamini1/DFTB/dftbplus-24.1.x86_64-linux/bin/dftb+"



atm = 'O33H66' # Define the atomic structure
# Load the atomic coordinates from the 'finalstep_InTheBox.xyz' file
coordina_initial = np.loadtxt('initXYZ.dat', usecols=(2,3,4))
velocity_initial = np.loadtxt('initPxPyPz.dat', usecols=(0,1,2))
coordinates = np.array(coordina_initial) * bhr #a.u.

atoms = Atoms(atm,  positions =   coordinates /bhr) # atomic structure
mass = atoms.get_masses()  * 1822.8884
box = params.box
atoms.set_cell([box, box, box])
atoms.set_pbc(True)

coordinates = atoms.get_positions(wrap=False) * bhr
atoms.set_positions(coordinates/bhr)

write('WaterMD_Cavity.xyz',atoms,format='xyz', append = True)


rxj, ryj, rzj = coordinates[:,0], coordinates[:,1], coordinates[:,2]
rj = np.concatenate((rxj, ryj, rzj))
Rcom = np.sum(rj[:natoms] * mass) / np.sum(mass)

fj, charges = getForcesCharges(rj, natoms, atm, box)
µj = np.sum(charges * rxj)

np.savetxt('dmu.dat', [µj], fmt='%.6f')