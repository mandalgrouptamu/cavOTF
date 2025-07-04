import os
import numpy as np
#from ase.build import molecule
from ase.calculators.dftb import Dftb
from ase import Atoms

def caldftb(atm, coordinates, box, force=True, charge=True):
    # proform single point DFTB calculation in periodic box and return forces and charges
    atoms = Atoms(atm,  positions =   coordinates) # atomic structure
    cell = np.array([[box, 0.0, 0.0],  # X-direction
                    [0.0, box, 0.0],  # Y-direction
                    [0.0, 0.0, box]]) # Z-direction

    # Set the cell and enable periodic boundary conditions
    atoms.set_cell(cell)  # Set the unit cell
    atoms.set_pbc(True)   # Enable periodic boundary conditions

    calc = Dftb(label = atm,
            Hamiltonian_SCC='Yes',
            Hamiltonian_SCCTolerance=1e-005,
            Hamiltonian_MaxSCCIterations=500,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_C='p',
            Hamiltonian_MaxAngularMomentum_O='p',
            Hamiltonian_MaxAngularMomentum_H='s',
            kpts=(3,3,3))

    atoms.calc = calc
    if force:
        force = atoms.get_forces()
        
    if charge:
        force = atoms.get_forces()
        charge = atoms.get_charges()

    return force, charge


def getForcesCharges(rj, natoms, atm, box):
    bhr = 1.8897259886
    evdivA = 27.2114 * bhr
    rxj = rj[:natoms]
    ryj = rj[natoms: 2*natoms]
    rzj = rj[2*natoms: 3*natoms]
    coordinates = np.column_stack((rxj, ryj, rzj))

    forces, charges = caldftb(atm, coordinates/bhr, box)
    fj = np.concatenate((forces[:,0], forces[:,1], forces[:,2]))/evdivA
    return fj, charges

def getCharges(rj, natoms, atm, box):
    bhr = 1.8897259886
    rxj = rj[:natoms]
    ryj = rj[natoms:2*natoms]
    rzj = rj[2*natoms:3*natoms]
    coordinates = np.column_stack((rxj, ryj, rzj))
    _, charges = caldftb(atm,coordinates/bhr,box, True, True)
    return charges

def getdµ(natoms, rj, μj, atm, box, dr=0.0001):
        """Forward finite difference to calculate the dipole derivative
          If dr is False it will return charges instead of dipole derivative.
        """
        if not dr:
            return getCharges(rj, natoms, atm, box)
        else:
            dµ = np.zeros(natoms) # initialize the dipole derivative
            for j in range(natoms):
                    rj2 = rj * 1.0
                    rj2[j] += dr
                    charges2 = getCharges(rj2, natoms, atm, box)
                    μ2 = np.sum(charges2 * (rj2[:natoms]))
                    dµ[j] = (μ2-μj)/dr
        return dµ