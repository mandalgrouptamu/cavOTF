import os
import numpy as np
#from ase.build import molecule
from ase.calculators.dftb import Dftb
from ase import Atoms


# Default calculator options (can be overridden via cavotf [dftb] config section)
DEFAULT_CALCULATOR_OPTIONS = {
    "label": "O33H66",
    "Hamiltonian_SCC": "Yes",
    "Hamiltonian_SCCTolerance": 1e-5,
    "Hamiltonian_MaxSCCIterations": 400,
    "Hamiltonian_Mixer": "Anderson { MixingParameter = 0.026 }",
    "Hamiltonian_ConvergentSccOnly": "No",
    "Hamiltonian_MaxAngularMomentum_": "",
    "Hamiltonian_MaxAngularMomentum_O": "p",
    "Hamiltonian_MaxAngularMomentum_H": "s",
    "Hamiltonian_Charge": 0.0,
    "Hamiltonian_ReadInitialCharges": "Yes",
    "Hamiltonian_SpinConstants": "{O = { -0.035 -0.030 -0.030 -0.028 } H = { -0.072 }}",
    "Options_WriteDetailedXml": "No",
    "Options_WriteEigenvectors": "No",
    "Options_WriteResultsTag": "Yes",
    "kpts": (3, 3, 3),
}

_CALCULATOR_OPTIONS = dict(DEFAULT_CALCULATOR_OPTIONS)


def set_calculator_options(overrides: dict[str, object]):
    """Update calculator options from configuration overrides."""

    global _CALCULATOR_OPTIONS
    merged = dict(DEFAULT_CALCULATOR_OPTIONS)
    for key, value in overrides.items():
        if isinstance(value, str):
            value = " ".join(value.replace("\\n", " ").replace("\\\n", " ").split())
        if value is None:
            merged.pop(key, None)
        else:
            merged[key] = value
    _CALCULATOR_OPTIONS = merged

def caldftb(atm, coordinates, box, force=True, charge=True):
    # proform single point DFTB calculation in periodic box and return forces and charges
    atoms = Atoms(atm,  positions =   coordinates) # atomic structure
    cell = np.array([[box, 0.0, 0.0],  # X-direction
                    [0.0, box, 0.0],  # Y-direction
                    [0.0, 0.0, box]]) # Z-direction

    # Set the cell and enable periodic boundary conditions
    atoms.set_cell(cell)  # Set the unit cell
    atoms.set_pbc(True)   # Enable periodic boundary conditions

    options = dict(_CALCULATOR_OPTIONS)
    options.setdefault("label", atm)
    calc = Dftb(**options)

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
    f = open('check.out', 'w')
    print(forces,file =f)
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

def getdµ(natoms, rj, μj, atm, box, dr=0.01):
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
