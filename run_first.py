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

os.system("rm -rf *.txt")  # Clean up previous run files
os.system(f"rm -rf {foldprefix}*")

location = "/scratch/user/u.sw216206/vsc-fluxside/project1/cavityDFTB/with_dmu/"
folders = [f for f in os.listdir(location) if os.path.isdir(os.path.join(location, f))]

random.shuffle(folders)

for i in range(N):
    
    
    os.system(f"mkdir -p {foldprefix}{i}")  # Create a directory for each client
    os.system(f"cp -r {cleanFolder}/* {foldprefix}{i}/")  # Copy the clean folder contents to each client's directory
    os.system(f"cp run_client.sbatch {foldprefix}{i}/")  # Copy the run_client.sbatch file to each client's directory

    os.system(f"cp /scratch/user/u.sw216206/vsc-fluxside/project1/cavityDFTB/with_dmu/{folders[i]}/thermaliz_water__InTheBox.dat {foldprefix}{i}/initXYZ.dat")  # Copy the initial coordinates
    os.system(f"cp /scratch/user/u.sw216206/vsc-fluxside/project1/cavityDFTB/with_dmu/{folders[i]}/thermaliz_water_vel__InTheBox.dat {foldprefix}{i}/initPxPyPz.dat")
    
    os.chdir(foldprefix + str(i))  # Change to the client's directory
    os.system(f"sbatch muDFTB.slurm ")  # Start each client with its index
    print(f"{GREEN}mu computation {i} started.{RESET}")
    os.chdir("..")  # Change back to the parent directory


