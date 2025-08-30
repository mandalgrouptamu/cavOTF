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

os.system(f"sbatch run_server.sbatch {N}")  # Start the server
# green color
print(f"{GREEN}Server started with {N} clients.{RESET}")

time.sleep(30)  # Wait for the server to start


for i in range(N):
    # os.system(f"cp {cleanFolder}/client_DFTB_extend.py {foldprefix}{i}/")
    os.system(f"cp run_client_extend.sbatch {foldprefix}{i}/") 
    
    os.chdir(foldprefix + str(i))  # Change to the client's directory

    os.system(f"sbatch run_client_extend.sbatch {i}")  # Start each client with its index
    print(f"{RED}Client {i} started.{RESET}")
    os.chdir("..")  # Change back to the parent directory




