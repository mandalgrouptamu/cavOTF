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

os.system(f"sbatch run_server.sbatch {N}")  # Start the server
# green color
print(f"{GREEN}Server started with {N} clients.{RESET}")

time.sleep(30)  # Wait for the server to start


for i in range(N):
    
    os.chdir(foldprefix + str(i))  # Change to the client's directory

    for fname in ["qt.out", "WaterMD_Cavity.xyz"]:
        try:
            os.remove(fname)
        except FileNotFoundError:
            pass  # Ignore if file doesn't exist


    os.system(f"sbatch run_client.sbatch {i}")  # Start each client with its index
    print(f"{RED}Client {i} started.{RESET}")
    os.chdir("..")  # Change back to the parent directory




