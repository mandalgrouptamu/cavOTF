import os 
import sys
import numpy as np
import time
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

os.system(f"sbatch run_server.sbatch {N}")  # Start the server
# green color
print(f"{GREEN}Server started with {N} clients.{RESET}")

time.sleep(30)  # Wait for the server to start




for i in range(N):
    
    
    os.system(f"mkdir -p {foldprefix}{i}")  # Create a directory for each client
    os.system(f"cp -r {cleanFolder}/* {foldprefix}{i}/")  # Copy the clean folder contents to each client's directory
    os.system(f"cp run_client.sbatch {foldprefix}{i}/")  # Copy the run_client.sbatch file to each client's directory
    
    os.chdir(foldprefix + str(i))  # Change to the client's directory
    os.system(f"sbatch run_client.sbatch {i}")  # Start each client with its index
    print(f"{RED}Client {i} started.{RESET}")
    os.chdir("..")  # Change back to the parent directory
