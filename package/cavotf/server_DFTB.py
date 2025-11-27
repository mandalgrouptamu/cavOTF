# server.py
#!/sw/eb/sw/Anaconda3/2024.02-1/bin/python
import socketserver, json, threading
import numpy as np
import sys, time, os
from operator import attrgetter
from numpy.random import normal as gran
from scipy.fft import fft, dct, ifft


sys.path.append(os.popen("pwd").read().replace("\n","")+"/DFTB_clean")
from funcLM import param

param = param()  # Initialize parameters

lock = threading.Lock()


def vvlKC(q, p, param): #only for 1 cavity
    dt,  ω0, ωk , λ = attrgetter("dt",  "ω0", "ωk", "λ")(param)
    
    ndof = len(q)  # number of degrees of freedom
    β  = param.β
    dt2 = dt/1
    λ = λ * 0.0#/ param.m
    
    σ = (2.0 * λ/(β * param.m )) ** 0.5
    ξ = gran(0, 1, ndof)  #np.array([0.5 * gran() for i in range(len(x))])
    θ = gran(0, 1, ndof) #np.array([gran() * 0.28867513459  for i in range(len(x))])
    c = 0.28867513459

    qk, pk = jtok(q[:], p[:], ω0, ωk)

    # f1 = -ωk * qk

    # A = (0.5 * dt2**2) * (f1 - λ * p) + (σ * dt2**(3.0/2.0)) * (0.5 * ξ + c * θ) 
    # #---- q update -----------
    # qk += (pk * dt2 + A) 
    # #-------------------------
    # f2 = -ωk * qk
    # #---- p update ----------- 
    # pk += ( 0.5 * dt2 * (f1+f2) - dt2 * λ * pk +  σ * (dt**0.5) * ξ - A * λ ) 
    # #-------------------------

    
    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2)/ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0

    q[:], p[:] = ktoj(qk, pk, ω0, ωk)

    return q, p



def jtok(qj, pj, ω, ωk):
    an = np.sqrt(ω/2) * (qj + 1j*pj/ω)
    ak = fft(an, norm='ortho')
    akd =  np.conj(ak)
    qk = (ak + akd)/np.sqrt(2*ωk)
    pk = -1j*(ak - akd)*np.sqrt(ωk/2)
     
    return qk.real, pk.real


def ktoj(qk, pk, ω, ωk):
    
    ak = np.sqrt(ωk/2) * (qk + 1j*pk/ωk)

    
    aj = ifft(ak, norm='ortho')
    ajd =  np.conj(aj)
    qj = (aj + ajd)/np.sqrt(2*ω)
    pj = -1j*(aj - ajd)*np.sqrt(ω/2)
    return qj.real, pj.real

# ------------ SERVER CODE -------------------
def server(q, p, param):
    print("Server got q and p:", q, p)
    q, p = vvlKC(q, p, param)
    print("Server updated q and p:", q, p)
    # q +=1
    # p += 0.1
    return q , p
#----------------------------------------------




def logfile(msg):
        with open("logfile.txt", "a") as f:
            f.write(f"{msg}\n")

def resultsfile(msg):
        with open("results.txt", "a") as f:
            f.write(f"{msg}\n")

class Handler(socketserver.BaseRequestHandler):
    def handle(self):
        with lock:
            # --- DATA -----------------------------
            data = self.request.recv(4096).strip()
            data = json.loads(data)
            data['idx'] = int(data['idx'])  # ensure idx is an integer
            srv.stepdata[data['idx']] = data['step']
            # update the step based on the minimum step of all clients
            if srv.loglevel >= 3:
                logfile(f"\n{data}")
            if srv.loglevel >= 2:
                print(srv.stepdata, "pings from:", data['idx'], "| Server step:", srv.step)
                
            if srv.stepdata[data['idx']] > srv.step and srv.update[data['idx']]:
                srv.qs[data['idx']] = data['q']
                srv.ps[data['idx']] = data['p']    
                srv.update[data['idx']] = False  # Mark this client as updated
            
            if all(srv.stepdata == srv.stepdata[data['idx']]) and srv.stepdata[data['idx']] > srv.step:
                if srv.loglevel >= 2:
                    print("Ding ding ding! All clients at the same step. Proceeding with calculations.")
                    logfile(f"All clients at step {srv.stepdata[data['idx']]}. Proceeding with calculations.")   
                srv.qs[:] , srv.ps[:] = server(srv.qs[:], srv.ps[:], param)  # Call the server function to update qs and ps
                
                
                resultsfile(f"{srv.step + 1} {' '.join(srv.qs.astype(str))}")

                srv.step += 1 
                if srv.loglevel >= 1:
                    logfile(f"Step {srv.step} | {time.time() - srv.t0:.2f} s")
                    srv.t0 = time.time()  # reset timer
                    logfile("---------------------")
                srv.update = [True for i in range(srv.N)]  # Reset update flags for all clients
                
            #------------------ REPLY -----------------
            reply = {'N': srv.N,
                    'ids': list(srv.ids), 
                    'step': srv.step, 
                    'q': srv.qs[data['idx']],
                    'p': srv.ps[data['idx']]}
            self.request.sendall(json.dumps(reply).encode())
            
            
            #------------------ KILL SERVER -------------------------------------------------------
            thisKilled = eval(data['killed']) if isinstance(data['killed'], str) else data['killed']
            if thisKilled:
                srv.killed[data['idx']] = True
            
            if all(srv.killed):
                if srv.loglevel >= 1:
                    logfile("All clients killed. Exiting server.")
                self.request.close()
                srv.shutdown()

            
if __name__ == "__main__":
    port = 16012 + np.random.randint(100)  # Random port for each run
    with open("server_hostname.txt", "a") as f:
            f.write(f"{port}\n")
            
    



    
    srv = socketserver.ThreadingTCPServer(("", port), Handler)
    srv.N      = int(sys.argv[1]) 
    srv.ids    = range(srv.N)  
    srv.update = [True for i in range(srv.N)]
    #-----------------------------------------------------------
    srv.qs = np.zeros(srv.N)  # Placeholder for any calculations
    srv.ps = np.zeros(srv.N)  # Placeholder for any calculations
    #-----------------------------------------------------------

    srv.step = -1
    srv.stepdata = np.zeros(srv.N) -1  # Placeholder for step data
    
    srv.loglevel = 2  # 0, 1, 2 (0: no logging, 1: basic logging, 2: detailed logging)
    srv.killed = [False for i in range(srv.N)]
    srv.t0 = time.time()
    srv.serve_forever()
    
    srv.server_close()
    sys.exit(0)

