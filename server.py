# server.py
#!/sw/eb/sw/Anaconda3/2024.02-1/bin/python
import socketserver, json, threading
import numpy as np
import sys, time

lock = threading.Lock()

# ------------ SERVER CODE -------------------
def server(q, p):
    # NOTHING 
    return q, p
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
            # switch client off
            srv.runs[data['idx']] = False 
            if srv.loglevel >= 2:
                logfile(f"\n{data}")
                logfile(f"{srv.runs}")
                
            srv.qs[data['idx']] = data['q']  # update the q value for the client
            srv.ps[data['idx']] = data['p']
            # ---------------------------------------

            # --- UPDATE CLIENTS -------------
            srv.activeids.append(data['idx'])
            srv.activeids = list(set(srv.activeids)) # unique ids
            #--------------------------------------
            #is string
            srv.runs = [eval(i) if isinstance(i, str) else i for i in srv.runs]  # ensure runs is a list of booleans
            # ------ All client dead -----------
            if all([not i for i in srv.runs]):
                # CALCULATION---------------------------------------
                srv.qs[:], srv.ps[:] = server(srv.qs[:], srv.ps[:])
                if srv.loglevel >= 1:
                    logfile(f"Step {srv.step} | {time.time() - srv.t0:.2f} s")
                    srv.t0 = time.time()  # reset timer
                    logfile("---------------------")
                resultsfile(f"{srv.step} {' '.join(srv.qs.astype(str))} {' '.join(srv.ps.astype(str))}")
                #------------------------------------------
                srv.revive = True 
                srv.step += 1
            #-----------------------------------
            
            #---- Revive clients ----------------
            if srv.revive:
                if data['idx'] not in srv.reviveids:
                    if srv.loglevel >= 2:
                        logfile(f"Reviving client {data['idx']} ...")
                    srv.runs[data['idx']] = True
                    srv.reviveids.append(data['idx'])
                    srv.reviveids = list(set(srv.reviveids))  # unique ids
                    
                else:
                    if srv.loglevel >= 2:
                        logfile(f"Client {data['idx']} already revived")
            #-----------------------------------
            
            # --- end revive clients -------------
            if len(srv.reviveids) == srv.N:
                if srv.loglevel >= 2:
                    logfile(f"All clients revived (Step {srv.step})")
                srv.revive = False
                srv.reviveids = []
            #-----------------------------------
            reply = {'N': srv.N, 'ids': list(srv.ids), 'runs': srv.runs[data['idx']], 'q': srv.qs[data['idx']], 'p': srv.ps[data['idx']]}
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
    srv.activeids = []
    srv.revive = False
    srv.reviveids = []
    srv.runs   = [True for i in range(srv.N)]
    srv.qs = np.zeros(srv.N)  # Placeholder for any calculations
    srv.ps = np.zeros(srv.N)  # Placeholder for any calculations
    srv.step = 0
    
    srv.loglevel = 1  # 0, 1, 2 (0: no logging, 1: basic logging, 2: detailed logging)
    srv.killed = [False for i in range(srv.N)]
    srv.t0 = time.time()
    srv.serve_forever()
    
    srv.server_close()
    sys.exit(0)
