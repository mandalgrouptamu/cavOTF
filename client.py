# client.py
#!/sw/eb/sw/Anaconda3/2024.02-1/bin/python
import socket, sys, time, json
import numpy as np

fob = open("server_hostname.txt", "r")
HOST = fob.readline().strip()
PORT = int(fob.readline().strip())

idx = sys.argv[1] if len(sys.argv) > 1 else "0"

def comm(msg):
    msg  = json.dumps(msg) 
    with socket.socket() as cli:
        cli.connect((HOST, PORT))
        cli.sendall(msg.encode())
        reply = cli.recv(4096)
        # print(f"[CLIENT] Reply: {reply.decode().strip()}")
    return reply.decode().strip()

# actual code
def client(q, param):
    return q + 1


# Actual code
q = np.random.rand()
steps = 10
sleeptime = 1.5
param = {}

for i in range(steps): # time loop
    dat = {'q': q, 'idx': idx, 'killed': False}
    reply = comm(dat)
    reply = json.loads(reply)
    runs = reply['runs']
    
    pause = True
    while pause:
        time.sleep(sleeptime)
        if runs:
            pause = False
            # Client CALCULATION
            q = client(q, param) # np.random.rand() # one time step propagation
            break
        reply = comm(dat)
        reply = json.loads(reply)
        
        
        #---- update from server----
        runs = reply['runs']
        q = float(reply['q'])

dat['killed'] = True
comm(dat)