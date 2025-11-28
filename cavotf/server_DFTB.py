"""Shared server script for CAVOTF dynamics.

This adapts the legacy server_DFTB.py to operate from the project root
without copying into each run directory. Configuration overrides from
``input.txt`` are applied to keep cavity dimensions consistent with the
number of active clients.
"""
from __future__ import annotations

import argparse
import json
import os
import socketserver
import sys
import threading
import time
from operator import attrgetter
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.random import normal as gran
from scipy.fft import fft, ifft

from cavotf.resources import default_template_dir

sys.path.append(str(default_template_dir()))

from funcLM import param  # noqa: E402

try:
    from cavotf.config import OutputConfig, load_config
    from cavotf.dynamics import _recompute_mode_grid
except Exception:  # noqa: BLE001
    load_config = None
    _recompute_mode_grid = None
    OutputConfig = None

lock = threading.Lock()
OUTPUT_CONFIG = None


def _default_output_config():
    if OutputConfig:
        return OutputConfig()
    return SimpleNamespace(
        write_logfile=True,
        write_results=True,
        record_k_space=True,
        print_k_space=False,
    )


def jtok(qj, pj, ω, ωk):
    an = np.sqrt(ω / 2) * (qj + 1j * pj / ω)
    ak = fft(an, norm="ortho")
    akd = np.conj(ak)
    qk = (ak + akd) / np.sqrt(2 * ωk)
    pk = -1j * (ak - akd) * np.sqrt(ωk / 2)
    return qk.real, pk.real


def ktoj(qk, pk, ω, ωk):
    ak = np.sqrt(ωk / 2) * (qk + 1j * pk / ωk)
    aj = ifft(ak, norm="ortho")
    ajd = np.conj(aj)
    qj = (aj + ajd) / np.sqrt(2 * ω)
    pj = -1j * (aj - ajd) * np.sqrt(ω / 2)
    return qj.real, pj.real


def vvlKC(q, p, param_obj, i, output_config):  # only for 1 cavity
    dt, ω0, ωk, λ = attrgetter("dt", "ω0", "ωk", "λ")(param_obj)
    ndof = len(q)
    β = param_obj.β
    dt2 = dt / 1

    σ = (2.0 * λ / (β * param_obj.m)) ** 0.5
    ξ = gran(0, 1, ndof)
    θ = gran(0, 1, ndof)
    c = 0.28867513459

    qk, pk = jtok(q[:], p[:], ω0, ωk)

    if output_config.record_k_space:
        with open("qkpk.dat", "a") as f:
            f.write("\t".join(map(str, qk)) + "\t" + "\t".join(map(str, pk)) + "\n")
        if output_config.print_k_space:
            print(f"Step {i}: qk={qk}, pk={pk}")

    qk1 = qk * np.cos(ωk * dt2) + pk * np.sin(ωk * dt2) / ωk
    pk1 = pk * np.cos(ωk * dt2) - ωk * qk * np.sin(ωk * dt2)
    qk, pk = qk1 * 1.0, pk1 * 1.0

    q[:], p[:] = ktoj(qk, pk, ω0, ωk)

    return q, p


# ------------ SERVER CODE -------------------
def server(q, p, param_obj, i, output_config):
    q, p = vvlKC(q, p, param_obj, i, output_config)
    return q, p


def logfile(msg):
    if OUTPUT_CONFIG and OUTPUT_CONFIG.write_logfile:
        with open("logfile.txt", "a") as f:
            f.write(f"{msg}\n")


def resultsfile(msg):
    if OUTPUT_CONFIG and OUTPUT_CONFIG.write_results:
        with open("results.txt", "a") as f:
            f.write(f"{msg}\n")


class Handler(socketserver.BaseRequestHandler):
    def handle(self):
        with lock:
            data = self.request.recv(4096).strip()
            data = json.loads(data)
            data["idx"] = int(data["idx"])
            srv.stepdata[data["idx"]] = data["step"]
            if srv.loglevel >= 2:
                print(srv.stepdata, "pings from:", data["idx"], "| Server step:", srv.step)

            if srv.stepdata[data["idx"]] > srv.step and srv.update[data["idx"]]:
                srv.qs[data["idx"]] = data["q"]
                srv.ps[data["idx"]] = data["p"]
                srv.update[data["idx"]] = False

            if all(srv.stepdata == srv.stepdata[data["idx"]]) and srv.stepdata[data["idx"]] > srv.step:
                srv.qs[:], srv.ps[:] = server(srv.qs[:], srv.ps[:], param_obj, srv.step, srv.output_config)
                resultsfile(f"{srv.step + 1} {' '.join(srv.qs.astype(str))}")

                srv.step += 1
                if srv.loglevel >= 1:
                    logfile(f"Step {srv.step} | {time.time() - srv.t0:.2f} s")
                    srv.t0 = time.time()
                    logfile("---------------------")
                srv.update = [True for _ in range(srv.N)]

            reply = {
                "N": srv.N,
                "ids": list(srv.ids),
                "step": srv.step,
                "q": srv.qs[data["idx"]],
                "p": srv.ps[data["idx"]],
            }
            self.request.sendall(json.dumps(reply).encode())

            thisKilled = eval(data["killed"]) if isinstance(data["killed"], str) else data["killed"]
            if thisKilled:
                srv.killed[data["idx"]] = True

            if all(srv.killed):
                if srv.loglevel >= 1:
                    logfile("All clients killed. Exiting server.")
                self.request.close()
                srv.shutdown()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int, help="Number of clients")
    parser.add_argument("--config", type=str, default=None, help="Path to cavotf input.txt")
    args = parser.parse_args()

    param_obj = param()
    OUTPUT_CONFIG = _default_output_config()
    config_path = Path(args.config) if args.config else None
    if config_path and load_config and _recompute_mode_grid:
        try:
            cfg = load_config(config_path)
            OUTPUT_CONFIG = cfg.outputs
            overrides = {
                "nk": cfg.physics.nk,
                "β": cfg.physics.beta,
                "λ": cfg.physics.lambda_,
                "ωc": cfg.physics.omega_c,
                "ηb": cfg.physics.eta_b,
                "thermal_steps": cfg.physics.thermal_steps,
            }
            for key, value in overrides.items():
                if hasattr(param_obj, key):
                    setattr(param_obj, key, value)
            if param_obj.nk != args.N:
                param_obj.nk = args.N
            _recompute_mode_grid(param_obj)
        except Exception as exc:  # noqa: BLE001
            print(f"Warning: failed to apply config overrides: {exc}")
    else:
        param_obj.nk = args.N

    os.system('echo "$(hostname)" ')
    os.system('echo "$(hostname)" > server_hostname.txt')
    port = 16012 + np.random.randint(100)
    with open("server_hostname.txt", "a") as f:
        f.write(f"{port}\n")

    print(f" Server starting on localhost:{port}")

    srv = socketserver.ThreadingTCPServer(("", port), Handler)
    srv.N = int(args.N)
    srv.ids = range(srv.N)
    srv.output_config = OUTPUT_CONFIG
    srv.update = [True for _ in range(srv.N)]
    srv.qs = np.zeros(param_obj.nk)
    srv.ps = np.zeros(param_obj.nk)
    srv.step = -1
    srv.stepdata = np.zeros(srv.N) - 1
    srv.loglevel = 2
    srv.killed = [False for _ in range(srv.N)]
    srv.t0 = time.time()

    srv.serve_forever()
    srv.server_close()
    sys.exit(0)
