from __future__ import annotations

import random
import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

from .config import CavOTFConfig
from .dynamics import init, vvl
from .parameters import load_param


class CavOTFWorkflow:
    def __init__(self, config: CavOTFConfig):
        self.config = config
        self.param, self.param_module = load_param(config.clean_template_dir)
        self._apply_physics_overrides()

    def run_all(self) -> None:
        self.stage_clients()
        if self.config.hpc.run_get_mu:
            self.wait_for_mu_results()
        if self.config.hpc.run_dynamics:
            self.prepare_initial_conditions()
            self.launch_server_and_clients()

    def stage_clients(self) -> None:
        working = self.config.working_dir
        working.mkdir(parents=True, exist_ok=True)

        self._clean_previous_runs()

        source_dirs = self._gather_source_dirs()
        if len(source_dirs) < self.param.nk:
            raise RuntimeError(
                f"Not enough source folders in {self.config.data_source_dir} to stage {self.param.nk} clients"
            )

        random.shuffle(source_dirs)

        for idx in range(self.param.nk):
            run_dir = working / f"{self.config.run_prefix}{idx}"
            self._prepare_client_dir(run_dir, source_dirs[idx])
            if self.config.hpc.run_get_mu:
                self._submit_mu_job(run_dir)

    def wait_for_mu_results(self) -> None:
        required = [self._run_dir(i) / self.config.mu_results_filename for i in range(self.param.nk)]
        deadline = None
        if self.config.mu_timeout_seconds > 0:
            deadline = time.time() + self.config.mu_timeout_seconds

        while True:
            missing = [path for path in required if not path.exists()]
            if not missing:
                return

            if deadline and time.time() > deadline:
                missing_str = ", ".join(str(p) for p in missing)
                raise TimeoutError(f"Timed out waiting for mu results: {missing_str}")

            time.sleep(self.config.mu_poll_interval)

    def prepare_initial_conditions(self) -> None:
        mu_values = np.zeros(self.param.nk)
        for idx in range(self.param.nk):
            run_dir = self._run_dir(idx)
            mu_values[idx] = np.loadtxt(run_dir / self.config.mu_results_filename)

        q, p = init(mu_values, self.param)
        for _ in range(self.param.thermal_steps * self.config.thermalization_multiplier):
            q, p = vvl(q, p, mu_values, self.param)

        for idx in range(self.param.nk):
            run_dir = self._run_dir(idx)
            np.savetxt(run_dir / "initial.dat", np.c_[q[idx], p[idx]])

    def launch_server_and_clients(self) -> None:
        server_script = self._ensure_script_available(self.config.server_sbatch, self.config.working_dir)
        subprocess.run(
            ["sbatch", *self._sbatch_options(), *self._export_environment(), server_script.name, str(self.param.nk)],
            cwd=self.config.working_dir,
            check=True,
        )
        time.sleep(self.config.server_wait_seconds)

        for idx in range(self.param.nk):
            run_dir = self._run_dir(idx)
            for fname in ["qt.out", "WaterMD_Cavity.xyz"]:
                path = run_dir / fname
                if path.exists():
                    path.unlink()

            client_script = self._ensure_script_available(self.config.client_sbatch, run_dir)
            subprocess.run(
                [
                    "sbatch",
                    *self._sbatch_options(),
                    *self._export_environment(),
                    client_script.name,
                    str(idx),
                ],
                cwd=run_dir,
                check=True,
            )

    def _prepare_client_dir(self, run_dir: Path, source_dir: Path) -> None:
        run_dir.mkdir(parents=True, exist_ok=True)
        shutil.copytree(self.config.clean_template_dir, run_dir, dirs_exist_ok=True)
        self._ensure_script_available(self.config.client_sbatch, run_dir)

        for source_name, target_name in self._initial_file_pairs():
            source_path = source_dir / source_name
            if not source_path.exists():
                raise FileNotFoundError(f"Missing source file {source_path}")
            shutil.copy2(source_path, run_dir / target_name)

    def _submit_mu_job(self, run_dir: Path) -> None:
        script_path = self._ensure_script_available(self.config.mu_submission_script, run_dir)
        subprocess.run(
            ["sbatch", *self._sbatch_options(), *self._export_environment(), Path(script_path).name],
            cwd=run_dir,
            check=True,
        )

    def _ensure_script_available(self, script: Path | str, destination: Path) -> Path:
        script_path = Path(script)
        if not script_path.is_absolute():
            candidate = destination / script_path.name
            if candidate.exists():
                return candidate
            script_path = Path(self.config.clean_template_dir / script_path)

        if script_path.exists():
            target = destination / script_path.name
            if target.resolve() != script_path.resolve():
                shutil.copy2(script_path, target)
            return target

        raise FileNotFoundError(f"Unable to locate submission script: {script_path}")

    def _clean_previous_runs(self) -> None:
        for txt_file in self.config.working_dir.glob("*.txt"):
            txt_file.unlink()

        for folder in self.config.working_dir.glob(f"{self.config.run_prefix}*"):
            if folder.is_dir():
                shutil.rmtree(folder)

    def _gather_source_dirs(self) -> list[Path]:
        return [p for p in self.config.data_source_dir.iterdir() if p.is_dir()]

    def _run_dir(self, idx: int) -> Path:
        return self.config.working_dir / f"{self.config.run_prefix}{idx}"

    def _initial_file_pairs(self) -> Iterable[Tuple[str, str]]:
        return (
            (self.config.initial_positions_name, self.config.initial_position_target),
            (self.config.initial_velocities_name, self.config.initial_velocity_target),
        )

    def _sbatch_options(self) -> list[str]:
        opts: list[str] = []
        hpc = self.config.hpc
        if hpc.partition:
            opts.extend(["--partition", hpc.partition])
        if hpc.account:
            opts.extend(["--account", hpc.account])
        if hpc.cpus_per_job > 0:
            opts.extend(["--cpus-per-task", str(hpc.cpus_per_job)])
        return opts

    def _export_environment(self) -> list[str]:
        exports = []
        env_pairs = []
        if self.config.hpc.dftb_prefix:
            env_pairs.append(f"DFTB_PREFIX={shlex.quote(self.config.hpc.dftb_prefix)}")
        if self.config.hpc.dftb_command:
            env_pairs.append(f"DFTB_COMMAND={shlex.quote(self.config.hpc.dftb_command)}")
        if env_pairs:
            exports.extend(["--export", "ALL," + ",".join(env_pairs)])
        return exports

    def _apply_physics_overrides(self) -> None:
        physics = self.config.physics
        param = self.param
        if physics.use_param_nk:
            param.nk = physics.nk
        else:
            param.nk = self.config.general.n_trajectories

        param.β = physics.beta
        param.λ = physics.lambda_
        param.ωc = physics.omega_c
        param.ω0 = physics.omega_c
        param.ηb = physics.eta_b

        param.ky = np.fft.fftfreq(param.nk) * param.nk * 2 * np.pi / (param.dL * param.nk)
        param.ωk = np.sqrt(param.ωc**2 + (param.c * param.ky) ** 2)
