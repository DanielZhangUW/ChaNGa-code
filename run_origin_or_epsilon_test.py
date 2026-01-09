#!/usr/bin/env python3
"""Run the origin test with epsilon fallback.

Usage:
  python3 task2/run_origin_or_epsilon_test.py
"""

from __future__ import annotations

import glob
from pathlib import Path
import subprocess
import sys

import numpy as np
import pynbody as pb

DKPC_UNIT = 4.848e-9  # 1 code length unit = 4.848e-9 kpc
MASS_THRESHOLD = 1.0e-6
KPC_IN_KM = 3.085677581e16
YEAR_SECONDS = 365.25 * 86400.0
LENGTH_UNIT_KM = DKPC_UNIT * KPC_IN_KM
TIME_UNIT_S = YEAR_SECONDS / (2.0 * np.pi)
VELOCITY_UNIT_KM_S = LENGTH_UNIT_KM / TIME_UNIT_S


def select_physical_particle(snap: pb.SimSnap):
    masses = snap["mass"].in_units("Msol")
    mask = masses > MASS_THRESHOLD
    if mask.sum() == 0:
        raise RuntimeError("No physical particle found in snapshot")
    if mask.sum() > 1:
        idx = np.argmax(masses[mask])
        return snap[mask][idx]
    return snap[mask][0]


REPO_ROOT = Path(__file__).resolve().parents[1]
CHANGA_DIR = REPO_ROOT / "third_party" / "ChaNGa"
PARAM_PATH = CHANGA_DIR / "single_particle_hill.param"


def print_first_snapshot_state():
    snap = pb.load(str(CHANGA_DIR / "single_particle_hill.000001"))
    snap.physical_units()
    particle = select_physical_particle(snap)
    pos_kpc = np.asarray(particle["pos"].in_units("kpc"), dtype=float)[0]
    vel_kms = np.asarray(particle["vel"].in_units("km s**-1"), dtype=float)[0]
    x = pos_kpc[0] / DKPC_UNIT
    y = pos_kpc[1] / DKPC_UNIT
    vx = vel_kms[0] / VELOCITY_UNIT_KM_S
    vy = vel_kms[1] / VELOCITY_UNIT_KM_S
    print(
        "First snapshot (code units): "
        f"x={x:.6e}, y={y:.6e}, vx={vx:.6e}, vy={vy:.6e}"
    )


def run_changa():
    subprocess.run(
        ["./ChaNGa", str(PARAM_PATH), "+p1"],
        cwd=str(CHANGA_DIR),
        check=True,
    )


def cleanup_snapshots():
    patterns = [
        str(CHANGA_DIR / "single_particle_hill.0000*"),
        str(CHANGA_DIR / "single_particle_hill.chk*"),
    ]
    for pattern in patterns:
        for path in glob.glob(pattern):
            try:
                Path(path).unlink(missing_ok=True)
            except PermissionError:
                print(f"[warn] Could not delete {path} (permission denied)")


def ensure_snapshots_exist():
    snap = CHANGA_DIR / "single_particle_hill.000001"
    if not snap.exists():
        raise RuntimeError("ChaNGa did not produce single_particle_hill.000001")


def run_pipeline(x0: float, output_path: str):
    cleanup_snapshots()
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "tools" / "make_single_particle_tipsy.py"),
            "--x0",
            str(x0),
            "--y0",
            "0.0",
            "--vx0",
            "0.0",
            "--vy0",
            "0.0",
            "--add-guard",
        ],
        check=True,
    )
    run_changa()
    ensure_snapshots_exist()
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "task2" / "plot_time_series.py"),
            "--x0",
            str(x0),
            "--y0",
            "0.0",
            "--vx0",
            "0.0",
            "--vy0",
            "0.0",
            "--output",
            output_path,
        ],
        check=True,
    )


def main() -> None:
    eps_values = [0.0, 1.0e-12, 1.0e-11, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6]
    retries = 0
    for eps in eps_values:
        try:
            if eps == 0.0:
                output = "task2/single_particle_timeseries_origin.png"
            else:
                output = f"task2/single_particle_timeseries_eps{eps:.0e}.png"
            run_pipeline(eps, output)
            print(f"SUCCESS using x0={eps:.0e}")
            print(f"Figure saved to: {output}")
            print_first_snapshot_state()
            print(f"Retries needed: {retries}")
            return
        except subprocess.CalledProcessError:
            retries += 1
            print(f"[warn] ChaNGa failed for x0={eps:.0e}; retrying with larger epsilon")
    raise RuntimeError("ChaNGa failed for all epsilon values up to 1e-6")


if __name__ == "__main__":
    main()
