#!/usr/bin/env python3
"""Plot ChaNGa vs Quinn Hill time series for x, y, vx, vy.

Defaults match the epicycle IC from task2/plot_task2.py and the ChaNGa run card.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pynbody as pb

from plot_task2 import IC as DEFAULT_IC
from plot_task2 import run_quinn_states_for_ic

DKPC_UNIT = 4.848e-9  # 1 code length unit = 4.848e-9 kpc (≈1 AU)
MASS_THRESHOLD = 1.0e-6
KPC_IN_KM = 3.085677581e16
YEAR_SECONDS = 365.25 * 86400.0
LENGTH_UNIT_KM = DKPC_UNIT * KPC_IN_KM
TIME_UNIT_S = YEAR_SECONDS / (2.0 * np.pi)  # 2π code time = 1 year
VELOCITY_UNIT_KM_S = LENGTH_UNIT_KM / TIME_UNIT_S


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot single-particle time series.")
    parser.add_argument("--x0", type=float, default=DEFAULT_IC["x"])
    parser.add_argument("--y0", type=float, default=DEFAULT_IC["y"])
    parser.add_argument("--vx0", type=float, default=DEFAULT_IC["vx"])
    parser.add_argument("--vy0", type=float, default=DEFAULT_IC["vy"])
    parser.add_argument(
        "--snap-pattern",
        type=str,
        default="third_party/ChaNGa/single_particle_hill.0000*",
        help="Glob pattern for ChaNGa snapshots (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("single_particle_timeseries.png"),
        help="Where to save the plot (default: %(default)s)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Skip displaying the plot window",
    )
    return parser.parse_args()


def select_physical_particle(snap: pb.SimSnap):
    masses = snap["mass"].in_units("Msol")
    mask = masses > MASS_THRESHOLD
    if mask.sum() == 0:
        raise RuntimeError("No physical particle (mass>threshold) found in snapshot")
    if mask.sum() > 1:
        idx = np.argmax(masses[mask])
        return snap[mask][idx]
    return snap[mask][0]


def load_changa_series(pattern: str):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No snapshots match pattern {pattern}")

    xs = []
    ys = []
    vxs = []
    vys = []

    # Add t=0 from the input std (code units).
    std = pb.load("third_party/ChaNGa/single_particle.std", fmt="tipsy")
    particle_std = select_physical_particle(std)
    pos0 = np.asarray(particle_std["pos"], dtype=float)[0]
    vel0 = np.asarray(particle_std["vel"], dtype=float)[0]
    xs.append(pos0[0])
    ys.append(pos0[1])
    vxs.append(vel0[0])
    vys.append(vel0[1])

    for fname in paths:
        snap = pb.load(fname)
        snap.physical_units()
        particle = select_physical_particle(snap)
        pos_vec = np.asarray(particle["pos"].in_units("kpc"), dtype=float)[0]
        vel_vec = np.asarray(particle["vel"].in_units("km s**-1"), dtype=float)[0]
        xs.append(pos_vec[0] / DKPC_UNIT)
        ys.append(pos_vec[1] / DKPC_UNIT)
        vxs.append(vel_vec[0] / VELOCITY_UNIT_KM_S)
        vys.append(vel_vec[1] / VELOCITY_UNIT_KM_S)

    xs = np.array(xs)
    ys = np.array(ys)
    vxs = np.array(vxs)
    vys = np.array(vys)
    nsteps = len(paths)
    dt = 2.0 * np.pi / nsteps
    times = dt * np.arange(len(xs))  # include t=0 input + snapshots 1..N
    first_state = {"pos": (xs[0], ys[0], 0.0), "vel": (vxs[0], vys[0], 0.0)}
    return times, xs, ys, vxs, vys, dt, first_state


def read_param_timing(param_path: Path):
    """Return (nsteps, dt) from a ChaNGa param file if present."""
    nsteps = None
    dt = None
    try:
        with param_path.open("r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "#" in line:
                    line = line.split("#", 1)[0].strip()
                    if not line:
                        continue
                if line.startswith("nSteps"):
                    nsteps = int(line.split("=")[1].strip())
                elif line.startswith("dDelta"):
                    dt = float(line.split("=")[1].strip())
    except FileNotFoundError:
        return None, None
    return nsteps, dt


def hill_time_series(x0, y0, vx0, vy0, nsteps, dt):
    ic = {"x": x0, "y": y0, "z": 0.0, "vx": vx0, "vy": vy0, "vz": 0.0}
    x, y, vx, vy = run_quinn_states_for_ic(ic, dt=dt, nsteps=nsteps)
    # run_quinn_states_for_ic returns nsteps+1 samples including t=0
    return (
        np.array(x),
        np.array(y),
        np.array(vx),
        np.array(vy),
    )


def plot_time_series(times, changa_series, quinn_series, output: Path, show: bool):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    labels = [("x(t)", "x"), ("y(t)", "y"), ("vx(t)", "vx"), ("vy(t)", "vy")]
    data_pairs = [
        (quinn_series[0], changa_series[0]),
        (quinn_series[1], changa_series[1]),
        (quinn_series[2], changa_series[2]),
        (quinn_series[3], changa_series[3]),
    ]

    for ax, (title, ylabel), (q_data, c_data) in zip(axes.flat, labels, data_pairs):
        ax.plot(times, q_data, label="Quinn (Python)", color="tab:blue")
        ax.plot(times, c_data, label="ChaNGa", color="tab:red", linestyle="--")
        ax.set_ylabel(f"{ylabel} [code units]")
        ax.set_title(title)
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useOffset=False)
        ax.grid(True, alpha=0.3)

    axes[-1, 0].set_xlabel("time [code units]")
    axes[-1, 1].set_xlabel("time [code units]")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output, dpi=200)
    print(f"Saved time-series plot to {output}")
    if not show:
        return
    plt.show()


def main():
    args = parse_args()
    times, xs, ys, vxs, vys, dt, first_state = load_changa_series(args.snap_pattern)
    param_path = Path("third_party/ChaNGa/single_particle_hill.param")
    param_nsteps, param_dt = read_param_timing(param_path)
    if param_dt is not None:
        dt = param_dt
        times = dt * np.arange(len(times))
    if param_nsteps is not None:
        nsteps = min(len(times), param_nsteps + 1)
    else:
        nsteps = len(times)
    print(f"Using nsteps={nsteps}, dt={dt:.6f} (code units) for Quinn overlay")
    print(
        "Quinn requested IC (code units): "
        f"x={args.x0:.6e}, y={args.y0:.6e}, vx={args.vx0:.6e}, vy={args.vy0:.6e}"
    )
    print(
        "ChaNGa t=0 input (code units): "
        f"x={first_state['pos'][0]:.6e}, y={first_state['pos'][1]:.6e}, "
        f"vx={first_state['vel'][0]:.6e}, vy={first_state['vel'][1]:.6e}"
    )
    tol = 1.0e-6
    dx = abs(first_state["pos"][0] - args.x0)
    dy = abs(first_state["pos"][1] - args.y0)
    dvx = abs(first_state["vel"][0] - args.vx0)
    dvy = abs(first_state["vel"][1] - args.vy0)
    if any(val > tol for val in (dx, dy, dvx, dvy)):
        print("WARNING: ChaNGa snapshot IC does not match requested IC!", file=sys.stderr)
        print(
            f"  Requested (x0,y0,vx0,vy0) = ({args.x0}, {args.y0}, {args.vx0}, {args.vy0})",
            file=sys.stderr,
        )
        print(
            "  Snapshot (x,y,vx,vy)       = "
            f"({first_state['pos'][0]}, {first_state['pos'][1]}, "
            f"{first_state['vel'][0]}, {first_state['vel'][1]})",
            file=sys.stderr,
        )
        print("  Did you forget to regenerate Tipsy and rerun ChaNGa?", file=sys.stderr)
    qx, qy, qvx, qvy = hill_time_series(args.x0, args.y0, args.vx0, args.vy0, nsteps - 1, dt)
    changa_series = (xs[:nsteps], ys[:nsteps], vxs[:nsteps], vys[:nsteps])
    times = times[:nsteps]
    quinn_series = (qx[:nsteps], qy[:nsteps], qvx[:nsteps], qvy[:nsteps])
    plot_time_series(times, changa_series, quinn_series, args.output, not args.no_show)


if __name__ == "__main__":
    main()
