#!/usr/bin/env python3
"""Plotter for the 2-particle Hill encounter test.

3-step workflow (from repo root):
1) make IC
   python3 micro-test-encounter/make_encounter_tipsy.py --mode y-far --b-factor 1.0 --y-span 10 --T-run 6.283185 --out micro-test-encounter/encounter.std
2) run ChaNGa
   cd third_party/ChaNGa && ./ChaNGa ../../micro-test-encounter/encounter_hill.param +p1 && cd ../..
3) make plots
   MPLBACKEND=Agg python3 micro-test-encounter/plot_time_series.py --no-show
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pynbody as pb

DKPC_UNIT = 4.848e-9
KPC_IN_KM = 3.085677581e16
YEAR_SECONDS = 365.25 * 86400.0
LENGTH_UNIT_KM = DKPC_UNIT * KPC_IN_KM
TIME_UNIT_S = YEAR_SECONDS / (2.0 * np.pi)
V_UNIT = LENGTH_UNIT_KM / TIME_UNIT_S


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="2-particle Hill encounter diagnostics")
    p.add_argument("--snap-pattern", default="third_party/ChaNGa/encounter_hill.000*")
    p.add_argument("--ic-file", type=Path, default=Path("third_party/ChaNGa/single_particle.std"))
    p.add_argument("--param-file", type=Path, default=Path("micro-test-encounter/encounter_hill.param"))
    p.add_argument("--mass-threshold", type=float, default=1.0e-20)
    p.add_argument("--m1", type=float, default=1.0e-3)
    p.add_argument("--dOrbdist", type=float, default=1.0)
    p.add_argument("--dCentMass", type=float, default=1.0)
    p.add_argument("--q", type=float, default=1.5)
    p.add_argument("--Omega", type=float, default=1.0)
    p.add_argument("--output-dir", type=Path, default=Path("micro-test-encounter"))
    p.add_argument("--prefix", default="encounter")
    p.add_argument("--no-show", action="store_true")
    return p.parse_args()


def discover_snapshots(pattern: str) -> tuple[np.ndarray, list[str]]:
    paths = glob.glob(pattern)
    parsed = []
    for path in paths:
        m = re.search(r"\.(\d+)$", Path(path).name)
        if m:
            parsed.append((int(m.group(1)), path))
    if not parsed:
        raise FileNotFoundError(f"No snapshots found with pattern: {pattern}")
    parsed.sort(key=lambda x: x[0])
    return np.array([p[0] for p in parsed], dtype=int), [p[1] for p in parsed]


def read_param_dt_nsteps(param_file: Path) -> tuple[float | None, int | None]:
    dt = None
    nsteps = None
    if not param_file.exists():
        return dt, nsteps
    for raw in param_file.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("dDelta"):
            dt = float(line.split("=", 1)[1].strip())
        elif line.startswith("nSteps"):
            nsteps = int(line.split("=", 1)[1].strip())
    return dt, nsteps


def to_code_units(snap: pb.SimSnap) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pos = np.asarray(snap["pos"].in_units("kpc"), dtype=float) / DKPC_UNIT
    vel = np.asarray(snap["vel"].in_units("km s**-1"), dtype=float) / V_UNIT
    mass = np.asarray(snap["mass"], dtype=float)
    return pos, vel, mass


def pick_two_physical_by_mass(mass: np.ndarray, threshold: float) -> tuple[int, int]:
    cand = np.where(mass > threshold)[0]
    if cand.size < 2:
        raise RuntimeError(
            f"Need at least 2 particles above threshold={threshold:g}; got {cand.size}."
        )
    order = cand[np.argsort(mass[cand])[::-1]]
    return int(order[0]), int(order[1])


def load_series(args: argparse.Namespace):
    step_ids_snap, paths = discover_snapshots(args.snap_pattern)

    std = pb.load(str(args.ic_file), fmt="tipsy")
    std_pos = np.asarray(std["pos"], dtype=float)
    std_vel = np.asarray(std["vel"], dtype=float)
    std_mass = np.asarray(std["mass"], dtype=float)
    i1_0, i2_0 = pick_two_physical_by_mass(std_mass, args.mass_threshold)

    x1 = [std_pos[i1_0, 0]]
    y1 = [std_pos[i1_0, 1]]
    vx1 = [std_vel[i1_0, 0]]
    vy1 = [std_vel[i1_0, 1]]

    x2 = [std_pos[i2_0, 0]]
    y2 = [std_pos[i2_0, 1]]
    vx2 = [std_vel[i2_0, 0]]
    vy2 = [std_vel[i2_0, 1]]

    for path in paths:
        s = pb.load(path)
        s.physical_units()
        pos, vel, mass = to_code_units(s)
        i1, i2 = pick_two_physical_by_mass(mass, args.mass_threshold)

        x1.append(pos[i1, 0])
        y1.append(pos[i1, 1])
        vx1.append(vel[i1, 0])
        vy1.append(vel[i1, 1])

        x2.append(pos[i2, 0])
        y2.append(pos[i2, 1])
        vx2.append(vel[i2, 0])
        vy2.append(vel[i2, 1])

    dt_param, nsteps_param = read_param_dt_nsteps(args.param_file)
    n = len(x2)
    if nsteps_param is not None:
        n = min(n, nsteps_param + 1)

    step_ids = np.concatenate([[0], step_ids_snap])[:n]
    dt = dt_param if dt_param is not None else 2.0 * np.pi / max(n - 1, 1)
    t = dt * np.arange(n)

    return {
        "step_ids": step_ids,
        "t": t,
        "x1": np.asarray(x1[:n]),
        "y1": np.asarray(y1[:n]),
        "vx1": np.asarray(vx1[:n]),
        "vy1": np.asarray(vy1[:n]),
        "x2": np.asarray(x2[:n]),
        "y2": np.asarray(y2[:n]),
        "vx2": np.asarray(vx2[:n]),
        "vy2": np.asarray(vy2[:n]),
    }


def save_timeseries(data: dict, out: Path) -> None:
    t = data["t"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)

    ax[0, 0].plot(t, data["x2"], color="tab:red", lw=1.6, label="x2 (test)")
    ax[0, 0].plot(t, data["x1"], color="tab:gray", lw=1.0, ls="--", label="x1 (massive)")
    ax[0, 0].set_ylabel("x")
    ax[0, 0].legend(fontsize=8)

    ax[0, 1].plot(t, data["y2"], color="tab:blue", lw=1.6)
    ax[0, 1].set_ylabel("y2")

    ax[1, 0].plot(t, data["vx2"], color="tab:green", lw=1.6)
    ax[1, 0].set_xlabel("time")
    ax[1, 0].set_ylabel("vx2")

    ax[1, 1].plot(t, data["vy2"], color="tab:purple", lw=1.6)
    ax[1, 1].set_xlabel("time")
    ax[1, 1].set_ylabel("vy2")

    for a in ax.flat:
        a.grid(True, alpha=0.3)

    fig.suptitle("2-particle Hill encounter: time series (test particle)")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    print(f"Saved {out}")


def save_xy(data: dict, out: Path) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.plot(data["x2"], data["y2"], color="tab:red", lw=1.3)
    ax.scatter(data["x2"][0], data["y2"][0], marker="*", s=90, color="black")
    ax.text(data["x2"][0], data["y2"][0], "START", fontsize=9)
    ax.scatter(data["x2"][-1], data["y2"][-1], marker="o", s=55, color="black")
    ax.text(data["x2"][-1], data["y2"][-1], "END", fontsize=9)
    r12 = np.sqrt((data["x2"] - data["x1"]) ** 2 + (data["y2"] - data["y1"]) ** 2)
    i_closest = int(np.argmin(r12))
    ax.scatter(data["x2"][i_closest], data["y2"][i_closest], marker="x", s=90, color="tab:blue")
    ax.text(data["x2"][i_closest], data["y2"][i_closest], "CLOSEST", fontsize=9, color="tab:blue")
    ax.set_xlabel("x2")
    ax.set_ylabel("y2")
    ax.set_title("Encounter trajectory: x2 vs y2")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    print(f"Saved {out}")


def find_y_zero_crossing_time(t: np.ndarray, y: np.ndarray) -> float | None:
    for i in range(1, len(y)):
        y0 = y[i - 1]
        y1 = y[i]
        if y0 == 0.0:
            return float(t[i - 1])
        if (y0 < 0.0 and y1 > 0.0) or (y0 > 0.0 and y1 < 0.0):
            frac = -y0 / (y1 - y0)
            return float(t[i - 1] + frac * (t[i] - t[i - 1]))
    return None


def save_xscaled(data: dict, rh: float, out: Path, metrics: dict) -> None:
    t = data["t"]
    x_scaled = data["x2"] / rh
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.8))
    ax.plot(t, x_scaled, color="tab:red", lw=1.5)
    ax.axhline(metrics["x_start"], color="tab:green", ls="--", lw=1.2, alpha=0.9, label="x_start/R_H")
    ax.axvline(metrics["t_min"], color="tab:blue", ls="--", lw=1.2, alpha=0.9, label="t_min")
    i_peak = int(np.argmax(np.abs(x_scaled - metrics["x_start"])))
    peak_dev = float(np.abs(x_scaled[i_peak] - metrics["x_start"]))
    ax.scatter([t[i_peak]], [x_scaled[i_peak]], color="tab:orange", s=36, zorder=4)
    ax.text(t[i_peak], x_scaled[i_peak], f"peak |Δx|={peak_dev:.3f}", fontsize=8, color="tab:orange")
    ax.set_xlabel("time")
    ax.set_ylabel("x2 / R_H")
    ax.set_title("Scaled encounter coordinate x2/R_H")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    print(f"Saved {out}")


def save_yscaled(data: dict, rh: float, out: Path, t_y_zero: float | None) -> None:
    t = data["t"]
    y_scaled = data["y2"] / rh
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.8))
    ax.plot(t, y_scaled, color="tab:blue", lw=1.5)
    ax.axhline(-10.0, color="black", ls="--", lw=1.0, alpha=0.7)
    ax.axhline(0.0, color="black", ls=":", lw=1.0, alpha=0.8)
    ax.axhline(+10.0, color="black", ls="--", lw=1.0, alpha=0.7)
    if t_y_zero is not None:
        ax.axvline(t_y_zero, color="tab:orange", ls="--", lw=1.2, alpha=0.9)
        ax.text(t_y_zero, 0.0, "y=0", fontsize=9, color="tab:orange", ha="left", va="bottom")
    ax.set_xlabel("time")
    ax.set_ylabel("y2 / R_H")
    ax.set_title("Scaled encounter coordinate y2/R_H")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    print(f"Saved {out}")


def encounter_metrics(data: dict, rh: float) -> dict:
    t = data["t"]
    x_scaled = data["x2"] / rh
    xmin = float(np.min(x_scaled))
    xmax = float(np.max(x_scaled))
    y_scaled = data["y2"] / rh
    ymin = float(np.min(y_scaled))
    ymax = float(np.max(y_scaled))
    x_start = float(x_scaled[0])
    x_end = float(x_scaled[-1])
    y_start = float(y_scaled[0])
    y_end = float(y_scaled[-1])
    delta_x_end = x_end - x_start
    max_abs_dx = float(np.max(np.abs(x_scaled - x_start)))
    r12 = np.sqrt((data["x2"] - data["x1"]) ** 2 + (data["y2"] - data["y1"]) ** 2)
    i_closest = int(np.argmin(r12))
    r_min = float(r12[i_closest])
    r_min_scaled = r_min / rh
    t_closest = float(t[i_closest])
    t_y_zero = find_y_zero_crossing_time(t, data["y2"])
    return {
        "xmin": xmin,
        "xmax": xmax,
        "ymin": ymin,
        "ymax": ymax,
        "x_start": x_start,
        "x_end": x_end,
        "y_start": y_start,
        "y_end": y_end,
        "delta_x_end": delta_x_end,
        "max_abs_dx": max_abs_dx,
        "r_min_scaled": r_min_scaled,
        "t_min": t_closest,
        "t_y_zero": t_y_zero,
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rh = args.dOrbdist * (args.m1 / (3.0 * args.dCentMass)) ** (1.0 / 3.0)
    print(f"Computed R_H = {rh:.9e} from m1={args.m1:.3e}, dOrbdist={args.dOrbdist:.3f}, dCentMass={args.dCentMass:.3f}")

    data = load_series(args)
    print(
        f"Loaded {len(data['t'])-1} snapshots (+t0), "
        f"time span [{data['t'][0]:.6f}, {data['t'][-1]:.6f}]"
    )

    out_ts = args.output_dir / f"{args.prefix}_timeseries.png"
    out_xy = args.output_dir / f"{args.prefix}_xy.png"
    out_xs = args.output_dir / f"{args.prefix}_x_scaled.png"
    out_ys = args.output_dir / f"{args.prefix}_y_scaled.png"

    metrics = encounter_metrics(data, rh)
    save_timeseries(data, out_ts)
    save_xy(data, out_xy)
    save_xscaled(data, rh, out_xs, metrics)
    save_yscaled(data, rh, out_ys, metrics["t_y_zero"])

    print("Encounter summary:")
    print(f"  min(x2/R_H) = {metrics['xmin']:+.6f}")
    print(f"  max(x2/R_H) = {metrics['xmax']:+.6f}")
    print(f"  min(y2/R_H) = {metrics['ymin']:+.6f}")
    print(f"  max(y2/R_H) = {metrics['ymax']:+.6f}")
    print(f"  r_min/R_H = {metrics['r_min_scaled']:+.6f}")
    print(f"  t_min (closest approach) = {metrics['t_min']:.6f}")
    print(
        f"  x_start/R_H = {metrics['x_start']:+.6f}, "
        f"x_end/R_H = {metrics['x_end']:+.6f}, "
        f"delta_x_end = {metrics['delta_x_end']:+.6f}"
    )
    print(
        f"  y_start/R_H = {metrics['y_start']:+.6f}, "
        f"y_end/R_H = {metrics['y_end']:+.6f}"
    )
    if metrics["t_y_zero"] is None:
        print("  y2=0 crossing time = none (no sign change detected)")
    else:
        print(f"  y2=0 crossing time = {metrics['t_y_zero']:.6f}")

    if metrics["r_min_scaled"] > 2.0 or metrics["max_abs_dx"] < 0.2:
        print("  Encounter weak: try --b-factor 0.5 then 0.2")
    if metrics["r_min_scaled"] < 0.3:
        print("  Encounter very strong: try --b-factor 2.0")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
