"""Compare ChaNGa's single-particle Hill trajectory with the Python Quinn integrator."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pynbody as pb

from plot_task2 import IC, run_quinn_points_for_ic

TWOPI = 2.0 * np.pi
N_STEPS_ORBIT = 100
DT_ORBIT = TWOPI / N_STEPS_ORBIT
T_FINAL = TWOPI
DKPC_UNIT = 4.848e-9  # Must match dKpcUnit in single_particle_hill.param
PHYS_MIN_MASS = 1e-6  # physical Hill particle: m ~ 1; guard: m << 1
DEBUG_GUARD_MAX_MASS = 1e-6  # anything below this is considered a guard particle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare ChaNGa and Quinn orbits")
    parser.add_argument(
        "--snapshot-base",
        type=str,
        default="single_particle_hill",
        help="Base name of ChaNGa snapshots (default: single_particle_hill)",
    )
    parser.add_argument(
        "--snapshot-step",
        type=int,
        default=100,
        help="Final snapshot index to report diagnostics for (default: 100)",
    )
    parser.add_argument("--dt-fraction", type=float, default=1.0 / N_STEPS_ORBIT)
    parser.add_argument("--nsteps", type=int, default=N_STEPS_ORBIT)
    parser.add_argument(
        "--output-figure",
        type=str,
        help="Optional path to save the comparison plot",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Skip showing the matplotlib window",
    )
    return parser.parse_args()


def select_physical_particle(sim: pb.SimSnap) -> pb.SimSnap | None:
    """Return a sub-sim containing exactly the physical particle (mass > PHYS_MIN_MASS)."""
    m = sim["mass"]
    mask_phys = m > PHYS_MIN_MASS
    if mask_phys.sum() == 0:
        print(f"[warn] no particles with m > {PHYS_MIN_MASS}")
        return None
    if mask_phys.sum() == 1:
        return sim[mask_phys]

    # Choose the physical particle closest to the origin (smallest r^2 in kpc)
    x_kpc_all = sim["x"].in_units("kpc")
    y_kpc_all = sim["y"].in_units("kpc")
    x_kpc = x_kpc_all[mask_phys]
    y_kpc = y_kpc_all[mask_phys]
    r2 = x_kpc**2 + y_kpc**2
    idx_min = int(np.argmin(r2))
    phys_indices = np.where(mask_phys)[0]
    chosen_index = phys_indices[idx_min]
    full_mask = np.zeros(len(sim), dtype=bool)
    full_mask[chosen_index] = True
    return sim[full_mask]


def select_main_particle(snap: pb.SimSnap, guard_max_mass: float = DEBUG_GUARD_MAX_MASS):
    """Return the main Hill particle, ignoring guard particles by mass."""
    masses = snap["mass"]
    mask_phys = masses > guard_max_mass

    if mask_phys.sum() != 1:
        print(f"[warn] expected 1 main particle, found {mask_phys.sum()}")

    # If there is at least one physical particle, take the first one
    if mask_phys.sum() >= 1:
        return snap[mask_phys][0]

    # Fallback: no particle above guard_max_mass, just take the most massive
    idx = masses.argmax()
    return snap[idx]


def debug_changa_initial_states():
    """Print ChaNGa input and first-step states for debugging."""
    std = pb.load("third_party/ChaNGa/single_particle.std")
    p_in = select_main_particle(std)
    mass_in = float(p_in["mass"])
    print("ChaNGa INPUT single_particle.std (raw units):")
    print(f"  mass = {mass_in:.6e}")
    print(f"  pos  = {p_in['pos']}")
    print(f"  vel  = {p_in['vel']}")

    snap1 = pb.load("third_party/ChaNGa/single_particle_hill.000001")
    p1 = select_main_particle(snap1)
    mass_out = float(p1["mass"].in_units("Msol"))
    print("ChaNGa OUTPUT single_particle_hill.000001:")
    print(f"  mass (Msol) = {mass_out:.6e}")
    print(f"  pos (kpc) = {p1['pos'].in_units('kpc')}")
    print(f"  vel (km/s) = {p1['vel'].in_units('km s**-1')}")


def load_changa_trajectory(snapshot_base: Path, n_steps: int):
    xs, ys = [], []
    for i in range(1, n_steps + 1):
        fname = snapshot_base.with_name(f"{snapshot_base.name}.{i:06d}")
        if not fname.exists():
            print(f"[warn] snapshot not found: {fname}")
            continue
        sim = pb.load(str(fname))
        sim.physical_units()
        sim_phys = select_physical_particle(sim)
        if sim_phys is None:
            continue
        x_kpc = sim_phys["x"].in_units("kpc")
        y_kpc = sim_phys["y"].in_units("kpc")
        xs.append(float(x_kpc[0]))
        ys.append(float(y_kpc[0]))
    return np.array(xs), np.array(ys)


def run_quinn_orbit(dt: float, nsteps: int) -> Tuple[Iterable[float], Iterable[float]]:
    X, Y = run_quinn_points_for_ic(IC, dt=dt, nsteps=nsteps)
    return np.asarray(X, dtype=float), np.asarray(Y, dtype=float)


def main() -> None:
    args = parse_args()
    snap_base = Path("third_party/ChaNGa") / args.snapshot_base
    final_snapshot = snap_base.with_name(f"{snap_base.name}.{args.snapshot_step:06d}")
    if not final_snapshot.exists():
        print(f"[warn] Final snapshot {final_snapshot} missing; trajectory loader will handle available steps")

    x_ch_arr, y_ch_arr = load_changa_trajectory(snap_base, args.nsteps)
    if x_ch_arr.size == 0:
        raise RuntimeError("No ChaNGa trajectory points found; aborting")

    nsteps = args.nsteps
    dt = args.dt_fraction * TWOPI
    total_time = nsteps * dt
    X_py_code, Y_py_code = run_quinn_orbit(dt=dt, nsteps=nsteps)
    X_py_kpc = X_py_code * DKPC_UNIT
    Y_py_kpc = Y_py_code * DKPC_UNIT

    step_count = min(len(x_ch_arr), len(X_py_kpc) - 1)
    if step_count == 0:
        raise RuntimeError("No overlapping steps between ChaNGa and Quinn")

    quinn_step_x = X_py_kpc[1 : step_count + 1]
    quinn_step_y = Y_py_kpc[1 : step_count + 1]
    changa_x = x_ch_arr[:step_count]
    changa_y = y_ch_arr[:step_count]

    print("Per-step diagnostics (kpc, radians):")
    for i in range(step_count):
        t = i * dt
        x_q, y_q = quinn_step_x[i], quinn_step_y[i]
        x_c, y_c = changa_x[i], changa_y[i]
        r_q = np.hypot(x_q, y_q)
        phi_q = np.arctan2(y_q, x_q)
        r_c = np.hypot(x_c, y_c)
        phi_c = np.arctan2(y_c, x_c)
        dx = x_c - x_q
        dy = y_c - y_q
        dr = r_c - r_q
        raw_dphi = phi_c - phi_q
        dphi = (raw_dphi + np.pi) % (2.0 * np.pi) - np.pi
        dphi_deg = np.degrees(dphi)
        print(f"Step {i}, t={t:.6f}:")
        print(f"  Quinn : x={x_q:.6e}, y={y_q:.6e}, r={r_q:.6e}, phi={phi_q:.6e}")
        print(f"  ChaNGa: x={x_c:.6e}, y={y_c:.6e}, r={r_c:.6e}, phi={phi_c:.6e}")
        print(
            f"  delta : dx={dx:.6e}, dy={dy:.6e}, dr={dr:.6e}, "
            f"raw_dphi={raw_dphi:.6e}, wrapped_dphi={dphi:.6e} (deg={dphi_deg:.3f})"
        )

    dx_steps = changa_x - quinn_step_x
    dy_steps = changa_y - quinn_step_y
    final_dx = dx_steps[-1]
    final_dy = dy_steps[-1]

    print("ChaNGa trajectory:")
    print(f"  final position (x, y) = ({changa_x[-1]:.8f}, {changa_y[-1]:.8f})")
    print("Python Quinn integrator:")
    print(
        f"  nsteps, dt           = {nsteps}, {dt:.6f} (total T = {total_time:.6f}, target {T_FINAL:.6f})"
    )
    print(f"  final position (x, y) = ({quinn_step_x[-1]:.8f}, {quinn_step_y[-1]:.8f})")
    print("Difference at final step (ChaNGa - Quinn):")
    print(f"  dx = {final_dx:.3e}, dy = {final_dy:.3e}, |Δr| = {np.hypot(final_dx, final_dy):.3e}")
    print("Max absolute error over all steps:")
    print(f"  max|dx| = {np.max(np.abs(dx_steps)):.3e} kpc, max|dy| = {np.max(np.abs(dy_steps)):.3e} kpc")

    x_ch_final = changa_x[-1]
    y_ch_final = changa_y[-1]
    x_q_final = quinn_step_x[-1]
    y_q_final = quinn_step_y[-1]
    r_c = np.hypot(x_ch_final, y_ch_final)
    phi_c = np.arctan2(y_ch_final, x_ch_final)
    r_q = np.hypot(x_q_final, y_q_final)
    phi_q = np.arctan2(y_q_final, x_q_final)
    d_r = r_c - r_q
    d_phi = (phi_c - phi_q + np.pi) % (2.0 * np.pi) - np.pi
    d_phi_deg = np.degrees(d_phi)

    print(f"=== Single-particle Hill test ({nsteps} steps per orbit, dt = {dt:.6f}) ===")
    print(f"ChaNGa radius : {r_c:.6e} kpc")
    print(f"Quinn radius  : {r_q:.6e} kpc")
    print(f"Δr            : {d_r:.6e} kpc")
    print(f"ChaNGa phase  : {phi_c:.6e} rad")
    print(f"Quinn phase   : {phi_q:.6e} rad")
    print(f"Δφ            : {d_phi:.6e} rad  ({d_phi_deg:.3f} deg)")

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(X_py_kpc, Y_py_kpc, "o-", label=f"Python Quinn orbit ({nsteps} steps)")
    ax.plot(changa_x, changa_y, "o-", color="red", label=f"ChaNGa orbit ({nsteps} steps)")
    ax.scatter(changa_x[-1], changa_y[-1], color="red", s=70, zorder=5)
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    if args.output_figure:
        fig.savefig(args.output_figure)
        print(f"Saved plot to {args.output_figure}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    debug_changa_initial_states()
    main()
