#!/usr/bin/env python3
"""Plot ChaNGa time series for shearing-sheet tests.

One-command workflow (from repo root):
  python3 tools/make_single_particle_tipsy.py --epicycle-boundary-test --D 1.0 --q 1.5 --Omega 1.0 --out third_party/ChaNGa/single_particle.std
  cd third_party/ChaNGa && ./ChaNGa single_particle_hill.param +p1 && cd ../..
  MPLBACKEND=Agg python3 plot_time_series.py --Lx 1.4 --q 1.5 --Omega 1.0 --epicycle --output epi_boundary_test.png --no-show
"""

from __future__ import annotations

import argparse
import glob
import re
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


def mark_start_end(ax, xvals: np.ndarray, yvals: np.ndarray, color):
    """Highlight first and last points on phase plots."""
    if len(xvals) == 0:
        return
    ax.scatter([xvals[0]], [yvals[0]], s=95, marker="*", color=color, edgecolors="black", linewidths=0.6, zorder=7)
    ax.text(xvals[0], yvals[0], "S", fontsize=9, color="black", ha="left", va="bottom")
    ax.scatter([xvals[-1]], [yvals[-1]], s=60, marker="o", color=color, edgecolors="black", linewidths=0.6, zorder=7)
    ax.text(xvals[-1], yvals[-1], "E", fontsize=9, color="black", ha="left", va="bottom")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot ChaNGa time series and epicycle diagnostics.")
    parser.add_argument("--x0", type=float, default=DEFAULT_IC["x"])
    parser.add_argument("--y0", type=float, default=DEFAULT_IC["y"])
    parser.add_argument("--vx0", type=float, default=DEFAULT_IC["vx"])
    parser.add_argument("--vy0", type=float, default=DEFAULT_IC["vy"])
    parser.add_argument(
        "--snap-pattern",
        type=str,
        default="third_party/ChaNGa/single_particle_hill.000*",
        help="Glob pattern for ChaNGa snapshots",
    )
    parser.add_argument(
        "--ic-file",
        type=Path,
        default=Path("third_party/ChaNGa/single_particle.std"),
        help="Tipsy IC file used to seed t=0",
    )
    parser.add_argument(
        "--param-file",
        type=Path,
        default=Path("third_party/ChaNGa/single_particle_hill.param"),
        help="ChaNGa param file (for dt/nSteps/dyPeriod parsing)",
    )
    parser.add_argument(
        "--mass-threshold",
        type=float,
        default=MASS_THRESHOLD,
        help="Mass threshold for selecting physical particles",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("single_particle_timeseries.png"),
        help="Wrapped/time-series output figure path",
    )
    parser.add_argument("--Lx", type=float, default=3.8, help="x-period for wrap/unwrap")
    parser.add_argument("--q", type=float, default=1.5, help="Shear parameter q")
    parser.add_argument("--Omega", type=float, default=1.0, help="Angular frequency Omega")
    parser.add_argument(
        "--D",
        type=float,
        default=None,
        help="Expected epicycle scale D for printed diagnostics (default: infer as Lx/1.4)",
    )
    parser.add_argument("--epicycle", action="store_true", help="Also produce epicycle-space diagnostics")
    parser.add_argument(
        "--epicycle-output",
        type=Path,
        default=None,
        help="Optional output path for epicycle figure (default: <output>_epicycle.png)",
    )
    parser.add_argument(
        "--boundary-output",
        type=Path,
        default=None,
        help="Optional output path for boundary-view figure (default: <output>_boundary_view.png)",
    )
    parser.add_argument(
        "--epicycle-jump-plot",
        action="store_true",
        help="Generate ChaNGa-only wrapped epicycle jump plot for the primary physical particle",
    )
    parser.add_argument(
        "--jump-style",
        type=str,
        choices=["break", "dashed"],
        default="dashed",
        help="How crossing jumps are drawn in jump plot",
    )
    parser.add_argument(
        "--jump-output",
        type=Path,
        default=None,
        help="Optional output path for jump plot (default: <output>_epicycle_jump.png)",
    )
    parser.add_argument("--no-show", action="store_true", help="Skip plt.show()")
    parser.add_argument(
        "--debug-load",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Print snapshot ingestion/selection debug tables",
    )
    return parser.parse_args()


def discover_snapshot_paths(pattern: str, debug: bool = True):
    """Find and numerically sort snapshot files by dump index."""
    paths = glob.glob(pattern)
    tried = [pattern]
    if not paths:
        p = Path(pattern)
        alt1 = p.name
        if alt1 != pattern:
            tried.append(alt1)
            paths = glob.glob(alt1)
    if not paths:
        p = Path(pattern)
        alt2 = str(Path("third_party/ChaNGa") / p.name)
        if alt2 not in tried:
            tried.append(alt2)
            paths = glob.glob(alt2)

    parsed = []
    for path in paths:
        m = re.search(r"\.(\d+)$", Path(path).name)
        if m:
            parsed.append((int(m.group(1)), path))
    if not parsed:
        raise FileNotFoundError(
            "No snapshots found. "
            f"Tried patterns: {tried}. Current cwd: {Path.cwd()}. "
            "Run from repo root, or pass --snap-pattern explicitly."
        )

    parsed.sort(key=lambda it: it[0])
    step_ids = np.array([p[0] for p in parsed], dtype=int)
    sorted_paths = [p[1] for p in parsed]

    if debug:
        print(f"Snapshot discovery: found {len(sorted_paths)} files")
        print(f"  first: {Path(sorted_paths[0]).name}")
        print(f"  last : {Path(sorted_paths[-1]).name}")
        print(f"  sorted indices: {step_ids.tolist()}")

    if len(sorted_paths) < 10:
        raise RuntimeError(
            f"WARNING: only {len(sorted_paths)} snapshots found (<10). "
            "Stop: series is too short for boundary/epicycle diagnostics."
        )

    return step_ids, sorted_paths


def _candidate_indices(snap: pb.SimSnap, mass_threshold: float):
    masses = np.asarray(snap["mass"], dtype=float)
    cand = np.where(masses > mass_threshold)[0]
    reason = f"mass > {mass_threshold:g}"
    if cand.size == 0:
        cand = np.where(masses > 0.0)[0]
        reason = "fallback mass > 0"
    return cand, masses, reason


def _to_code_units(snap: pb.SimSnap):
    pos_code = np.asarray(snap["pos"].in_units("kpc"), dtype=float) / DKPC_UNIT
    vel_code = np.asarray(snap["vel"].in_units("km s**-1"), dtype=float) / VELOCITY_UNIT_KM_S
    return pos_code, vel_code


def _match_by_nearest(prev_pos: np.ndarray, curr_pos: np.ndarray) -> np.ndarray:
    n = prev_pos.shape[0]
    used = set()
    order = []
    for i in range(n):
        d2 = np.sum((curr_pos - prev_pos[i]) ** 2, axis=1)
        for j in np.argsort(d2):
            if j not in used:
                used.add(j)
                order.append(j)
                break
    return np.array(order, dtype=int)


def load_changa_series(
    snap_pattern: str,
    ic_file: Path,
    mass_threshold: float,
    debug: bool,
):
    step_ids_snap, paths = discover_snapshot_paths(snap_pattern, debug=debug)

    std = pb.load(str(ic_file), fmt="tipsy")
    pos0_all = np.asarray(std["pos"], dtype=float)
    vel0_all = np.asarray(std["vel"], dtype=float)
    ids0 = np.asarray(std["iord"], dtype=int) if "iord" in std.loadable_keys() else None
    cand0, masses0, reason0 = _candidate_indices(std, mass_threshold)

    if cand0.size == 0:
        raise RuntimeError("No selectable particles in IC file")

    single_target = cand0.size == 1
    selected_iord = None

    if single_target:
        idx0 = int(cand0[0])
        if ids0 is not None:
            selected_iord = int(ids0[idx0])
        if debug:
            print(
                f"Particle selection: single physical particle detected in IC ({reason0}); "
                f"selected index={idx0}, mass={masses0[idx0]:.6e}, iord={selected_iord}"
            )

        xs = [np.array([pos0_all[idx0, 0]], dtype=float)]
        ys = [np.array([pos0_all[idx0, 1]], dtype=float)]
        vxs = [np.array([vel0_all[idx0, 0]], dtype=float)]
        vys = [np.array([vel0_all[idx0, 1]], dtype=float)]
        selected_reason = ["IC single physical particle"]

        for path in paths:
            snap = pb.load(path)
            snap.physical_units()
            pos_code, vel_code = _to_code_units(snap)
            ids = np.asarray(snap["iord"], dtype=int) if "iord" in snap.loadable_keys() else None
            cand, masses, reason = _candidate_indices(snap, mass_threshold)
            if cand.size == 0:
                raise RuntimeError(f"No selectable particles in snapshot {path}")

            chosen_reason = "max-mass among candidates"
            if selected_iord is not None and ids is not None:
                match = cand[ids[cand] == selected_iord]
                if match.size > 0:
                    idx = int(match[0])
                    chosen_reason = f"matched iord={selected_iord}"
                else:
                    idx = int(cand[np.argmax(masses[cand])])
                    chosen_reason = f"iord missing; fallback {chosen_reason}"
            else:
                idx = int(cand[np.argmax(masses[cand])])
                chosen_reason = f"{chosen_reason} ({reason})"

            xs.append(np.array([pos_code[idx, 0]], dtype=float))
            ys.append(np.array([pos_code[idx, 1]], dtype=float))
            vxs.append(np.array([vel_code[idx, 0]], dtype=float))
            vys.append(np.array([vel_code[idx, 1]], dtype=float))
            selected_reason.append(chosen_reason)

        xs = np.vstack(xs)
        ys = np.vstack(ys)
        vxs = np.vstack(vxs)
        vys = np.vstack(vys)
        step_ids = np.concatenate([[0], step_ids_snap])

    else:
        if debug:
            print(
                f"Particle selection: {cand0.size} IC particles selected by {reason0}; "
                "tracking all selected particles"
            )
        if ids0 is not None:
            order0 = np.argsort(ids0[cand0])
        else:
            order0 = np.argsort(pos0_all[cand0, 0])

        idxs0 = cand0[order0]
        xs = [pos0_all[idxs0, 0]]
        ys = [pos0_all[idxs0, 1]]
        vxs = [vel0_all[idxs0, 0]]
        vys = [vel0_all[idxs0, 1]]

        prev_pos = pos0_all[idxs0, :2]
        prev_ids = ids0[idxs0] if ids0 is not None else None

        for path in paths:
            snap = pb.load(path)
            snap.physical_units()
            pos_code, vel_code = _to_code_units(snap)
            ids = np.asarray(snap["iord"], dtype=int) if "iord" in snap.loadable_keys() else None
            cand, _, _ = _candidate_indices(snap, mass_threshold)
            if cand.size == 0:
                raise RuntimeError(f"No selectable particles in snapshot {path}")

            curr_pos = pos_code[cand, :2]
            curr_vel = vel_code[cand, :]
            if ids is not None and prev_ids is not None:
                curr_ids = ids[cand]
                id_to_row = {int(pid): i for i, pid in enumerate(curr_ids)}
                if all(int(pid) in id_to_row for pid in prev_ids):
                    order = np.array([id_to_row[int(pid)] for pid in prev_ids], dtype=int)
                else:
                    order = _match_by_nearest(prev_pos, curr_pos)
            else:
                order = _match_by_nearest(prev_pos, curr_pos)

            curr_pos = curr_pos[order]
            curr_vel = curr_vel[order]

            xs.append(curr_pos[:, 0])
            ys.append(curr_pos[:, 1])
            vxs.append(curr_vel[:, 0])
            vys.append(curr_vel[:, 1])

            prev_pos = curr_pos
            if ids is not None and prev_ids is not None:
                prev_ids = prev_ids

        xs = np.vstack(xs)
        ys = np.vstack(ys)
        vxs = np.vstack(vxs)
        vys = np.vstack(vys)
        step_ids = np.concatenate([[0], step_ids_snap])
        selected_reason = None

    return {
        "step_ids": step_ids,
        "x_raw": xs,
        "y_raw": ys,
        "vx_raw": vxs,
        "vy_raw": vys,
        "single_target": single_target,
        "selected_reason": selected_reason,
        "primary_col": 0,
    }


def read_param_values(param_path: Path):
    nsteps = None
    dt = None
    dy_period = None
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
                elif line.startswith("dyPeriod"):
                    dy_period = float(line.split("=")[1].strip())
    except FileNotFoundError:
        return None, None, None
    return nsteps, dt, dy_period


def wrap_to_box(values: np.ndarray, period: float) -> np.ndarray:
    return ((values + period / 2.0) % period) - period / 2.0


def unwrap_x_and_crossings(x_wrap: np.ndarray, Lx: float) -> tuple[np.ndarray, np.ndarray]:
    x_unwrap = np.copy(x_wrap)
    crossing = np.zeros_like(x_wrap, dtype=bool)
    for p in range(x_wrap.shape[1]):
        nshift = 0
        for i in range(1, x_wrap.shape[0]):
            dx = x_wrap[i, p] - x_wrap[i - 1, p]
            if dx < -Lx / 2.0:
                nshift += 1
                crossing[i, p] = True
            elif dx > Lx / 2.0:
                nshift -= 1
                crossing[i, p] = True
            x_unwrap[i, p] = x_wrap[i, p] + nshift * Lx
    return x_unwrap, crossing


def build_y_wrap_from_x_crossings(
    y_raw: np.ndarray,
    x_wrap: np.ndarray,
    times: np.ndarray,
    q: float,
    omega: float,
    Lx: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply shear remap offsets to y at x-wrap events to produce y_wrap."""
    y_wrap = np.copy(y_raw)
    y_remap = np.zeros_like(y_raw)

    for p in range(y_raw.shape[1]):
        shift = 0.0
        y_wrap[0, p] = y_raw[0, p]
        for i in range(1, y_raw.shape[0]):
            dx = x_wrap[i, p] - x_wrap[i - 1, p]
            if abs(dx) > Lx / 2.0:
                direction = 1.0 if dx < -Lx / 2.0 else -1.0
                remap = direction * q * omega * Lx * times[i]
                shift += remap
                y_remap[i, p] = remap
            y_wrap[i, p] = y_raw[i, p] + shift

    return y_wrap, y_remap


def unwrap_periodic_series(
    values: np.ndarray,
    dy_period: float | None,
) -> np.ndarray:
    """Optionally unwrap periodic y-series; keeps raw values if no period provided."""
    y_unwrap = np.copy(values)
    if dy_period is None or dy_period <= 0:
        return y_unwrap

    half = dy_period / 2.0
    for p in range(y_unwrap.shape[1]):
        nshift = 0
        for i in range(1, y_unwrap.shape[0]):
            dy = y_unwrap[i, p] - y_unwrap[i - 1, p]
            if dy < -half:
                nshift += 1
            elif dy > half:
                nshift -= 1
            y_unwrap[i, p] = y_unwrap[i, p] + nshift * dy_period
    return y_unwrap


def print_raw_samples(
    step_ids: np.ndarray,
    times: np.ndarray,
    x_raw: np.ndarray,
    y_raw: np.ndarray,
    vx_raw: np.ndarray,
    vy_raw: np.ndarray,
    debug: bool,
):
    if not debug:
        return
    n = len(step_ids)
    first = list(range(min(5, n)))
    last = list(range(max(0, n - 2), n))
    idxs = first + [i for i in last if i not in first]

    print("Raw selected-particle samples (no wrap/unwrap):")
    print("  step_id      t             x_raw           y_raw           vx_raw          vy_raw")
    for i in idxs:
        print(
            f"  {step_ids[i]:>6d}  {times[i]:.6f}  {x_raw[i,0]:+.6e}  {y_raw[i,0]:+.6e}  "
            f"{vx_raw[i,0]:+.6e}  {vy_raw[i,0]:+.6e}"
        )


def print_shear_slope_table(
    times: np.ndarray,
    x_ref: np.ndarray,
    y_unwrap: np.ndarray,
    q: float,
    omega: float,
    nfit: int = 10,
):
    nfit = min(nfit, len(times))
    if nfit < 2:
        return
    t_fit = times[:nfit]
    print(f"Initial shear-slope check (first {nfit} samples):")
    print("  pid      x_ref             expected            measured            diff")
    for p in range(x_ref.shape[0]):
        expected = -q * omega * x_ref[p]
        measured = float(np.polyfit(t_fit, y_unwrap[:nfit, p], 1)[0])
        diff = measured - expected
        print(f"  p{p:<2d}  {x_ref[p]:+.6e}  {expected:+.6e}  {measured:+.6e}  {diff:+.6e}")


def print_crossing_table(
    step_ids: np.ndarray,
    times: np.ndarray,
    x_wrap: np.ndarray,
    crossing: np.ndarray,
    y_remap: np.ndarray,
):
    print("Boundary crossing events (from x_wrap jumps):")
    print("  event  pid  step_id      t           dx_wrap        x_before       x_after        dy_remap")
    event = 0
    for p in range(x_wrap.shape[1]):
        idx = np.where(crossing[:, p])[0]
        for i in idx:
            event += 1
            dx = x_wrap[i, p] - x_wrap[i - 1, p]
            print(
                f"  {event:>5d}  p{p:<2d} {step_ids[i]:>7d}  {times[i]:.6f}  "
                f"{dx:+.6e}  {x_wrap[i-1,p]:+.6e}  {x_wrap[i,p]:+.6e}  {y_remap[i,p]:+.6e}"
            )
    if event == 0:
        print("  (none)")


def default_epicycle_output(output: Path) -> Path:
    if output.suffix:
        return output.with_name(f"{output.stem}_epicycle{output.suffix}")
    return output.with_name(f"{output.name}_epicycle.png")


def default_boundary_output(output: Path) -> Path:
    if output.suffix:
        return output.with_name(f"{output.stem}_boundary_view{output.suffix}")
    return output.with_name(f"{output.name}_boundary_view.png")


def default_jump_output(output: Path) -> Path:
    if output.suffix:
        return output.with_name(f"{output.stem}_epicycle_jump{output.suffix}")
    return output.with_name(f"{output.name}_epicycle_jump.png")


def plot_wrapped_figure(
    times: np.ndarray,
    x_raw: np.ndarray,
    x_wrap: np.ndarray,
    y_raw: np.ndarray,
    y_wrap: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    crossing: np.ndarray,
    output: Path,
    show: bool,
):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex="col")
    colors = plt.cm.tab10(np.linspace(0, 1, x_raw.shape[1]))
    crossing_times = sorted(set(times[np.where(crossing)[0]].tolist()))

    for p in range(x_raw.shape[1]):
        c = colors[p]
        axes[0, 0].plot(times, x_raw[:, p], color=c, lw=1.2, label=f"p{p}")
        axes[0, 1].plot(times, x_wrap[:, p], color=c, lw=1.2)
        axes[0, 2].plot(times, y_raw[:, p], color=c, lw=1.2)
        axes[0, 2].plot(times, y_wrap[:, p], color=c, lw=1.0, ls="--", alpha=0.8)
        axes[1, 0].plot(times, vx[:, p], color=c, lw=1.2)
        axes[1, 1].plot(times, vy[:, p], color=c, lw=1.2)
        axes[1, 2].plot(x_raw[:, p], y_raw[:, p], color=c, lw=1.2)
        mark_start_end(axes[1, 2], x_raw[:, p], y_raw[:, p], c)

        idx = np.where(crossing[:, p])[0]
        if idx.size:
            axes[1, 2].scatter(x_raw[idx, p], y_raw[idx, p], color=c, s=16)

    for ax in [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1]]:
        for tc in crossing_times:
            ax.axvline(tc, color="gray", lw=0.8, ls=":", alpha=0.4)

    for ax in axes.flat:
        ax.grid(True, alpha=0.3)

    axes[0, 0].set_title("x_raw(t)")
    axes[0, 1].set_title("x_wrap(t)")
    axes[0, 2].set_title("y_raw(t) and y_wrap(t)")
    axes[1, 0].set_title("vx(t)")
    axes[1, 1].set_title("vy(t)")
    axes[1, 2].set_title("x_raw vs y_raw")

    axes[0, 0].set_ylabel("x_raw")
    axes[0, 1].set_ylabel("x_wrap")
    axes[0, 2].set_ylabel("y_raw")
    axes[1, 0].set_ylabel("vx")
    axes[1, 1].set_ylabel("vy")
    axes[1, 2].set_ylabel("y_raw")

    axes[1, 0].set_xlabel("time [code]")
    axes[1, 1].set_xlabel("time [code]")
    axes[1, 2].set_xlabel("x_raw")
    axes[0, 0].legend(loc="best", fontsize=8)

    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print(f"Saved time-series plot to {output}")
    if show:
        plt.show()


def plot_boundary_view_figure(
    x_raw: np.ndarray,
    y_raw: np.ndarray,
    x_wrap: np.ndarray,
    y_wrap: np.ndarray,
    crossing: np.ndarray,
    Lx: float,
    output: Path,
    show: bool,
):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, x_raw.shape[1]))

    # A) Raw epicycle ellipse
    for p in range(x_raw.shape[1]):
        axes[0].plot(x_raw[:, p], y_raw[:, p], color=colors[p], lw=1.2, label=f"p{p}")
        mark_start_end(axes[0], x_raw[:, p], y_raw[:, p], colors[p])
    axes[0].axvline(-Lx / 2.0, color="k", ls="--", lw=1.0, alpha=0.6)
    axes[0].axvline(+Lx / 2.0, color="k", ls="--", lw=1.0, alpha=0.6)
    axes[0].set_title("Raw epicycle: x_raw vs y_raw")
    axes[0].set_xlabel("x_raw")
    axes[0].set_ylabel("y_raw")
    axes[0].grid(True, alpha=0.3)

    # B) Wrapped patch trajectory with crossing markers
    for p in range(x_wrap.shape[1]):
        axes[1].plot(x_wrap[:, p], y_wrap[:, p], color=colors[p], lw=1.2, label=f"p{p}")
        mark_start_end(axes[1], x_wrap[:, p], y_wrap[:, p], colors[p])
        idx = np.where(crossing[:, p])[0]
        for i in idx:
            axes[1].scatter(x_wrap[i - 1, p], y_wrap[i - 1, p], color="black", s=18, zorder=4)
            axes[1].scatter(x_wrap[i, p], y_wrap[i, p], color="black", s=22, marker="x", zorder=5)
            axes[1].annotate(
                "",
                xy=(x_wrap[i, p], y_wrap[i, p]),
                xytext=(x_wrap[i - 1, p], y_wrap[i - 1, p]),
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8, alpha=0.7),
            )
    axes[1].axvline(-Lx / 2.0, color="k", ls="--", lw=1.0, alpha=0.6)
    axes[1].axvline(+Lx / 2.0, color="k", ls="--", lw=1.0, alpha=0.6)
    axes[1].set_xlim(-Lx / 2.0, Lx / 2.0)
    axes[1].set_title("Wrapped patch: x_wrap vs y_wrap")
    axes[1].set_xlabel("x_wrap")
    axes[1].set_ylabel("y_wrap")
    axes[1].grid(True, alpha=0.3)

    # C) Tiled patch view
    for p in range(x_wrap.shape[1]):
        yv = y_wrap[:, p]
        xv = x_wrap[:, p]
        axes[2].plot(xv - Lx, yv, color=colors[p], lw=1.0, alpha=0.6)
        axes[2].plot(xv, yv, color=colors[p], lw=1.2)
        axes[2].plot(xv + Lx, yv, color=colors[p], lw=1.0, alpha=0.6)
        mark_start_end(axes[2], xv, yv, colors[p])
    for xline in (-1.5 * Lx, -0.5 * Lx, 0.5 * Lx, 1.5 * Lx):
        axes[2].axvline(xline, color="k", ls="--", lw=1.0, alpha=0.5)
    axes[2].set_xlim(-1.6 * Lx, 1.6 * Lx)
    axes[2].set_title("Tiled patch view")
    axes[2].set_xlabel("x_tile")
    axes[2].set_ylabel("y_wrap")
    axes[2].grid(True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print(f"Saved boundary-view plot to {output}")
    if show:
        plt.show()


def print_jump_crossing_table(
    step_ids: np.ndarray,
    times: np.ndarray,
    x_wrap_1d: np.ndarray,
    y_raw_1d: np.ndarray,
    y_wrap_1d: np.ndarray,
    crossing_1d: np.ndarray,
    y_remap_1d: np.ndarray,
):
    print("Epicycle-jump crossing events (primary particle):")
    print("  event   i-1    i  step_id      t           x_before       x_after        dx_wrap        y_raw(i-1)     y_raw(i)       dy_remap       y_wrap(i)")
    idx = np.where(crossing_1d)[0]
    if idx.size == 0:
        print("  (none)")
        return
    for k, i in enumerate(idx, start=1):
        dx = x_wrap_1d[i] - x_wrap_1d[i - 1]
        print(
            f"  {k:>5d}  {i-1:>4d} {i:>4d} {step_ids[i]:>7d}  {times[i]:.6f}  "
            f"{x_wrap_1d[i-1]:+.6e}  {x_wrap_1d[i]:+.6e}  {dx:+.6e}  "
            f"{y_raw_1d[i-1]:+.6e}  {y_raw_1d[i]:+.6e}  {y_remap_1d[i]:+.6e}  {y_wrap_1d[i]:+.6e}"
        )


def plot_epicycle_jump_figure(
    step_ids: np.ndarray,
    x_wrap_1d: np.ndarray,
    y_raw_1d: np.ndarray,
    y_wrap_1d: np.ndarray,
    crossing_1d: np.ndarray,
    Lx: float,
    jump_style: str,
    output: Path,
    show: bool,
):
    start = 1 if len(step_ids) > 1 else 0  # snapshots only (exclude IC t=0 row)
    step_plot = step_ids[start:]
    x_plot = x_wrap_1d[start:]
    y_raw_plot = y_raw_1d[start:]
    y_wrap_plot = y_wrap_1d[start:]
    cross_plot = crossing_1d[start:]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True)
    for ax, y_plot, title, ylabel in [
        (axes[0], y_raw_plot, "A) x_wrap vs y_raw", "y_raw"),
        (axes[1], y_wrap_plot, "B) x_wrap vs y_wrap", "y_wrap"),
    ]:
        ax.scatter(x_plot, y_plot, s=20, color="tab:blue", alpha=0.9, zorder=3)
        event_num = 0
        for i in range(1, len(x_plot)):
            x0, x1 = x_plot[i - 1], x_plot[i]
            y0, y1 = y_plot[i - 1], y_plot[i]
            if cross_plot[i]:
                event_num += 1
                if jump_style == "dashed":
                    ax.plot([x0, x1], [y0, y1], color="tab:orange", lw=1.3, ls="--", zorder=2)
                ax.scatter([x1], [y1], marker="x", s=95, color="black", zorder=4)
                ax.text(x1, y1, f"E{event_num}", fontsize=9, color="black", ha="left", va="bottom")
            else:
                ax.plot([x0, x1], [y0, y1], color="tab:blue", lw=1.2, zorder=1)
        mark_start_end(ax, x_plot, y_plot, "tab:blue")

        ax.axvline(-Lx / 2.0, color="k", ls="--", lw=1.0, alpha=0.7)
        ax.axvline(+Lx / 2.0, color="k", ls="--", lw=1.0, alpha=0.7)
        ax.set_xlim(-Lx / 2.0, Lx / 2.0)
        ax.set_xlabel("x_wrap")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"ChaNGa epicycle boundary jumps (N={len(step_plot)})", fontsize=12)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    print(f"Saved epicycle jump plot to {output}")
    if show:
        plt.show()


def plot_epicycle_figure(
    times: np.ndarray,
    x_unwrap: np.ndarray,
    y_unwrap_raw: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    q: float,
    omega: float,
    output: Path,
    show: bool,
):
    x_ref = x_unwrap[0, :]
    y_epi = y_unwrap_raw + (q * omega * x_ref)[None, :] * times[:, None]
    vy_epi = vy + (q * omega * x_ref)[None, :]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, x_unwrap.shape[1]))
    for p in range(x_unwrap.shape[1]):
        axes[0].plot(x_unwrap[:, p], y_epi[:, p], color=colors[p], lw=1.2, label=f"p{p}")
        axes[1].plot(vx[:, p], vy_epi[:, p], color=colors[p], lw=1.2, label=f"p{p}")
        mark_start_end(axes[0], x_unwrap[:, p], y_epi[:, p], colors[p])
        mark_start_end(axes[1], vx[:, p], vy_epi[:, p], colors[p])
        y_range = np.ptp(y_epi[:, p])
        print(f"Epicycle p{p}: x0={x_ref[p]:+.6e}, y_epi_range={y_range:.6e}")

    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y_epi")
    axes[0].set_title(f"x vs y_epi (q={q}, Omega={omega})")
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("vx")
    axes[1].set_ylabel("vy_epi")
    axes[1].set_title(f"vx vs vy_epi (q={q}, Omega={omega})")
    axes[1].grid(True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    print(f"Saved epicycle plot to {output}")
    if show:
        plt.show()


def hill_time_series(x0, y0, vx0, vy0, nsteps, dt):
    ic = {"x": x0, "y": y0, "z": 0.0, "vx": vx0, "vy": vy0, "vz": 0.0}
    x, y, vx, vy = run_quinn_states_for_ic(ic, dt=dt, nsteps=nsteps)
    return np.array(x), np.array(y), np.array(vx), np.array(vy)


def main():
    args = parse_args()

    data = load_changa_series(
        snap_pattern=args.snap_pattern,
        ic_file=args.ic_file,
        mass_threshold=args.mass_threshold,
        debug=args.debug_load,
    )

    step_ids = data["step_ids"]
    x_raw = data["x_raw"]
    y_raw = data["y_raw"]
    vx = data["vx_raw"]
    vy = data["vy_raw"]

    param_nsteps, param_dt, dy_period = read_param_values(args.param_file)
    nsteps_from_files = len(step_ids) - 1
    dt = param_dt if param_dt is not None else (2.0 * np.pi / max(nsteps_from_files, 1))

    nsteps = len(step_ids)
    if param_nsteps is not None:
        nsteps = min(nsteps, param_nsteps + 1)

    step_ids = step_ids[:nsteps]
    x_raw = x_raw[:nsteps]
    y_raw = y_raw[:nsteps]
    vx = vx[:nsteps]
    vy = vy[:nsteps]
    times = dt * np.arange(nsteps)

    print(f"Using nsteps={nsteps}, dt={dt:.6f} (code units)")

    if data["single_target"]:
        print_raw_samples(step_ids, times, x_raw, y_raw, vx, vy, debug=args.debug_load)

    x_wrap = wrap_to_box(x_raw, args.Lx)
    x_unwrap, crossing = unwrap_x_and_crossings(x_wrap, args.Lx)
    y_wrap, y_remap = build_y_wrap_from_x_crossings(
        y_raw=y_raw,
        x_wrap=x_wrap,
        times=times,
        q=args.q,
        omega=args.Omega,
        Lx=args.Lx,
    )
    y_unwrap_raw = unwrap_periodic_series(y_raw, dy_period)

    x_ref = x_unwrap[0, :]
    print_crossing_table(step_ids, times, x_wrap, crossing, y_remap)
    print_shear_slope_table(times, x_ref, y_unwrap_raw, args.q, args.Omega, nfit=10)

    d_use = args.D if args.D is not None else args.Lx / 1.4
    print(
        f"Expected ellipse scales (q={args.q}, Omega={args.Omega}): "
        f"x amplitude ~ {d_use:.6f}, y_epi amplitude ~ {2.0 * d_use:.6f}"
    )

    if x_raw.shape[1] == 1:
        # Explicit mismatch check: stale snapshots or wrong run/config are common.
        dx0 = abs(x_raw[1, 0] - x_raw[0, 0]) if nsteps > 1 else 0.0
        dy0 = abs(y_raw[1, 0] - y_raw[0, 0]) if nsteps > 1 else 0.0
        if dx0 > 0.25 * args.Lx or dy0 > 5.0:
            print(
                "WARNING: first snapshot differs strongly from IC t=0. "
                "This often means stale outputs or a run/config mismatch."
            )

        cha_range = float(np.ptp(x_raw[:, 0]))
        cha_range_after_step0 = float(np.ptp(x_raw[1:, 0])) if nsteps > 2 else cha_range
        if cha_range < 1e-10 or cha_range_after_step0 < 1e-10:
            print("ChaNGa raw x is constant -> issue is in simulation output or run config, not unwrap.")

    plot_wrapped_figure(
        times=times,
        x_raw=x_raw,
        x_wrap=x_wrap,
        y_raw=y_raw,
        y_wrap=y_wrap,
        vx=vx,
        vy=vy,
        crossing=crossing,
        output=args.output,
        show=not args.no_show,
    )

    boundary_out = args.boundary_output if args.boundary_output is not None else default_boundary_output(args.output)
    plot_boundary_view_figure(
        x_raw=x_raw,
        y_raw=y_raw,
        x_wrap=x_wrap,
        y_wrap=y_wrap,
        crossing=crossing,
        Lx=args.Lx,
        output=boundary_out,
        show=not args.no_show,
    )

    if args.epicycle:
        epi_out = args.epicycle_output if args.epicycle_output is not None else default_epicycle_output(args.output)
        plot_epicycle_figure(
            times=times,
            x_unwrap=x_unwrap,
            y_unwrap_raw=y_unwrap_raw,
            vx=vx,
            vy=vy,
            q=args.q,
            omega=args.Omega,
            output=epi_out,
            show=not args.no_show,
        )

    if args.epicycle_jump_plot:
        primary_col = int(data.get("primary_col", 0))
        xw = x_wrap[:, primary_col]
        yr = y_raw[:, primary_col]
        yw = y_wrap[:, primary_col]
        cross = crossing[:, primary_col]
        yrm = y_remap[:, primary_col]
        print_jump_crossing_table(step_ids, times, xw, yr, yw, cross, yrm)
        jump_out = args.jump_output if args.jump_output is not None else default_jump_output(args.output)
        plot_epicycle_jump_figure(
            step_ids=step_ids,
            x_wrap_1d=xw,
            y_raw_1d=yr,
            y_wrap_1d=yw,
            crossing_1d=cross,
            Lx=args.Lx,
            jump_style=args.jump_style,
            output=jump_out,
            show=not args.no_show,
        )


if __name__ == "__main__":
    main()
