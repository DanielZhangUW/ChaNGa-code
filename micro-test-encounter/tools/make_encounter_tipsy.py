#!/usr/bin/env python3
"""Generate 2-particle Hill encounter ICs (Tipsy binary) for ChaNGa.

This is an isolated script for the encounter study and does not modify
existing single-particle epicycle/boundary workflows.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Create 2-particle Hill encounter IC")
    p.add_argument("--encounter-test", action="store_true", help="Compatibility flag (no longer required)")
    p.add_argument("--mode", choices=["x-far", "y-far"], default="y-far")
    p.add_argument("--m1", type=float, default=1.0e-3)
    p.add_argument("--m2", type=float, default=1.0e-10)
    p.add_argument("--dOrbdist", type=float, default=1.0)
    p.add_argument("--dCentMass", type=float, default=1.0)
    p.add_argument("--q", type=float, default=1.5)
    p.add_argument("--Omega", type=float, default=1.0)
    p.add_argument("--b-factor", type=float, default=1.0, help="Impact parameter factor in units of R_H")
    p.add_argument("--y-span", type=float, default=10.0, help="Start/end span in y for y-far mode (units of R_H)")
    p.add_argument(
        "--T-run",
        type=float,
        default=None,
        help="Encounter run duration used for default vy2; default is 2*pi/Omega",
    )
    p.add_argument("--far-factor", type=float, default=10.0, help="Far-start factor in units of R_H")
    p.add_argument("--x0-factor", type=float, default=None, help="Deprecated alias for x-far mode only")
    p.add_argument("--vx2", type=float, default=0.0)
    p.add_argument("--vy2-mode", choices=["shear", "manual"], default="shear")
    p.add_argument("--vy2", type=float, default=None, help="Override vy2 when provided")
    p.add_argument("--eps1", type=float, default=1.0e-4)
    p.add_argument("--eps2", type=float, default=1.0e-4)
    p.add_argument("--phi", type=float, default=0.0)
    p.add_argument("--time", type=float, default=0.0)
    p.add_argument(
        "--out",
        "--output",
        dest="output",
        type=Path,
        default=Path("third_party/ChaNGa/single_particle.std"),
    )
    p.add_argument("--add-guard", action="store_true", help="Optional tiny guard for DD stability")
    p.add_argument("--guard-mass", type=float, default=1.0e-30)
    p.add_argument("--guard-x", type=float, default=100.0)
    p.add_argument("--guard-y", type=float, default=100.0)
    p.add_argument("--guard-eps", type=float, default=1.0e-4)
    return p


def pack_header(time: float, nbodies: int) -> bytes:
    ndim = 3
    nsph = 0
    ndark = nbodies
    nstar = 0
    pad = 0
    return struct.pack(">d6i", time, nbodies, ndim, nsph, ndark, nstar, pad)


def pack_dark(mass: float, pos: tuple[float, float, float], vel: tuple[float, float, float], eps: float, phi: float) -> bytes:
    buf = bytearray()
    buf.extend(struct.pack(">f", float(mass)))
    buf.extend(struct.pack(">3f", float(pos[0]), float(pos[1]), float(pos[2])))
    buf.extend(struct.pack(">3f", float(vel[0]), float(vel[1]), float(vel[2])))
    buf.extend(struct.pack(">f", float(eps)))
    buf.extend(struct.pack(">f", float(phi)))
    return bytes(buf)


def main() -> None:
    args = build_parser().parse_args()

    rh = args.dOrbdist * (args.m1 / (3.0 * args.dCentMass)) ** (1.0 / 3.0)
    y_span = args.y_span
    if args.far_factor != 10.0:
        # Backward compatibility: allow old flag to override y-span.
        y_span = args.far_factor
    t_run = args.T_run if args.T_run is not None else (2.0 * 3.141592653589793 / args.Omega)
    if args.mode == "x-far":
        # Legacy geometry: start far on x-, with small y impact parameter.
        x_factor = args.x0_factor if args.x0_factor is not None else -args.far_factor
        x2 = x_factor * rh
        y2 = +args.b_factor * rh
    else:
        # Requested Mode Y-far: start at bottom (large negative y).
        y2 = -y_span * rh
        x2 = -args.b_factor * rh

    vx2 = args.vx2
    if args.vy2 is not None:
        vy2 = args.vy2
    else:
        # Default baseline: move from y=-y_span*R_H to +y_span*R_H over T_run.
        vy2 = (2.0 * y_span * rh) / t_run

    py2 = vy2 + 2.0 * args.Omega * x2
    y_end_pred = y2 + vy2 * t_run

    particles = [
        # Massive body at origin, zero velocity.
        (args.m1, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), args.eps1, args.phi),
        # Light encounter/test particle.
        (args.m2, (x2, y2, 0.0), (vx2, vy2, 0.0), args.eps2, args.phi),
    ]
    if args.add_guard:
        particles.append(
            (
                args.guard_mass,
                (args.guard_x, args.guard_y, 0.0),
                (0.0, 0.0, 0.0),
                args.guard_eps,
                0.0,
            )
        )

    out = args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as fh:
        fh.write(pack_header(args.time, len(particles)))
        for mass, pos, vel, eps, phi in particles:
            fh.write(pack_dark(mass, pos, vel, eps, phi))

    print("Encounter IC sanity (code units):")
    print(
        f"  mode = {args.mode}, b_factor={args.b_factor:.3f}, y_span={y_span:.3f}, "
        f"T_run={t_run:.6f}"
    )
    print(f"  R_H = {rh:.9e}  (m1={args.m1:.3e}, dOrbdist={args.dOrbdist:.3f}, dCentMass={args.dCentMass:.3f})")
    print("  Body 1 (massive):")
    print(f"    m1={args.m1:.3e}, x1=0, y1=0, vx1=0, vy1=0")
    print("  Body 2 (test):")
    print(f"    m2={args.m2:.3e}, x2={x2:+.9e}, y2={y2:+.9e}, vx2={vx2:+.9e}, vy2={vy2:+.9e}")
    print(f"    x0/R_H={x2/rh:+.6f}")
    print(f"    y0/R_H={y2/rh:+.6f}")
    print(f"    Py2={py2:+.9e}")
    print(f"    predicted y_end/R_H (drift model) = {y_end_pred/rh:+.6f} (target +{y_span:.6f})")
    if args.vy2 is None:
        print("    vy2 source = drift-target default (2*y_span*R_H/T_run)")
    else:
        print("    vy2 source = user override (--vy2)")
    if args.mode == "y-far":
        sign_ok = vy2 > 0.0
        print(f"    y-far sign check: vy2 > 0 ? {sign_ok} (vy2={vy2:+.6e})")
    if args.add_guard:
        print(
            "  Guard added: "
            f"m={args.guard_mass:.3e}, pos=({args.guard_x:.3f},{args.guard_y:.3f},0)"
        )

    print(f"Wrote {out} with {len(particles)} dark particles.")
    print("Tuning hint: if no strong encounter, try --b-factor 0.5 or 0.2; if too strong/capture, try --b-factor 2.0")


if __name__ == "__main__":
    main()
