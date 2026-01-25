#!/usr/bin/env python3
"""Generate a single-particle Tipsy binary file for ChaNGa Hill tests.

Defaults match the epicycle IC used in task2/plot_task2.py: (x0,y0,vx0,vy0)=(1,0,0,-2).
"""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a big-endian Tipsy binary with one collisionless particle.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("third_party/ChaNGa/single_particle.std"),
        help="Destination file (default: %(default)s)",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=0.0,
        help="Snapshot time stored in the Tipsy header (default: %(default)s)",
    )
    parser.add_argument(
        "--x0",
        type=float,
        default=1.0,
        help="Initial x position for the physical particle (code units, default: %(default)s)",
    )
    parser.add_argument(
        "--y0",
        type=float,
        default=0.0,
        help="Initial y position for the physical particle (default: %(default)s)",
    )
    parser.add_argument(
        "--z",
        type=float,
        default=0.0,
        help="Initial z position (default: %(default)s)",
    )
    parser.add_argument(
        "--vx0",
        type=float,
        default=0.0,
        help="Initial vx for the physical particle (code units, default: %(default)s)",
    )
    parser.add_argument(
        "--vy0",
        type=float,
        default=-2.0,
        help="Initial vy for the physical particle (default: %(default)s)",
    )
    parser.add_argument(
        "--vz",
        type=float,
        default=0.0,
        help="Initial z-velocity (default: %(default)s)",
    )
    parser.add_argument(
        "--omega",
        type=float,
        default=1.0,
        help="Hill angular frequency Ω used when vy is omitted (default: %(default)s)",
    )
    parser.add_argument("--mass", type=float, default=1.0, help="Particle mass")
    parser.add_argument(
        "--eps",
        type=float,
        default=1.0e-4,
        help="Softening length for the dark particle record (default: %(default)s)",
    )
    parser.add_argument(
        "--phi",
        type=float,
        default=0.0,
        help="Potential placeholder stored in the Tipsy record (default: %(default)s)",
    )
    parser.add_argument(
        "--add-guard",
        action="store_true",
        help="Include a tiny-mass guard particle to keep domain decomposition happy",
    )
    parser.add_argument(
        "--guard-mass",
        type=float,
        default=1.0e-12,
        help="Mass of the optional guard particle (default: %(default)s)",
    )
    parser.add_argument(
        "--guard-eps",
        type=float,
        default=1.0e-4,
        help="Softening length for the guard particle (default: %(default)s)",
    )
    return parser


def _pack_header(time: float, nbodies: int, nsph: int, ndark: int, nstar: int) -> bytes:
    """Return the standard Tipsy header (big-endian, padded)."""
    ndim = 3
    pad = 0
    return struct.pack(">d6i", time, nbodies, ndim, nsph, ndark, nstar, pad)


def _pack_dark_particle(
    mass: float,
    pos: tuple[float, float, float],
    vel: tuple[float, float, float],
    eps: float,
    phi: float,
) -> bytes:
    """Pack a single collisionless (dark) particle in big-endian Tipsy order."""
    buf = bytearray()
    buf.extend(struct.pack(">f", float(mass)))
    buf.extend(struct.pack(">3f", float(pos[0]), float(pos[1]), float(pos[2])))
    buf.extend(struct.pack(">3f", float(vel[0]), float(vel[1]), float(vel[2])))
    buf.extend(struct.pack(">f", float(eps)))
    buf.extend(struct.pack(">f", float(phi)))
    return bytes(buf)


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    nbodies = 1 + (1 if args.add_guard else 0)
    header = _pack_header(args.time, nbodies, 0, nbodies, 0)

    output_path = args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as fh:
        fh.write(header)
        fh.write(
            _pack_dark_particle(
                args.mass,
                (args.x0, args.y0, args.z),
                (args.vx0, args.vy0, args.vz),
                args.eps,
                args.phi,
            )
        )
        if args.add_guard:
            fh.write(
                _pack_dark_particle(
                    args.guard_mass,
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    args.guard_eps,
                    0.0,
                )
            )

    print("Physical particle:")
    print(
        f"  mass={args.mass:.6e}, "
        f"pos=({args.x0:.6f}, {args.y0:.6f}, {args.z:.6f}), "
        f"vel=({args.vx0:.6f}, {args.vy0:.6f}, {args.vz:.6f})"
    )
    if args.add_guard:
        print("Guard particle:")
        print(
            f"  mass={args.guard_mass:.6e}, "
            "pos=(0.000000, 0.000000, 0.000000), "
            "vel=(0.000000, 0.000000, 0.000000)"
        )
        print(f"Wrote {output_path} with one physical particle + guard.")
    else:
        print(f"Wrote {output_path} with one physical particle.")


if __name__ == "__main__":
    main()
