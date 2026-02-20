#!/usr/bin/env python3
"""Wrapper entrypoint for encounter IC generation.

Example:
  python3 micro-test-encounter/make_encounter_tipsy.py --mode y-far --b-factor 1.0 --far-factor 10 --out third_party/ChaNGa/single_particle.std
"""

from __future__ import annotations

import runpy
from pathlib import Path


def main() -> None:
    target = Path(__file__).resolve().parent / "tools" / "make_encounter_tipsy.py"
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()
