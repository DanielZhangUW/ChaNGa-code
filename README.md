Single-Particle Hill Test (ChaNGa vs Quinn)
===========================================

Purpose
-------
This repo contains a Python reference integrator for Hill's equations and a
ChaNGa sliding-patch run card. The main goal is to compare a single-particle
orbit step-by-step in code units (Omega=1, 2*pi code time per orbit).

Key Files
---------
- third_party/ChaNGa/single_particle_hill.param
  ChaNGa run card (sliding patch, 100 steps per orbit).
- tools/make_single_particle_tipsy.py
  Generates third_party/ChaNGa/single_particle.std from a chosen IC.
- task2/plot_task2.py
  Quinn/Hill reference integrator (Python).
- task2/plot_time_series.py
  Plots x(t), y(t), vx(t), vy(t) for ChaNGa vs Quinn.

Quick Start (Epicycle IC)
-------------------------
From the repo root:

  # Defaults are epicycle (x0=1, y0=0, vx0=0, vy0=-2).
  python3 tools/make_single_particle_tipsy.py --add-guard
  cd third_party/ChaNGa
  ./ChaNGa single_particle_hill.param +p1
  cd ..
  python3 task2/plot_time_series.py --output task2/single_particle_timeseries.png

Expected Outputs
----------------
- third_party/ChaNGa/single_particle.std
  Tipsy input with one physical particle + guard particle.
- third_party/ChaNGa/single_particle_hill.000001 ... .000100
  ChaNGa snapshots (one per step).
- task2/single_particle_timeseries.png
  4-panel plot comparing Quinn vs ChaNGa time series.

Units and Time Step
-------------------
- Code units are chosen so Omega = 1 and one orbit is 2*pi in code time.
- Length unit: dKpcUnit = 4.848e-9 so 1 code length ~ 1 AU.
- Time step: dDelta = 2*pi/100 (100 steps per orbit).

Consistency Checks
------------------
- The IC defaults are consistent across:
  - tools/make_single_particle_tipsy.py
  - task2/plot_task2.py
  - task2/plot_time_series.py
- plot_time_series.py prints:
  - Quinn requested IC (code units)
  - ChaNGa first snapshot (code units)
  This is useful to confirm the run used the expected IC.

Common Issues
-------------
- "No snapshots found": ChaNGa did not run or outputs were deleted.
  Re-run ./ChaNGa single_particle_hill.param +p1.
- "ChaNGa step-1 snapshot does not match Quinn step-1":
  This indicates a physics mismatch between ChaNGa sliding patch and the
  Quinn integrator for the same IC and time step.
  The message is expected if the two trajectories diverge.

Notes
-----
- Units: dKpcUnit = 4.848e-9 so 1 code length unit ~ 1 AU.
- ChaNGa reads third_party/ChaNGa/single_particle.std as input.
- Plot script prints both Quinn requested IC and ChaNGa first snapshot IC.
