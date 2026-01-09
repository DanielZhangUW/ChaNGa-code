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
  ChaNGa run card (sliding patch, 10 steps per orbit).
- tools/make_single_particle_tipsy.py
  Generates third_party/ChaNGa/single_particle.std from a chosen IC.
- task2/plot_task2.py
  Quinn/Hill reference integrator (Python).
- task2/plot_time_series.py
  Plots x(t), y(t), vx(t), vy(t) for ChaNGa vs Quinn.
- task2/run_origin_or_epsilon_test.py
  Runs the origin test with epsilon fallback if ChaNGa crashes at x0=0.

Quick Start (Epicycle IC)
-------------------------
From the repo root:

  python3 tools/make_single_particle_tipsy.py --x0 1 --y0 0 --vx0 0 --vy0 -2 --add-guard
  cd third_party/ChaNGa
  ./ChaNGa single_particle_hill.param +p1
  cd ..
  python3 task2/plot_time_series.py --x0 1 --y0 0 --vx0 0 --vy0 -2 \
    --output task2/single_particle_timeseries.png

Origin Test with Epsilon Fallback
---------------------------------
ChaNGa may segfault for (x0,y0,vx0,vy0)=(0,0,0,0). Use the epsilon driver:

  python3 task2/run_origin_or_epsilon_test.py

This retries x0 = 0, 1e-12, 1e-11, ... up to 1e-6, and saves a figure like:
  task2/single_particle_timeseries_eps1e-12.png

Notes
-----
- Units: dKpcUnit = 4.848e-9 so 1 code length unit ~ 1 AU.
- ChaNGa reads third_party/ChaNGa/single_particle.std as input.
- Plot script prints both Quinn requested IC and ChaNGa first snapshot IC.
