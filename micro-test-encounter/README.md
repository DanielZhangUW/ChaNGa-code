# Micro Test Encounter

2-particle Hill encounter test (isolated from legacy single-particle workflows).

## 1) Generate encounter IC
```bash
python3 micro-test-encounter/make_encounter_tipsy.py \
  --mode y-far \
  --b-factor 1.0 \
  --y-span 10 \
  --T-run 6.283185 \
  --m1 1e-3 --m2 1e-10 \
  --dOrbdist 1.0 --dCentMass 1.0 \
  --q 1.5 --Omega 1.0 \
  --vx2 0.0 --vy2-mode shear \
  --out third_party/ChaNGa/single_particle.std
```

Optional guard (if ChaNGa crashes with tiny-N DD):
```bash
python3 micro-test-encounter/make_encounter_tipsy.py --mode y-far --add-guard --out third_party/ChaNGa/single_particle.std
```

## 2) Run ChaNGa
```bash
cd third_party/ChaNGa
./ChaNGa ../../micro-test-encounter/encounter_hill.param +p1
cd ../..
```

## 3) Plot diagnostics
```bash
MPLBACKEND=Agg python3 micro-test-encounter/plot_time_series.py --no-show
```

Outputs are written to `micro-test-encounter/`:
- `encounter_timeseries.png`
- `encounter_xy.png`
- `encounter_x_scaled.png`
- `encounter_y_scaled.png`
