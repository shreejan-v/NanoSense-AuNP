# NanoSense-AuNP

Computational design of a **gold-nanoparticle (AuNP) colorimetric / LSPR biosensor** for low-cost point-of-care diagnostics.
Mini-project for 21NT0404T – Nanotechnology in Societal Development, SRMIST.

## What it simulates
| Experiment | Model | Key result |
|---|---|---|
| Size vs plasmon peak | Full Mie theory, Etchegoin (2006) gold permittivity | 10 nm → 526.5 nm, 100 nm → 573 nm |
| Biomolecule binding | Quasi-static core–shell (protein shell, n = 1.45) | 10 nm shell → +7.5 nm red-shift |
| Refractive-index sensitivity | Mie, n = 1.333–1.41 | ≈ 77.5 nm/RIU (20 nm AuNP) |
| Aggregation (colour change) | Mie on Maxwell-Garnett effective clusters | peak 528 → 663 nm, A620/A520 0.04 → 1.67 |

## Run
```bash
pip install numpy scipy matplotlib
python nanosense.py
```
Outputs: `figures/fig1–fig4*.png`, `results/results.csv`.

## Limitations
Simulation only (no wet-lab synthesis). The quasi-static shell model ignores retardation, and the Maxwell-Garnett clusters approximate real fractal aggregates, so trends are reliable but absolute values are approximate.
