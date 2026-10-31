# NanoSense-AuNP

**Computational design of a gold-nanoparticle (AuNP) plasmonic colorimetric biosensor for low-cost point-of-care diagnostics**

Mini-project for **21NT0404T – Nanotechnology in Societal Development**, SRM Institute of Science and Technology, Kattankulathur
Author: **Shreejan V** (RA2311070010025), B.Tech Biotechnology · Guide: **Dr. Karthika P C**, Department of Physics and Nanotechnology

![Python](https://img.shields.io/badge/Python-3.x-blue) ![NumPy](https://img.shields.io/badge/NumPy-SciPy-orange) ![Type](https://img.shields.io/badge/Study-Simulation-green)

---

## Table of contents
1. [Background](#background)
2. [Problem statement](#problem-statement)
3. [How the sensor works](#how-the-sensor-works)
4. [Physics and models](#physics-and-models)
5. [Repository structure](#repository-structure)
6. [Installation and usage](#installation-and-usage)
7. [Results](#results)
8. [Proposed sensor design](#proposed-sensor-design)
9. [Societal impact, safety and ethics](#societal-impact-safety-and-ethics)
10. [Limitations](#limitations)
11. [Future scope](#future-scope)
12. [References](#references)

---

## Background
Bulk gold is yellow, but gold nanoparticles of 10–100 nm in water look **ruby red**. Light drives the particle's conduction electrons into a collective oscillation called the **localized surface plasmon resonance (LSPR)**, which strongly absorbs green light around 520–530 nm. The position of this resonance depends on:

- **particle size and shape**
- **the refractive index right at the particle surface**, which changes when biomolecules bind
- **the distance between particles**: when particles aggregate, the plasmons couple and the colour turns **purple-blue**

This is the principle behind the gold-nanoparticle line in home pregnancy tests and COVID-19 rapid antigen tests, among the most widely used nano-enabled products.

## Problem statement
Laboratory diagnostics such as ELISA, PCR and culture are accurate but need trained staff, instruments, electricity and hours to days. In rural and low-resource settings this delays treatment. A sensor that **changes colour visibly within minutes** would decentralise testing. This project uses physics-based simulation to find which particle size and which readout (binding shift or aggregation) give the strongest, most reliable signal, before any costly synthesis.

**SDG relevance:** SDG 3 (Good Health), SDG 9 (Innovation), SDG 10 (Reduced Inequalities).

## How the sensor works
```
 Antibody/aptamer-coated AuNPs (red, peak ~528 nm)
              │  + target analyte (antigen / DNA)
              ▼
 Analyte cross-links particles → clusters form
              │  plasmon coupling
              ▼
 Peak red-shifts to 600–660 nm → solution turns purple/blue
              │
              ▼
 Read by eye, or by smartphone camera (A620/A520 ratio)
```

## Physics and models
| Component | Model | Why |
|---|---|---|
| Optical response of gold | Analytical Drude + two critical-point dielectric function (Etchegoin et al., 2006) | Reproduces gold's interband transitions in the visible; no lookup tables needed |
| Single-particle spectra | Full **Mie theory** (Bohren & Huffman), Riccati–Bessel functions via SciPy | Exact for spheres of any size; captures size red-shift and broadening |
| Biomolecule binding | Quasi-static **core–shell polarizability** (Au core + protein shell, n = 1.45) | Models the local refractive-index change from an adsorbed antibody layer |
| Bulk sensitivity | Mie theory with medium index n = 1.333–1.41 | Standard sensor figure of merit (nm/RIU) |
| Aggregation | **Maxwell-Garnett** effective medium (Au fill fraction 0.4) + Mie for the cluster | Approximates plasmon coupling and retardation in growing clusters |

Key equation (Mie extinction efficiency):

$$Q_{ext} = \frac{2}{x^2}\sum_{n=1}^{\infty}(2n+1)\,\mathrm{Re}(a_n+b_n), \qquad x = \frac{\pi d\, n_m}{\lambda}$$

Aggregation readout: **A620 / A520**, the absorbance ratio used in experimental AuNP assays.

## Repository structure
```
NanoSense-AuNP/
├── nanosense.py          # all models + the four virtual experiments
├── extended_analysis.py  # dielectric check, absorption vs scattering, sensitivity vs size, gold per test
├── figures/
│   ├── fig1_size_spectra.png     # Mie spectra, 10–100 nm
│   ├── fig2_binding_shift.png    # LSPR shift vs protein-shell thickness
│   ├── fig3_ri_sensitivity.png   # peak vs refractive index
│   └── fig4_aggregation.png      # aggregation-induced red-shift
├── results/
│   └── results.csv       # every number reported below
└── report/
    └── report.html       # full project report (SRMIST format)
```

## Installation and usage
```bash
git clone https://github.com/shreejan-v/NanoSense-AuNP.git
cd NanoSense-AuNP
pip install numpy scipy matplotlib
python nanosense.py
```
It runs in under a minute on a laptop, regenerates all four figures in `figures/` and writes `results/results.csv`.

To try your own particle, call the functions directly:
```python
from nanosense import mie_qext, WL
q = mie_qext(50, WL)            # 50 nm gold sphere in water
print(WL[q.argmax()])           # plasmon peak (nm)
```

## Results

### 1. Size dependence
![size](figures/fig1_size_spectra.png)

| Diameter (nm) | 10 | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|---|
| Peak (nm) | 526.5 | 528.0 | 532.5 | 542.0 | 555.5 | 573.0 |
| Max Q<sub>ext</sub> | 0.74 | 1.53 | 3.48 | 5.83 | 7.60 | 7.81 |

The peak barely moves up to 40 nm, then red-shifts and broadens. **20–40 nm particles** give a sharp band with strong extinction, the best choice for colorimetric sensing.

### 2. Biomolecule binding
![binding](figures/fig2_binding_shift.png)

| Protein shell (nm) | 0 | 1 | 2 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|
| Red-shift (nm) | 0 | 2.0 | 3.5 | 5.0 | 6.5 | 7.0 | 7.5 |

The shift saturates because the plasmon near-field decays within ~10 nm. A spectrometer can detect it, but the eye cannot.

### 3. Refractive-index sensitivity
![ri](figures/fig3_ri_sensitivity.png)

The peak moves from 528.0 to 534.0 nm for n = 1.333 → 1.41: **≈ 77.5 nm/RIU** for a 20 nm sphere, within the range reported for small spheres.

### 4. Aggregation (the visible colour change)
![agg](figures/fig4_aggregation.png)

| Cluster size (nm) | 20 (single) | 60 | 100 | 150 | 200 |
|---|---|---|---|---|---|
| ≈ particles | 1 | 11 | 50 | 169 | 400 |
| Peak (nm) | 528.0 | 532.5 | 543.0 | 616.0 | 663.0 |
| A620/A520 | 0.04 | 0.07 | 0.26 | 1.67 | 1.39 |

Aggregation gives a **>130 nm shift**, more than 15× the binding shift, and a **~40-fold rise in A620/A520**. That makes it the preferred mode for a by-eye or smartphone test.

### Summary
| Metric | Value |
|---|---|
| Optimal particle size | 20–40 nm |
| Max binding shift (10 nm shell) | 7.5 nm |
| RI sensitivity | ≈ 77.5 nm/RIU |
| Aggregation shift | 528 → 663 nm |
| A620/A520 change | 0.04 → 1.67 |

## Proposed sensor design
- **Probe:** 20–40 nm citrate-capped AuNPs (Turkevich synthesis) functionalised with antibodies or aptamers.
- **Format:** single-tube assay or paper strip.
- **Signal:** red → purple/blue within minutes when the target is present.
- **Readout:** visual, or smartphone RGB analysis calibrated against A620/A520.
- **Validation needed:** UV-Vis, TEM, DLS, zeta potential, and tests in real samples.

## Societal impact, safety and ethics
- **Benefits:** fast, cheap (micrograms of gold per test), no electricity needed for a visual read, usable by community health workers.
- **Nanosafety:** AuNPs have low acute toxicity but are biopersistent. Here they stay inside the device, so exposure is low (*Risk = Hazard × Exposure*). Used tests go to biomedical waste.
- **Ethics:** report sensitivity and false results honestly, with no "nano-washing". Smartphone data falls under the DPDP Act 2023. Affordable, open designs help close the nano-divide.
- **Regulation (India):** in-vitro diagnostics are regulated by CDSCO under the Medical Devices Rules 2017. Characterise following ISO/TC 229 and OECD guidance, with safe-by-design choices.

## Limitations
- **Simulation only:** no wet-lab synthesis.
- The quasi-static shell model ignores retardation.
- Maxwell-Garnett clusters approximate real fractal aggregates.
- Real samples cause non-specific aggregation, which is not modelled.

The trends are reliable, but the absolute numbers are approximate.

## Future scope
- Wet-lab synthesis and experimental validation
- Nanorods and nanostars for higher sensitivity
- DDA/FDTD simulation of realistic aggregates
- ML-based smartphone colour-reader app
- Multiplexed strips for several diseases

## References
1. G. Mie, *Annalen der Physik* 330, 377 (1908).
2. C. F. Bohren & D. R. Huffman, *Absorption and Scattering of Light by Small Particles*, Wiley (1983).
3. P. B. Johnson & R. W. Christy, *Phys. Rev. B* 6, 4370 (1972).
4. P. G. Etchegoin, E. C. Le Ru & M. Meyer, *J. Chem. Phys.* 125, 164705 (2006).
5. C. A. Mirkin et al., *Nature* 382, 607 (1996).
6. R. Elghanian et al., *Science* 277, 1078 (1997).
7. S. Link & M. A. El-Sayed, *J. Phys. Chem. B* 103, 4212 (1999).
8. W. Haiss et al., *Anal. Chem.* 79, 4215 (2007).
9. K. A. Willets & R. P. Van Duyne, *Annu. Rev. Phys. Chem.* 58, 267 (2007).
10. P. K. Jain et al., *J. Phys. Chem. B* 110, 7238 (2006).
11. G. A. Posthuma-Trumpie et al., *Anal. Bioanal. Chem.* 393, 569 (2009).
12. J. Turkevich et al., *Discuss. Faraday Soc.* 11, 55 (1951).

---
*Academic mini-project. All results in this README come from running `nanosense.py`.*
