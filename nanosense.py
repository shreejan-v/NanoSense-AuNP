"""NanoSense: simulation of a gold-nanoparticle (AuNP) colorimetric / LSPR biosensor.

Physics used
- Gold permittivity: analytical Drude + two critical-point model (Etchegoin et al., J. Chem. Phys. 125, 164705, 2006).
- Size-dependent extinction of a single sphere: full Mie theory (Bohren & Huffman, 1983).
- Biomolecule binding: quasi-static coated-sphere polarizability (core = Au, shell = protein layer, n = 1.45).
- Aggregation: Maxwell-Garnett effective medium for AuNP clusters of increasing gold volume fraction.

Run:  python nanosense.py      -> writes figures/ and results/results.csv
"""
import os, csv
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WL = np.linspace(400, 800, 801)          # wavelength, nm
N_WATER = 1.333
os.makedirs("figures", exist_ok=True); os.makedirs("results", exist_ok=True)


def eps_gold(wl):
    """Etchegoin 2006 analytical dielectric function of gold (wl in nm)."""
    e_inf, lp, gp = 1.53, 145.0, 17000.0
    eps = e_inf - 1.0 / (lp**2 * (1.0 / wl**2 + 1j / (gp * wl)))
    for A, phi, li, gi in [(0.94, -np.pi / 4, 468.0, 2300.0), (1.36, -np.pi / 4, 331.0, 940.0)]:
        eps += A / li * (np.exp(1j * phi) / (1 / li - 1 / wl - 1j / gi) +
                         np.exp(-1j * phi) / (1 / li + 1 / wl + 1j / gi))
    return eps


def mie_qext(d, wl, n_med=N_WATER, eps_fn=None):
    """Extinction efficiency of a sphere (diameter d, nm) by Mie theory; default material = gold."""
    eps_fn = eps_fn or eps_gold
    q = np.zeros_like(wl)
    for k, lam in enumerate(wl):
        m = np.sqrt(eps_fn(lam)) / n_med
        x = np.pi * d * n_med / lam
        nmax = int(x + 4 * x ** (1 / 3) + 2)
        n = np.arange(1, nmax + 1)
        mx = m * x
        # Riccati-Bessel functions and derivatives
        jx, jxp = spherical_jn(n, x), spherical_jn(n, x, derivative=True)
        yx, yxp = spherical_yn(n, x), spherical_yn(n, x, derivative=True)
        jm, jmp = spherical_jn(n, mx), spherical_jn(n, mx, derivative=True)
        psi_x, dpsi_x = x * jx, jx + x * jxp
        xi_x, dxi_x = x * (jx + 1j * yx), (jx + 1j * yx) + x * (jxp + 1j * yxp)
        psi_m, dpsi_m = mx * jm, jm + mx * jmp
        a = (m * psi_m * dpsi_x - psi_x * dpsi_m) / (m * psi_m * dxi_x - xi_x * dpsi_m)
        b = (psi_m * dpsi_x - m * psi_x * dpsi_m) / (psi_m * dxi_x - m * xi_x * dpsi_m)
        q[k] = 2 / x**2 * np.sum((2 * n + 1) * (a + b).real)
    return q


def coated_cext(d, t, wl, n_shell=1.45, n_med=N_WATER):
    """Quasi-static extinction cross-section (nm^2) of Au core + protein shell of thickness t."""
    e1, e2, em = eps_gold(wl), n_shell**2, n_med**2
    r1, r2 = d / 2, d / 2 + t
    f = (r1 / r2) ** 3
    alpha = 4 * np.pi * r2**3 * ((e2 - em) * (e1 + 2 * e2) + f * (e1 - e2) * (em + 2 * e2)) / \
            ((e2 + 2 * em) * (e1 + 2 * e2) + f * (2 * e2 - 2 * em) * (e1 - e2))
    k = 2 * np.pi * n_med / wl
    return k * alpha.imag


def eps_mg(fill, n_med=N_WATER):
    """Maxwell-Garnett effective permittivity of an AuNP cluster at gold volume fraction `fill`."""
    def f(wl):
        e1, em = eps_gold(wl), n_med**2
        beta = (e1 - em) / (e1 + 2 * em)
        return em * (1 + 2 * fill * beta) / (1 - fill * beta)
    return f


def peak(y):
    return WL[np.argmax(y)]


def ratio(y):                                            # standard aggregation readout
    return np.interp(620, WL, y) / np.interp(520, WL, y)


rows = []
# 1. size dependence
plt.figure(figsize=(4.5, 3))
for d in [10, 20, 40, 60, 80, 100]:
    qe = mie_qext(d, WL)
    rows.append(["size", d, round(peak(qe), 1), round(qe.max(), 3), ""])
    plt.plot(WL, qe / qe.max(), label=f"{d} nm")
plt.xlabel("Wavelength (nm)"); plt.ylabel("Normalised extinction"); plt.title("Mie extinction of AuNPs in water")
plt.legend(); plt.tight_layout(); plt.savefig("figures/fig1_size_spectra.png", dpi=300); plt.close()

# 2. binding sensing (protein shell growth)
plt.figure(figsize=(4.5, 3)); shifts = []
base = peak(coated_cext(20, 0, WL))
for t in [0, 1, 2, 4, 6, 8, 10]:
    c = coated_cext(20, t, WL)
    shifts.append(peak(c) - base)
    rows.append(["binding", t, round(peak(c), 1), "", round(peak(c) - base, 2)])
    plt.plot(WL, c, label=f"shell {t} nm")
plt.xlim(450, 650); plt.xlabel("Wavelength (nm)"); plt.ylabel("C_ext (nm$^2$)")
plt.title("20 nm AuNP: LSPR shift on biomolecule binding"); plt.legend(fontsize=8); plt.tight_layout()
plt.savefig("figures/fig2_binding_shift.png", dpi=300); plt.close()

# 3. bulk refractive-index sensitivity
ns = np.array([1.333, 1.35, 1.37, 1.39, 1.41])
pk = np.array([peak(mie_qext(20, WL, n)) for n in ns])
S = np.polyfit(ns, pk, 1)[0]
for n, p in zip(ns, pk):
    rows.append(["RI", n, round(p, 1), "", ""])
plt.figure(figsize=(3.75, 3)); plt.plot(ns, pk, "o-")
plt.xlabel("Medium refractive index"); plt.ylabel("LSPR peak (nm)")
plt.title(f"Sensitivity = {S:.0f} nm/RIU"); plt.tight_layout(); plt.savefig("figures/fig3_ri_sensitivity.png", dpi=300); plt.close()

# 4. aggregation colorimetric readout: clusters of 20 nm AuNPs (fill 0.4) growing in size, Mie of effective sphere
plt.figure(figsize=(4.5, 3))
for D in [20, 60, 100, 150, 200]:
    c = mie_qext(D, WL) if D == 20 else mie_qext(D, WL, eps_fn=eps_mg(0.4))
    nparts = round(0.4 * (D / 20) ** 3) if D > 20 else 1
    rows.append(["aggregation", D, round(peak(c), 1), nparts, round(ratio(c), 3)])
    plt.plot(WL, c / c.max(), label=f"cluster {D} nm (~{nparts} NPs)")
plt.xlabel("Wavelength (nm)"); plt.ylabel("Normalised extinction")
plt.title("Aggregation-induced red-shift (red -> blue/purple)", fontsize=10); plt.ylim(0, 1.45); plt.legend(fontsize=7, ncol=2, loc="upper center"); plt.tight_layout()
plt.savefig("figures/fig4_aggregation.png", dpi=300); plt.close()

with open("results/results.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["experiment", "parameter", "peak_nm", "max_Qext_or_NPs", "shift_or_A620/A520"]); w.writerows(rows)
print(f"RI sensitivity: {S:.1f} nm/RIU")
for r in rows: print(r)
