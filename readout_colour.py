"""NanoSense readout analysis: quasi-static vs Mie check, aggregation readout sweep,
simulated solution colour (CIE 1931 -> sRGB / CIELAB) and a simulated dose-response.
Run after nanosense.py (it imports the models from there).

Colour method
- CIE 1931 2-degree colour-matching functions: multi-lobe analytic fit of
  Wyman, Sloan & Shirley, J. Comput. Graph. Tech. 2(2), 1-11 (2013).
- Illuminant: 6504 K blackbody, white point scaled to D65 (XYZ scaling).
- XYZ -> linear sRGB (IEC 61966-2-1 matrix) -> gamma-encoded sRGB; CIELAB relative to D65.
- Colloid: 1 cm path, 20 nm dispersed AuNPs at the concentration giving peak absorbance 1.0.
  All other samples hold the same total gold content (gold is conserved on aggregation).
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from nanosense import eps_gold, mie_qext, eps_mg, WL, N_WATER

rows = []


def peak(wl, y):
    return wl[np.argmax(y)]


def r620(wl, y):
    return np.interp(620, wl, y) / np.interp(520, wl, y)


# ---------------------------------------------------------------- 1. quasi-static vs Mie
def qs_qext(d, wl, n_med=N_WATER):
    """Dipole (quasi-static) extinction efficiency: absorption + dipole scattering."""
    e, em = eps_gold(wl), n_med**2
    r = d / 2
    alpha = 4 * np.pi * r**3 * (e - em) / (e + 2 * em)
    k = 2 * np.pi * n_med / wl
    return (k * alpha.imag + k**4 * np.abs(alpha) ** 2 / (6 * np.pi)) / (np.pi * r**2)


sizes = [10, 20, 30, 40, 50, 60, 80, 100]
pk_qs, pk_mie = [], []
for d in sizes:
    a, b = qs_qext(d, WL), mie_qext(d, WL)
    pk_qs.append(peak(WL, a)); pk_mie.append(peak(WL, b))
    rows.append(["qs_vs_mie", d, pk_qs[-1], pk_mie[-1], round(pk_mie[-1] - pk_qs[-1], 1)])
    print(f"d={d}: QS peak {pk_qs[-1]}, Mie peak {pk_mie[-1]}, diff {pk_mie[-1]-pk_qs[-1]:.1f}")
fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(sizes, pk_mie, "o-", label="Full Mie theory")
ax.plot(sizes, pk_qs, "s--", label="Quasi-static (dipole)")
ax.set_xlabel("AuNP diameter (nm)"); ax.set_ylabel("LSPR peak (nm)")
ax.set_title("Quasi-static approximation vs Mie theory"); ax.legend(); fig.tight_layout()
fig.savefig("figures/fig8_quasistatic_vs_mie.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- 2. aggregation readout sweep
Ds = np.arange(40, 255, 5)
fills = [0.2, 0.3, 0.4, 0.5]
single = mie_qext(20, WL)
r_single, p_single = r620(WL, single), peak(WL, single)
sweep = {}
for f in fills:
    R, P = [], []
    for D in Ds:
        q = mie_qext(D, WL, eps_fn=eps_mg(f))
        R.append(r620(WL, q)); P.append(peak(WL, q))
    R, P = np.array(R), np.array(P)
    sweep[f] = (R, P)
    i_max = int(np.argmax(R))
    cross = Ds[np.argmax(R >= 1)] if (R >= 1).any() else None
    print(f"fill {f}: max A620/A520 {R[i_max]:.2f} at D={Ds[i_max]} nm (peak {P[i_max]} nm); "
          f"first D with ratio>=1: {cross}; peak at D=200: {P[Ds == 200][0]}")
    rows.append(["agg_sweep_summary", f, round(R[i_max], 3), int(Ds[i_max]), cross])
    for D, r_, p_ in zip(Ds, R, P):
        rows.append(["agg_sweep", f"{f}/{D}", p_, round(f * (D / 20) ** 3, 1), round(r_, 3)])
fig, axs = plt.subplots(1, 2, figsize=(7.15, 2.535))
for f in fills:
    axs[0].plot(Ds, sweep[f][0], label=f"fill {f}")
    P = sweep[f][1].astype(float).copy(); P[P >= 799.5] = np.nan       # peak beyond the 800 nm window
    axs[1].plot(Ds, P, ".", ms=5, label=f"fill {f}")
axs[0].axhline(r_single, ls=":", c="gray", lw=1, label="single 20 nm AuNP")
axs[0].axhline(1.0, ls="--", c="gray", lw=0.8)
axs[1].axhline(p_single, ls=":", c="gray", lw=1, label="single 20 nm AuNP")
axs[0].set_xlabel("Cluster diameter (nm)"); axs[0].set_ylabel("A620 / A520"); axs[0].set_title("(a) Colorimetric ratio")
axs[1].set_xlabel("Cluster diameter (nm)"); axs[1].set_ylabel("Extinction peak (nm)"); axs[1].set_title("(b) Peak position")
axs[0].legend(fontsize=7); axs[1].legend(fontsize=7); fig.tight_layout()
fig.savefig("figures/fig9_aggregation_readout.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- 3. colour of the colloid
WLc = np.arange(380.0, 781.0, 1.0)


def g(x, mu, s1, s2):
    s = np.where(x < mu, s1, s2)
    return np.exp(-0.5 * ((x - mu) / s) ** 2)


xbar = 1.056 * g(WLc, 599.8, 37.9, 31.0) + 0.362 * g(WLc, 442.0, 16.0, 26.7) - 0.065 * g(WLc, 501.1, 20.4, 26.2)
ybar = 0.821 * g(WLc, 568.8, 46.9, 40.5) + 0.286 * g(WLc, 530.9, 16.3, 31.1)
zbar = 1.217 * g(WLc, 437.0, 11.8, 36.0) + 0.681 * g(WLc, 459.0, 26.0, 13.8)
h, c, kB = 6.626e-34, 2.998e8, 1.381e-23
lam = WLc * 1e-9
S = 1 / (lam**5 * (np.exp(h * c / (lam * kB * 6504.0)) - 1))      # blackbody, arbitrary units
XYZ_w = np.array([np.sum(S * xbar), np.sum(S * ybar), np.sum(S * zbar)])
D65 = np.array([0.95047, 1.0, 1.08883])
M = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])


def colour(T):
    xyz = np.array([np.sum(S * T * xbar), np.sum(S * T * ybar), np.sum(S * T * zbar)]) / XYZ_w * D65
    rgb = M @ xyz
    rgb = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * np.clip(rgb, 0, None) ** (1 / 2.4) - 0.055)
    rgb = np.clip(rgb, 0, 1)
    t = xyz / D65
    ft = np.where(t > (6 / 29) ** 3, np.cbrt(t), t / (3 * (6 / 29) ** 2) + 4 / 29)
    L, a, b = 116 * ft[1] - 16, 500 * (ft[0] - ft[1]), 200 * (ft[1] - ft[2])
    return rgb, np.array([L, a, b])


def hexrgb(rgb):
    return "#" + "".join(f"{int(round(v * 255)):02X}" for v in rgb)


C20 = mie_qext(20, WLc) * np.pi * 10**2 * 1e-14                   # cm^2
N0 = np.log(10) / C20.max()                                        # particles per mL for peak A = 1 (1 cm)
print(f"20 nm AuNP number density for peak absorbance 1.0: {N0:.3e} per mL")
A_disp = N0 * C20 / np.log(10)


def absorbance_sphere(d):
    N = N0 * (20 / d) ** 3                                        # same total gold
    return N * mie_qext(d, WLc) * np.pi * (d / 2) ** 2 * 1e-14 / np.log(10)


def absorbance_cluster(D, f=0.4):
    n_c = f * (D / 20) ** 3
    N = N0 / n_c
    return N * mie_qext(D, WLc, eps_fn=eps_mg(f)) * np.pi * (D / 2) ** 2 * 1e-14 / np.log(10)


samples = [("size", f"{d} nm", absorbance_sphere(d)) for d in [10, 20, 40, 60, 80, 100]]
samples += [("aggregation", "20 nm single", A_disp)] + \
           [("aggregation", f"{D} nm cluster", absorbance_cluster(D)) for D in [60, 100, 150, 200]]
lab_ref = colour(10 ** (-A_disp))[1]
swatch = []
for grp, name, A in samples:
    rgb, lab = colour(10 ** (-A))
    hue = np.degrees(np.arctan2(lab[2], lab[1])) % 360
    dE = np.linalg.norm(lab - lab_ref)
    swatch.append((grp, name, rgb, lab, hue))
    print(f"{grp:12s} {name:16s} peakA={A.max():.2f} at {peak(WLc, A):.0f} nm  sRGB {hexrgb(rgb)}  "
          f"L*={lab[0]:.1f} a*={lab[1]:.1f} b*={lab[2]:.1f} hue={hue:.0f} dE76_vs_dispersed={dE:.1f}")
    rows.append(["colour", f"{grp}:{name}", round(peak(WLc, A), 1), hexrgb(rgb),
                 f"L={lab[0]:.1f};a={lab[1]:.1f};b={lab[2]:.1f};h={hue:.0f};dE={dE:.1f}"])

fig, ax = plt.subplots(figsize=(9, 3.6))
for row, grp, title in [(1, "size", "Dispersed AuNPs of different size (same gold content)"),
                        (0, "aggregation", "Aggregation of 20 nm AuNPs (fill fraction 0.4, same gold content)")]:
    items = [s for s in swatch if s[0] == grp]
    ax.text(0, row * 1.6 + 1.22, title, fontsize=10, fontweight="bold")
    for i, (_, name, rgb, lab, hue) in enumerate(items):
        ax.add_patch(Rectangle((i * 1.45, row * 1.6 + 0.3), 1.25, 0.82, color=rgb, ec="k", lw=0.6))
        ax.text(i * 1.45 + 0.625, row * 1.6 + 0.1, name, ha="center", va="center", fontsize=8)
ax.set_xlim(-0.1, 8.8); ax.set_ylim(-0.1, 3.2); ax.axis("off"); fig.tight_layout()
fig.savefig("figures/fig10_colour_swatches.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- 4. simulated dose-response
A_agg = absorbance_cluster(150)
phis = np.linspace(0, 1, 101)
R_phi, dE_phi, rgb_phi = [], [], []
for p in phis:
    A = (1 - p) * A_disp + p * A_agg
    rgb, lab = colour(10 ** (-A))
    R_phi.append(r620(WLc, A)); dE_phi.append(np.linalg.norm(lab - lab_ref)); rgb_phi.append(rgb)
R_phi, dE_phi = np.array(R_phi), np.array(dE_phi)
phi_jnd = np.interp(2.3, dE_phi, phis)
phi_10 = np.interp(10, dE_phi, phis)
print(f"dose-response: A620/A520 {R_phi[0]:.3f} -> {R_phi[-1]:.3f}; dE at phi=1 {dE_phi[-1]:.1f}")
print(f"phi at dE=2.3 (JND): {phi_jnd:.3f}, ratio there {np.interp(phi_jnd, phis, R_phi):.3f}")
print(f"phi at dE=10: {phi_10:.3f}, ratio there {np.interp(phi_10, phis, R_phi):.3f}")
for p in [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]:
    i = int(round(p * 100))
    print(f"  phi={p:.2f}: ratio {R_phi[i]:.3f}, dE {dE_phi[i]:.1f}, sRGB {hexrgb(rgb_phi[i])}")
    rows.append(["dose_response", p, round(R_phi[i], 3), hexrgb(rgb_phi[i]), round(dE_phi[i], 1)])
rows.append(["dose_response_threshold", "dE=2.3", round(phi_jnd, 3), "", round(np.interp(phi_jnd, phis, R_phi), 3)])
rows.append(["dose_response_threshold", "dE=10", round(phi_10, 3), "", round(np.interp(phi_10, phis, R_phi), 3)])

fig = plt.figure(figsize=(6.5, 3.38))
ax1 = fig.add_axes([0.11, 0.36, 0.36, 0.56]); ax2 = fig.add_axes([0.61, 0.36, 0.36, 0.56])
ax3 = fig.add_axes([0.07, 0.05, 0.90, 0.16])
ax1.plot(phis * 100, R_phi); ax1.set_xlabel("Gold in aggregates (%)"); ax1.set_ylabel("A620 / A520")
ax1.set_title("(a) Ratio readout")
ax2.plot(phis * 100, dE_phi)
ax2.axhline(2.3, ls="--", c="gray", lw=0.8); ax2.axhline(10, ls=":", c="gray", lw=0.8)
ax2.text(100, 2.3, "ΔE = 2.3 (just noticeable)", ha="right", va="bottom", fontsize=7)
ax2.text(100, 10, "ΔE = 10", ha="right", va="bottom", fontsize=7)
ax2.set_xlabel("Gold in aggregates (%)"); ax2.set_ylabel("Colour difference ΔE*ab"); ax2.set_title("(b) Visual colour change")
for i, p in enumerate(np.linspace(0, 1, 11)):
    ax3.add_patch(Rectangle((i, 0), 0.9, 1, color=rgb_phi[int(round(p * 100))], ec="k", lw=0.5))
    ax3.text(i + 0.45, -0.25, f"{int(p * 100)}%", ha="center", va="top", fontsize=7)
ax3.set_xlim(-0.05, 11); ax3.set_ylim(-0.7, 1.05); ax3.axis("off")
fig.savefig("figures/fig11_dose_response.png", dpi=300); plt.close(fig)

with open("results/readout_colour.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["analysis", "parameter", "value1", "value2", "value3"]); w.writerows(rows)
