"""Extended analysis for NanoSense: dielectric model check, absorption vs scattering,
size-dependent sensitivity, and gold usage per test.  Run after nanosense.py."""
import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nanosense import eps_gold, mie_qext, WL, N_WATER


def mie_abs_sca(d, wl, n_med=N_WATER):
    qe, qs = np.zeros_like(wl), np.zeros_like(wl)
    for k, lam in enumerate(wl):
        m = np.sqrt(eps_gold(lam)) / n_med
        x = np.pi * d * n_med / lam
        n = np.arange(1, int(x + 4 * x ** (1 / 3) + 2) + 1)
        mx = m * x
        jx, jxp = spherical_jn(n, x), spherical_jn(n, x, derivative=True)
        yx, yxp = spherical_yn(n, x), spherical_yn(n, x, derivative=True)
        jm, jmp = spherical_jn(n, mx), spherical_jn(n, mx, derivative=True)
        px, dpx = x * jx, jx + x * jxp
        hx, dhx = x * (jx + 1j * yx), (jx + 1j * yx) + x * (jxp + 1j * yxp)
        pm, dpm = mx * jm, jm + mx * jmp
        a = (m * pm * dpx - px * dpm) / (m * pm * dhx - hx * dpm)
        b = (pm * dpx - m * px * dpm) / (pm * dhx - m * hx * dpm)
        qe[k] = 2 / x**2 * np.sum((2 * n + 1) * (a + b).real)
        qs[k] = 2 / x**2 * np.sum((2 * n + 1) * (abs(a) ** 2 + abs(b) ** 2))
    return qe - qs, qs


# Fig 5: dielectric function
e = eps_gold(WL)
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(WL, e.real, label="Re(ε)"); ax.plot(WL, e.imag, label="Im(ε)")
ax.axhline(-2 * N_WATER**2, ls="--", c="gray", lw=0.8, label="Fröhlich condition (−2ε$_m$)")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Permittivity"); ax.legend(); ax.set_title("Gold dielectric function (Etchegoin 2006)")
fig.tight_layout(); fig.savefig("figures/fig5_dielectric.png", dpi=300)
froh = WL[np.argmin(abs(e.real + 2 * N_WATER**2))]
print(f"Frohlich wavelength (Re eps = -2 eps_m): {froh:.1f} nm")

# Fig 6: absorption vs scattering
fig, axs = plt.subplots(1, 3, figsize=(7.7, 2.52))
for ax, d in zip(axs, [20, 50, 80]):
    qa, qs = mie_abs_sca(d, WL)
    ax.plot(WL, qa, label="absorption"); ax.plot(WL, qs, label="scattering"); ax.set_title(f"{d} nm")
    ax.set_xlabel("Wavelength (nm)")
    i = np.argmax(qa + qs)
    print(f"d={d}: scattering fraction at peak = {qs[i] / (qa[i] + qs[i]):.3f}")
axs[0].set_ylabel("Efficiency"); axs[0].legend(); fig.tight_layout(); fig.savefig("figures/fig6_abs_sca.png", dpi=300)

# Fig 7: RI sensitivity vs size
ns = np.array([1.333, 1.35, 1.37, 1.39, 1.41]); sizes = [20, 40, 60, 80]; S = []
for d in sizes:
    pk = [WL[np.argmax(mie_qext(d, WL, n))] for n in ns]
    S.append(np.polyfit(ns, pk, 1)[0]); print(f"d={d}: sensitivity {S[-1]:.1f} nm/RIU")
fig, ax = plt.subplots(figsize=(3.75, 3)); ax.plot(sizes, S, "o-")
ax.set_xlabel("AuNP diameter (nm)"); ax.set_ylabel("Sensitivity (nm/RIU)"); ax.set_title("Refractive-index sensitivity vs size")
fig.tight_layout(); fig.savefig("figures/fig7_sensitivity_vs_size.png", dpi=300)

# Gold used per test (100 uL of 20 nm AuNP at ~7.2e11 particles/mL, typical OD-1 colloid)
vol = 4 / 3 * np.pi * (10e-7) ** 3                 # cm^3 per particle
mass = 19.3 * vol * 7.2e11 * 0.1                    # g in 100 uL
print(f"Gold per 100 uL test: {mass * 1e6:.2f} ug")
