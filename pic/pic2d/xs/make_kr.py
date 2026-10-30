"""Krypton cross sections for the PIC, fitted to HallThruster.jl's krypton rate tables so that the PIC and the 1-D
reference (../../hallthruster) share the same collision rates.

HallThruster.jl tabulates rate coefficients k against the mean electron energy eps = 3/2 Te. A cross section
sigma(E) reproduces a table if its Maxwellian average at Te = 2 eps / 3 does:

    k(Te) = int sigma(E) v(E) f(E; Te) dE,   f = 2 sqrt(E / pi) Te^-3/2 exp(-E / Te),   v = sqrt(2 e E / m_e)

- Ionization (Kr -> Kr+): threshold law sigma = a (x - 1)^p / x^q, x = E / E_th.
- Lumped excitation and elastic (momentum transfer, scattered isotropically by WarpX): a smooth non-negative
  piecewise-linear curve in log E (regularised non-negative least squares), zero below the excitation threshold.
  The threshold law misses the krypton excitation table by up to 37 % (Te 2-3 eV) and 10 % (Te 10-20 eV). The
  excitation penalty (lam = 10) is the weakest of 3e-3 ... 10 that leaves a single-peaked curve with no gaps; weaker
  ones fit the coarse table with spurious bumps. Its rate error is <= 4.2 % for Te 4-80 eV, 15 % at Te 2-3 eV.

Method check: the same fits on HallThruster.jl's xenon tables are compared with the measured xenon set that ships
with WarpX (warpx-data/MCC_cross_sections/Xe).

    python make_kr.py        -> Kr/{ionization,excitation,elastic}.dat, Kr/processes.json, kr_fit.png
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, nnls

HERE = Path(__file__).resolve().parent
HT = next((Path.home() / ".julia/packages/HallThruster").glob("*/reactions"))
WD = Path.home() / "tools/warpx-data/MCC_cross_sections"
Q_E, M_E = 1.602176634e-19, 9.1093837015e-31
E = np.logspace(-3, np.log10(5000.0), 6000)                    # integration grid [eV]
TE_FIT = (1.0, 100.0)                                           # Te range fitted [eV]
TE_REPORT = (2.0, 80.0)                                         # Te range reported (HET channel and plume)


def table(name):
    """HallThruster.jl table -> (threshold energy or None, Te [eV], k [m^3/s])."""
    lines = (HT / f"{name}.dat").read_text().splitlines()
    eth = float(lines[0].split(":")[1]) if ":" in lines[0] else None
    d = np.array([[float(v) for v in ln.split()] for ln in lines[1 + (eth is not None):] if ln.strip()])
    m = (d[:, 0] > 0) & (d[:, 1] > 0)
    return eth, 2.0 / 3.0 * d[m, 0], d[m, 1]


def kernel(te):
    """Rows: Maxwellian weights v f dE at each Te, on the grid E (trapezoid rule)."""
    w = np.gradient(E)
    f = 2.0 * np.sqrt(E[None, :] / np.pi) * te[:, None] ** -1.5 * np.exp(-E[None, :] / te[:, None])
    return f * np.sqrt(2.0 * Q_E * E / M_E)[None, :] * w[None, :]


def threshold_law(par, eth, e=E):
    a, p, q = np.exp(par[0]), par[1], par[2]
    x = np.maximum(e / eth, 1.0)
    return np.where(e > eth, a * (x - 1.0) ** p / x ** q, 0.0)


def fit_inelastic(eth, te, k):
    m = (te >= TE_FIT[0]) & (te <= TE_FIT[1])
    K = kernel(te[m])
    res = lambda par: np.log(K @ threshold_law(par, eth) / k[m])  # noqa: E731
    best = None
    for p0 in (0.5, 1.0, 1.5):
        for q0 in (1.0, 1.5, 2.0):
            s = least_squares(res, [np.log(3e-20), p0, q0 + p0], bounds=([-80, 0.05, 0.0], [-30, 4.0, 8.0]))
            if best is None or s.cost < best.cost:
                best = s
    return best.x


def fit_curve(te, k, eth=None, scale=1e-19, lam=3e-3):
    """Non-negative piecewise-linear sigma in log E, second differences penalised. With a threshold the first node
    sits at E_th with sigma = 0; without one the curve is flat below 0.01 eV. Flat above 1 keV."""
    m = (te >= TE_FIT[0]) & (te <= TE_FIT[1])
    nodes = np.logspace(np.log10(eth) if eth else -2, 3, 41)
    le, ln = np.log(np.maximum(E, nodes[0])), np.log(nodes)
    basis = np.array([np.interp(le, ln, np.eye(len(nodes))[i], right=(i == len(nodes) - 1))
                      for i in range(len(nodes))])                # hat functions in log E
    if eth:
        basis[:, E < eth] = 0.0
    A = kernel(te[m]) @ basis.T / k[m][:, None]
    D = np.diff(np.eye(len(nodes)), 2, axis=0)
    if eth:                                                     # sigma(E_th) = 0
        A, D = A[:, 1:], D[:, 1:]
    c, _ = nnls(np.vstack([A * scale, np.sqrt(lam) * D]), np.concatenate([np.ones(m.sum()), np.zeros(len(D))]))
    c = np.concatenate([[0.0], c]) * scale if eth else c * scale
    if eth:
        return lambda e: np.where(e >= eth, np.interp(np.log(np.maximum(e, eth)), ln, c), 0.0)
    return lambda e: np.interp(np.log(np.maximum(e, nodes[0])), ln, c)


def rate_error(sig, te, k):
    m = (te >= TE_REPORT[0]) & (te <= TE_REPORT[1])
    r = kernel(te[m]) @ sig(E) / k[m] - 1.0
    return float(np.max(np.abs(r))), float(np.sqrt(np.mean(r ** 2)))


def warpx_xs(path):
    d = np.loadtxt(path)
    return lambda e: np.interp(e, d[:, 0], d[:, 1], left=0.0, right=d[-1, 1])


def fit_gas(gas, ion_name):
    eth_i, te_i, k_i = table(f"ionization_{gas}_{ion_name}")
    eth_x, te_x, k_x = table(f"excitation_{gas}")
    _, te_e, k_e = table(f"elastic_{gas}")
    pi = fit_inelastic(eth_i, te_i, k_i)
    out = {"ionization": (eth_i, lambda e: threshold_law(pi, eth_i, e), te_i, k_i, pi),
           "excitation": (eth_x, fit_curve(te_x, k_x, eth=eth_x, scale=1e-20, lam=10.0), te_x, k_x, None),
           "elastic": (None, fit_curve(te_e, k_e), te_e, k_e, None)}
    return out


def write_warpx(path, sig, e0, e1=1000.0, de=0.01):
    """WarpX MCC tables: two columns, equally spaced energies [eV], cross section [m^2]."""
    n = int(round((e1 - e0) / de)) + 1
    e = e0 + de * np.arange(n)
    np.savetxt(path, np.column_stack([e, sig(e)]), fmt="%.8e")


def main():
    log = []
    say = lambda s="": (print(s), log.append(s))  # noqa: E731

    # ---------------------------------------------------------------- method check on xenon
    xe = fit_gas("Xe", "Xe+")
    meas = {"ionization": warpx_xs(WD / "Xe/ionization.dat"), "excitation": warpx_xs(WD / "Xe/excitation_1.dat"),
            "elastic": warpx_xs(WD / "Xe/electron_scattering.dat")}
    say("Method check on xenon (HallThruster.jl tables vs the measured set in warpx-data):")
    say(f"  {'process':11s} {'measured set vs table':>24s}   {'fit vs table':>15s}   fitted / measured cross section")
    probes = {"ionization": (20, 50, 100, 300), "excitation": (10, 20, 50, 100), "elastic": (1, 5, 10, 50)}
    for name, (eth, sig, te, k, _) in xe.items():
        em, er = rate_error(meas[name], te, k)
        fm, fr = rate_error(sig, te, k)
        ratios = ", ".join(f"{e:g} eV {sig(np.array([e]))[0] / meas[name](np.array([e]))[0]:.2f}" for e in probes[name])
        say(f"  {name:11s} max {em:6.1%}, rms {er:6.1%}   max {fm:5.1%}, rms {fr:5.1%}   {ratios}")

    # ---------------------------------------------------------------- krypton
    kr = fit_gas("Kr", "Kr+")
    say()
    say("Krypton fits (Maxwellian rate error over Te 2-80 eV; HallThruster.jl tables):")
    for name, (eth, sig, te, k, par) in kr.items():
        fm, fr = rate_error(sig, te, k)
        e = np.linspace(eth or 0.01, 1000, 200000)
        s = sig(e)
        law = (f"sigma = {np.exp(par[0]):.3e} (x - 1)^{par[1]:.3f} / x^{par[2]:.3f}, x = E / {eth:g} eV"
               if par is not None else "piecewise linear in log E")
        say(f"  {name:11s} max {fm:5.1%}, rms {fr:5.1%}; peak {s.max():.3e} m^2 at {e[np.argmax(s)]:.1f} eV; {law}")
    ar = {"ionization": warpx_xs(WD / "Ar/ionization.dat"), "elastic": warpx_xs(WD / "Ar/electron_scattering.dat")}
    say()
    say("Plausibility: krypton should sit between argon and xenon (warpx-data measured sets):")
    for name, es in (("ionization", (30, 70, 100, 200)), ("elastic", (2, 5, 10, 20, 50))):
        row = []
        for e in es:
            v = [ar[name](np.array([e]))[0], kr[name][1](np.array([e]))[0], meas[name](np.array([e]))[0]]
            row.append(f"{e:g} eV: Ar {v[0]:.2e}, Kr {v[1]:.2e}, Xe {v[2]:.2e}")
        say(f"  {name}: " + "; ".join(row))

    # ---------------------------------------------------------------- write the WarpX tables
    d = HERE / "Kr"
    d.mkdir(exist_ok=True)
    eth_i, eth_x = kr["ionization"][0], kr["excitation"][0]
    write_warpx(d / "ionization.dat", kr["ionization"][1], eth_i)
    write_warpx(d / "excitation.dat", kr["excitation"][1], eth_x)
    write_warpx(d / "elastic.dat", kr["elastic"][1], 0.0)
    proc = {"source": "fitted to HallThruster.jl krypton rate tables (Dragnea 2020 ionization and elastic, Bordage "
                      "2013 excitation) by xs/make_kr.py",
            "electron": {"elastic": ["elastic.dat", None], "excitation1": ["excitation.dat", eth_x],
                         "ionization": ["ionization.dat", eth_i]},
            "ion": {},
            "ion_note": "no ion-neutral collisions: HallThruster.jl has none, and the PIC is compared with it"}
    (d / "processes.json").write_text(json.dumps(proc, indent=2))
    (HERE / "kr_fit.txt").write_text("\n".join(log) + "\n")
    plot(xe, meas, kr, ar)
    print(f"wrote {d}/*.dat, processes.json, {HERE / 'kr_fit.png'}")


def plot(xe, meas, kr, ar):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))
    for j, name in enumerate(("ionization", "excitation", "elastic")):
        ax = axs[0, j]
        e = np.logspace(-1 if name == "elastic" else np.log10(kr[name][0] * 1.001), 3, 800)
        ax.loglog(e, kr[name][1](e), "C0", lw=2, label="Kr fitted (used by the PIC)")
        ax.loglog(e, xe[name][1](e), "C1--", label="Xe fitted (method check)")
        ax.loglog(e, meas[name](e), "C1", lw=0.8, label="Xe measured (warpx-data)")
        if name in ar:
            ax.loglog(e, ar[name](e), "0.5", lw=0.8, label="Ar measured (warpx-data)")
        ax.set_ylim(1e-22 if name != "elastic" else 1e-21, 1e-18)
        ax.set_xlabel("electron energy [eV]")
        ax.set_ylabel("cross section [m^2]")
        ax.set_title(f"{name} cross section", fontsize=10)
        ax.legend(fontsize=8)
        ax = axs[1, j]
        for gas, fits, c in (("Kr", kr, "C0"), ("Xe", xe, "C1")):
            eth, sig, te, k, _ = fits[name]
            ax.loglog(te, k, c + "o", ms=3, label=f"{gas} HallThruster.jl table")
            ax.loglog(te, kernel(te) @ sig(E), c, label=f"{gas} fitted cross section, Maxwellian")
        te, k = xe[name][2], xe[name][3]
        ax.loglog(te, kernel(te) @ meas[name](E), "C1:", label="Xe measured cross section, Maxwellian")
        ax.axvspan(*TE_REPORT, color="0.94", zorder=0)
        ax.set_xlim(0.5, 120)
        ax.set_ylim(1e-17 if name != "elastic" else 1e-14, 2e-12)
        ax.set_xlabel("Te [eV] (table energy = 3/2 Te)")
        ax.set_ylabel("rate coefficient [m^3/s]")
        ax.set_title(f"{name} rate (grey: Te range checked)", fontsize=10)
        ax.legend(fontsize=7)
    for ax in axs.ravel():
        ax.grid(alpha=0.3, which="both")
    fig.suptitle("Krypton cross sections for the PIC, fitted to the 1-D reference's rate tables", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE / "kr_fit.png", dpi=130)


if __name__ == "__main__":
    sys.exit(main())
