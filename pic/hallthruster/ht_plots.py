"""Figures for the 1-D reference: the breathing mode of one run (discharge current vs time and its spectrum,
plus time-averaged axial profiles) and the transport-coefficient scan.

    python ht_plots.py runs/base        -> results/breathing_<tag>.png, results/scan.png
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE / "results"


def breathing(run):
    run = Path(run)
    o = json.loads((run / "output.json").read_text())
    fr, inp = o["output"]["frames"], o["input"]
    t = np.array([f["t"] for f in fr])
    idis = np.array([f["discharge_current"] for f in fr])
    thr = np.array([f["thrust"] for f in fr]) * 1e3
    t0 = inp["postprocess"]["average_start_time"]
    m = t >= t0
    a = o["output"]["average"]
    z = np.array(a["z"]) * 1e3
    L = inp["config"]["thruster"]["geometry"]["channel_length"] * 1e3
    fig, axs = plt.subplots(2, 2, figsize=(14, 8.5))
    ax = axs[0, 0]
    ax.plot(t * 1e6, idis, lw=0.8)
    ax.axvspan(t0 * 1e6, t[-1] * 1e6, color="0.92", zorder=0, label="averaging window")
    ax.axhline(3.53, color="C3", ls=":", label="set point 3.53 A")
    ax.set_xlabel("t [us]")
    ax.set_ylabel("discharge current [A]")
    ax.legend(fontsize=8)
    c = inp["config"]["anom_model"]
    ax.set_title(f"HallThruster.jl, Kr 3.0 mg/s, 300 V, two-zone Bohm c1 = 1/{1 / c['c1']:.0f}, "
                 f"c2 = 1/{1 / c['c2']:.0f}", fontsize=10)
    ax = axs[0, 1]
    tt, ii = t[m], idis[m]
    dt = np.median(np.diff(tt))
    f = np.fft.rfftfreq(len(ii), dt) / 1e3
    p = np.abs(np.fft.rfft((ii - ii.mean()) * np.hanning(len(ii))))
    ax.plot(f[1:], p[1:])
    k = np.argmax(p[1:]) + 1
    ax.set_xlim(0, 100)
    ax.set_xlabel("frequency [kHz]")
    ax.set_ylabel("|FFT| of the discharge current")
    ax.set_title(f"spectrum over the averaging window: peak {f[k]:.1f} kHz; mean {ii.mean():.2f} A, "
                 f"std {ii.std():.2f} A; thrust {thr[m].mean():.1f} mN", fontsize=10)
    ax = axs[1, 0]
    ax.semilogy(z, a["neutrals"]["Kr"]["n"], label="neutral density")
    ax.semilogy(z, a["ne"], label="electron density")
    ax.set_ylabel("m^-3")
    ax2 = ax.twinx()
    ax2.plot(z, np.array(a["B"]) * 1e4, color="0.5", ls="--", label="B [G]")
    ax2.set_ylabel("B [G]", color="0.5")
    ax.legend(fontsize=8, loc="lower left")
    ax = axs[1, 1]
    ax.plot(z, a["potential"], label="potential [V]")
    ax.plot(z, np.array(a["Tev"]) * 10, label="Te x 10 [eV]")
    ax.plot(z, np.array(a["ions"]["Kr"][0]["u"]) / 100, label="Kr+ velocity [100 m/s]")
    ax.legend(fontsize=8)
    for ax in axs[1]:
        ax.axvline(L, color="k", lw=0.6, ls="--")
        ax.set_xlabel("z from the anode [mm] (exit plane dashed)")
        ax.set_title("time-averaged axial profiles", fontsize=10)
    for ax in axs.ravel():
        ax.grid(alpha=0.3)
    fig.tight_layout()
    RES.mkdir(exist_ok=True)
    fig.savefig(RES / f"breathing_{run.name}.png", dpi=140)
    plt.close(fig)


def scan():
    rows = list(csv.DictReader(open(HERE / "scan.csv")))
    c1s = sorted({int(r["c1_inv"]) for r in rows}, reverse=True)
    c2s = sorted({int(r["c2_inv"]) for r in rows}, reverse=True)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, key, lab, ref in ((axs[0], "I_d_mean", "discharge current [A]", 3.53), (axs[1], "thrust_mN", "thrust [mN]", 56),
                              (axs[2], "osc_rel", "breathing amplitude (std / mean)", None)):
        for c2 in c2s:
            pts = sorted(((int(r["c1_inv"]), float(r[key])) for r in rows if int(r["c2_inv"]) == c2), reverse=True)
            ax.plot([1 / p[0] for p in pts], [p[1] for p in pts], "o-", label=f"c2 = 1/{c2}")
        if ref:
            ax.axhline(ref, color="C3", ls=":", label="set point")
        ax.set_xscale("log")
        ax.set_xlabel("c1 (channel)")
        ax.set_title(lab, fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("1-D reference (HallThruster.jl), Kr 3.0 mg/s, 300 V: two-zone Bohm scan", fontsize=11)
    fig.tight_layout()
    fig.savefig(RES / "scan.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    for r in sys.argv[1:]:
        breathing(r)
    if (HERE / "scan.csv").exists():
        scan()
