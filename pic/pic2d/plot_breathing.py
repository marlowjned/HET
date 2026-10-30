"""Breathing-mode comparison: PIC (het_breathing.py, runs/<tag>/monitor.csv) against the 1-D reference
(../hallthruster/runs/<ref>).

    python plot_breathing.py <tag> [--ref ref] [--skip 20]   -> results/breathing_<tag>.png, printed summary

The PIC currents are tallied over --monitor-every steps (10 ns at the default dt); they are smoothed over 1 us for
the plots. Mean, std and the spectrum use the record after --skip microseconds (the start-up from the 1-D state).
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent


def pic_record(tag):
    rows = list(csv.DictReader(open(HERE / "runs" / tag / "monitor.csv")))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def ref_record(tag):
    o = json.loads((HERE.parent / "hallthruster" / "runs" / tag / "output.json").read_text())
    fr = o["output"]["frames"]
    t = np.array([f["t"] for f in fr]) * 1e6
    return {"t_us": t, "I_d": np.array([f["discharge_current"] for f in fr]),
            "t0_us": o["input"]["postprocess"]["average_start_time"] * 1e6}


def smooth(y, n):
    if n <= 1 or len(y) < n:
        return y
    return np.convolve(y, np.ones(n) / n, mode="same")


def spectrum(t_us, y):
    """Uniformly resampled, mean removed, Hann window. Returns f [kHz], |FFT|^2 normalised to its peak."""
    tt = np.linspace(t_us[0], t_us[-1], len(t_us))
    yy = np.interp(tt, t_us, y)
    x = (yy - yy.mean()) * np.hanning(len(yy))
    f = np.fft.rfftfreq(len(x), (tt[1] - tt[0]) * 1e-6) / 1e3
    p = np.abs(np.fft.rfft(x)) ** 2
    p[0] = 0.0
    return f, p / p.max()


def summary(t_us, i_d):
    f, p = spectrum(t_us, i_d)
    k = int(np.argmax(p))
    span = t_us[-1] - t_us[0]
    return {"mean": float(i_d.mean()), "std": float(i_d.std()), "osc": float(i_d.std() / i_d.mean()),
            "f_kHz": float(f[k]), "span_us": float(span), "periods": float(span * 1e-3 * f[k]) if f[k] else 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--ref", default="ref")
    ap.add_argument("--skip", type=float, default=20.0, help="start-up to leave out [us]")
    a = ap.parse_args()
    d = pic_record(a.tag)
    t = d["t_us"]
    dt_rec = float(np.median(np.diff(t)))
    n1 = max(1, int(round(1.0 / dt_rec)))                      # 1 us smoothing
    m = t >= min(a.skip, 0.5 * t[-1])
    i_d = smooth(d["I_discharge"], n1)
    s_pic = summary(t[m], i_d[m])
    r = ref_record(a.ref)
    mr = r["t_us"] >= r["t0_us"]
    s_ref = summary(r["t_us"][mr], r["I_d"][mr])

    fig, axs = plt.subplots(2, 2, figsize=(15, 8.8))
    ax = axs[0, 0]
    ax.plot(t, d["I_discharge"], color="0.8", lw=0.4, label="PIC discharge current (10 ns tallies)")
    ax.plot(t, i_d, "C0", lw=1.2, label="PIC discharge current (1 us mean)")
    ax.plot(t, smooth(d["I_i_beam"], n1), "C1", lw=1.0, label="PIC ion beam current at the cathode plane")
    ax.plot(t, smooth(d["I_e_wall"], n1), "C2", lw=0.8, label="PIC electrons lost to the walls (re-emitted)")
    ax.axhline(s_ref["mean"], color="C3", ls="--", lw=0.8, label=f"1-D reference mean {s_ref['mean']:.2f} A")
    ax.axhline(3.53, color="k", ls=":", lw=0.8, label="set point 3.53 A")
    ax.axvspan(t[0], t[m][0], color="0.93", zorder=0)
    ax.set_ylim(0, max(10.0, 1.2 * float(np.percentile(i_d, 99.5))))
    ax.set_xlabel("t [us] (grey: start-up left out)")
    ax.set_ylabel("current [A]")
    ax.set_title(f"PIC {a.tag}: mean {s_pic['mean']:.2f} A, std {s_pic['std']:.2f} A ({s_pic['osc']:.0%}), "
                 f"spectral peak {s_pic['f_kHz']:.1f} kHz over {s_pic['periods']:.1f} periods", fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    ax = axs[0, 1]
    ax.plot(t, d["neutral_mass_ug"], "C4")
    ax.set_xlabel("t [us]")
    ax.set_ylabel("neutral mass in the domain [ug]")
    ax.set_title("PIC neutral inventory (falls as the discharge burns it, refills from the anode)", fontsize=10)
    ax = axs[1, 0]
    w = r["t_us"] >= r["t_us"][-1] - max(300.0, t[-1])
    ax.plot(r["t_us"][w], r["I_d"][w], "C3", lw=0.9)
    ax.axhline(s_ref["mean"], color="C3", ls="--", lw=0.8)
    ax.set_xlabel("t [us]")
    ax.set_ylabel("discharge current [A]")
    ax.set_title(f"1-D reference ({a.ref}): mean {s_ref['mean']:.2f} A, std {s_ref['std']:.2f} A "
                 f"({s_ref['osc']:.0%}), peak {s_ref['f_kHz']:.1f} kHz", fontsize=10)
    ax = axs[1, 1]
    f, p = spectrum(t[m], i_d[m])
    ax.plot(f, p, "C0", label=f"PIC (after {t[m][0]:.0f} us)")
    f, p = spectrum(r["t_us"][mr], r["I_d"][mr])
    ax.plot(f, p, "C3", label="1-D reference (averaging window)")
    ax.set_xlim(0, 60)
    ax.set_xlabel("frequency [kHz]")
    ax.set_ylabel("power of the discharge current (normalised)")
    ax.set_title("spectra; a PIC record shorter than ~3 periods resolves only ~1/span", fontsize=10)
    ax.legend(fontsize=8)
    for ax in axs.ravel():
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    fig.savefig(out / f"breathing_{a.tag}.png", dpi=130)
    print(f"PIC {a.tag}: {t[-1]:.1f} us recorded; after {t[m][0]:.1f} us: I_d {s_pic['mean']:.3f} A, std "
          f"{s_pic['std']:.3f} A ({s_pic['osc']:.1%}), peak {s_pic['f_kHz']:.2f} kHz ({s_pic['periods']:.1f} periods)")
    print(f"1-D {a.ref}: I_d {s_ref['mean']:.3f} A, std {s_ref['std']:.3f} A ({s_ref['osc']:.1%}), "
          f"peak {s_ref['f_kHz']:.2f} kHz")
    print(f"wrote {out / f'breathing_{a.tag}.png'}")


if __name__ == "__main__":
    main()
