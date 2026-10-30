"""Summaries of HallThruster.jl runs: time-averaged performance and the breathing mode (FFT of the discharge
current over the averaging window).

    python ht_analyse.py runs/<tag> [...]
"""
import json
import sys
from pathlib import Path

import numpy as np


def summary(run):
    o = json.loads((Path(run) / "output.json").read_text())
    out, inp = o["output"], o["input"]
    fr = out["frames"]
    t = np.array([f["t"] for f in fr])
    idis = np.array([f["discharge_current"] for f in fr])
    thr = np.array([f["thrust"] for f in fr])
    t0 = inp["postprocess"]["average_start_time"]
    m = t >= t0
    tt, ii = t[m], idis[m]
    s = {"run": Path(run).name, "retcode": out["retcode"], "t_end_us": float(t[-1] * 1e6),
         "I_d_mean": float(ii.mean()), "I_d_std": float(ii.std()), "I_d_pp": float(ii.max() - ii.min()),
         "thrust_mN": float(thr[m].mean() * 1e3)}
    a = out["average"]
    for k in ("anode_eff", "mass_eff", "current_eff", "voltage_eff", "divergence_eff"):
        v = a.get(k)
        s[k] = float(v[0] if isinstance(v, list) else v) if v is not None else None
    if len(tt) > 16:
        dt = np.median(np.diff(tt))
        x = ii - ii.mean()
        f = np.fft.rfftfreq(len(x), dt)
        p = np.abs(np.fft.rfft(x * np.hanning(len(x)))) ** 2
        k = np.argmax(p[1:]) + 1
        s["f_peak_kHz"] = float(f[k] / 1e3)
        s["osc_rel"] = float(ii.std() / ii.mean())
    return s


if __name__ == "__main__":
    for r in sys.argv[1:]:
        s = summary(r)
        print(f"{s['run']:14s} I_d {s['I_d_mean']:.3f} A (std {s['I_d_std']:.3f}, p-p {s['I_d_pp']:.3f}), "
              f"T {s['thrust_mN']:.1f} mN, eta_a {s['anode_eff']:.3f}, eta_m {s['mass_eff']:.3f}, "
              f"f_peak {s.get('f_peak_kHz', float('nan')):.1f} kHz, osc {s.get('osc_rel', float('nan')):.1%}, "
              f"{s['retcode']}")
