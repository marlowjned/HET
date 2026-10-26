"""
shieldless.py

Best configuration WITHOUT the inner and outer screens (user, 2026-10-21).
Fixed: flat front plate (no chamfer), pole cap radius 20.5, everything else as
built; nothing past the exit plane. Free: the pole cap chamfer -- pole cap
edge face length (cap_t) and pole cap chamfer width (cap_w). Objective: RMS
misfit to the reference profile (target_sigma.py). Names follow
figures/glossary.png.

    python shieldless.py grid|search|figures|all
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import target_sigma as ts  # noqa: E402
import peak_forward as pf  # noqa: E402
import stock_design as sd  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import Y_EXIT  # noqa: E402

FIG = HERE / "figures" / "shieldless"
RES = HERE / "results" / "shieldless"
X0 = dict(sd.X0, no_is=1, no_os=1)
CAP_T = [1.5, 2.5, 3.5, 5.0, 7.0, 9.525]
CAP_W = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0]
VARS = ["cap_t", "cap_w"]


def label(x):
    return "flat pole cap" if x["cap_t"] >= 9.525 else f"cap {x['cap_t']:g}/{x['cap_w']:g}"


def grid():
    xs = [dict(X0, cap_t=t, cap_w=w) for t, w in itertools.product(CAP_T[:-1], CAP_W)] + [dict(X0, cap_t=9.525)]
    rows = pf.run_geoms([(label(x), x) for x in xs], "shieldless grid")
    pf.RES = RES
    pf.write("grid", rows)
    pf.show(rows, "shieldless grid (pole cap edge face length / chamfer width)")
    return rows


def search():
    pf.RES = RES
    b = pf.read("grid")[0]
    x = dict(X0, cap_t=b["cap_t"], cap_w=b["cap_w"])
    best = pf.run_geoms([("start", x)], "shieldless search")[0]
    trace = [best]
    for st in (0.5, 0.25):
        while True:
            cands = [(label(c), c) for v in VARS for s in (-1, 1)
                     if pf.feasible(c := dict(x, **{v: round(min(x[v] + s * st, 9.525), 3)})) is None and c != x]
            rows = pf.run_geoms(cands, "shieldless search")
            if rows and rows[0]["ref_rms_pct"] < best["ref_rms_pct"] - 1e-3:
                best = rows[0]
                x = dict(X0, cap_t=best["cap_t"], cap_w=best["cap_w"])
                trace.append(best)
                print(f"step {st}: misfit {best['ref_rms_pct']:.3f}  {label(x)}", flush=True)
            else:
                break
    pf.write("search_trace", trace)
    (RES / "optimum.json").write_text(json.dumps(x, indent=1))
    pf.show(trace, "search trace")
    return x


def figures():
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    xo = json.loads((RES / "optimum.json").read_text())
    picks = [("shieldless optimum", xo, fd.PALETTE[3])]
    so = sd.RES / "optimum.json"
    if so.exists():
        picks.insert(0, ("with stock screens (stock_design optimum)", json.loads(so.read_text()), fd.PALETTE[1]))
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    for ax in axs:
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45,
                label="reference (ideal)")
    for lab, x, c in picks:
        m, (y, Br, By) = ts.score(x)
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.7,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.2f}%, anode {m['anode_Br_pct']:.1f}%, "
                          f"peak {m['peak_offset_mm']:+.1f} mm")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.clip(Br[sel], 0.1, None), color=c, lw=1.7)
    axs[0].set_ylim(0, 185)
    axs[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axs[1].set_yscale("log")
    axs[1].set_ylim(1, 250)
    axs[0].set_title("B_r along the channel (165 G peak)", loc="left", fontsize=10)
    axs[1].set_title("Log scale (anode end)", loc="left", fontsize=10)
    for ax in axs:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.3", lw=0.7)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("channel-averaged B_r (G)")
    fig.savefig(FIG / "profiles.png", dpi=140)
    plt.close(fig)
    fd.FIG = FIG
    fd.drawing(xo, fname="dimensioned_shieldless",
               note="Shieldless: no inner or outer screen, flat front plate, optimised pole cap chamfer.")


if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in what:
        if w in ("grid", "all"):
            grid()
        if w in ("search", "all"):
            search()
        if w in ("figures", "all"):
            figures()
