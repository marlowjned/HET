"""
stock_design.py

Re-optimisation with the purchasable parts (user, 2026-10-21):
  - front plate: FLAT, no chamfer (bore face length = full 9.525 plate)
  - inner screen: 1-1/2" OD x 0.120" wall C1020 DOM  -> r 16.00-19.05
  - outer screen: 3-1/2" OD x 0.188" wall C1020 DOM  -> r 39.67-44.45
  - rule: nothing past the exit plane
Free: inner screen tip gap (is_gap), outer screen tip gap (os_gap) -- i.e. the
two screen lengths -- and the pole cap edge face length (cap_t) and pole cap
chamfer width (cap_w). Objective: RMS misfit to the reference profile
(target_sigma.py: sigma 11.66 mm, 165 G). Names follow figures/glossary.png.

    python stock_design.py grid|search|check|figures|all
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import target_sigma as ts  # noqa: E402
import peak_forward as pf  # noqa: E402
import permeability as pm  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import PLATE_T, Y_COIL_TOP, Y_EXIT, Y_FP_IN  # noqa: E402

FIG = HERE / "figures" / "stock"
RES = HERE / "results" / "stock"
X0 = dict(pm.STOCK, fp_t=PLATE_T, fp_w=14.0)          # flat front plate (fp_w unused when flat)
GRID = dict(is_gap=[12.0, 15.0, 18.0, 21.0, 24.0, 27.0], os_gap=[1.5, 3.5, 5.5, 8.0, 12.0],
            cap_t=[1.5, 2.5, 3.5, 5.0], cap_w=[1.0, 2.0, 3.0, 4.0, 6.0])
VARS = list(GRID)
SCREEN_LEN = Y_COIL_TOP - Y_FP_IN                      # back plate face to the tip-gap reference plane


def label(x):
    return (f"is_gap {x['is_gap']:g} os_gap {x['os_gap']:g} cap {x['cap_t']:g}/{x['cap_w']:g}")


def grid():
    import itertools
    xs = [dict(X0, **dict(zip(VARS, v))) for v in itertools.product(*GRID.values())]
    rows = pf.run_geoms([(label(x), x) for x in xs], "stock grid")
    pf.RES = RES
    pf.write("grid", rows)
    pf.show(rows, "grid, best 20", 20)
    return rows


def search():
    pf.RES = RES
    best = pf.read("grid")[0]
    x = dict(X0, **{v: best[v] for v in VARS})
    best = pf.run_geoms([("start", x)], "stock search")[0]
    trace = [best]
    it = 0
    for st in (1.0, 0.5, 0.25):
        while it < 80:
            it += 1
            cands = []
            for v in VARS:
                for sgn in (-1, 1):
                    c = dict(x, **{v: round(x[v] + sgn * st, 3)})
                    if pf.feasible(c) is None:
                        cands.append((label(c), c))
            rows = pf.run_geoms(cands, "stock search")
            if rows and rows[0]["ref_rms_pct"] < best["ref_rms_pct"] - 1e-3:
                best = dict(rows[0], case=f"it{it} step{st:g}: {rows[0]['case']}")
                x = dict(X0, **{v: best[v] for v in VARS})
                trace.append(best)
                print(f"it {it} step {st}: misfit {best['ref_rms_pct']:.3f}  {label(x)}", flush=True)
            else:
                break
    pf.write("search_trace", trace)
    (RES / "optimum.json").write_text(json.dumps(x, indent=1))
    pf.show(trace, "search trace")
    return x


def check(x=None):
    """Optimum: +-0.5 mm tolerance on each free knob, and the low-permeability screens."""
    x = x or json.loads((RES / "optimum.json").read_text())
    cases = [("optimum", x)]
    cases += [(f"{v} {s * 0.5:+g}", dict(x, **{v: x[v] + s * 0.5})) for v in VARS for s in (-1, 1)]
    cases += [(f"screens mu {ms:g}, rest 5000", dict(x, mu=5000.0, mu_screen=ms)) for ms in (100.0, 200.0, 500.0)]
    cases += [("all iron mu 100", dict(x, mu=100.0))]
    rows = pf.run_geoms(cases, "stock check")
    order = {n: i for i, (n, _) in enumerate(cases)}
    rows.sort(key=lambda r: order[r["case"]])
    pf.RES = RES
    pf.write("check", rows)
    pf.show(rows, "tolerance (+-0.5 mm) and permeability at the optimum")
    return rows


def figures():
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    xo = json.loads((RES / "optimum.json").read_text())
    xp = json.loads((HERE / "results" / "peak" / "best_chamfers.json").read_text())
    picks = [("previous: chamfered front plate, modelled screens", xp, fd.PALETTE[3]),
             ("this: flat front plate, stock tubes, re-optimised", xo, fd.PALETTE[1])]
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
    fd.drawing(xo, fname="dimensioned_stock",
               note="Buildable design: flat front plate, stock C1020 DOM tubes (1-1/2\" x 0.120\", 3-1/2\" x 0.188\"), "
                    "re-optimised pole cap chamfer and screen tip gaps.")


if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in what:
        if w in ("grid", "all"):
            grid()
        if w in ("search", "all"):
            search()
        if w in ("check", "all"):
            check()
        if w in ("figures", "all"):
            figures()
