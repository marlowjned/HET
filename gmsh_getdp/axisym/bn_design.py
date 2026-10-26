"""
bn_design.py  --  run with AXISYM_CONFIG=bn

Screened design re-optimised for the boron nitride channel part (drawing Boron_Nitride_Channel,
2026-10-20; see geometry.py "bn" config). The screens sit directly against the BN part:
  inner screen inside the 1.395 in bore (outer radius 17.62, 0.1 mm slip clearance)
  outer screen over the 3.5 in body (inner radius 44.45)
Pole cap radius 17.6 (inside the bore); front plate flat with a 3.0 in bore; posts at 2.5 in.
The channel is 41.43 mm long, so the reference (sigma = 0.353 L_ch) has sigma 14.6 mm.

Free: inner and outer screen tip gaps (screen lengths), pole cap edge face length and chamfer
width. Screen walls: stock tube and exact-fit variants.
Scored on the reference misfit plus field shape (noshield.objective: reversal / axial area).

    AXISYM_CONFIG=bn python bn_design.py check|grid|search|variants|figures
"""
import csv
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np

assert os.environ.get("AXISYM_CONFIG") == "bn", "run with AXISYM_CONFIG=bn"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import target_sigma as ts  # noqa: E402
import noshield as ns  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import BN, CHANNEL, DEFAULTS, PLATE_T, Y_EXIT  # noqa: E402

FIG = HERE / "figures" / "bn"
RES = HERE / "results" / "bn"
IN = 25.4
# screens
IS_EXACT = dict(is_r0=DEFAULTS["icoil_ro"] + 0.4, is_t=round(opt.R_KEEP_IN - DEFAULTS["icoil_ro"] - 0.4, 3))
IS_STOCK = dict(is_r0=round((1.375 / 2 - 0.095) * IN, 3), is_t=round(0.095 * IN, 3))      # 1-3/8 x .095
OS_STOCK = dict(os_r=BN["od"] / 2, os_t=round(0.125 * IN, 3))                                # 3-3/4 x .125
OS_THICK = dict(os_r=BN["od"] / 2, os_t=round(opt.R_KEEP_OS - BN["od"] / 2, 3))              # fills to the coils
X0 = dict(fp_t=PLATE_T, fp_w=14.0, cap_r=17.6, cap_t=1.5, cap_w=4.0, is_gap=20.0, is_flare=0.0, os_gap=6.0,
          **IS_EXACT, **OS_STOCK)
VARS = ["is_gap", "os_gap", "cap_t", "cap_w"]
COLS = ns.COLS + ["is_gap", "os_gap", "is_r0", "is_t", "os_r", "os_t"]


def run(cases, name, top=None, sort=True):
    cases = [(n, x) for n, x in cases if opt.feasible(x) is None]
    opt.solve_all([x for _, x in cases])
    rows = []
    for n, x in cases:
        m, _ = ts.score(x)
        m.update(ns.field_shape(x)[0])
        m.update(case=n, objective=ns.objective(m))
        m.update({k: x.get(k, "") for k in COLS if k in x})
        rows.append(m)
    if sort:
        rows.sort(key=lambda r: r["objective"])
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / f"{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    ns.show(rows[:top] if top else rows, name)
    return rows


def label(x):
    return f"is_gap {x['is_gap']:g} os_gap {x['os_gap']:g} cap {x['cap_t']:g}/{x['cap_w']:g}"


def check():
    """Geometry drawing + a few reference cases."""
    FIG.mkdir(parents=True, exist_ok=True)
    fd.FIG = FIG
    fd.drawing(X0, fname="geometry_check",
               note=f"BN channel config: channel r {CHANNEL['r_in']:.2f}-{CHANNEL['r_out']:.2f}, "
                    f"anode {CHANNEL['y_anode'] - Y_EXIT:.2f} mm upstream of the exit plane.")
    print("channel", CHANNEL, "L_ch", ts.L_CH, "sigma_ref", ts.SIG_REF)
    print("screens: inner exact", IS_EXACT, "inner stock", IS_STOCK, "outer stock", OS_STOCK, "outer thick", OS_THICK)
    print("limits: R_KEEP_IN", opt.R_KEEP_IN, "R_CAP_MAX", opt.R_CAP_MAX, "OS_R_MIN", opt.OS_R_MIN,
          "R_KEEP_OS", opt.R_KEEP_OS)
    no = dict(X0, no_is=1, no_os=1)
    return run([("start (screens, gaps 20/6, cap 1.5/4)", X0),
                ("no screens, cap 1.5/4", no),
                ("no screens, flat pole cap", dict(no, cap_t=9.525))], "check", sort=False)


def grid():
    g = dict(is_gap=[10.0, 15.0, 20.0, 25.0, 30.0], os_gap=[2.0, 6.0, 10.0, 15.0, 20.0],
             caps=[(t, w) for t in (1.5, 3.0, 5.0) for w in (2.0, 4.0, 6.0, 8.0)] + [(9.525, 4.0)])
    xs = [dict(X0, is_gap=a, os_gap=b, cap_t=t, cap_w=w) for a, b, (t, w) in itertools.product(g["is_gap"], g["os_gap"], g["caps"])]
    return run([(label(x), x) for x in xs], "grid", top=20)


def search(name="search", x=None):
    if x is None:
        b = next(csv.DictReader(open(RES / "grid.csv")))
        x = dict(X0, **{v: float(b[v]) for v in VARS})
    best = run([("start", x)], f"{name}_start")[0]
    trace = [best]
    for st in (2.0, 1.0, 0.5, 0.25):
        while True:
            cands = []
            for v in VARS:
                for s in (-1, 1):
                    val = round(x[v] + s * st, 3)
                    if v == "cap_t":
                        val = min(val, 9.525)
                    c = dict(x, **{v: val})
                    if c != x and opt.feasible(c) is None:
                        cands.append((label(c), c))
            rows = run(cands, f"{name}_step")
            if rows and rows[0]["objective"] < best["objective"] - 1e-3:
                best = rows[0]
                x = dict(x, **{v: float(best[v]) for v in VARS})
                trace.append(best)
                print(f"step {st}: obj {best['objective']:.3f} misfit {best['ref_rms_pct']:.3f}  {label(x)}", flush=True)
            else:
                break
    (RES / f"{name}.json").write_text(json.dumps(x, indent=1))
    ns.show(trace, f"{name} trace")
    return x


def variants():
    """Screen wall variants at the optimum, each re-searched."""
    xo = json.loads((RES / "search.json").read_text())
    out = {}
    for nm, sc in (("inner stock 1-3/8x.095, outer stock 3-3/4x.125", dict(**IS_STOCK, **OS_STOCK)),
                   ("inner exact, outer thick", dict(**IS_EXACT, **OS_THICK)),
                   ("inner stock, outer thick", dict(**IS_STOCK, **OS_THICK))):
        tag = nm.split(",")[0].split()[1] + "_" + nm.split(",")[1].split()[1]
        out[nm] = search(f"search_{tag}", dict(xo, **sc))
    return out


VARIANTS = [("inner exact, outer stock", "search"), ("inner stock, outer stock", "search_stock_stock"),
            ("inner exact, outer thick", "search_exact_thick"), ("inner stock, outer thick", "search_stock_thick")]


def summary():
    cases = [(lab, json.loads((RES / f"{f}.json").read_text())) for lab, f in VARIANTS]
    return run(cases, "summary", sort=False)


def robust():
    x = json.loads((RES / "search_stock_stock.json").read_text())
    cases = [("stock tubes optimum", x)]
    cases += [(f"{v} {d:+g}", dict(x, **{v: x[v] + d})) for v in VARS for d in (-0.5, 0.5)]
    cases += [(f"screens mu {m:g}, rest 5000", dict(x, mu=5000.0, mu_screen=m)) for m in (100.0, 200.0, 500.0)]
    cases += [(f"all iron mu {m:g}", dict(x, mu=m)) for m in (200.0, 500.0, 5000.0)]
    cases += [(f"joint gap {g:g} mm", dict(x, core_gap=g)) for g in (0.05, 0.1)]
    return run(cases, "robust", sort=False)


def figures():
    import shaping
    from solve import solve_case
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    x = json.loads((RES / "search_stock_stock.json").read_text())
    picks = [("BN part, stock screens (re-optimised)", x, fd.PALETTE[1]),
             ("BN part, start point (gaps 20/6, cap 1.5/4)", X0, "0.55"),
             ("BN part, no screens", dict(x, no_is=1, no_os=1), fd.PALETTE[3])]
    opt.solve_all([xx for _, xx, _ in picks])
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    for ax in axs:
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45,
                label=f"reference (sigma {ts.SIG_REF:.1f} mm on the {ts.L_CH:.1f} mm channel)")
    for lab, xx, c in picks:
        m, (y, Br, By) = ts.score(xx)
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.7,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.2f}%, anode {m['anode_Br_pct']:+.1f}%, {m['P165_W']:.2f} W")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.abs(Br[sel]), color=c, lw=1.7)
    axs[0].set_ylim(-10, 185)
    axs[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axs[1].set_yscale("log")
    axs[1].set_ylim(0.3, 250)
    axs[0].set_title("Channel-averaged B_r (165 G peak)", loc="left", fontsize=10)
    axs[1].set_title("|B_r|, log (anode end)", loc="left", fontsize=10)
    for ax in axs:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.3", lw=0.7)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("B_r (G)")
    fig.savefig(FIG / "profiles.png", dpi=140)
    plt.close(fig)
    opt.solve_all([xx for _, xx, _ in picks[:1] + picks[2:]], map_=True)
    items = []
    for lab, xx, _ in (picks[0], picks[2]):
        m, _ = ts.score(xx)
        m = dict(m, NI_inner=m["NI_inner"] * ts.K, NI_outer=m["NI_outer"] * ts.K, P_W=m["P165_W"])
        items.append((lab, solve_case(opt.key(xx), opt.to_params(xx), map_=True), m, opt.to_params(xx)))
    shaping.FIG = FIG
    shaping.map_fig(items, "bn", ncol=2)
    fd.FIG = FIG
    fd.drawing(x, fname="dimensioned_bn_stock",
               note="BN channel part. Screens directly against it: inner 1-3/8\" x 0.095\" DOM in the 1.395\" bore, "
                    "outer 3-3/4\" x 0.125\" DOM over the 3.5\" body. Flat front plate, 3.0\" bore.")


if __name__ == "__main__":
    for w in sys.argv[1:]:
        {"check": check, "grid": grid, "search": search, "variants": variants, "summary": summary,
         "robust": robust, "figures": figures}[w]()
