"""
report_figs.py  --  run with AXISYM_CONFIG=bn

Figures and the no-shim robustness set for the design review report (2026-10-22):
  Design A = screened, stock DOM tubes (results/bn/search_stock_stock.json)
  Design B = shieldless, coils forward, no post shim (results/bn_noshield/best.json minus ring_gap)
Writes figures/report/*.png and results/bn_noshield/robust_noshim.csv.

    AXISYM_CONFIG=bn python report_figs.py robust|figures
"""
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

assert os.environ.get("AXISYM_CONFIG") == "bn"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import target_sigma as ts  # noqa: E402
import noshield as ns  # noqa: E402
import bn_design as bd  # noqa: E402
import bn_noshield as bns  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import CHANNEL, Y_EXIT  # noqa: E402

FIG = HERE / "figures" / "report"
A = json.loads((HERE / "results" / "a_reopt" / "A_shape.json").read_text())   # A re-optimised for field shape
B = {k: v for k, v in bns.best_design().items() if k != "ring_gap"}
COL_A, COL_B, COL_REF = "#2f6db5", "#d0612b", "#9a9a9a"


def robust():
    x = B
    cases = [("nominal", x)]
    cases += [(f"all iron mu {m:g}", dict(x, mu=m)) for m in (200.0, 500.0, 1000.0, 5000.0)]
    cases += [(f"inner core joint gap {g:g} mm", dict(x, core_gap=g)) for g in (0.02, 0.05, 0.1)]
    cases += [(f"inner coil {d:g} mm off the pole cap", dict(x, in_off=d)) for d in (0.5, 1.0)]
    cases += [(f"outer coils {d:g} mm off the front plate", dict(x, out_off=d)) for d in (0.5, 1.0)]
    cases += [(f"pole cap edge face {d:+g} mm", dict(x, cap_t=x["cap_t"] + d)) for d in (0.5,)]
    cases += [(f"pole cap chamfer {d:+g} mm", dict(x, cap_w=x["cap_w"] + d)) for d in (-0.5, 0.5)]
    cases += [("worst low: mu 500 + 0.05 mm joint gap", dict(x, mu=500.0, core_gap=0.05)),
              ("high: mu 5000", dict(x, mu=5000.0))]
    return bns.run(cases, "robust_noshim", sort=False)


def robustA():
    """Robustness set for the re-optimised A, named to match the B set."""
    import a_reopt as ar
    x = A
    cases = [("nominal", x)]
    cases += [(f"all iron mu {m:g}", dict(x, mu=m)) for m in (200.0, 500.0, 5000.0)]
    cases += [(f"screens mu {m:g}, rest 5000", dict(x, mu=5000.0, mu_screen=m)) for m in (500.0, 200.0)]
    cases += [(f"inner core joint gap {g:g} mm", dict(x, core_gap=g)) for g in (0.05, 0.1)]
    cases += [(f"{v} {d:+g}", dict(x, **{v: x[v] + d})) for v in ("is_gap", "os_gap", "cap_t", "cap_w") for d in (-0.5, 0.5)]
    cases += [("worst low: mu 500 + 0.05 mm joint gap", dict(x, mu=500.0, core_gap=0.05))]
    return ar.run(cases, "robust_A_report", sort=False)


def score(x):
    m, prof = ts.score(x)
    sh, fld = ns.field_shape(x)
    m.update(sh)
    if "n_in" in x:
        m["P165"] = bns.power_fixed_turns(x, m["NI_inner"] * ts.K, m["NI_outer"] * ts.K)
    else:
        m["P165"] = m["P165_W"]
    return m, prof, fld


def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shaping
    from solve import solve_case
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    FIG.mkdir(parents=True, exist_ok=True)
    opt.solve_all([A, B])
    mA, pA, fA = score(A)
    mB, pB, fB = score(B)
    d = np.linspace(0, ts.L_CH, 300)
    xr, br = 1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d)

    # 1. profiles: linear + log
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    for ax in axs:
        ax.plot(xr, br, color=COL_REF, lw=7, alpha=0.4, solid_capstyle="butt",
                label=f"reference: half-Gaussian, sigma {ts.SIG_REF:.1f} mm, 165 G")
        for (lab, m, (y, Br, By), c) in (("A  screened", mA, pA, COL_A), ("B  shieldless, coils forward", mB, pB, COL_B)):
            dd = y - Y_EXIT
            sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
            ax.plot(1 - dd[sel] / ts.L_CH, np.abs(Br[sel]) if ax is axs[1] else Br[sel], color=c, lw=2,
                    label=f"{lab}: {m['ref_rms_pct']:.2f} % RMS, anode {m['anode_Br_pct']:.1f} %")
        ax.axvline(1, color="0.35", lw=0.8)
        ax.set_xlim(0, 1.08)
        ax.set_xlabel("x / L_ch   (0 = anode, 1 = exit plane)")
        ax.set_ylabel("channel-averaged B_r (G)")
        ax.grid(alpha=0.25)
    axs[0].text(1.005, 5, "exit\nplane", fontsize=8, color="0.35")
    axs[0].set_ylim(0, 185)
    axs[0].legend(frameon=False, fontsize=8.5, loc="upper left")
    axs[0].set_title("Both designs track the reference", loc="left", fontsize=11)
    axs[1].set_yscale("log")
    axs[1].set_ylim(1, 300)
    axs[1].set_title("Log scale: the anode end", loc="left", fontsize=11)
    fig.savefig(FIG / "fig_profiles.png", dpi=170)
    plt.close(fig)

    # 2. field-line angle, one panel per design
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.0), constrained_layout=True, sharey=True)
    rmid = 0.5 * (CHANNEL["r_in"] + CHANNEL["r_out"])
    for ax, (lab, (r, yy, bR, bZ), c) in zip(axs, (("A  screened", fA, COL_A), ("B  shieldless, coils forward", fB, COL_B))):
        for rw, ls, nm in ((CHANNEL["r_in"], ":", "inner wall"), (rmid, "-", "mid-channel"), (CHANNEL["r_out"], "--", "outer wall")):
            k = np.abs(r - rw) == np.min(np.abs(r - rw))
            dd = yy[k] - Y_EXIT
            o = np.argsort(dd)
            sel = (dd[o] >= 0) & (dd[o] <= ts.L_CH)
            ax.plot(1 - dd[o][sel] / ts.L_CH, np.degrees(np.arctan2(bZ[k][o][sel], bR[k][o][sel])), color=c, ls=ls,
                    lw=2, label=nm)
        for a_ in (-90, 90):
            ax.axhline(a_, color="0.6", lw=0.8, ls=":")
        ax.axhline(0, color="0.35", lw=0.8)
        ax.text(0.01, 92, "axial", fontsize=8, color="0.45")
        ax.text(0.01, 3, "radial", fontsize=8, color="0.45")
        ax.set_xlim(0, 1)
        ax.set_ylim(-100, 100)
        ax.set_yticks([-90, -45, 0, 45, 90])
        ax.set_xlabel("x / L_ch   (0 = anode, 1 = exit plane)")
        ax.set_title(lab, loc="left", fontsize=11, color=c)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        ax.grid(alpha=0.25)
    axs[0].set_ylabel("field direction from radial (deg)")
    fig.suptitle("Field-line angle near the anode: A is straighter mid-channel, B is gentler at the walls", x=0.01, ha="left",
                 fontsize=11)
    fig.savefig(FIG / "fig_angle.png", dpi=170)
    plt.close(fig)

    # 3. field maps at 165 G
    opt.solve_all([A, B], map_=True)
    items = []
    for lab, x, m in (("A  screened (stock DOM tubes)", A, mA), ("B  shieldless, coils forward", B, mB)):
        mm = dict(m, NI_inner=m["NI_inner"] * ts.K, NI_outer=m["NI_outer"] * ts.K, P_W=m["P165"])
        items.append((lab, solve_case(opt.key(x), opt.to_params(x), map_=True), mm, opt.to_params(x)))
    shaping.FIG = FIG
    shaping.map_fig(items, "fig", ncol=2)            # writes fig_maps.png

    # 4. dimensioned drawings
    fd.FIG = FIG
    fd.drawing(A, fname="fig_design_A",
               note='Design A (re-optimised for field shape): screens directly against the BN part. Inner 1-3/8" x '
                    '0.065" DOM in the 1.395" bore, outer 3-3/4" x 0.125" DOM over the 3.5" body; spec coils.')
    bns.FIG = FIG
    bns.figures(B, "tmp_B")
    os.replace(FIG / "tmp_B" / "dimensioned.png", FIG / "fig_design_B.png")

    # 5. robustness: anode B_r and worst reversal, both designs
    def rows_of(path):
        return list(csv.DictReader(open(path)))
    ra = {r["case"]: r for r in rows_of(HERE / "results" / "a_reopt" / "robust_A_report.csv")}
    rb = {r["case"]: r for r in rows_of(bns.RES / "robust_noshim.csv")}
    pairs = [("nominal", "nominal", "nominal"),
             ("all iron mu 5000", "all iron mu 5000", "all iron mu 5000"),
             ("all iron mu 500", "all iron mu 500", "all iron mu 500"),
             ("all iron mu 200", "all iron mu 200", "all iron mu 200"),
             ("screens mu 200 (A only)", "screens mu 200, rest 5000", None),
             ("inner core joint gap 0.05 mm", "inner core joint gap 0.05 mm", "inner core joint gap 0.05 mm"),
             ("inner core joint gap 0.1 mm", "inner core joint gap 0.1 mm", "inner core joint gap 0.1 mm"),
             ("inner coil 1 mm loose (B only)", None, "inner coil 1 mm off the pole cap"),
             ("worst low: mu 500 + 0.05 mm gap", "worst low: mu 500 + 0.05 mm joint gap",
              "worst low: mu 500 + 0.05 mm joint gap")]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True, sharey=True)
    ylab = [p[0] for p in pairs]
    yy = np.arange(len(pairs))[::-1]
    for ax, key, src, c, name in ((axs[0], 1, ra, COL_A, "A  screened"), (axs[1], 2, rb, COL_B, "B  shieldless")):
        for yi, p in zip(yy, pairs):
            r = src.get(p[key]) if p[key] else None
            if r is None:
                ax.text(0.2, yi, "n/a", va="center", fontsize=8, color="0.6")
                continue
            an, bm = float(r["anode_Br_pct"]), float(r["br_min"])
            ax.plot([min(an, bm), max(an, bm)], [yi, yi], color="0.8", lw=1.2, zorder=1)
            ax.scatter([an], [yi], color=c, s=40, zorder=3, label="B_r at anode" if yi == yy[0] else None)
            ax.scatter([bm], [yi], facecolor="white", edgecolor=c, s=40, lw=1.5, zorder=3,
                       label="worst B_r in channel (< 0 = reversal)" if yi == yy[0] else None)
        ax.axvline(0, color="0.35", lw=0.8)
        ax.axvline(1.8, color=COL_REF, lw=1, ls="--")
        ax.text(1.85, yy[0] + 0.45, "reference anode 1.8 %", fontsize=8, color="0.45")
        ax.set_xlim(-4, 6)
        ax.set_xlabel("% of peak B_r")
        ax.set_title(name, loc="left", fontsize=11, color=c)
        ax.grid(alpha=0.25, axis="x")
    axs[0].legend(frameon=False, fontsize=8, loc="center left")
    axs[0].set_yticks(yy, ylab)
    fig.suptitle("Tolerance cases: anode field and worst reversal", x=0.01, ha="left", fontsize=11)
    fig.savefig(FIG / "fig_robust.png", dpi=170)
    plt.close(fig)
    for lab, m in (("A", mA), ("B", mB)):
        print(lab, {k: round(float(m[k]), 3) for k in ("ref_rms_pct", "anode_Br_pct", "br_min", "axial_back",
                                                        "peak_offset_mm", "P165", "I165_inner_A", "I165_outer_A",
                                                        "sat165", "sigma_mm", "exit_Br_pct")})


if __name__ == "__main__":
    for w in sys.argv[1:]:
        if w != "inner_wall":
            {"robust": robust, "figures": figures, "robustA": robustA}[w]()


def inner_wall():
    """Design A with the inner screen walls Speedy Metals stocks in 1-3/8 in OD (.065, .049)."""
    IN = 25.4

    def inner(wall):
        return dict(is_r0=round((1.375 / 2 - wall) * IN, 3), is_t=round(wall * IN, 3))
    cases = [("A nominal (1-3/8 x .095 inner)", A),
             ("A inner 1-3/8 x .065", dict(A, **inner(0.065))),
             ("A inner 1-3/8 x .049", dict(A, **inner(0.049))),
             ("A inner .065, screens mu 500", dict(A, **inner(0.065), mu=5000.0, mu_screen=500.0)),
             ("A inner .065, screens mu 200", dict(A, **inner(0.065), mu=5000.0, mu_screen=200.0))]
    return bd.run(cases, "inner_wall_check", sort=False)


if __name__ == "__main__" and sys.argv[1:] == ["inner_wall"]:
    inner_wall()
