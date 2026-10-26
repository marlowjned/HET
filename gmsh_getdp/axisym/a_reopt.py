"""
a_reopt.py  --  run with AXISYM_CONFIG=bn

Screened design re-optimised with field shape weighted heavily, and the hybrid (screens + coils wound
forward), user 2026-10-22. Screens = the stocked tubes:
  inner 1-3/8 in OD x 0.065 in wall DOM (R 15.81-17.46) in the BN bore
  outer 3-3/4 in OD x 0.125 in wall DOM (R 44.45-47.63) over the BN body

Field-shape objective (lower is better), all in % of peak or degrees:
  J = misfit + 0.15 * axial_back + 0.03 * ang_mid + 0.01 * ang_wall + 2 * max(0, -br_min)
  axial_back  % of back-half points with |B_z| > |B_r| and |B_z| > 5 % of peak
  ang_mid     max |angle of B from radial| mid-channel, within 15 % of L_ch of the anode
  ang_wall    max |angle| at either wall in the same zone
  br_min      most negative pointwise B_r in the channel (reversal)
Design B scores about 3.3 on it, the 10-21 screened design about 4.9.

    AXISYM_CONFIG=bn python a_reopt.py refs|gridA|gridH|search|figures
"""
import csv
import itertools
import json
import math
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
from geometry import CHANNEL, DEFAULTS, Y_EXIT  # noqa: E402

RES = HERE / "results" / "a_reopt"
IN = 25.4
SCR = dict(is_r0=round((1.375 / 2 - 0.065) * IN, 3), is_t=round(0.065 * IN, 3),
           os_r=3.5 / 2 * IN, os_t=round(0.125 * IN, 3))
A0 = dict(json.loads((bd.RES / "search_stock_stock.json").read_text()), **SCR)
B0 = {k: v for k, v in bns.best_design().items() if k != "ring_gap"}
R_COIL_MAX_SCR = SCR["is_r0"] - 0.3                      # inner coil OD inside the inner screen
COLS = ["case", "J", "ref_rms_pct", "anode_Br_pct", "axial_back", "ang_mid", "ang_wall", "br_min", "b_min",
        "peak_offset_mm", "exit_Br_pct", "P165", "I165_inner_A", "I165_outer_A", "sat165",
        "is_gap", "os_gap", "cap_t", "cap_w", "n_in", "n_out"]


def angles(fld):
    r, yy, br, bz = fld
    d = yy - Y_EXIT
    band = (d >= 0.85 * ts.L_CH) & (d <= ts.L_CH)
    ang = np.abs(np.degrees(np.arctan2(bz, br)))
    rmid = 0.5 * (CHANNEL["r_in"] + CHANNEL["r_out"])
    rs = np.unique(np.round(r, 6))
    r_mid = rs[np.argmin(np.abs(rs - rmid))]
    mid = band & (np.abs(r - r_mid) < 1e-6)
    wall = band & ((np.abs(r - rs.min()) < 1e-6) | (np.abs(r - rs.max()) < 1e-6))
    return float(ang[mid].max()), float(ang[wall].max())


def objective(m):
    return (m["ref_rms_pct"] + 0.15 * m["axial_back"] + 0.03 * m["ang_mid"] + 0.01 * m["ang_wall"]
            + 2.0 * max(0.0, -m["br_min"]))


def feasible(x):
    if "n_in" in x and not x.get("no_is"):
        li, bi = opt.layered(int(x["n_in"]), opt.N_TURNS_IN)
        if DEFAULTS["core_r"] + bi > R_COIL_MAX_SCR + 1e-9:
            return "inner coil hits the inner screen"
    if "n_out" in x and not x.get("no_os"):
        # outer coils reach in to R 63.5 - post r - build at the 4 points nearest the axis; where that is inside
        # the outer screen's OD, the screen tip must stop behind the coils
        lo, bo = opt.layered(int(x["n_out"]), opt.N_TURNS_OUT)
        if 2.5 * IN - DEFAULTS["ocore_r"] - bo < x["os_r"] + x["os_t"] + 0.5 and x["os_gap"] < lo + 0.5:
            return "outer coils clash with the outer screen"
    if "n_in" in x or "n_out" in x:
        why = bns.feasible(dict(x, n_in=x.get("n_in", 46), n_out=x.get("n_out", 46)))
        if why:
            return why
    return opt.feasible(x)


def score(x):
    m, prof = ts.score(x)
    sh, fld = ns.field_shape(x)
    m.update(sh)
    m["ang_mid"], m["ang_wall"] = angles(fld)
    m["P165"] = bns.power_fixed_turns(dict(x, n_in=x.get("n_in", 46), n_out=x.get("n_out", 46)),
                                      m["NI_inner"] * ts.K, m["NI_outer"] * ts.K)
    m["J"] = objective(m)
    return m, prof, fld


def run(cases, name, top=None, sort=True):
    cases = [(n, x) for n, x in cases if feasible(x) is None]
    opt.solve_all([x for _, x in cases])
    rows = []
    for n, x in cases:
        m, _, _ = score(x)
        m.update(case=n, **{k: x.get(k, "") for k in ("is_gap", "os_gap", "cap_t", "cap_w", "n_in", "n_out")})
        rows.append(m)
    if sort:
        rows.sort(key=lambda r: r["J"])
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / f"{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    show(rows[:top] if top else rows, name)
    return rows


def show(rows, title):
    print(f"\n{title}")
    print(f"{'case':>50s} {'J':>5s} {'misfit':>6s} {'anBr%':>6s} {'axl%':>5s} {'angM':>5s} {'angW':>5s} "
          f"{'Brmin':>6s} {'pk':>5s} {'P165':>5s} {'sat':>5s}")
    for r in rows:
        print(f"{r['case']:>50s} {r['J']:5.2f} {r['ref_rms_pct']:6.2f} {r['anode_Br_pct']:6.2f} "
              f"{r['axial_back']:5.1f} {r['ang_mid']:5.0f} {r['ang_wall']:5.0f} {r['br_min']:6.2f} "
              f"{r['peak_offset_mm']:5.1f} {r['P165']:5.2f} {r['sat165']:5.2f}")


def label(x):
    s = f"is {x['is_gap']:g} os {x['os_gap']:g} cap {x['cap_t']:g}/{x['cap_w']:g}"
    if "n_in" in x or "n_out" in x:
        s += (f" in {int(x.get('n_in', 46))}x{math.ceil(230 / x.get('n_in', 46))}"
              f" out {int(x.get('n_out', 46))}x{math.ceil(184 / x.get('n_out', 46))}")
    return s


CAPS = [(1.5, 4.0), (1.5, 6.0), (2.5, 4.5), (2.5, 6.5), (3.5, 6.0), (1.5, 8.0)]


def refs():
    return run([("A as reported (stocked tubes)", A0), ("B shieldless no shim", B0)], "refs", sort=False)


def gridA():
    xs = [dict(A0, is_gap=a, os_gap=b, cap_t=t, cap_w=w) for a, b, (t, w) in
          itertools.product((10.0, 14.0, 18.0, 22.0, 26.0), (6.0, 12.0, 17.0, 21.0, 25.0, 30.0), CAPS)]
    return run([(label(x), x) for x in xs], "gridA", top=20)


def gridH():
    xs = []
    for n_in, n_out in itertools.product((46, 39), (46, 31, 23, 21)):
        for a, b, (t, w) in itertools.product((14.0, 18.0, 22.0, 26.0), (6.0, 12.0, 18.0, 21.0, 25.0),
                                              ((1.5, 6.0), (2.5, 6.5), (2.5, 4.5))):
            xs.append(dict(A0, is_gap=a, os_gap=b, cap_t=t, cap_w=w, n_in=n_in, n_out=n_out))
    return run([(label(x), x) for x in xs], "gridH", top=25)


def search_from(name, x, vars_):
    best = run([("start", x)], f"{name}_start")[0]
    trace = [best]
    for st in (2.0, 1.0, 0.5):
        while True:
            cands = []
            for v in vars_:
                for s in (-1, 1):
                    val = x[v] + s * st
                    if v in ("n_in", "n_out"):
                        continue
                    c = dict(x, **{v: round(min(val, 9.525) if v == "cap_t" else val, 3)})
                    if c != x and feasible(c) is None:
                        cands.append((label(c), c))
            rows = run(cands, f"{name}_step")
            if rows and rows[0]["J"] < best["J"] - 1e-3:
                best = rows[0]
                x = dict(x, **{v: float(best[v]) for v in vars_ if best.get(v, "") != ""})
                trace.append(best)
                print(f"{name} step {st}: J {best['J']:.3f} misfit {best['ref_rms_pct']:.2f}  {label(x)}", flush=True)
            else:
                break
    (RES / f"{name}.json").write_text(json.dumps(x, indent=1))
    show(trace, f"{name} trace")
    return x


def best_row(name):
    r = next(csv.DictReader(open(RES / f"{name}.csv")))
    x = dict(A0, is_gap=float(r["is_gap"]), os_gap=float(r["os_gap"]), cap_t=float(r["cap_t"]), cap_w=float(r["cap_w"]))
    if r.get("n_in"):
        x.update(n_in=int(float(r["n_in"])), n_out=int(float(r["n_out"])))
    return x


def search():
    xh = search_from("H_shape", best_row("gridH"), ["is_gap", "os_gap", "cap_t", "cap_w"])
    return xh


if __name__ == "__main__":
    for w in sys.argv[1:]:
        if w in ("refs", "gridA", "gridH", "search"):
            {"refs": refs, "gridA": gridA, "gridH": gridH, "search": search}[w]()


def robust():
    out = []
    for name in ("A_shape", "H_shape"):
        x = json.loads((RES / f"{name}.json").read_text())
        cases = [(f"{name} nominal", x)]
        cases += [(f"{name} screens mu {m:g}", dict(x, mu=5000.0, mu_screen=m)) for m in (500.0, 200.0)]
        cases += [(f"{name} all iron mu {m:g}", dict(x, mu=m)) for m in (500.0, 5000.0)]
        cases += [(f"{name} joint gap 0.05", dict(x, core_gap=0.05))]
        cases += [(f"{name} {v} {d:+g}", dict(x, **{v: x[v] + d})) for v in ("is_gap", "os_gap") for d in (-0.5, 0.5)]
        out += run(cases, f"robust_{name}", sort=False)
    return out


def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shaping
    from solve import solve_case
    fig_dir = HERE / "figures" / "a_reopt"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    xa = json.loads((RES / "A_shape.json").read_text())
    xh = json.loads((RES / "H_shape.json").read_text())
    picks = [("A as reported", A0, "#9ab7dc", "--"), ("A re-optimised for field shape", xa, "#2f6db5", "-"),
             ("Hybrid: screens + coils forward", xh, "#2a9d6f", "-"), ("B shieldless, coils forward", B0, "#d0612b", "-")]
    opt.solve_all([x for _, x, _, _ in picks])
    res = [(lab, *score(x), c, ls, x) for lab, x, c, ls in picks]
    d = np.linspace(0, ts.L_CH, 300)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    for ax in axs:
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#9a9a9a", lw=7, alpha=0.4, solid_capstyle="butt",
                label="reference")
    for lab, m, (y, Br, By), fld, c, ls, x in res:
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, ls=ls, lw=1.8,
                    label=f"{lab}: {m['ref_rms_pct']:.2f} %, anode {m['anode_Br_pct']:.1f} %, {m['P165']:.2f} W")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.abs(Br[sel]), color=c, ls=ls, lw=1.8)
    axs[0].set_ylim(0, 185)
    axs[0].legend(frameon=False, fontsize=8, loc="upper left")
    axs[1].set_yscale("log")
    axs[1].set_ylim(0.5, 300)
    for ax in axs:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.35", lw=0.8)
        ax.set_xlabel("x / L_ch   (0 = anode, 1 = exit plane)")
        ax.set_ylabel("channel-averaged B_r (G)")
        ax.grid(alpha=0.25)
    axs[0].set_title("B_r profiles (165 G peak)", loc="left", fontsize=11)
    axs[1].set_title("Log scale: the anode end", loc="left", fontsize=11)
    fig.savefig(fig_dir / "profiles.png", dpi=170)
    plt.close(fig)
    fig, axs = plt.subplots(1, 4, figsize=(17, 3.9), constrained_layout=True, sharey=True)
    rmid = 0.5 * (CHANNEL["r_in"] + CHANNEL["r_out"])
    for ax, (lab, m, prof, (r, yy, bR, bZ), c, ls, x) in zip(axs, res):
        for rw, lst, nm in ((CHANNEL["r_in"], ":", "inner wall"), (rmid, "-", "mid-channel"), (CHANNEL["r_out"], "--", "outer wall")):
            k = np.abs(r - rw) == np.min(np.abs(r - rw))
            dd = yy[k] - Y_EXIT
            o = np.argsort(dd)
            sel = (dd[o] >= 0) & (dd[o] <= ts.L_CH)
            ax.plot(1 - dd[o][sel] / ts.L_CH, np.degrees(np.arctan2(bZ[k][o][sel], bR[k][o][sel])), color=c, ls=lst,
                    lw=1.8, label=nm)
        for a_ in (-90, 90):
            ax.axhline(a_, color="0.6", lw=0.8, ls=":")
        ax.axhline(0, color="0.35", lw=0.8)
        ax.set_ylim(-100, 100)
        ax.set_yticks([-90, -45, 0, 45, 90])
        ax.set_xlim(0, 1)
        ax.set_xlabel("x / L_ch")
        ax.set_title(f"{lab}\nscore {m['J']:.2f}", loc="left", fontsize=10, color=c)
        ax.grid(alpha=0.25)
    axs[0].set_ylabel("field direction from radial (deg)")
    axs[-1].legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(fig_dir / "angles.png", dpi=170)
    plt.close(fig)
    opt.solve_all([x for _, x, _, _ in picks[1:]], map_=True)
    items = []
    for lab, m, prof, fld, c, ls, x in res[1:]:
        mm = dict(m, NI_inner=m["NI_inner"] * ts.K, NI_outer=m["NI_outer"] * ts.K, P_W=m["P165"])
        items.append((lab, solve_case(opt.key(x), opt.to_params(x), map_=True), mm, opt.to_params(x)))
    shaping.FIG = fig_dir
    shaping.map_fig(items, "fig", ncol=3)
    import final_design as fd
    fd.FIG = fig_dir
    fd.drawing(xa, fname="dimensioned_A_shape",
               note="A re-optimised for field shape: stocked tubes (1-3/8 x 0.065 inner, 3-3/4 x 0.125 outer), spec coils.")
    fd.drawing(xh, fname="dimensioned_hybrid",
               note="Hybrid: same screens; inner coil 46 x 5 layers, outer coils 31 x 6 layers, wound forward.")


if __name__ == "__main__" and sys.argv[1:2] in (["robust"], ["figures"]):
    {"robust": robust, "figures": figures}[sys.argv[1]]()
