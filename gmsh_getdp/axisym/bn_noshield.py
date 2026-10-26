"""
bn_noshield.py  --  run with AXISYM_CONFIG=bn

Shieldless magnetic circuit for the boron nitride channel part (user, 2026-10-22: rerun normal and
forward coils, then keep sweeping to a solution worth committing to; consider manufacturability).

Fixed: BN part (geometry.py "bn"), no screens, flat front plate with 3.0 in bore, pole cap R 17.6
(inside the 1.395 in bore), outer posts at 2.5 in, nothing past the exit plane.

Knobs (names as figures/glossary.png):
  inner coil   n_in turns per layer, layers = ceil(230 / n_in) of AWG 20 at 0.89 mm, seated against the
               pole cap back face; OD must clear the BN bore (<= R 17.4)
  outer coils  n_out turns per layer, layers = ceil(184 / n_out), seated against the front plate; build
               must clear the BN body (<= 9.0 mm, posts at R 63.5)
  core_r       inner core radius (stock bar: 1/2, 5/8, 11/16, 3/4 in)
  cap_T        pole cap thickness; cap_t / cap_w edge face length / chamfer width
  ring_gap     optional non-magnetic shim under the outer posts (trim)
Power is computed for the fixed turns (resistance scales with mean turn length).

    AXISYM_CONFIG=bn python bn_noshield.py refs|grid|capT|chamfer|shim|robust|figures
"""
import csv
import itertools
import json
import math
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
import bn_design as bd  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import BN, DEFAULTS, Y_COIL_TOP, Y_EXIT, Y_FP_IN  # noqa: E402
from shaping import K_IN, K_OUT  # noqa: E402

RES = HERE / "results" / "bn_noshield"
FIG = HERE / "figures" / "bn_noshield"
IN = 25.4
CORES = {"3/4": 0.375 * IN, "11/16": 0.34375 * IN, "5/8": 0.3125 * IN, "1/2": 0.25 * IN}
R_COIL_MAX = BN["bore"] / 2 - 0.32          # inner coil OD inside the BN bore (insulation / clearance)
B_OUT_MAX = 2.5 * IN - DEFAULTS["ocore_r"] - BN["od"] / 2 - 0.5
BASE = dict(bd.X0, no_is=1, no_os=1, cap_t=2.5, cap_w=6.5)
EXTRA = ["core_r", "n_in", "n_out", "cap_T", "ring_gap", "L_in", "B_in", "L_out", "B_out", "P165_fixed",
         "core_flux_T"]
COLS = ns.COLS + EXTRA


def coil_dims(x):
    li, bi = opt.layered(int(x.get("n_in", 46)), opt.N_TURNS_IN)
    lo, bo = opt.layered(int(x.get("n_out", 46)), opt.N_TURNS_OUT)
    return li, bi, lo, bo


def feasible(x):
    li, bi, lo, bo = coil_dims(x)
    rc = x.get("core_r", DEFAULTS["core_r"])
    T = x.get("cap_T", 9.525)
    if rc + bi > R_COIL_MAX + 1e-9:
        return "inner coil OD past the BN bore"
    if bo > B_OUT_MAX + 1e-9:
        return "outer coil hits the BN body"
    if Y_EXIT + T + li > Y_COIL_TOP + 1e-9 or Y_FP_IN + lo > Y_COIL_TOP + 1e-9:
        return "coil longer than its core"
    if x["cap_t"] > T + 1e-9 or x["cap_w"] > x["cap_r"] - rc + 1e-9:
        return "pole cap chamfer out of range"
    return opt.feasible(x)


def power_fixed_turns(x, ni_in, ni_out):
    """Hot power for the fixed turns: R scales with mean turn length from the spec windings."""
    li, bi, lo, bo = coil_dims(x)
    rc = x.get("core_r", DEFAULTS["core_r"])
    mlt_in = (rc + bi / 2) / (DEFAULTS["core_r"] + (DEFAULTS["icoil_ro"] - DEFAULTS["core_r"]) / 2)
    mlt_out = (DEFAULTS["ocore_r"] + bo / 2) / (DEFAULTS["ocore_r"] + DEFAULTS["ocoil_build"] / 2)
    return ni_in ** 2 * K_IN * mlt_in + 4 * ni_out ** 2 * K_OUT * mlt_out


def objective(m):
    """Field fit + field shape (ns.objective) + a mild power term (1 point per W at 165 G)."""
    return ns.objective(m) + 1.0 * m["P165_fixed"]


def label(x):
    core = next((k for k, v in CORES.items() if abs(v - x.get("core_r", DEFAULTS["core_r"])) < 1e-3), "?")
    li, bi, lo, bo = coil_dims(x)
    s = f"core {core} in {int(x['n_in'])}x{math.ceil(230 / x['n_in'])} out {int(x['n_out'])}x{math.ceil(184 / x['n_out'])}"
    s += f" cap {x['cap_t']:g}/{x['cap_w']:g}"
    if x.get("cap_T", 9.525) != 9.525:
        s += f" T{x['cap_T']:g}"
    if x.get("ring_gap"):
        s += f" shim {x['ring_gap']:g}"
    return s


def run(cases, name, top=None, sort=True):
    cases = [(n, x) for n, x in cases if feasible(x) is None]
    opt.solve_all([x for _, x in cases])
    rows = []
    for n, x in cases:
        m, _ = ts.score(x)
        m.update(ns.field_shape(x)[0])
        li, bi, lo, bo = coil_dims(x)
        m.update(L_in=li, B_in=bi, L_out=lo, B_out=bo,
                 P165_fixed=power_fixed_turns(x, m["NI_inner"] * ts.K, m["NI_outer"] * ts.K),
                 core_flux_T=m.get("inner_core_Bflux_max_T", float("nan")) * ts.K)
        m.update(case=n, objective=objective(m))
        m.update({k: x[k] for k in ("core_r", "n_in", "n_out", "cap_T", "ring_gap", "cap_t", "cap_w") if k in x})
        rows.append(m)
    if sort:
        rows.sort(key=lambda r: r["objective"])
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
    print(f"{'case':>52s} {'obj':>5s} {'misfit':>6s} {'anBr%':>6s} {'axl%':>5s} {'Brmin':>6s} {'Bmin':>5s} "
          f"{'pk':>5s} {'P165':>5s} {'coreT':>5s} {'sat':>5s}")
    for r in rows:
        print(f"{r['case']:>52s} {r['objective']:5.2f} {r['ref_rms_pct']:6.2f} {r['anode_Br_pct']:6.2f} "
              f"{r['axial_back']:5.1f} {r['br_min']:6.2f} {r['b_min']:5.2f} {r['peak_offset_mm']:5.1f} "
              f"{r['P165_fixed']:5.2f} {r['core_flux_T']:5.2f} {r['sat165']:5.2f}")


def mk(core="3/4", n_in=46, n_out=46, **kw):
    return dict(BASE, core_r=CORES[core], n_in=n_in, n_out=n_out, **kw)


def n_for_layers(m, turns):
    return math.ceil(turns / m)


# --------------------------------------------------------------------------
def refs():
    cases = [("normal coils, cap 2.5/6.5", mk()),
             ("normal coils, cap 1.5/4", mk(cap_t=1.5, cap_w=4.0)),
             ("normal coils, flat pole cap", mk(cap_t=9.525)),
             ("coils forward (max for 3/4 core), cap 2.5/6.5", mk(n_in=n_for_layers(8, 230), n_out=n_for_layers(10, 184))),
             ("coils forward (max for 3/4 core), cap 1.5/4", mk(n_in=n_for_layers(8, 230), n_out=n_for_layers(10, 184),
                                                               cap_t=1.5, cap_w=4.0))]
    return run(cases, "refs", sort=False)


def grid():
    cases = []
    for core, rc in CORES.items():
        m_max = int((R_COIL_MAX - rc) // opt.WIRE)
        for m_in in range(5, m_max + 1):
            for m_out in (4, 6, 8, 10):
                x = mk(core, n_for_layers(m_in, 230), n_for_layers(m_out, 184))
                cases.append((label(x), x))
    return run(cases, "grid", top=30)


def best_from(name, k=1):
    rows = list(csv.DictReader(open(RES / f"{name}.csv")))[:k]
    out = []
    for r in rows:
        x = dict(BASE, core_r=float(r["core_r"]), n_in=int(float(r["n_in"])), n_out=int(float(r["n_out"])),
                 cap_t=float(r["cap_t"]), cap_w=float(r["cap_w"]))
        for key in ("cap_T", "ring_gap"):
            if r.get(key) not in (None, ""):
                x[key] = float(r[key])
        out.append(x)
    return out


def capT():
    cases = []
    for x in best_from("grid", 6):
        rc = x["core_r"]
        for T in (9.525, 7.9375, 6.35, 4.7625):
            # with a thinner cap the inner coil can start further forward; keep its layers, re-check fit
            xx = dict(x, cap_T=T, cap_t=min(x["cap_t"], T), cap_w=min(x["cap_w"], x["cap_r"] - rc))
            cases.append((label(xx), xx))
    return run(cases, "capT", top=20)


def chamfer():
    cases = []
    for x in best_from("capT", 3):
        T = x.get("cap_T", 9.525)
        for t, w in itertools.product((1.5, 2.5, 3.5, 5.0), (3.0, 4.5, 6.0, 7.5, 9.0)):
            if t <= T and w <= x["cap_r"] - x["core_r"]:
                xx = dict(x, cap_t=t, cap_w=w)
                cases.append((label(xx), xx))
        cases.append((label(dict(x, cap_t=T)), dict(x, cap_t=T)))
    return run(cases, "chamfer", top=20)


def shim():
    x = best_from("chamfer", 1)[0]
    cases = [(label(dict(x, ring_gap=g)) if g else label(x), dict(x, ring_gap=g) if g else x)
             for g in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)]
    return run(cases, "shim", sort=False)


BEST = None


def best_design():
    return json.loads((RES / "best.json").read_text())


def robust():
    x = best_design()
    cases = [("chosen design", x)]
    cases += [(f"all iron mu {m:g}", dict(x, mu=m)) for m in (200.0, 500.0, 1000.0, 5000.0)]
    cases += [(f"inner core joint gap {g:g} mm", dict(x, core_gap=g)) for g in (0.02, 0.05, 0.1)]
    cases += [(f"post shim {g:g} mm", dict(x, ring_gap=g)) for g in (0.0, 0.75, 1.25, 1.5)]
    cases += [(f"inner coil {d:g} mm off the pole cap", dict(x, in_off=d)) for d in (0.5, 1.0, 2.0)]
    cases += [(f"outer coils {d:g} mm off the front plate", dict(x, out_off=d)) for d in (0.5, 1.0, 2.0)]
    cases += [(f"cap_t {d:+g}", dict(x, cap_t=x["cap_t"] + d)) for d in (0.5,)]
    cases += [(f"cap_w {d:+g}", dict(x, cap_w=x["cap_w"] + d)) for d in (-0.5, 0.5)]
    cases += [("worst: mu 500 + joint gap 0.05", dict(x, mu=500.0, core_gap=0.05)),
              ("worst: mu 500 + joint gap 0.05 + inner coil 1 mm off", dict(x, mu=500.0, core_gap=0.05, in_off=1.0)),
              ("high side: mu 5000, shim 0", dict(x, mu=5000.0, ring_gap=0.0))]
    return run(cases, "robust", sort=False)


def figures(x=None, sub=""):
    import shaping
    from matplotlib.patches import Rectangle
    from geometry import CHANNEL, PLATE_T
    from solve import solve_case
    plt = fd._mpl()
    fig_dir = FIG / sub if sub else FIG
    fig_dir.mkdir(parents=True, exist_ok=True)
    x = x or best_design()
    g = x.get("ring_gap", 0.0)
    screened = json.loads((bd.RES / "search_stock_stock.json").read_text())
    picks = [(f"shieldless, coils forward, {'%g mm post shim' % g if g else 'no shim'}", x, fd.PALETTE[1]),
             ("shieldless, normal coils", mk(), fd.PALETTE[3]),
             ("screened (stock tubes), normal coils", screened, fd.PALETTE[0])]
    opt.solve_all([xx for _, xx, _ in picks])
    # ---- profiles + field-line angle
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    for ax in axs[:2]:
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45,
                label=f"reference (sigma {ts.SIG_REF:.1f} mm)")
    rows = []
    for lab, xx, c in picks:
        m, (y, Br, By) = ts.score(xx)
        sh, fld = ns.field_shape(xx)
        m.update(sh)
        P = power_fixed_turns(xx, m["NI_inner"] * ts.K, m["NI_outer"] * ts.K) if "n_in" in xx else m["P165_W"]
        rows.append((lab, m, c, fld, P))
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.7,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.2f}%, anode {m['anode_Br_pct']:+.1f}%, {P:.2f} W")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.abs(Br[sel]), color=c, lw=1.7)
    axs[0].set_ylim(-10, 185)
    axs[0].legend(frameon=False, fontsize=7.3, loc="upper left")
    axs[0].set_title("Channel-averaged B_r, BN part (165 G peak)", loc="left", fontsize=10)
    axs[1].set_yscale("log")
    axs[1].set_ylim(0.3, 250)
    axs[1].set_title("|B_r|, log (anode end)", loc="left", fontsize=10)
    for ax in axs[:2]:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.3", lw=0.7)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("B_r (G)")
    ax = axs[2]
    rmid = 0.5 * (CHANNEL["r_in"] + CHANNEL["r_out"])
    for lab, m, c, (r, yy, br, bz), _ in rows:
        for rw, ls in ((CHANNEL["r_in"], ":"), (rmid, "-"), (CHANNEL["r_out"], "--")):
            k = np.abs(r - rw) == np.min(np.abs(r - rw))
            dd = yy[k] - Y_EXIT
            o = np.argsort(dd)
            sel = (dd[o] >= 0) & (dd[o] <= ts.L_CH)
            ax.plot(1 - dd[o][sel] / ts.L_CH, np.degrees(np.arctan2(bz[k][o][sel], br[k][o][sel])), color=c, ls=ls,
                    lw=1.3)
    for a_ in (-90, 0, 90):
        ax.axhline(a_, color="0.6", lw=0.6, ls=":" if a_ else "-")
    ax.set_ylim(-180, 180)
    ax.set_yticks([-180, -90, -45, 0, 45, 90, 180])
    ax.set_xlim(0, 1)
    ax.set_xlabel("x / L_ch")
    ax.set_ylabel("field direction from radial (deg)")
    ax.set_title("Field-line angle: solid mid-channel, dotted inner wall, dashed outer wall", loc="left", fontsize=9)
    fig.savefig(fig_dir / "profiles.png", dpi=140)
    plt.close(fig)
    # ---- maps at 165 G
    opt.solve_all([xx for _, xx, _ in picks], map_=True)
    items = []
    for lab, xx, _ in picks:
        m, _ = ts.score(xx)
        P = power_fixed_turns(xx, m["NI_inner"] * ts.K, m["NI_outer"] * ts.K) if "n_in" in xx else m["P165_W"]
        m = dict(m, NI_inner=m["NI_inner"] * ts.K, NI_outer=m["NI_outer"] * ts.K, P_W=P)
        items.append((lab, solve_case(opt.key(xx), opt.to_params(xx), map_=True), m, opt.to_params(xx)))
    shaping.FIG = fig_dir
    shaping.map_fig(items, "bn_noshield", ncol=3)
    # ---- dimensioned drawing with the retention hardware
    li, bi, lo, bo = coil_dims(x)
    rc = x["core_r"]
    z = fd.z
    zc0, zo0, zb = z(Y_EXIT + 9.525), z(Y_FP_IN), z(Y_COIL_TOP)
    R0 = DEFAULTS["ring_R"] - DEFAULTS["ocore_r"]
    fig, ax = plt.subplots(figsize=(12, 8.8), constrained_layout=True)
    fd._section(ax, x)
    SLV = dict(fc="#dfe7ef", ec="#5b7083", lw=1.0, hatch="..")
    BNS = dict(fill=False, ec="#7a4fa0", lw=1.2, ls="--")
    # BN part outline (non-magnetic)
    bn_back = zc0 - 9.525 + BN["lip_t"] + BN["body_t"]
    ax.add_patch(Rectangle((0, BN["bore"] / 2), BN["lip_t"], BN["lip_od"] / 2 - BN["bore"] / 2, **BNS))
    ax.add_patch(Rectangle((BN["lip_t"], BN["bore"] / 2), BN["body_t"], BN["od"] / 2 - BN["bore"] / 2, **BNS))
    ax.text(bn_back - 1, BN["od"] / 2 - 1.2, "BN channel part", color="#7a4fa0", fontsize=7.5, ha="right", va="top")
    # spacer sleeves behind the coils (non-magnetic)
    ax.add_patch(Rectangle((zc0 + li + 0.3, rc + 0.3), zb - (zc0 + li) - 1.3, bi - 0.6, **SLV))
    ax.add_patch(Rectangle((zo0 + lo + 0.3, R0 - bo + 0.3), zb - (zo0 + lo) - 1.3, bo - 0.6, **SLV))
    ax.text((zc0 + li + zb) / 2, rc + bi / 2, "spacer sleeve\n+ wave spring", fontsize=7, ha="center", va="center",
            color="#2f4252")
    ax.text((zo0 + lo + zb) / 2, R0 - bo / 2, "spacer sleeve + wave spring\n(each post)", fontsize=7,
            ha="center", va="center", color="#2f4252")
    # shim under the posts
    if g:
        r0p, r1p = DEFAULTS["ring_R"] - DEFAULTS["ocore_r"], DEFAULTS["ring_R"] + DEFAULTS["ocore_r"]
        ax.add_patch(Rectangle((zb - g, r0p), g, r1p - r0p, fc="#f3c7c7", ec="#b03030", lw=1.0))
        ax.annotate(f"{g:g} mm non-magnetic shim\nunder each post (field trim)", xy=(zb - g / 2, r0p + 4),
                    xytext=(zb - 22, 70), fontsize=7.5, color="#b03030",
                    arrowprops=dict(arrowstyle="-|>", color="#b03030", lw=0.9, mutation_scale=8))
    # dimensions: coils
    fd._ext(ax, [(zc0, rc + bi + 0.3), (zc0, rc + bi + 3.4)])
    fd._ext(ax, [(zc0 + li, rc + bi + 0.3), (zc0 + li, rc + bi + 3.4)])
    fd._dim(ax, (zc0, rc + bi + 2.7), (zc0 + li, rc + bi + 2.7),
            f"inner coil {li:.1f} long: 29 turns x 8 layers", off=(0, 1.2), fs=7.3, color="#8a5a00")
    fd._ext(ax, [(zc0 + li + 0.3, rc + bi), (zc0 + li + 9, rc + bi)])
    fd._dim(ax, (zc0 + li + 8, 0), (zc0 + li + 8, rc + bi), f"inner coil OD R {rc + bi:.2f}", rot=90, tpos=0.25,
            fs=7.3, color="#8a5a00")
    fd._ext(ax, [(zo0, R0 - bo - 0.3), (zo0, R0 - bo - 3.4)])
    fd._ext(ax, [(zo0 + lo, R0 - bo - 0.3), (zo0 + lo, R0 - bo - 3.4)])
    fd._dim(ax, (zo0, R0 - bo - 2.7), (zo0 + lo, R0 - bo - 2.7),
            f"outer coils {lo:.1f} long: 21 turns x 9 layers", off=(0, -1.3), fs=7.3, color="#8a5a00")
    ax.text(zo0 + lo + 2, R0 - bo - 6.8, f"build {bo:.2f}: clears the BN body by "
            f"{2.5 * IN - DEFAULTS['ocore_r'] - bo - BN['od'] / 2:.2f} (posts at R 63.5)", fontsize=7, color="#8a5a00")
    # pole cap
    R, t, w = x["cap_r"], x["cap_t"], x["cap_w"]
    fd._ext(ax, [(-1, R), (-6.5, R)])
    fd._dim(ax, (-5.5, 0), (-5.5, R), f"pole cap R {R:g}", rot=90, tpos=0.45, fs=7.3)
    ax.text(-13, R + 6.0, f"pole cap edge face length {t:g},\nchamfer width {w:g}", fontsize=7.3, ha="left")
    ax.text(PLATE_T / 2, DEFAULTS["fp_bore_r"] + 2, "flat front plate", fontsize=7, ha="center", color=fd.CHG)
    ax.text(11, 3, "inner core 3/4 in: press-fit into a counterbore\nin the back plate (or one piece)",
            fontsize=7, ha="left", color="#6e6e6e")
    ax.set_xlim(-14, zb + PLATE_T + 6)
    ax.set_ylim(-1.5, 76)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("axial distance upstream of the exit plane (mm)")
    ax.set_ylabel("radius (mm)")
    ax.set_title("Shieldless, BN channel part: coils forward. No screens, no inserts in the flux path. "
                 "AWG 20, 230 / 184 turns (spec).\nBlue = A36 parts that changed; orange = windings; dotted grey = "
                 "non-magnetic spacer sleeves; dashed purple = BN part" + ("; red = post shim." if g else "."),
                 loc="left", fontsize=9)
    for ext in ("png", "svg"):
        fig.savefig(fig_dir / f"dimensioned.{ext}", dpi=160)
    plt.close(fig)
    for lab, m, c, _, P in rows:
        print(f"{lab:48s} misfit {m['ref_rms_pct']:.2f} anode {m['anode_Br_pct']:+.2f} axial {m['axial_back']:.1f} "
              f"Brmin {m['br_min']:.2f} Bmin {m['b_min']:.2f} pk {m['peak_offset_mm']:.1f} P165 {P:.2f} "
              f"I {m['I165_inner_A']:.2f}/{m['I165_outer_A']:.2f} sat {m['sat165']:.2f}")


def finish():
    """3/4 in core (existing part): outer layering for clearance, pole cap variants, post shim trim."""
    cases = []
    for n_out in (19, 21, 23):
        for t, w in ((1.5, 6.0), (2.5, 6.0), (2.5, 4.5)):
            for g in (0.0, 0.5, 1.0, 2.0):
                x = mk("3/4", 29, n_out, cap_t=t, cap_w=w)
                if g:
                    x["ring_gap"] = g
                cases.append((label(x), x))
    return run(cases, "finish", top=24)


if __name__ == "__main__":
    for w in sys.argv[1:]:
        {"refs": refs, "grid": grid, "capT": capT, "chamfer": chamfer, "shim": shim, "finish": finish, "robust": robust, "figures": figures,
         "noshim": lambda: figures({k: v for k, v in best_design().items() if k != "ring_gap"}, "no_shim")}[w]()
