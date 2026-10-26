"""
noshield.py

Shieldless configurations, scored on the reference profile AND on field shape
(user, 2026-10-21: the 2 mm disc insert bends the field lines axial near the
anode -- acceptable only if the axial field is weak).

Fixed: no inner or outer screen, flat front plate, pole cap radius 20.5,
nothing past the exit plane. Knobs (names as figures/glossary.png):
  pole cap edge face length / chamfer width   cap_t, cap_w
  inner coil wound over the front fraction f_in of the inner core (same turns,
      build grows as 1/f_in up to the channel keep-out), outer coils f_out
  rear disc insert: inner core + inner coil shortened by core_gap, disc as wide
      as the pole cap (coil_short=1) -- the 10-21 disc design
  split insert: non-magnetic layer of thickness split_g cutting the inner core
      split_z mm upstream of the pole cap back face (radius split_r)

Field-shape metrics (pointwise on the channel grid, at 165 G peak, % of peak):
  bz_back_p95     95th percentile |B_z| over the back half (x/L < 0.5)
  bz_back_max     max |B_z| over the back half
  axial_back      % of back-half points where |B_z| > |B_r| AND |B_z| > 5 % of peak
                  (field lines running anode-to-exit with real strength)
  br_min          most negative pointwise B_r in the channel (reversal if < 0)
  b_min           smallest pointwise |B| in the channel (null if ~0)

    python noshield.py refs|coils|inserts|combos|chamfer|figures
"""
import csv
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import target_sigma as ts  # noqa: E402
import peak_forward as pf  # noqa: E402
import spacer as sp  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import CHANNEL, Y_EXIT  # noqa: E402
from metrics import combine, shape_metrics  # noqa: E402
from solve import solve_case  # noqa: E402

FIG = HERE / "figures" / "noshield"
RES = HERE / "results" / "noshield"
BASE = {k: v for k, v in sp.DISC.items() if k not in ("core_gap", "coil_short")}   # cap 2.5/4, no insert
STOCK = json.loads((HERE / "results" / "stock" / "optimum.json").read_text())
SHAPE = ["bz_back_p95", "bz_back_max", "axial_back", "br_min", "b_min", "bz_anode"]
COLS = ["case", "objective", "ref_rms_pct", "anode_Br_pct", "anode_B_pct", "peak_offset_mm", "exit_Br_pct",
        "sigma_mm", "P165_W", "I165_inner_A", "I165_outer_A", "sat165"] + SHAPE + \
       ["cap_t", "cap_w", "f_in", "f_out", "core_gap", "coil_short", "split_z", "split_g", "split_r"]


def field_shape(x):
    """Pointwise channel field at 165 G peak -> shape metrics (% of peak)."""
    res = solve_case(opt.key(x), opt.to_params(x))
    y, Br, By = combine(res, 1.0, opt.Q_SPEC)
    _, s = shape_metrics(y, Br, By)                       # signed scale to 300 G
    s *= ts.K                                              # -> 165 G
    ci, co = res["channel"]["inner"], res["channel"]["outer"]
    r, yy = ci[:, 0] * 1e3, ci[:, 1] * 1e3
    br = s * (ci[:, 2] + opt.Q_SPEC * co[:, 2]) * 1e4
    bz = s * (ci[:, 3] + opt.Q_SPEC * co[:, 3]) * 1e4
    d = yy - Y_EXIT
    inside = (d >= 0) & (d <= ts.L_CH)
    back = inside & (d >= 0.5 * ts.L_CH)                  # x/L < 0.5
    pk = ts.B_REF
    bmag = np.hypot(br, bz)
    anode = np.abs(d - ts.L_CH) < 0.13
    return dict(bz_back_p95=float(np.percentile(np.abs(bz[back]), 95) / pk * 100),
                bz_back_max=float(np.abs(bz[back]).max() / pk * 100),
                axial_back=float(np.mean((np.abs(bz[back]) > np.abs(br[back])) & (np.abs(bz[back]) > 0.05 * pk)) * 100),
                br_min=float(br[inside].min() / pk * 100),
                b_min=float(bmag[inside].min() / pk * 100),
                bz_anode=float(np.mean(bz[anode]) / pk * 100)), (r, yy, br, bz)


# The screened design (accepted) has bz_back_p95 15.2 %, axial_back 10 %, br_min -1.6 %. The axial MAGNITUDE
# is intrinsic: curl B = 0 gives dB_z/dr = dB_r/dz, so a B_r decaying like the reference forces about
# +-20 % of peak B_z at the walls in the back half. What makes field lines turn axial is B_r going through
# zero (reversal / null), so the shape terms penalise reversal and axial-dominated area beyond the screened
# design, not B_z itself.
SCREENED_AXIAL, SCREENED_BRMIN = 10.0, -1.6


def objective(m):
    return (m["ref_rms_pct"] + 2.0 * max(0.0, SCREENED_BRMIN - m["br_min"])
            + 0.25 * max(0.0, m["axial_back"] - SCREENED_AXIAL))


def run(cases, name, top=None):
    cases = [(n, x) for n, x in cases if pf.feasible(x) is None and feasible(x) is None]
    opt.solve_all([x for _, x in cases])
    rows = []
    for n, x in cases:
        m, _ = ts.score(x)
        m.update(field_shape(x)[0])
        m.update(case=n, objective=objective(m))
        m.update({k: x.get(k, "") for k in ("f_in", "f_out", "core_gap", "coil_short", "split_z", "split_g", "split_r")})
        rows.append(m)
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
    print(f"{'case':>40s} {'obj':>5s} {'misfit':>6s} {'anBr%':>6s} {'Bz95':>5s} {'Bzmx':>5s} {'axl%':>5s} "
          f"{'Brmin':>6s} {'Bmin':>5s} {'pk':>5s} {'P165':>5s}")
    for r in rows:
        print(f"{r['case']:>40s} {r['objective']:5.2f} {r['ref_rms_pct']:6.2f} {r['anode_Br_pct']:6.2f} "
              f"{r['bz_back_p95']:5.1f} {r['bz_back_max']:5.1f} {r['axial_back']:5.1f} {r['br_min']:6.2f} "
              f"{r['b_min']:5.2f} {r['peak_offset_mm']:5.1f} {r['P165_W']:5.2f}")


def feasible(x):
    f = x.get("f_in", 1.0)
    if f < 1.0 and 9.525 + 4.7625 / f > 20.2 + 1e-9:
        return "inner coil build past the channel keep-out"
    return None


# --------------------------------------------------------------------------
def refs():
    return run([("stock screened design", STOCK), ("shieldless, no insert", BASE),
                ("shieldless, 2 mm disc insert", sp.DISC),
                ("shieldless, 1.5 mm disc insert", dict(sp.DISC, core_gap=1.5)),
                ("shieldless, 1 mm disc insert", dict(sp.DISC, core_gap=1.0))], "refs")


def coils():
    cases = [(f"f_in {fi:g} f_out {fo:g}", dict(BASE, f_in=fi, f_out=fo))
             for fi in (1.0, 0.85, 0.7, 0.6, 0.5, 0.45) for fo in (1.0, 0.7, 0.5)]
    return run(cases, "coils")


def inserts():
    cases = [(f"split z{zz:g} g{g:g}", dict(BASE, split_z=zz, split_g=g))
             for zz in (1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0) for g in (0.5, 1.0, 2.0)]
    cases += [(f"disc {g:g}", dict(sp.DISC, core_gap=g)) for g in (0.5, 0.75, 1.0, 1.25, 1.5)]
    return run(cases, "inserts", top=25)


def coils2():
    """Outer coils pushed further forward (build grows as 1/f_out; inner sheet stays clear of the channel)."""
    cases = [(f"f_in {fi:g} f_out {fo:g}", dict(BASE, f_in=fi, f_out=fo))
             for fi in (0.45, 0.5) for fo in (0.4, 0.3, 0.25)]
    return run(cases, "coils2")


def combos():
    cases = []
    for fi, fo in ((0.45, 0.5), (0.45, 0.3), (0.5, 0.5), (0.6, 0.5)):
        tag = f"f_in {fi:g} f_out {fo:g}"
        for g in (0.1, 0.2, 0.3, 0.5, 0.75):
            cases.append((f"{tag} + disc {g:g}", dict(BASE, f_in=fi, f_out=fo, core_gap=g, coil_short=1)))
            cases.append((f"{tag} + core gap {g:g}", dict(BASE, f_in=fi, f_out=fo, core_gap=g)))
            for zz in (10.0, 25.0):
                cases.append((f"{tag} + split z{zz:g} g{g:g}", dict(BASE, f_in=fi, f_out=fo, split_z=zz, split_g=g)))
    return run(cases, "combos", top=30)


def fine():
    cases = [(f"f_in {fi:g} f_out {fo:g} cap 1.5/{w:g}", dict(BASE, f_in=fi, f_out=fo, cap_t=1.5, cap_w=w))
             for fi in (0.447, 0.47, 0.5) for fo in (0.22, 0.25, 0.28) for w in (3.5, 4.0, 4.5, 5.0)]
    return run(cases, "fine", top=12)


def robust(x):
    """Permeability, an accidental contact gap at the inner core / back plate joint, and +-0.5 mm on the cap."""
    cases = [("best", x)]
    cases += [(f"all iron mu {m:g}", dict(x, mu=m)) for m in (200.0, 500.0, 3000.0, 5000.0)]
    cases += [(f"joint gap {g:g} mm", dict(x, core_gap=g)) for g in (0.02, 0.05, 0.1)]
    cases += [(f"cap_w {x['cap_w'] + d:+g}", dict(x, cap_w=x["cap_w"] + d)) for d in (-0.5, 0.5)]
    cases += [("cap_t +0.5", dict(x, cap_t=x["cap_t"] + 0.5))]
    rows = run(cases, "robust")
    return rows


def figures():
    import shaping
    from geometry import DEFAULTS, PLATE_T, Y_COIL_TOP, Y_FP_IN, ring_geometry
    from shaping import B_IN0, B_OUT0, L_CORE
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    best = json.loads((RES / "best.json").read_text())
    picks = [("screened design (stock tubes)", STOCK, fd.PALETTE[0]),
             ("shieldless + 2 mm disc insert", sp.DISC, fd.PALETTE[3]),
             ("shieldless, coils forward (this study)", best, fd.PALETTE[1])]
    # ---- profiles: B_r and channel-averaged B_z
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    axs[0].plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45, label="reference")
    axs[1].plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45)
    rows = []
    for lab, x, c in picks:
        m, (y, Br, By) = ts.score(x)
        sh, (r, yy, br, bz) = field_shape(x)
        m.update(sh)
        rows.append((lab, m, c, (r, yy, br, bz)))
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.7,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.2f}%, anode {m['anode_Br_pct']:+.1f}%, {m['P165_W']:.2f} W")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.abs(Br[sel]), color=c, lw=1.7)
    axs[0].set_ylim(-10, 185)
    axs[0].axhline(0, color="0.3", lw=0.6)
    axs[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axs[0].set_title("Channel-averaged B_r (165 G peak)", loc="left", fontsize=10)
    axs[1].set_yscale("log")
    axs[1].set_ylim(0.3, 250)
    axs[1].set_title("|B_r|, log (dips = reversal)", loc="left", fontsize=10)
    for ax in axs[:2]:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.3", lw=0.7)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("B_r (G)")
    # field-line angle from radial along the channel mid-radius and the walls, back half
    ax = axs[2]
    rmid = 0.5 * (CHANNEL["r_in"] + CHANNEL["r_out"])
    for lab, m, c, (r, yy, br, bz) in rows:
        for rw, ls in ((CHANNEL["r_in"], ":"), (rmid, "-"), (CHANNEL["r_out"], "--")):
            k = np.abs(r - rw) == np.min(np.abs(r - rw))
            dd = yy[k] - Y_EXIT
            o = np.argsort(dd)
            sel = (dd[o] >= 0) & (dd[o] <= ts.L_CH)
            ang = np.degrees(np.arctan2(bz[k][o][sel], br[k][o][sel]))
            ax.plot(1 - dd[o][sel] / ts.L_CH, ang, color=c, ls=ls, lw=1.3)
    ax.axhline(0, color="0.3", lw=0.6)
    for a_ in (-90, 90):
        ax.axhline(a_, color="0.6", lw=0.6, ls=":")
    ax.set_ylim(-180, 180)
    ax.set_yticks([-180, -90, -45, 0, 45, 90, 180])
    ax.set_xlim(0, 1)
    ax.set_xlabel("x / L_ch")
    ax.set_ylabel("field direction from radial (deg)")
    ax.set_title("Field-line angle: solid mid-channel, dotted inner wall, dashed outer wall\n"
                 "(0 = radial, +-90 = axial, +-180 = reversed)", loc="left", fontsize=9)
    fig.savefig(FIG / "compare_profiles.png", dpi=140)
    plt.close(fig)
    # ---- field maps at 165 G
    opt.solve_all([x for _, x, _ in picks], map_=True)
    items = []
    for lab, x, _ in picks:
        m, _ = ts.score(x)
        m = dict(m, NI_inner=m["NI_inner"] * ts.K, NI_outer=m["NI_outer"] * ts.K, P_W=m["P165_W"])
        items.append((lab, solve_case(opt.key(x), opt.to_params(x), map_=True), m, opt.to_params(x)))
    shaping.FIG = FIG
    shaping.map_fig(items, "compare", ncol=3)
    # ---- dimensioned drawing of the coil-forward design
    fig, ax = plt.subplots(figsize=(11, 8.3), constrained_layout=True)
    fd._section(ax, best)
    z = fd.z
    zf = z(Y_FP_IN)
    li, lo = best["f_in"] * L_CORE, best["f_out"] * L_CORE
    bi, bo = B_IN0 / best["f_in"], B_OUT0 / best["f_out"]
    rc = DEFAULTS["core_r"]
    R0 = DEFAULTS["ring_R"] - DEFAULTS["ocore_r"]
    R, t, w = best["cap_r"], best["cap_t"], best["cap_w"]
    fd._ext(ax, [(-1, R), (-6.5, R)])
    fd._dim(ax, (-5.5, 0), (-5.5, R), f"pole cap R {R:g}", rot=90, tpos=0.45)
    fd._ext(ax, [(t, R + 0.3), (t, R + 5.5)])
    fd._ext(ax, [(0, R + 0.3), (0, R + 5.5)])
    ax.annotate("", xy=(0, R + 4.5), xytext=(-3.2, R + 4.5),
                arrowprops=dict(arrowstyle="-|>", color=fd.INK, lw=0.7, mutation_scale=7))
    ax.annotate("", xy=(t, R + 4.5), xytext=(t + 3.2, R + 4.5),
                arrowprops=dict(arrowstyle="-|>", color=fd.INK, lw=0.7, mutation_scale=7))
    ax.text(t + 3.5, R + 4.5, f"pole cap edge face length {t:g}", fontsize=7.5, va="center")
    fd._ext(ax, [(PLATE_T - 0.3, R - w), (PLATE_T - 3.4, R - w)])
    fd._dim(ax, (PLATE_T - 2.8, R - w), (PLATE_T - 2.8, R), f"chamfer width {w:g}", off=(-1.7, -0.6), rot=90, fs=7)
    ax.text(PLATE_T / 2, DEFAULTS["fp_bore_r"] - 2.0, "flat front plate\n(no chamfer)", fontsize=7.2,
            ha="center", va="top", color=fd.CHG)
    # inner coil: length and outer radius
    fd._ext(ax, [(zf + li, rc + bi + 0.3), (zf + li, rc + bi + 3.5)])
    fd._ext(ax, [(zf, rc + bi + 0.3), (zf, rc + bi + 3.5)])
    fd._dim(ax, (zf, rc + bi + 2.8), (zf + li, rc + bi + 2.8), f"inner coil length {li:.1f}", off=(4, -1.2), fs=7.5,
            color="#8a5a00")
    fd._ext(ax, [(zf + li + 0.3, rc + bi), (zf + li + 6, rc + bi)])
    fd._dim(ax, (zf + li + 5, 0), (zf + li + 5, rc + bi), f"inner coil OD R {rc + bi:.2f}", rot=90, tpos=0.3, fs=7.5,
            color="#8a5a00")
    # outer coil: length and build
    fd._ext(ax, [(zf + lo, R0 - bo - 0.3), (zf + lo, R0 - bo - 3.5)])
    fd._ext(ax, [(zf, R0 - bo - 0.3), (zf, R0 - bo - 3.5)])
    fd._dim(ax, (zf, R0 - bo - 2.8), (zf + lo, R0 - bo - 2.8), f"outer coil length {lo:.1f}", off=(0, -1.3), fs=7.5,
            color="#8a5a00")
    fd._ext(ax, [(zf + lo + 0.3, R0 - bo), (zf + lo + 4.5, R0 - bo)])
    fd._ext(ax, [(zf + lo + 0.3, R0), (zf + lo + 4.5, R0)])
    fd._dim(ax, (zf + lo + 3.5, R0 - bo), (zf + lo + 3.5, R0), f"build {bo:.1f}", rot=90, fs=7.5, color="#8a5a00")
    ax.text(zf + lo + 6, R0 - bo / 2, f"each of the 4 outer coils:\nwinding OD {2 * (rc + bo):.1f} around its post",
            fontsize=7.2, va="center", color="#8a5a00")
    ax.text(zf + li + 8, (rc + bi) / 2 + 2, "bare inner core\n(no winding)", fontsize=7.2, color="#6e6e6e")
    ax.text(z(Y_COIL_TOP) - 1, R0 - 6, "bare posts", fontsize=7.2, color="#6e6e6e", ha="right")
    ax.set_xlim(-14, z(Y_COIL_TOP) + PLATE_T + 8)
    ax.set_ylim(-1.5, 76)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("axial distance upstream of the exit plane (mm)")
    ax.set_ylabel("radius (mm)")
    ax.set_title("Shieldless, coils forward: r-z half-section (mm). No screens, no insert. Same turns as the spec "
                 "windings,\nwound shorter and deeper against the front plate / pole cap. "
                 f"Inner coil over the front {best['f_in'] * 100:.0f} % of the core, outer coils over the front "
                 f"{best['f_out'] * 100:.0f} % of the posts.", loc="left", fontsize=9)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"dimensioned_coils_forward.{ext}", dpi=160)
    plt.close(fig)
    for lab, m, c, _ in rows:
        print(f"{lab:40s} misfit {m['ref_rms_pct']:.2f} anode {m['anode_Br_pct']:+.2f} Bz95 {m['bz_back_p95']:.1f} "
              f"axial {m['axial_back']:.1f} Brmin {m['br_min']:.2f} Bmin {m['b_min']:.2f} pk {m['peak_offset_mm']} "
              f"sigma {m['sigma_mm']} P165 {m['P165_W']:.2f} I {m['I165_inner_A']:.2f}/{m['I165_outer_A']:.2f} "
              f"sat {m['sat165']:.2f}")


def thick():
    """Coils-forward design: how much of the rear iron (behind the coils) does the circuit need?"""
    from shaping import L_CORE
    x0 = json.loads((RES / "best.json").read_text())
    bare_in = L_CORE * (1 - x0["f_in"])          # bare inner core behind the inner coil
    bare_out = L_CORE * (1 - x0["f_out"])        # bare post behind the outer coils
    cases = [("coils forward (all iron)", x0)]
    for g in (0.5, 1.0, 2.0, 5.0, 10.0, round(bare_in - 0.05, 2)):
        cases.append((f"inner core insert {g:g} mm", dict(x0, core_gap=g)))
    for g in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, round(bare_out - 0.05, 2)):
        cases.append((f"post insert {g:g} mm", dict(x0, ring_gap=g)))
    for g in (2.0, 10.0):
        cases.append((f"both inserts {g:g} mm", dict(x0, core_gap=g, ring_gap=g)))
    cases.append(("no back plate", dict(x0, no_bp=1)))
    cases.append(("no rear iron (no back plate, bare core and posts removed)",
                  dict(x0, no_bp=1, core_gap=round(bare_in - 0.05, 2), ring_gap=round(bare_out - 0.05, 2))))
    rows = run(cases, "thick")
    order = {n: i for i, (n, _) in enumerate(cases)}
    rows.sort(key=lambda r: order[r["case"]])
    show(rows, "thick inserts / rear iron removed (sorted as run)")
    return rows


def chamfer(name, x0):
    xs = [dict(x0, cap_t=t, cap_w=w) for t, w in itertools.product((1.5, 2.5, 3.5, 5.0), (2.0, 3.0, 4.0, 6.0, 8.0))]
    xs.append(dict(x0, cap_t=9.525))
    return run([(f"{name} cap {x['cap_t']:g}/{x['cap_w']:g}", x) for x in xs], f"chamfer_{name}", top=10)


if __name__ == "__main__":
    for w in sys.argv[1:]:
        if w == "refs":
            refs()
        if w == "coils":
            coils()
        if w == "inserts":
            inserts()
        if w == "combos":
            combos()
        if w == "coils2":
            coils2()
        if w == "thick":
            thick()
        if w == "figures":
            figures()
        if w == "fine":
            fine()
        if w == "robust":
            robust(json.loads((RES / "best.json").read_text()))
        if w == "chamfers":
            chamfer("f0.45-0.3", dict(BASE, f_in=0.45, f_out=0.3))
            chamfer("f0.45-0.25", dict(BASE, f_in=0.45, f_out=0.25))
