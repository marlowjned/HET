"""
spacer.py

Does lengthening the magnetic circuit help the SHIELDLESS design? (user,
2026-10-21). The back plate moves L upstream; the channel, anode, front plate
and pole cap stay put. Three ways to fill the extra length:
  1 taller        cores, outer posts and both coils all longer by L (same
                  ampere-turns spread over a longer winding: lower power)
  2 iron spacer   cores and posts longer by L, coils unchanged (bare iron behind
                  the coils)
  3 inert spacer  cores, posts and coils unchanged; non-magnetic spacer of
                  length L between their rear ends and the back plate
Start: shieldless optimum (flat front plate, pole cap edge face length 4.5,
chamfer width 9.525, no screens). Scored against the reference profile.

    python spacer.py sweep|chamfer|figures
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
import shieldless as sl  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import Y_EXIT  # noqa: E402

FIG = HERE / "figures" / "spacer"
RES = HERE / "results" / "spacer"
X0 = json.loads((sl.RES / "optimum.json").read_text())
MODES = {1: "taller", 2: "iron spacer", 3: "inert spacer", 4: "inert spacer, inner core only"}
LS = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0]


def sweep():
    cases = [("no spacer", X0)]
    cases += [(f"{MODES[m]} L{L:g}", dict(X0, spacer=L, spacer_mode=m)) for m in MODES for L in LS]
    rows = pf.run_geoms(cases, "spacer")
    order = {n: i for i, (n, _) in enumerate(cases)}
    rows.sort(key=lambda r: order[r["case"]])
    pf.RES = RES
    pf.write("sweep", rows)
    pf.show(rows, "spacer sweep (shieldless optimum chamfer)")
    return rows


def small():
    """Small inert spacers, where the anode field crosses zero."""
    cases = [(f"{MODES[m]} L{L:g}", dict(X0, spacer=L, spacer_mode=m)) for m in (3, 4)
             for L in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)]
    rows = pf.run_geoms(cases, "spacer small")
    order = {n: i for i, (n, _) in enumerate(cases)}
    rows.sort(key=lambda r: order[r["case"]])
    pf.RES = RES
    pf.write("small", rows)
    pf.show(rows, "small inert spacers")
    return rows


def chamfer(mode, L):
    """Re-sweep the pole cap chamfer at one spacer setting."""
    xs = [dict(X0, spacer=L, spacer_mode=mode, cap_t=t, cap_w=w)
          for t, w in itertools.product(sl.CAP_T[:-1], sl.CAP_W)] + [dict(X0, spacer=L, spacer_mode=mode, cap_t=9.525)]
    rows = pf.run_geoms([(f"{MODES[mode]} L{L:g} {sl.label(x)}", x) for x in xs], "spacer chamfer")
    pf.RES = RES
    pf.write(f"chamfer_m{mode}_L{L:g}", rows)
    pf.show(rows, f"pole cap chamfer at {MODES[mode]} L{L:g}, best 10", 10)
    return rows


def figures(picks):
    """picks: [(label, x)]."""
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    for ax in axs:
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45,
                label="reference (ideal)")
    for (lab, x), c in zip(picks, fd.PALETTE):
        m, (y, Br, By) = ts.score(x)
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.6,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.1f}%, anode {m['anode_Br_pct']:.1f}%, "
                          f"{m['P165_W']:.2f} W")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.clip(Br[sel], 0.1, None), color=c, lw=1.6)
    axs[0].set_ylim(0, 185)
    axs[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axs[1].set_yscale("log")
    axs[1].set_ylim(1, 250)
    for ax in axs:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.3", lw=0.7)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("channel-averaged B_r (G)")
    axs[0].set_title("Shieldless with a spacer (165 G peak)", loc="left", fontsize=10)
    axs[1].set_title("Log scale (anode end)", loc="left", fontsize=10)
    fig.savefig(FIG / "profiles.png", dpi=140)
    plt.close(fig)


BEST = dict(X0, spacer_mode=4, cap_t=2.5, cap_w=4.0, spacer=2.0)


def best():
    """Spacer-thickness tolerance at the best shieldless + spacer design, and the comparison figure."""
    import stock_design as sd
    cases = [(f"inner-core spacer {L:g} mm", dict(BEST, spacer=L)) for L in (1.5, 1.75, 2.0, 2.25, 2.5)]
    rows = pf.run_geoms(cases, "spacer tolerance")
    rows.sort(key=lambda r: r["case"])
    pf.RES = RES
    pf.write("best_tolerance", rows)
    for r in rows:
        print(f"{r['case']:26s} misfit {r['ref_rms_pct']:.2f}  anode {r['anode_Br_pct']:+.2f}%  P165 {r['P165_W']:.2f} W  "
              f"I165 {r['I165_inner_A']:.2f}/{r['I165_outer_A']:.2f} A  sat165 {r['sat165']:.2f}  "
              f"pk {r['peak_offset_mm']:+.1f}  sigma {r['sigma_mm']:.1f}")
    (RES / "best.json").write_text(json.dumps(BEST, indent=1))
    figures([("shieldless, no spacer", X0),
             ("shieldless + 2 mm inert spacer behind the inner core, pole cap 2.5/4", BEST),
             ("with stock screens (buildable design)", json.loads((sd.RES / "optimum.json").read_text()))])


def short_core():
    """Inner core 2 mm shorter, insert in its place; nothing else changes (user, 2026-10-21)."""
    base = {k: v for k, v in BEST.items() if k not in ("spacer", "spacer_mode")}
    cases = [(f"shortened core, insert {g:g} mm", dict(base, core_gap=g)) for g in (1.5, 1.75, 2.0, 2.25, 2.5)]
    cases += [(f"back plate moved, insert {g:g} mm", dict(BEST, spacer=g)) for g in (1.5, 1.75, 2.0, 2.25, 2.5)]
    rows = pf.run_geoms(cases, "short core")
    order = {n: i for i, (n, _) in enumerate(cases)}
    rows.sort(key=lambda r: order[r["case"]])
    pf.RES = RES
    pf.write("short_core", rows)
    for r in rows:
        print(f"{r['case']:34s} misfit {r['ref_rms_pct']:.2f}  anode {r['anode_Br_pct']:+.2f}%  P165 {r['P165_W']:.2f} W  "
              f"I165 {r['I165_inner_A']:.2f}/{r['I165_outer_A']:.2f} A  sat165 {r['sat165']:.2f}  pk {r['peak_offset_mm']:+.1f}  "
              f"sigma {r['sigma_mm']:.1f}")


DISC = {**{k: v for k, v in BEST.items() if k not in ("spacer", "spacer_mode")}, "core_gap": 2.0, "coil_short": 1}


def disc():
    """Insert as wide as the pole cap (R 20.5) behind the inner core AND the inner coil; both 2 mm shorter
    (user, 2026-10-21). Back plate, posts, outer coils unchanged."""
    cases = [(f"disc insert {g:g} mm", dict(DISC, core_gap=g)) for g in (1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0)]
    cases += [(f"core-only insert {g:g} mm", {**DISC, "core_gap": g, "coil_short": 0}) for g in (1.75, 2.0)]
    rows = pf.run_geoms(cases, "disc insert")
    order = {n: i for i, (n, _) in enumerate(cases)}
    rows.sort(key=lambda r: order[r["case"]])
    pf.RES = RES
    pf.write("disc_insert", rows)
    for r in rows:
        print(f"{r['case']:26s} misfit {r['ref_rms_pct']:.2f}  anode {r['anode_Br_pct']:+.2f}%  P165 {r['P165_W']:.2f} W  "
              f"I165 {r['I165_inner_A']:.2f}/{r['I165_outer_A']:.2f} A  sat165 {r['sat165']:.2f}  pk {r['peak_offset_mm']:+.1f}  "
              f"sigma {r['sigma_mm']:.1f}")


def disc_plots():
    """Results for the disc insert (pole cap diameter, inner core + inner coil shortened)."""
    import optimise as opt
    import shaping
    from solve import solve_case
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    gaps = [1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0]
    base = {k: v for k, v in DISC.items() if k not in ("core_gap", "coil_short")}   # same pole cap, no insert
    cases = [("no insert", base)] + [(f"insert {g:g} mm", dict(DISC, core_gap=g)) for g in gaps]
    opt.solve_all([x for _, x in cases])
    out = [(lab, *ts.score(x), x) for lab, x in cases]
    cmap = plt.get_cmap("viridis")
    col = {lab: ("0.45" if lab == "no insert" else cmap(i / (len(gaps) - 1)))
           for i, (lab, *_rest) in enumerate(out[1:])} | {"no insert": "0.45"}
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    for ax in axs[:2]:
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="#888888", lw=6, alpha=0.45,
                label="reference (ideal)")
    for lab, m, (y, Br, By), x in out:
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        lw = 2.4 if lab == "insert 2 mm" else 1.3
        ls = "--" if lab == "no insert" else "-"
        axs[0].plot(1 - dd[sel] / ts.L_CH, Br[sel], color=col[lab], lw=lw, ls=ls,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.1f}%, anode {m['anode_Br_pct']:+.1f}%")
        axs[1].plot(1 - dd[sel] / ts.L_CH, np.abs(Br[sel]), color=col[lab], lw=lw, ls=ls)
    axs[0].set_ylim(-10, 185)
    axs[0].axhline(0, color="0.3", lw=0.6)
    axs[0].legend(frameon=False, fontsize=7.2, loc="upper left")
    axs[0].set_title("B_r along the channel (165 G peak)", loc="left", fontsize=10)
    axs[1].set_yscale("log")
    axs[1].set_ylim(0.3, 250)
    axs[1].set_title("|B_r|, log scale (anode end; dips = B_r crossing zero)", loc="left", fontsize=10)
    for ax in axs[:2]:
        ax.set_xlim(0, 1.08)
        ax.axvline(1, color="0.3", lw=0.7)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("channel-averaged B_r (G)")
    ax = axs[2]
    ms = [m for lab, m, *_ in out[1:]]
    ax.plot(gaps, [m["ref_rms_pct"] for m in ms], "o-", color=fd.PALETTE[0], label="misfit to reference (%)")
    ax.plot(gaps, [m["anode_Br_pct"] for m in ms], "s-", color=fd.PALETTE[1], label="B_r at anode (% of peak)")
    ax.axhline(1.8, color=fd.PALETTE[1], lw=0.8, ls=":")
    ax.text(gaps[0], 2.0, "reference anode 1.8 %", fontsize=7, color=fd.PALETTE[1])
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_xlabel("insert thickness (mm)")
    ax.set_ylabel("%")
    a2 = ax.twinx()
    a2.plot(gaps, [m["P165_W"] for m in ms], "^--", color=fd.PALETTE[3], label="coil power at 165 G (W)")
    a2.set_ylabel("coil power at 165 G (W)")
    a2.set_ylim(0, None)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=7.5, loc="upper right")
    m0 = out[0][1]
    ax.set_title(f"Against insert thickness (no insert: misfit {m0['ref_rms_pct']:.1f} %, anode "
                 f"{m0['anode_Br_pct']:.1f} %, {m0['P165_W']:.2f} W)", loc="left", fontsize=9)
    fig.suptitle("Shieldless design with the pole-cap-diameter insert (inner core and inner coil shortened); "
                 "pole cap edge face length 2.5, chamfer width 4", x=0.01, ha="left", fontsize=10)
    fig.savefig(FIG / "disc_results.png", dpi=140)
    plt.close(fig)
    # field maps: no insert vs 2 mm insert
    picks = [("shieldless, no insert (same pole cap)", base), ("shieldless, 2 mm disc insert", DISC)]
    opt.solve_all([x for _, x in picks], map_=True)
    items = []
    for lab, x in picks:
        m, _ = ts.score(x)
        # maps at the 165 G target (ts.score's drive is at 300 G)
        m = dict(m, NI_inner=m["NI_inner"] * ts.K, NI_outer=m["NI_outer"] * ts.K, P_W=m["P165_W"])
        res = solve_case(opt.key(x), opt.to_params(x), map_=True)
        items.append((lab, res, m, opt.to_params(x)))
    shaping.FIG = FIG
    shaping.map_fig(items, "disc_maps", ncol=2)


def mu_check():
    cases = [(f"all iron mu {m:g}", dict(BEST, mu=m)) for m in (200.0, 500.0, 1000.0, 1500.0, 3000.0, 5000.0)]
    rows = pf.run_geoms(cases, "spacer mu")
    rows.sort(key=lambda r: float(r["case"].split()[-1]))
    pf.RES = RES
    pf.write("best_mu", rows)
    for r in rows:
        print(f"{r['case']:20s} misfit {r['ref_rms_pct']:.2f}  anode {r['anode_Br_pct']:+.2f}%  P165 {r['P165_W']:.2f} W")


def insert_drawing(x=None):
    """Dimensioned half-section of the shieldless design with the inert insert behind the inner core,
    plus a zoomed detail of the insert. Names follow figures/glossary.png."""
    from matplotlib.patches import Rectangle
    from geometry import DEFAULTS, PLATE_T, Y_COIL_TOP, Y_FP_IN, ring_geometry
    x = dict(x or BEST)
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    short = bool(x.get("core_gap"))   # True: inner core shortened, insert in its place, nothing else moves
    L = x["core_gap"] if short else x["spacer"]
    z = fd.z
    if short:
        zb = z(Y_COIL_TOP)             # back plate inner face: unchanged
        zc = zb - L                    # rear end of the shortened inner core
    else:
        zc = z(Y_COIL_TOP)             # rear end of the inner core (and of both coils): unchanged
        zb = zc + L                    # back plate inner face, moved upstream by L
    rc, rh = DEFAULTS["core_r"], DEFAULTS["bp_hole_r"]
    disc = bool(x.get("coil_short"))  # insert as wide as the pole cap, under core and inner coil
    ri = x["cap_r"] if disc else rc   # insert radius
    r0, r1, _ = ring_geometry(DEFAULTS)
    INS = dict(fc="#f3c7c7", ec="#b03030", lw=1.4, hatch="////")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15, 7.6), gridspec_kw=dict(width_ratios=[1.55, 1]),
                                 constrained_layout=True)
    for a in (ax, bx):
        fd._section(a, x)
        a.add_patch(Rectangle((zc, 0), L, ri, **INS))
        a.grid(False)
        a.set_aspect("equal")
    # --- full view: pole cap dimensions and the insert callout
    R, t, w = x["cap_r"], x["cap_t"], x["cap_w"]
    zt = PLATE_T
    fd._ext(ax, [(-1, R), (-6.5, R)])
    fd._dim(ax, (-5.5, 0), (-5.5, R), f"pole cap R {R:g}", rot=90, tpos=0.45)
    fd._ext(ax, [(t, R + 0.3), (t, R + 5.5)])
    fd._ext(ax, [(0, R + 0.3), (0, R + 5.5)])
    ax.annotate("", xy=(0, R + 4.5), xytext=(-3.2, R + 4.5),
                arrowprops=dict(arrowstyle="-|>", color=fd.INK, lw=0.7, mutation_scale=7))
    ax.annotate("", xy=(t, R + 4.5), xytext=(t + 3.2, R + 4.5),
                arrowprops=dict(arrowstyle="-|>", color=fd.INK, lw=0.7, mutation_scale=7))
    ax.text(t + 3.5, R + 4.5, f"pole cap edge face length {t:g}", fontsize=7.5, va="center")
    fd._ext(ax, [(zt + 0.3, R - w), (zt + 3.5, R - w)])
    fd._ext(ax, [(t + 0.3, R), (zt + 3.5, R)])
    fd._dim(ax, (zt + 2.8, R - w), (zt + 2.8, R), f"chamfer width {w:g}", off=(0.6, 0), ha="left", fs=7.2)
    ax.text(PLATE_T / 2, DEFAULTS["fp_bore_r"] - 2.0, "flat front plate\n(no chamfer)", fontsize=7.2,
            ha="center", va="top", color=fd.CHG)
    ax.annotate(f"NON-MAGNETIC INSERT\n{L:g} mm thick, Ø {2 * ri:g} (see detail)", xy=(zc + L / 2, rc * 0.5),
                xytext=(zc - 24, -9), fontsize=8.5, color="#b03030", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="#b03030", lw=1.0, mutation_scale=9))
    ax.text((z(Y_FP_IN) + zc) / 2, (r0 + r1) / 2 + 3,
            "outer posts (unchanged)" if short else f"outer posts: {L:g} mm longer (iron)", fontsize=7.5,
            ha="center", color="#6e6e6e")
    ax.text(zb + PLATE_T / 2, 70, "back\nplate", fontsize=7, ha="center", va="center", color="#6e6e6e")
    ax.text(PLATE_T / 2, 70, "front\nplate", fontsize=7, ha="center", va="center", color="#6e6e6e")
    ax.text((z(Y_FP_IN) + zc) / 2, (rc + DEFAULTS["icoil_ro"]) / 2,
            f"inner coil ({L:g} mm shorter, same turns)" if disc else "inner coil (unchanged)", fontsize=7,
            ha="center", va="center", color="#8a5a00")
    ax.text((z(Y_FP_IN) + zc) / 2, rc / 2, f"inner core ({L:g} mm shorter)" if short else "inner core (unchanged)",
            fontsize=7, ha="center", va="center",
            color="#6e6e6e")
    ax.add_patch(Rectangle((zc - 6, -2), zb - zc + PLATE_T + 12, ri + 6, fill=False, ec="#b03030", lw=0.9, ls="--"))
    ax.text(zb + PLATE_T + 0.5, ri + 4.4, "detail", fontsize=7.5, color="#b03030")
    ax.set_xlim(-14, zb + PLATE_T + 8)
    ax.set_ylim(-12, 76)
    ax.set_xlabel("axial distance upstream of the exit plane (mm)")
    ax.set_ylabel("radius (mm)")
    ax.set_title("Shieldless design with the inner-core insert: r-z half-section (mm).\nNo inner or outer screen. "
                 "Blue = changed A36 parts; red hatched = non-magnetic insert."
                 + ((f"\nInner core and inner coil {L:g} mm shorter; back plate, outer posts and outer coils unchanged."
                     if disc else f"\nInner core {L:g} mm shorter; back plate, outer posts and coils unchanged.")
                    if short else ""),
                 loc="left", fontsize=9.5)
    # --- detail
    yd = max(DEFAULTS["icoil_ro"], ri) + 2.2        # dimension line above the coil / insert
    fd._ext(bx, [(zc, rc + 0.3), (zc, yd + 0.8)])
    fd._ext(bx, [(zb, rc + 0.3), (zb, yd + 0.8)])
    bx.annotate("", xy=(zc, yd), xytext=(zc - 3, yd),
                arrowprops=dict(arrowstyle="-|>", color=fd.INK, lw=0.8, mutation_scale=9))
    bx.annotate("", xy=(zb, yd), xytext=(zb + 3, yd),
                arrowprops=dict(arrowstyle="-|>", color=fd.INK, lw=0.8, mutation_scale=9))
    bx.text(zb + 3.3, yd, f"insert thickness {L:.2f} ± 0.10", fontsize=9, va="center")
    fd._ext(bx, [(zb + PLATE_T + 0.4, ri), (zb + PLATE_T + 7.5, ri)])
    fd._dim(bx, (zb + PLATE_T + 6.5, 0), (zb + PLATE_T + 6.5, ri),
            f"insert R {ri:g} (Ø {2 * ri:g})" + (" = pole cap" if disc else ""), rot=90, fs=8.5)
    fd._ext(bx, [(zb + PLATE_T + 0.4, rh), (zb + PLATE_T + 2.5, rh)])
    bx.text(zb + PLATE_T + 0.6, rh - 0.3, f"back plate\nhole R {rh:g}", fontsize=7, color="#6e6e6e", va="top")
    bx.text(zc + L / 2, -1.4, "insert", fontsize=8, ha="center", va="top", color="#b03030")
    bx.text(zc - 3.2, rc / 2, f"inner core\n(A36, {L:g} mm shorter)" if short else "inner core\n(A36, unchanged length)",
            fontsize=8, ha="right", va="center",
            color="#6e6e6e")
    bx.text(zb + PLATE_T / 2, (rh + rc) / 2 + 1.0, "back plate\n(unchanged)" if short else
            f"back plate\n(moved\n{L:g} mm\nupstream)", fontsize=7.5,
            ha="center", va="center", color="#6e6e6e")
    bx.text(zc - 3.2, (rc + DEFAULTS["icoil_ro"]) / 2, f"inner coil\n({L:g} mm shorter)" if disc else
            "inner coil\n(unchanged)", fontsize=8, ha="right",
            va="center", color="#8a5a00")
    bx.text(zc - 5, -4.0,
            f"Insert: disc Ø {2 * ri:g} x {L:g} mm, any non-magnetic material (aluminium, 300-series\n"
            f"stainless, titanium, ceramic); a centre hole is fine. Any bolt through it must also be\n"
            f"non-magnetic: a steel bolt bridges the gap and undoes it.\n"
            + (f"The inner coil is {L:g} mm shorter with the same ampere-turns; the outer coils are unchanged.\n"
               f"Tune with shims: B_r at the anode is ~0 at 1.6 mm ({L:g} mm: -1.9 %). " if disc else
               f"The inner coil is unchanged: its last {L:g} mm winds around the insert. Tune 1.75-2 mm\n"
               f"with shims (1.75 mm gives ~0 B_r at the anode). " if short else
               f"The {L:g} mm between the coil ends and the back plate is free space (coils unchanged).\n")
            + "Faces flat and parallel: extra air at either\nface adds to the gap "
              "(about 5 % of peak at the anode per mm of total gap).",
            fontsize=7.6, va="top", ha="left")
    bx.set_xlim(zc - 16, zb + PLATE_T + 9)
    bx.set_ylim(-12.5, max(DEFAULTS["icoil_ro"], ri) + 5)
    bx.set_xlabel("axial distance upstream of the exit plane (mm)")
    bx.set_title("Detail: the insert behind the inner core", loc="left", fontsize=9.5)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"insert_drawing{'_disc' if disc else '_short_core' if short else ''}.{ext}", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    if sys.argv[1:2] == ["sweep"]:
        sweep()
    elif sys.argv[1:2] == ["drawing"]:
        insert_drawing()
    elif sys.argv[1:2] == ["drawing_disc"]:
        insert_drawing(DISC)
    elif sys.argv[1:2] == ["drawing_short"]:
        insert_drawing({**{k: v for k, v in BEST.items() if k not in ("spacer", "spacer_mode")}, "core_gap": 2.0})
    elif sys.argv[1:2] == ["short"]:
        short_core()
    elif sys.argv[1:2] == ["disc"]:
        disc()
    elif sys.argv[1:2] == ["disc_plots"]:
        disc_plots()
    elif sys.argv[1:2] == ["mu"]:
        mu_check()
    elif sys.argv[1:2] == ["best"]:
        best()
    elif sys.argv[1:2] == ["small"]:
        small()
    elif sys.argv[1:2] == ["chamfer"]:
        chamfer(int(sys.argv[2]), float(sys.argv[3]))
