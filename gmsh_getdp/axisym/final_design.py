"""
final_design.py

Final field-shaping geometry (2026-10-21) and the sigma family.

FINAL is the 10-20 optimum (results/optimise/optimum.json) with the outer
screen wall cut from 7 mm (optimiser drift) to 4 mm (user, 2026-10-21). Sigma
is still open, so the inner-screen tip gap -- the sigma knob -- is swept and a
family of arrangements is picked at round sigma values. Every arrangement is
scored with optimise.evaluate_x (spec ratio, 300 G peak, saturation per part).

    python final_design.py sweep|family|drawing|figures|all

Outputs: results/final/*.csv|json, figures/final/*.png|svg.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
from geometry import CHANNEL, DEFAULTS, Y_COIL_TOP, Y_EXIT, Y_FP_IN, PLATE_T, layout, ring_geometry  # noqa: E402
from geometry import R_AIR, R_IRON, R_RING  # noqa: E402
from metrics import combine, shape_metrics  # noqa: E402
from shaping import _mpl, PALETTE  # noqa: E402
from solve import solve_case  # noqa: E402

FIG = HERE / "figures" / "final"
RES = HERE / "results" / "final"
FINAL = dict(json.loads((opt.RES / "optimum.json").read_text()), os_t=4.0)
GAPS = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0]
SIGMAS = [7.5, 8.0, 8.5, 9.0, 9.5]       # family targets (mm)
COLS = ["case", "is_gap", "sigma_mm", "exit_Br_pct", "peak_offset_mm", "anode_Br_pct", "anode_B_pct",
        "half_depth_mm", "gauss_rms_pct", "plume_B_10mm_pct", "P_W", "I_inner_A", "I_outer_A",
        "sat_worst", "sat_worst_part", "inner_core_Bflux_max_T", "outer_posts_Bflux_max_T",
        "inner_screen_Bp95_T", "outer_screen_Bp95_T", "objective"]


def write(name, rows):
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / f"{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


def show(rows):
    print(f"{'gap':>5s} {'sig':>5s} {'exit%':>6s} {'pk':>5s} {'an|B|%':>6s} {'gRMS':>5s} {'P_W':>5s} "
          f"{'sat':>5s} {'coreB':>5s}")
    for r in rows:
        print(f"{r['is_gap']:5.1f} {r['sigma_mm']:5.2f} {r['exit_Br_pct']:6.1f} {r['peak_offset_mm']:5.1f} "
              f"{r['anode_B_pct']:6.2f} {r['gauss_rms_pct']:5.2f} {r['P_W']:5.2f} {r['sat_worst']:5.2f} "
              f"{r['inner_core_Bflux_max_T']:5.2f}")


def run(xs, tag):
    opt.solve_all(xs)
    out = []
    for x in xs:
        m, prof, _ = opt.evaluate_x(x)
        m["case"] = f"{tag} gap {x['is_gap']:g}"
        out.append((m, prof))
    return out


def sweep():
    xs = [dict(FINAL, is_gap=g) for g in GAPS]
    rows = [m for m, _ in run(xs, "sweep")]
    write("gap_sweep", rows)
    show(rows)
    return rows


def family():
    """Gaps at the target sigmas (interpolated from the sweep, rounded to
    0.5 mm, the fabrication resolution that matters here), solved exactly."""
    rows = sweep()
    g = np.array([r["is_gap"] for r in rows])
    s = np.array([r["sigma_mm"] for r in rows])
    o = np.argsort(s)
    picks = sorted({float(np.round(np.interp(t, s[o], g[o]) * 2) / 2) for t in SIGMAS})
    out = run([dict(FINAL, is_gap=p) for p in picks], "family")
    fam = [m for m, _ in out]
    write("family", fam)
    print("\nfamily")
    show(fam)
    (RES / "final.json").write_text(json.dumps(dict(FINAL, family_gaps=picks), indent=1))
    return out


def as_built():
    """The as-built circuit (flat 9.525 mm poles, r 17.46 cap, no screens) at the spec ratio."""
    res = solve_case("final/as_built", {})
    y, Br, By = combine(res, 1.0, opt.Q_SPEC)
    m, s = shape_metrics(y, Br, By)
    return m, (y, Br * s, By * s)


# --------------------------------------------------------------------------
# drawing
INK, GREY, CHG = "#1f2328", "#b9b9b3", PALETTE[0]


def z(y):
    """Axial coordinate of the drawings: mm upstream of the exit plane."""
    return np.asarray(y) - Y_EXIT


def _poly(ax, pts, **kw):
    from matplotlib.patches import Polygon
    P = np.array(pts, float)
    ax.add_patch(Polygon(np.c_[z(P[:, 1]), P[:, 0]], closed=True, **kw))


def _section(ax, x, highlight=True, coils=True):
    """Half-section (axial horizontal, r vertical). Unchanged parts grey,
    changed parts in the accent colour."""
    p = opt.to_params(x)
    changed = [p["front_plate_poly"], p["pole_cap_poly"]] + p["extra_iron"]
    for reg, poly in layout(p):
        if any(poly is c for c in changed):
            continue
        if reg == R_AIR:               # non-magnetic insert
            _poly(ax, poly, fc="#f3c7c7", ec="#b03030", lw=1.2, hatch="////")
        elif reg in (R_IRON, R_RING):
            _poly(ax, poly, fc="#e4e3dd", ec="#8f8e88", lw=0.8)
        elif coils:
            _poly(ax, poly, fc="#f6dfb3", ec="#c98a1a", lw=0.6, hatch="xxx")
    for c in changed:
        _poly(ax, c, fc=("#cfe1f7" if highlight else "#e4e3dd"), ec=(CHG if highlight else "#8f8e88"),
              lw=1.4 if highlight else 0.8)
    _poly(ax, [(CHANNEL["r_in"], Y_EXIT), (CHANNEL["r_out"], Y_EXIT), (CHANNEL["r_out"], CHANNEL["y_anode"]),
               (CHANNEL["r_in"], CHANNEL["y_anode"])], fill=False, ec="#2a9d8f", lw=1.0, ls="--")
    ax.axhline(0, color=INK, lw=0.6, ls=(0, (8, 3, 2, 3)))


def _dim(ax, p0, p1, text, off=(0, 0), tpos=0.5, ha="center", va="center", fs=7.5, rot=0, color=INK):
    """Dimension line with arrows between p0 and p1 (drawing coords)."""
    (x0, y0), (x1, y1) = p0, p1
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="<|-|>", color=color, lw=0.7, shrinkA=0, shrinkB=0, mutation_scale=7))
    tx, ty = x0 + tpos * (x1 - x0) + off[0], y0 + tpos * (y1 - y0) + off[1]
    ax.text(tx, ty, text, fontsize=fs, ha=ha, va=va, rotation=rot, color=color,
            bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))


def _ext(ax, pts, **kw):
    P = np.array(pts, float)
    ax.plot(P[:, 0], P[:, 1], color="#6e6e6e", lw=0.45, **kw)


def drawing(x=None, fname="dimensioned", note="Inner-screen gap shown at 6.5 (nominal); see sigma family."):
    plt = _mpl()
    x = dict(x or FINAL)
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 8.3), constrained_layout=True)
    _section(ax, x)
    zfp, zbp = z(Y_FP_IN), z(Y_COIL_TOP)            # front-plate inner face / back-plate inner face
    zt = z(Y_EXIT + PLATE_T)                         # pole-cap back face plane (== zfp)
    B = DEFAULTS["fp_bore_r"]
    is0 = x.get("is_r0", opt.IS_R0)
    is1 = is0 + x["is_t"]
    os0, os1 = x["os_r"], x["os_r"] + x["os_t"]
    ztip_i, ztip_o = z(Y_FP_IN + x["is_gap"]), z(Y_FP_IN + x["os_gap"])

    # pole cap: radius, rim, taper width
    _ext(ax, [(-1, x["cap_r"]), (-6.5, x["cap_r"])])
    _dim(ax, (-5.5, 0), (-5.5, x["cap_r"]), f"pole cap R {x['cap_r']:g}", rot=90, tpos=0.45)
    if x["cap_t"] < PLATE_T - 1e-6:
        _ext(ax, [(x["cap_t"], x["cap_r"] + 0.3), (x["cap_t"], x["cap_r"] + 5.5)])
        _ext(ax, [(0, x["cap_r"] + 0.3), (0, x["cap_r"] + 5.5)])
        ax.annotate("", xy=(0, x["cap_r"] + 4.5), xytext=(-3.2, x["cap_r"] + 4.5),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=7))
        ax.annotate("", xy=(x["cap_t"], x["cap_r"] + 4.5), xytext=(x["cap_t"] + 3.2, x["cap_r"] + 4.5),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=7))
        ax.text(x["cap_t"] + 3.5, x["cap_r"] + 4.5, f"pole cap edge face length {x['cap_t']:g}", fontsize=7.5, va="center")
        # taper width: radial, chamfer runs from the rim to R cap_r - cap_w on the cap back face
        _ext(ax, [(zt - 0.3, x["cap_r"] - x["cap_w"]), (zt - 3.4, x["cap_r"] - x["cap_w"])])
        _ext(ax, [(x["cap_t"] + 0.3, x["cap_r"]), (zt - 2.6, x["cap_r"])])
        _dim(ax, (zt - 2.8, x["cap_r"] - x["cap_w"]), (zt - 2.8, x["cap_r"]), f"chamfer width {x['cap_w']:g}",
             off=(-1.7, -0.6), rot=90, fs=7)
    else:
        ax.text(PLATE_T / 2, x["cap_r"] * 0.62, f"flat cap\n{PLATE_T:g} thick", fontsize=7.2, ha="center",
                va="center", rotation=90, color=CHG)

    # front plate: bore face length, chamfer width (dims below the bore, in the plume side)
    if x["fp_t"] < PLATE_T - 1e-6:
        _ext(ax, [(x["fp_t"], B - 0.3), (x["fp_t"], B - 5.0)])
        _ext(ax, [(0, B - 0.3), (0, B - 5.0)])
        ax.annotate("", xy=(0, B - 4.0), xytext=(-3.2, B - 4.0),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=7))
        ax.annotate("", xy=(x["fp_t"], B - 4.0), xytext=(x["fp_t"] + 3.2, B - 4.0),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.7, mutation_scale=7))
        ax.text(-3.5, B - 4.0, f"bore face length {x['fp_t']:g}", fontsize=7.5, va="center", ha="right")
        _ext(ax, [(-0.5, B), (-6.5, B)])
        _ext(ax, [(-0.5, B + x["fp_w"]), (-6.5, B + x["fp_w"])])
        _dim(ax, (-5.5, B), (-5.5, B + x["fp_w"]), f"chamfer width {x['fp_w']:g}", rot=90,
             off=(0, 0 if x["fp_w"] >= 8 else x["fp_w"] / 2 + 6))
    else:
        ax.text(PLATE_T / 2, B - 2.0, "flat front plate\n(no chamfer)", fontsize=7.2, ha="center", va="top",
                color=CHG)
    _dim(ax, (-9.5, 0), (-9.5, B), f"front plate bore R {B:g}", rot=90, color="#6e6e6e", fs=7)

    # inner screen: radii, wall, tip gap + length (chained)
    zr = zbp + PLATE_T
    if not x.get("no_is"):
        _ext(ax, [(zr + 0.5, is0), (zr + 7, is0)])
        _ext(ax, [(zr + 0.5, is1), (zr + 7, is1)])
        _dim(ax, (zr + 3, 0), (zr + 3, is0), f"R {is0:.2f}", rot=90, tpos=0.55)
        _dim(ax, (zr + 6, 0), (zr + 6, is1), f"R {is1:.2f}", rot=90, tpos=0.55)
        ax.text(zr + 13.5, (is0 + is1) / 2, f"inner screen\nwall {x['is_t']:g}\n(clears inner coil OD "
                f"{DEFAULTS['icoil_ro']:.2f} by {is0 - DEFAULTS['icoil_ro']:.2f})", fontsize=7.2, va="center", color=CHG)
        yd = is1 + 1.5
        _ext(ax, [(ztip_i, is1 + 0.3), (ztip_i, yd + 0.8)])
        _ext(ax, [(zt, x["cap_r"] - x["cap_w"] + 0.3), (zt, yd + 0.8)])
        _dim(ax, (zt, yd), (ztip_i, yd), f"inner screen tip gap {x['is_gap']:g}", off=(0, 1.6), fs=7.2)
        _dim(ax, (ztip_i, yd), (zbp, yd), f"{zbp - ztip_i:.2f}", off=(0, 1.6), fs=7, color="#6e6e6e")

    # outer screen: radii, wall, tip gap + length (chained)
    if not x.get("no_os"):
        _ext(ax, [(zr + 0.5, os0), (zr + 10, os0)])
        _ext(ax, [(zr + 0.5, os1), (zr + 13, os1)])
        _dim(ax, (zr + 9, 0), (zr + 9, os0), f"R {os0:g}", rot=90, tpos=0.78)
        _dim(ax, (zr + 12, 0), (zr + 12, os1), f"R {os1:g}", rot=90, tpos=0.78)
        ax.text(zr + 13.5, (os0 + os1) / 2, f"outer screen\nwall {x['os_t']:g}", fontsize=7.2,
                va="center", color=CHG)
        yd = os0 - 1.5
        _ext(ax, [(ztip_o, os0 - 0.3), (ztip_o, yd - 0.8)])
        _ext(ax, [(zfp, os0 - 0.3), (zfp, yd - 0.8)])
        _dim(ax, (zfp, yd), (ztip_o, yd), f"outer screen tip gap {x['os_gap']:g}", off=(0, -1.7), fs=7.2)
        _dim(ax, (ztip_o, yd), (zbp, yd), f"{zbp - ztip_o:.2f}", off=(0, -1.7), fs=7, color="#6e6e6e")

    # labels
    ax.text(0.2, -3.0, "exit plane", fontsize=7, rotation=90, va="top", ha="left", color="#6e6e6e")
    ax.axvline(0, color="#6e6e6e", lw=0.5, ls=":")
    ax.text((z(CHANNEL["y_anode"])) / 2, (CHANNEL["r_in"] + CHANNEL["r_out"]) / 2, "channel",
            color="#2a9d8f", fontsize=8, ha="center", va="center")
    ax.text(z(CHANNEL["y_anode"]) + 0.3, CHANNEL["r_out"] + 0.6, "anode", color="#2a9d8f", fontsize=7)
    r0, r1, _ = ring_geometry(DEFAULTS)
    ax.text((zfp + zbp) / 2, DEFAULTS["ring_R"], "outer posts (4x, smeared ring)", fontsize=7,
            ha="center", va="center", color="#6e6e6e")
    ax.text((zfp + zbp) / 2, (DEFAULTS["core_r"] + DEFAULTS["icoil_ro"]) / 2, "inner coil", fontsize=7,
            ha="center", va="center", color="#8a5a00")
    ax.text((zfp + zbp) / 2, DEFAULTS["core_r"] / 2, "inner core", fontsize=7, ha="center", va="center",
            color="#6e6e6e")
    ax.text(zbp + PLATE_T / 2, 70, "back\nplate", fontsize=7, ha="center", va="center", color="#6e6e6e")
    ax.text(PLATE_T / 2, 70, "front\nplate", fontsize=7, ha="center", va="center", color="#6e6e6e")
    ax.text(-1.0, 1.0, "axis", fontsize=7, ha="right", color=INK)

    ax.set_xlim(-14, zbp + PLATE_T + 30)
    ax.set_ylim(-1.5, 76)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("axial distance upstream of the exit plane (mm)")
    ax.set_ylabel("radius (mm)")
    ax.set_title("BPL-1000 field-shaping geometry, r-z half-section (mm). Names as in figures/glossary.png. Blue = changed A36 parts, "
                 "dimensioned; grey = unchanged.\n" + note + "\n"
                 "Tip gaps are measured to the plane y = -28.32 (back face of the pole cap / upstream face of the front plate).",
                 loc="left", fontsize=9)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"{fname}.{ext}", dpi=160)
    plt.close(fig)


def family_fig(out, base):
    """Array: one zoomed half-section per sigma (gap dimensioned), then all the
    channel profiles together, then a table strip."""
    plt = _mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    n = len(out)
    fig = plt.figure(figsize=(3.1 * n, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(2, n, height_ratios=[1.0, 1.05])
    cols = [PALETTE[i] for i in range(n)]
    zfp = z(Y_FP_IN)
    for i, (m, _) in enumerate(out):
        x = {k: m[k] for k in opt.BASE}
        ax = fig.add_subplot(gs[0, i])
        _section(ax, x, coils=True)
        is1 = opt.IS_R0 + x["is_t"]
        ztip = z(Y_FP_IN + x["is_gap"])
        _ext(ax, [(ztip, is1 + 0.3), (ztip, is1 + 3.2)])
        _ext(ax, [(zfp, is1 + 0.3), (zfp, is1 + 3.2)])
        _dim(ax, (zfp, is1 + 2.6), (ztip, is1 + 2.6), f"{x['is_gap']:g}", off=(0, 1.6), fs=8, color=cols[i])
        ax.set_xlim(-3, 34)
        ax.set_ylim(0, 48)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(f"sigma {m['sigma_mm']:.1f} mm\ninner-screen gap {x['is_gap']:g} mm", fontsize=9,
                     loc="left", color=cols[i])
        if i == 0:
            ax.set_ylabel("r (mm)")
        ax.set_xlabel("mm upstream of exit")
    ax = fig.add_subplot(gs[1, : max(2, n - 2)])
    ya = CHANNEL["y_anode"]
    for (m, (y, Br, By)), c in zip(out, cols):
        d = z(y)
        ax.plot(d, Br, color=c, lw=1.8, label=f"gap {m['is_gap']:g}: sigma {m['sigma_mm']:.1f}")
        dd = np.linspace(0, z(ya), 100)
        ax.plot(dd, 300 * m["exit_Br_pct"] / 100 * np.exp(-dd ** 2 / (2 * m["sigma_mm"] ** 2)), color=c, lw=0.8,
                ls=":")
    m0, (y, Br, By) = base
    ax.plot(z(y), Br, color="0.55", lw=1.3, ls="--", label=f"as built (sigma {m0['sigma_mm']:.0f})")
    ax.axvline(0, color="0.3", lw=0.7)
    ax.axvline(z(ya), color="#2a9d8f", lw=0.8, ls="--")
    ax.text(z(ya) - 0.3, 290, "anode", ha="right", fontsize=8, color="#2a9d8f")
    ax.text(0.3, 290, "exit", fontsize=8)
    ax.set_xlim(-12, z(ya) + 4)
    ax.set_xlabel("mm upstream of the exit plane (negative = plume)")
    ax.set_ylabel("channel-averaged B_r (G), 300 G peak")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Channel profiles (solid) and half-Gaussian fits (dotted)", loc="left", fontsize=9)
    ax = fig.add_subplot(gs[1, max(2, n - 2):])
    ax.axis("off")
    rows = [["gap", "sigma", "exit %", "anode |B| %", "misfit %", "P (W)", "core B (T)"]]
    rows += [[f"{m['is_gap']:g}", f"{m['sigma_mm']:.1f}", f"{m['exit_Br_pct']:.1f}", f"{m['anode_B_pct']:.2f}",
              f"{m['gauss_rms_pct']:.1f}", f"{m['P_W']:.2f}", f"{m['inner_core_Bflux_max_T']:.2f}"]
             for m, _ in out]
    t = ax.table(cellText=rows, loc="center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1, 1.5)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#d0d0d0")
        if r == 0:
            cell.set_text_props(weight="bold")
        elif c == 0:
            cell.set_text_props(color=cols[r - 1], weight="bold")
    fig.suptitle("Sigma family: final geometry, only the inner-screen tip gap changes "
                 "(2-D axisym, spec ratio, 300 G peak)", x=0.01, ha="left", fontsize=10)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"sigma_family.{ext}", dpi=150)
    plt.close(fig)


def sweep_fig(rows, fam):
    plt = _mpl()
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.0), constrained_layout=True)
    g = [r["is_gap"] for r in rows]
    for ax, (k, lab) in zip(axs, [("sigma_mm", "sigma (mm)"), ("anode_B_pct", "|B| at anode (% of peak)"),
                                  ("exit_Br_pct", "B_r at exit (% of peak)"), ("P_W", "coil power (W)")]):
        ax.plot(g, [r[k] for r in rows], "o-", color=PALETTE[0], ms=3.5, lw=1.4)
        for i, m in enumerate(fam):
            ax.plot(m["is_gap"], m[k], "o", color=PALETTE[i], ms=7, mfc="none", mew=1.6)
        ax.set_xlabel("inner-screen tip gap (mm)")
        ax.set_ylabel(lab)
    fig.suptitle("Inner-screen gap sweep at the final geometry (circles = sigma family)", x=0.01, ha="left",
                 fontsize=10)
    fig.savefig(FIG / "gap_sweep.png", dpi=140)
    plt.close(fig)


def figures():
    FIG.mkdir(parents=True, exist_ok=True)
    out = family()
    rows = list(csv.DictReader(open(RES / "gap_sweep.csv")))
    rows = [{k: (float(v) if k not in ("case", "sat_worst_part") else v) for k, v in r.items()} for r in rows]
    sweep_fig(rows, [m for m, _ in out])
    base = as_built()
    print(f"\nas built: {base[0]}")
    (RES / "as_built.json").write_text(json.dumps(base[0], indent=1))
    family_fig(out, base)
    drawing()


if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in what:
        if w == "sweep":
            sweep()
        if w == "family":
            family()
        if w == "drawing":
            drawing()
        if w in ("figures", "all"):
            figures()
