"""
shaping.py

Field-shaping studies, 2026-10-18. Target: channel B_r a half-Gaussian with
its peak on the exit plane and as little field as possible at the anode.
Plates and cores are A36 (iron), the anode is 316L (air).

Theory (Ampere loop in the r-y plane from the inner core to the outer posts,
closed through the back plate): the magnetic potential difference across the
channel at axial station y is the ampere-turns located BEHIND y (toward the
back plate). With the coils filling the cores, ~43% of the turns sit behind
the anode plane, which is the ~20% anode floor. Levers, simplest first:

  coils     S1  axial placement of the windings (turns forward of the anode)
  rear      S2  core / post to back-plate air gaps and back-plate thickness
  front     S2b front plate (+ pole cap) thickness
  screens   S3  iron tubes on the back plate between coils and channel
  combo     S4  best of the above together

All studies use the best pole tips found on 2026-10-15 (tapered front plate
t1.5 w15 + r20 t1.5 w6 pole cap) unless named "flat". Each case is scored at
its minimum-power outer/inner ratio, peak scaled to 300 G.

    python shaping.py theory|coils|rear|front|screens|combo|summary|all
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from explore import Q_SCAN, front_plate, pole_cap  # noqa: E402
from geometry import (CHANNEL, DEFAULTS, R_IRON, R_RING, R_TRIM, R_WIN, R_WOI, R_WOO,  # noqa: E402
                      Y_COIL_TOP, Y_EXIT, Y_FP_IN, layout)
from metrics import combine, evaluate  # noqa: E402
from solve import solve_case  # noqa: E402

FIG = HERE / "figures" / "shaping"
RES = HERE / "results" / "shaping"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

POLES = dict(front_plate_poly=front_plate(t_tip=1.5, w=15.0), pole_cap_poly=pole_cap(r=20.0, t_tip=1.5, w=6.0))
L_CORE = Y_COIL_TOP - Y_FP_IN                      # 41.275 mm, front plate to back plate
B_IN0 = DEFAULTS["icoil_ro"] - DEFAULTS["core_r"]  # 4.7625 mm spec inner build
B_OUT0 = DEFAULTS["ocoil_build"]                   # 3.969 mm spec outer build
F_ANODE = (CHANNEL["y_anode"] - Y_FP_IN) / L_CORE  # 0.569: coil from the front plate to the anode plane
R_IN_MAX = 19.5                                    # inner coil OD limit (channel keep-out 20.62)

# spec windings, hot (EMAG_SPEC.md): inner 230 t 0.643 ohm, outer 184 t 0.494 ohm each
K_IN, K_OUT = 0.643 / 230 ** 2, 0.494 / 184 ** 2


# --------------------------------------------------------------------------
# geometry helpers (mm)
def coils(f_in=1.0, f_out=1.0, rear=False):
    """Coils occupying fraction f of the core length, against the front plate
    (or the back plate if rear), with the build scaled to keep the winding
    area, i.e. the same turns of the same wire."""
    def win(f):
        return (Y_COIL_TOP - f * L_CORE, Y_COIL_TOP) if rear else (Y_FP_IN, Y_FP_IN + f * L_CORE)
    p = {}
    if f_in != 1.0:
        p.update(icoil_y=win(f_in), icoil_ro=DEFAULTS["core_r"] + B_IN0 / f_in)
    if f_out != 1.0:
        p.update(ocoil_y=win(f_out), ocoil_build=B_OUT0 / f_out)
    return p


def tube(r0, t, y_tip):
    """Iron tube r0..r0+t from the back plate forward to y_tip."""
    return [(r0, Y_COIL_TOP), (r0 + t, Y_COIL_TOP), (r0 + t, y_tip), (r0, y_tip)]


def inner_screen(gap, ro_coil, t=1.5):
    """Screen just outside the inner coil, tip `gap` mm behind the pole cap."""
    r0 = ro_coil + 0.4
    assert r0 + t <= 20.2, "inner screen hits the channel keep-out"
    return tube(r0, t, Y_FP_IN + gap)


def outer_screen(gap, r0=40.0, t=3.0):
    """Screen between the channel and the outer coils (inner sheet from 49.3),
    tip `gap` mm behind the front plate's inner face."""
    return tube(r0, t, Y_FP_IN + gap)


def coil_power(params, ni_in, ni_out):
    """Hot coil power (W): spec windings scaled by mean turn length and
    inverse window area (same wire and packing)."""
    p = dict(DEFAULTS, **params)
    shapes = dict((reg, poly) for reg, poly in layout(p) if reg in (R_WIN, R_WOI))
    (r0, y0), (r1, _), (_, y1), _ = shapes[R_WIN]
    a_in = (r1 - r0) * (y1 - y0)
    a_in0 = B_IN0 * L_CORE
    p_in = ni_in ** 2 * K_IN * ((r0 + r1) / (DEFAULTS["core_r"] + DEFAULTS["icoil_ro"])) * (a_in0 / a_in)
    (_, y0), _, (_, y1), _ = shapes[R_WOI]
    b, rc = p["ocoil_build"], p["ocore_r"]
    p_out = 4 * ni_out ** 2 * K_OUT * ((rc + b / 2) / (rc + B_OUT0 / 2)) * (B_OUT0 * L_CORE / (b * (y1 - y0)))
    return p_in + p_out


def score(name, params):
    """Solve and score at the minimum-power ratio. Returns (metrics, profile, res)."""
    res = solve_case(name, params, map_=True)
    best = None
    for q in Q_SCAN:
        m, prof = evaluate(res, q)
        m["P_W"] = coil_power(params, m["NI_inner"], m["NI_outer"])
        if best is None or m["P_W"] < best[0]["P_W"]:
            best = (m, prof)
    m, prof = best
    m["case"] = name
    return m, prof, res


# --------------------------------------------------------------------------
# figures
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6})
    return plt


def profile_fig(curves, fname, title, ref=None):
    """curves: [(label, (y, Br, By), m)], at most 8. ref: optional (label, prof)
    drawn in grey as the reference."""
    plt = _mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    L = CHANNEL["y_anode"] - Y_EXIT
    if ref:
        y, Br, By = ref[1]
        axs[0].plot(y - Y_EXIT, Br, color="0.55", lw=2.5, alpha=0.6, label=ref[0])
        axs[1].plot(y - Y_EXIT, np.hypot(Br, By), color="0.55", lw=2.5, alpha=0.6)
    for i, (lab, (y, Br, By), m) in enumerate(curves[:8]):
        c = PALETTE[i]
        axs[0].plot(y - Y_EXIT, Br, color=c, lw=1.8,
                    label=f"{lab}: anode {m['anode_Br_pct']:+.0f}%, pk {m['peak_offset_mm']:+.1f} mm, "
                          f"{m['P_W']:.1f} W")
        axs[1].plot(y - Y_EXIT, np.hypot(Br, By), color=c, lw=1.8)
    for ax in axs:
        ax.axvspan(0, L, color="0.95", zorder=0)
        ax.axvline(0, color="0.3", lw=0.8)
        ax.axvline(L, color="0.3", lw=0.8)
        ax.axhline(0, color="0.3", lw=0.6)
        ax.set_xlabel("distance from exit plane (mm), <0 plume, >0 into channel")
        ax.set_xlim(-30, L + 3)
    axs[0].text(0.5, 3, "exit", fontsize=8, color="0.3")
    axs[0].text(L - 0.5, 3, "anode", ha="right", fontsize=8, color="0.3")
    axs[0].set_ylabel("B_r, channel-averaged (G), peak = 300 G")
    axs[0].set_title(f"{title}: radial field", loc="left")
    axs[0].legend(fontsize=7, frameon=False)
    axs[1].set_ylabel("|B| (G)")
    axs[1].set_title("field magnitude", loc="left")
    fig.savefig(FIG / f"{fname}_profiles.png", dpi=130)
    plt.close(fig)


def map_fig(items, fname, ncol=None):
    """items: [(label, res, m, params)]. |B| (log colour) and field lines
    (contours of r a_phi), with the iron and windings outlined."""
    plt = _mpl()
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Polygon
    import matplotlib.tri as mtri
    FIG.mkdir(parents=True, exist_ok=True)
    n = len(items)
    ncol = ncol or min(n, 4)
    nrow = -(-n // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 4.4 * nrow), constrained_layout=True, squeeze=False)
    for ax in axs.flat[n:]:
        ax.set_visible(False)
    for ax, (lab, res, m, params) in zip(axs.flat, items):
        ni = {"inner": m["NI_inner"], "outer": m["NI_outer"]}
        bi = res["map"]["inner"][0]
        x, y = bi[:, 0] * 1e3, bi[:, 1] * 1e3
        Bx = sum(ni[k] * res["map"][k][0][:, 2] for k in ni)
        By = sum(ni[k] * res["map"][k][0][:, 3] for k in ni)
        A = sum(ni[k] * res["map"][k][1][:, 2] for k in ni)
        B = np.hypot(Bx, By) * 1e4
        sel = (x < 85) & (y > -62) & (y < 30)
        tri = mtri.Triangulation(x[sel], y[sel])
        im = ax.tripcolor(tri, np.clip(B[sel], 5, None), cmap="magma", norm=LogNorm(vmin=10, vmax=10000),
                          shading="gouraud", rasterized=True)
        a = A[sel]
        ax.tricontour(tri, a, levels=np.linspace(np.percentile(a, 0.5), np.percentile(a, 99.5), 40),
                      colors="w", linewidths=0.45, alpha=0.8)
        for reg, poly in layout(params):
            ec = {R_IRON: "#c3c2b7", R_RING: "#c3c2b7", R_TRIM: "lime"}.get(reg, "#eda100")
            ax.add_patch(Polygon(poly, closed=True, fill=False, ec=ec, lw=0.9))
        ax.add_patch(Polygon([(CHANNEL["r_in"], Y_EXIT), (CHANNEL["r_out"], Y_EXIT),
                              (CHANNEL["r_out"], CHANNEL["y_anode"]), (CHANNEL["r_in"], CHANNEL["y_anode"])],
                             closed=True, fill=False, ec="cyan", lw=1.2, ls="--"))
        ax.set_xlim(0, 85)
        ax.set_ylim(-62, 30)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(f"{lab}\nanode {m['anode_Br_pct']:+.0f}%, pk {m['peak_offset_mm']:+.1f} mm, "
                     f"{m['P_W']:.1f} W", fontsize=8, loc="left")
        ax.set_xlabel("r (mm)")
        ax.set_ylabel("y (mm), plume toward -y")
    fig.colorbar(im, ax=axs.ravel().tolist(), label="|B| (G)", shrink=0.7)
    fig.savefig(FIG / f"{fname}_maps.png", dpi=125)
    plt.close(fig)


COLS = ["study", "case", "q", "NI_inner", "NI_outer", "P_W", "peak_offset_mm", "exit_Br_pct", "anode_Br_pct",
        "anode_B_pct", "half_depth_mm", "sigma_mm", "gauss_rms_pct", "iron_B_max_T"]


def write(study, rows):
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / f"{study}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 3) if isinstance(v, (float, np.floating)) else v)
                        for k, v in dict(r, study=study).items()})
    print(f"\n{study}")
    print(f"{'case':40s} {'q':>5s} {'P_W':>5s} {'pk':>5s} {'exit%':>6s} {'anBr%':>6s} {'an|B|%':>6s} "
          f"{'half':>5s} {'sig':>5s} {'gRMS':>5s} {'Bfe':>5s}")
    for r in rows:
        print(f"{r['case']:40s} {r['q']:5.2f} {r['P_W']:5.1f} {r['peak_offset_mm']:5.1f} {r['exit_Br_pct']:6.1f} "
              f"{r['anode_Br_pct']:6.1f} {r['anode_B_pct']:6.1f} {r['half_depth_mm']:5.1f} {r['sigma_mm']:5.1f} "
              f"{r['gauss_rms_pct']:5.1f} {r['iron_B_max_T']:5.2f}")


def run(study, cases, title, ref=None, maps=None):
    """cases: {name: params}. ref: name of the case drawn in grey; maps: names
    for the field-map figure (default all)."""
    out = {n: score(n, p) for n, p in cases.items()}
    rows = [m for m, _, _ in out.values()]
    write(study, rows)
    curves = [(n, prof, m) for n, (m, prof, _) in out.items() if n != ref]
    refc = (f"{ref} ({out[ref][0]['anode_Br_pct']:+.0f}% anode)", out[ref][1]) if ref else None
    for k in range(0, len(curves), 8):
        profile_fig(curves[k:k + 8], f"{study}{'' if k == 0 else k // 8 + 1}", title, refc)
    names = maps or list(cases)
    map_fig([(n, out[n][2], out[n][0], cases[n]) for n in names], study)
    return out


# --------------------------------------------------------------------------
# studies
REF = "poles (2026-10-15 best)"


def study_theory():
    """Cumulative-winding prediction vs FEM B_r for coils filling the core,
    pushed forward and pushed back."""
    plt = _mpl()
    cases = {"full-length coils": dict(POLES), "coils front 57%": dict(POLES, **coils(F_ANODE, F_ANODE)),
             "coils rear 57%": dict(POLES, **coils(F_ANODE, F_ANODE, rear=True))}
    fig, ax = plt.subplots(figsize=(7.5, 4.4), constrained_layout=True)
    names = {"full-length coils": REF}
    for i, (lab, p) in enumerate(cases.items()):
        res = solve_case(names.get(lab, f"theory {lab}"), p, map_=True)
        y, Br, _ = combine(res, 1.0, 0.2)
        Br = Br / np.interp(Y_EXIT, y, Br)
        pp = dict(DEFAULTS, **p)
        y0, y1 = pp["icoil_y"] or (Y_FP_IN, Y_COIL_TOP)
        behind = np.clip((y1 - np.clip(y, y0, y1)) / (y1 - y0), 0, 1)
        d = y - Y_EXIT
        ax.plot(d, Br, color=PALETTE[i], lw=1.8, label=f"{lab}: FEM")
        ax.plot(d, behind, color=PALETTE[i], lw=1.2, ls="--", label=f"{lab}: turns behind y (theory)")
    L = CHANNEL["y_anode"] - Y_EXIT
    ax.axvspan(0, L, color="0.95", zorder=0)
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_xlim(-10, L + 3)
    ax.set_xlabel("distance from exit plane (mm), >0 into channel")
    ax.set_ylabel("B_r / B_r(exit)   and   fraction of turns behind y")
    ax.set_title("Anode floor = ampere-turns behind the station (q = 0.2)", loc="left")
    ax.legend(fontsize=7, frameon=False)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "theory.png", dpi=130)
    plt.close(fig)


def study_coils():
    c = {REF: dict(POLES)}
    for f in (0.85, 0.70, F_ANODE, 0.48):
        c[f"coils front {f:.0%}"] = dict(POLES, **coils(f, f))
    c["coils rear 57% (check)"] = dict(POLES, **coils(F_ANODE, F_ANODE, rear=True))
    c["inner front 57% only"] = dict(POLES, **coils(F_ANODE, 1.0))
    c["outer front 57% only"] = dict(POLES, **coils(1.0, F_ANODE))
    c["inner 48%, outer 35%"] = dict(POLES, **coils(0.48, 0.35))
    c["flat poles"] = {}
    c["flat poles, coils front 57%"] = coils(F_ANODE, F_ANODE)
    return run("coils", c, "S1 coil axial placement", ref=REF,
               maps=[REF, "coils front 70%", "coils front 57%", "coils front 48%", "coils rear 57% (check)",
                     "inner 48%, outer 35%", "flat poles", "flat poles, coils front 57%"])


def study_rear():
    c = {REF: dict(POLES)}
    for g in (0.25, 0.5, 1.0, 2.0):
        c[f"core gap {g:g} mm"] = dict(POLES, core_gap=g)
    for g in (0.5, 1.0, 2.0):
        c[f"post gap {g:g} mm"] = dict(POLES, ring_gap=g)
    for t in (3.175, 6.35, 15.875):
        c[f"back plate {t:.2f} mm"] = dict(POLES, bp_t=t)
    return run("rear", c, "S2 rear air gaps and back-plate thickness", ref=REF,
               maps=[REF, "core gap 0.5 mm", "core gap 2 mm", "post gap 1 mm", "back plate 3.17 mm"])


def study_front():
    c = {"flat poles": {}}
    for t in (3.175, 4.7625, 6.35, 12.7):
        c[f"flat, front plate + cap {t:.2f} mm"] = dict(fp_t=t, pole_cap_t=t)
    return run("front", c, "S2b front plate + pole cap thickness (flat poles)", ref="flat poles",
               maps=["flat poles", "flat, front plate + cap 3.17 mm", "flat, front plate + cap 12.70 mm"])


def study_screens():
    ro = DEFAULTS["icoil_ro"]
    c = {REF: dict(POLES)}
    for g in (3.0, 6.0, 10.0, 15.0):
        c[f"screens gap {g:g}"] = dict(POLES, extra_iron=[inner_screen(g, ro), outer_screen(g)])
    c["inner screen only gap 6"] = dict(POLES, extra_iron=[inner_screen(6.0, ro)])
    c["outer screen only gap 6"] = dict(POLES, extra_iron=[outer_screen(6.0)])
    c["screens gap 6, outer r35"] = dict(POLES, extra_iron=[inner_screen(6.0, ro), outer_screen(6.0, r0=35.0)])
    c["screens gap 6, outer r45"] = dict(POLES, extra_iron=[inner_screen(6.0, ro), outer_screen(6.0, r0=45.0)])
    return run("screens", c, "S3 magnetic screens", ref=REF,
               maps=[REF, "screens gap 3", "screens gap 6", "screens gap 15", "inner screen only gap 6",
                     "outer screen only gap 6", "screens gap 6, outer r35", "screens gap 6, outer r45"])


# S4, chosen from S1-S3: the inner screen and the inner coil's placement are
# the levers that matter; outer screen, post gaps and back-plate thickness do
# nothing; a core gap works but costs power. Screens 3 mm thick for flux margin.
_RO = DEFAULTS["icoil_ro"]
_RO57 = DEFAULTS["core_r"] + B_IN0 / F_ANODE
_RO85 = DEFAULTS["core_r"] + B_IN0 / 0.85
_LIP = pole_cap(r=20.0, t_tip=1.5, w=6.0, ext=2.0)
COMBO = {
    "inner screen t3 gap 6": dict(POLES, extra_iron=[inner_screen(6.0, _RO, t=3.0)]),
    "inner screen t3 gap 3": dict(POLES, extra_iron=[inner_screen(3.0, _RO, t=3.0)]),
    "inner screen t3 gap 6 + outer screen": dict(POLES, extra_iron=[inner_screen(6.0, _RO, t=3.0),
                                                                     outer_screen(6.0)]),
    "inner coil 85% + inner screen t2.5 gap 6": dict(POLES, **coils(0.85, 1.0),
                                                       extra_iron=[inner_screen(6.0, _RO85, t=2.5)]),
    "inner coil 57% + inner screen t1.5 gap 6": dict(POLES, **coils(F_ANODE, 1.0),
                                                       extra_iron=[inner_screen(6.0, _RO57, t=1.5)]),
    "inner screen t3 gap 6 + 2 mm cap lip": dict(POLES, pole_cap_poly=_LIP,
                                                  extra_iron=[inner_screen(6.0, _RO, t=3.0)]),
    "both screens gap 6 + 2 mm cap lip": dict(POLES, pole_cap_poly=_LIP,
                                              extra_iron=[inner_screen(6.0, _RO, t=3.0), outer_screen(6.0)]),
    "inner coil 57% + screen + lip": dict(POLES, **coils(F_ANODE, 1.0), pole_cap_poly=_LIP,
                                          extra_iron=[inner_screen(6.0, _RO57, t=1.5)]),
}


def study_combo():
    c = {REF: dict(POLES)}
    c.update(COMBO)
    out = run("combo", c, "S4 combinations", ref=REF)
    screen_flux(out, c)
    return out


def screen_flux(out, cases):
    """Bulk |B| in the inner screen (median and 95th percentile over its
    elements, joints excluded) at each case's excitation, to check that the
    linear-iron result is not hiding saturation."""
    print(f"\n{'case':44s} {'screen B50':>10s} {'B95 (T)':>8s}")
    for n, (m, _, res) in out.items():
        polys = [p for p in cases[n].get("extra_iron", []) if p[0][0] < 21]
        if not polys:
            continue
        r0, r1, y1 = polys[0][0][0], polys[0][1][0], polys[0][2][1]
        bi, bo = res["iron"]["inner"], res["iron"]["outer"]
        x, y = bi[:, 0] * 1e3, bi[:, 1] * 1e3
        sel = (x > r0) & (x < r1) & (y > y1 + 1) & (y < Y_COIL_TOP - 1.5)
        B = np.linalg.norm(m["NI_inner"] * bi[sel, 2:4] + m["NI_outer"] * bo[sel, 2:4], axis=1)
        print(f"{n:44s} {np.median(B):10.2f} {np.percentile(B, 95):8.2f}")


def summary():
    """One figure of the key metrics for every case of every study."""
    plt = _mpl()
    studies = ["coils", "rear", "front", "screens", "combo"]
    rows = []
    for s in studies:
        f = RES / f"{s}.csv"
        if f.exists():
            rows += [r for r in csv.DictReader(open(f)) if not (s != "coils" and r["case"] == REF)]
    keys = [("anode_Br_pct", "B_r at anode (% of peak)"), ("peak_offset_mm", "peak position (mm, + = past exit)"),
            ("sigma_mm", "half-Gaussian sigma (mm)"), ("gauss_rms_pct", "Gaussian misfit RMS (% of peak)"),
            ("P_W", "coil power for 300 G (W)"), ("iron_B_max_T", "peak iron |B| (T)")]
    n = len(rows)
    fig, axs = plt.subplots(1, len(keys), figsize=(18, 0.24 * n + 1.6), sharey=True, constrained_layout=True)
    yy = np.arange(n)[::-1]
    col = {s: PALETTE[i] for i, s in enumerate(studies)}
    for ax, (k, lab) in zip(axs, keys):
        v = np.array([float(r[k]) for r in rows])
        ax.barh(yy, v, color=[col[r["study"]] for r in rows], height=0.7)
        ax.axvline(0, color="0.3", lw=0.6)
        ax.set_title(lab, fontsize=9, loc="left")
        ax.grid(axis="y", visible=False)
    axs[0].set_yticks(yy, [f"{r['study']}: {r['case']}" for r in rows], fontsize=7)
    fig.savefig(FIG / "summary.png", dpi=120)
    plt.close(fig)
    # power vs anode floor scatter, one colour per study
    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    for s in studies:
        sub = [r for r in rows if r["study"] == s]
        if sub:
            ax.scatter([float(r["P_W"]) for r in sub], [float(r["anode_Br_pct"]) for r in sub], s=40,
                       color=col[s], edgecolor="white", linewidth=1, label=s, zorder=3)
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_xlabel("coil power for 300 G peak (W)")
    ax.set_ylabel("B_r at anode (% of peak)")
    ax.set_title("Anode floor vs power, every case", loc="left")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "floor_vs_power.png", dpi=130)
    plt.close(fig)


def headline():
    """The as-built magnet, the 2026-10-15 pole tips, and the best case of
    each lever, on one profile plot; as-built vs winner field maps."""
    from metrics import TARGET_G  # noqa: F401
    picks = [("as built (flat poles)", "flat poles", {}),
             ("+ pole tips (10-15)", REF, dict(POLES)),
             ("+ inner coil front 48%", "coils front 48%", dict(POLES, **coils(0.48, 0.48))),
             ("+ 2 mm core gap", "core gap 2 mm", dict(POLES, core_gap=2.0)),
             ("+ inner screen", "inner screen t3 gap 6", COMBO["inner screen t3 gap 6"]),
             ("+ both screens + cap lip", "both screens gap 6 + 2 mm cap lip",
              COMBO["both screens gap 6 + 2 mm cap lip"])]
    out = [(lab, *score(n, p), p) for lab, n, p in picks]
    profile_fig([(lab, prof, m) for lab, m, prof, _, _ in out], "headline", "As built to best")
    map_fig([(lab, res, m, p) for lab, m, _, res, p in (out[0], out[1], out[5])], "headline", ncol=3)


STUDIES = {"headline": headline, "theory": study_theory, "coils": study_coils, "rear": study_rear, "front": study_front,
           "screens": study_screens, "combo": study_combo, "summary": summary}

if __name__ == "__main__":
    for s in (["theory", "coils", "rear", "front", "screens"] if sys.argv[1] == "all" else sys.argv[1:]):
        STUDIES[s]()
