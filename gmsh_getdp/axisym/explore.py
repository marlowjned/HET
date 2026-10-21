"""
explore.py

Field-shaping studies on the axisymmetric model (see geometry.py), scored by
metrics.py: peak B_r on the exit plane, minimum field at the anode.

    python explore.py ratio          # inner/outer ratio on the baseline
    python explore.py <study>        # see STUDIES

Each study writes figures/<study>_profiles.png (B_r along the axis, peak
scaled to 300 G), figures/<study>_fieldlines.png and results/<study>.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from geometry import CHANNEL, Y_EXIT  # noqa: E402
from metrics import evaluate  # noqa: E402
from solve import solve_case  # noqa: E402

FIG = HERE / "figures"
RES = HERE / "results"
Q_SPEC = 223.4 / 352.0
Q_SCAN = np.round(np.arange(0.0, 2.001, 0.05), 3)
COLS = ["case", "q", "NI_inner", "NI_outer", "peak_offset_mm", "exit_Br_pct", "anode_Br_pct",
        "anode_B_pct", "half_depth_mm", "sigma_mm", "gauss_rms_pct", "iron_B_max_T"]


def best_q(res, key="anode_B_pct"):
    """Ratio minimising key (default: |B| at the anode, % of peak)."""
    ms = [evaluate(res, q)[0] for q in Q_SCAN]
    return min(ms, key=lambda m: m[key])["q"]


def run_cases(cases, map_=True):
    """cases: {name: params}. Returns {name: res}."""
    return {n: solve_case(n, p, map_=map_) for n, p in cases.items()}


def table(rows, study):
    RES.mkdir(exist_ok=True)
    with open(RES / f"{study}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    hdr = (f"{'case':28s} {'q':>5s} {'NI_in':>6s} {'NI_out':>6s} {'pk_off':>6s} {'exit%':>6s} "
           f"{'anBr%':>6s} {'an|B|%':>6s} {'half':>5s} {'sigma':>5s} {'gRMS%':>5s} {'Bfe':>5s}")
    print(hdr + ("    P_W" if "P_W" in rows[0] else ""))
    for r in rows:
        print(f"{r['case']:28s} {r['q']:5.2f} {r['NI_inner']:6.0f} {r['NI_outer']:6.0f} "
              f"{r['peak_offset_mm']:6.1f} {r['exit_Br_pct']:6.1f} {r['anode_Br_pct']:6.1f} "
              f"{r['anode_B_pct']:6.1f} {r['half_depth_mm']:5.1f} {r['sigma_mm']:5.1f} "
              f"{r['gauss_rms_pct']:5.1f} {r['iron_B_max_T']:5.2f}" + (f" {r['P_W']:6.1f}" if "P_W" in r else ""))


def profile_fig(curves, study, title):
    """curves: list of (label, (y, Br, By), metrics)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    for i, (lab, (y, Br, By), m) in enumerate(curves):
        d = y - Y_EXIT   # mm from the exit, + = into the channel
        c = f"C{i % 10}"
        axs[0].plot(d, Br, color=c, label=f"{lab}: anode {m['anode_Br_pct']:.0f}%, "
                                          f"peak {m['peak_offset_mm']:+.1f} mm")
        axs[1].plot(d, np.hypot(Br, By), color=c, label=lab)
        axs[1].plot(d, np.abs(By), color=c, ls=":", lw=1)
    L = CHANNEL["y_anode"] - Y_EXIT
    s = np.linspace(0, L, 100)
    for ax in axs:
        ax.axvline(0, color="k", lw=0.8)
        ax.axvline(L, color="k", lw=0.8)
        ax.axvspan(0, L, color="0.93", zorder=0)
        ax.set_xlabel("distance from exit plane (mm): <0 downstream (plume), >0 into channel")
        ax.grid(alpha=0.3)
    axs[0].plot(s, 300 * np.exp(-s ** 2 / (2 * 8.0 ** 2)), "k--", lw=1, label="half-Gaussian, sigma 8 mm (reference)")
    axs[0].set_ylabel("B_r, channel-averaged (G), peak scaled to 300 G")
    axs[0].set_title(f"{title}: radial field")
    axs[0].legend(fontsize=7)
    axs[1].set_ylabel("|B| (solid), |B_axial| (dotted) (G)")
    axs[1].set_title("field magnitude and axial component")
    axs[1].text(L, 5, "anode ", ha="right", fontsize=8)
    axs[1].text(0, 5, " exit", fontsize=8)
    fig.savefig(FIG / f"{study}_profiles.png", dpi=130)
    plt.close(fig)


def fieldline_fig(items, study, extra_draw=None):
    """items: list of (label, res, NI_in, NI_out, params). Field lines are
    contours of r*a_phi; colour is |B| (log)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle
    import matplotlib.tri as mtri
    n = len(items)
    fig, axs = plt.subplots(1, n, figsize=(5.2 * n, 6.2), constrained_layout=True, squeeze=False)
    for ax, item in zip(axs[0], items):
        lab, res, ni_in, ni_out, params = item[:5]
        ni = {"inner": ni_in, "outer": ni_out, "trim": item[5] if len(item) > 5 else 0.0}
        bi = res["map"]["inner"][0]
        x, y = bi[:, 0] * 1e3, bi[:, 1] * 1e3
        Bx = sum(ni[k] * b[:, 2] for k, (b, _) in res["map"].items())
        By = sum(ni[k] * b[:, 3] for k, (b, _) in res["map"].items())
        B = np.hypot(Bx, By) * 1e4
        A = sum(ni[k] * r[:, 2] for k, (_, r) in res["map"].items())
        sel = (x < 100) & (y > -85) & (y < 35)
        tri = mtri.Triangulation(x[sel], y[sel])
        im = ax.tripcolor(tri, np.clip(B[sel], 5, None), cmap="magma", norm=LogNorm(vmin=10, vmax=10000),
                          shading="gouraud")
        ax.tricontour(tri, A[sel], levels=np.linspace(A[sel].min(), A[sel].max(), 60), colors="w",
                      linewidths=0.5)
        if params.get("trim"):
            r0, r1, y0, y1 = params["trim"]
            ax.add_patch(Rectangle((r0, y0), r1 - r0, y1 - y0, fill=False, ec="lime", lw=1.2))
        ax.add_patch(Rectangle((CHANNEL["r_in"], Y_EXIT), CHANNEL["r_out"] - CHANNEL["r_in"],
                               CHANNEL["y_anode"] - Y_EXIT, fill=False, ec="cyan", lw=1.4, ls="--"))
        ax.set_xlim(0, 95)
        ax.set_ylim(-75, 30)
        ax.set_aspect("equal")
        ax.set_title(lab, fontsize=9)
        ax.set_xlabel("r (mm)")
        ax.set_ylabel("y (mm), exit at -37.8, plume toward -y")
    fig.colorbar(im, ax=axs[0].tolist(), label="|B| (G)", shrink=0.8)
    fig.savefig(FIG / f"{study}_fieldlines.png", dpi=120)
    plt.close(fig)


# --------------------------------------------------------------------------
def study_ratio():
    res = solve_case("baseline", {}, map_=True)
    rows, curves = [], []
    for q in (0.0, 0.3, Q_SPEC, 1.0, 2.0):
        m, prof = evaluate(res, q)
        m["case"] = f"baseline q={q:.2f}"
        rows.append(m)
        curves.append((f"q={q:.2f}" + (" (spec)" if q == Q_SPEC else ""), prof, m))
    qb = best_q(res)
    m, prof = evaluate(res, qb)
    m["case"] = f"baseline q={qb:.2f} (min anode |B|)"
    rows.append(m)
    table(rows, "ratio")
    profile_fig(curves, "ratio", "Baseline, outer/inner A-turn ratio q")
    fieldline_fig([(f"baseline, q={q:.2f}", res, rows[i]["NI_inner"], rows[i]["NI_outer"], {})
                   for i, q in ((0, 0.0), (2, Q_SPEC), (4, 2.0))], "ratio")


# --------------------------------------------------------------------------
# Geometry helpers (mm). Keep-out: the channel annulus plus a 1.5 mm wall each
# side, r 20.62..33.10, anode to exit. CAD (chamber, anode) is ignored.
from geometry import DEFAULTS, Y_FP_IN, Y_COIL_TOP  # noqa: E402

KEEP_R = (CHANNEL["r_in"] - 1.5, CHANNEL["r_out"] + 1.5)


def pole_cap(t_tip=None, w=3.0, ext=0.0, r=None):
    """Inner pole cap polygon. t_tip: thickness of the rim face toward the
    channel (back-side chamfer of radial width w); ext: rim extended ext mm
    downstream of the exit plane; r: cap radius."""
    R = r or DEFAULTS["pole_cap_r"]
    T = DEFAULTS["pole_cap_t"]
    t = T if t_tip is None else t_tip
    y0 = Y_EXIT - ext
    return [(0, Y_EXIT), (R - w, Y_EXIT), (R - w, y0), (R, y0), (R, Y_EXIT + t),
            (R - w, Y_EXIT + T), (0, Y_EXIT + T)] if (t < T or ext > 0) else \
        [(0, Y_EXIT), (R, Y_EXIT), (R, Y_EXIT + T), (0, Y_EXIT + T)]


def front_plate(t_tip=None, w=6.0, ext=0.0, bore=None):
    B = bore or DEFAULTS["fp_bore_r"]
    P = DEFAULTS["plate_r"]
    T = Y_FP_IN - Y_EXIT
    t = T if t_tip is None else t_tip
    y0 = Y_EXIT - ext
    return [(B, y0), (B + w, y0), (B + w, Y_EXIT), (P, Y_EXIT), (P, Y_FP_IN), (B + w, Y_FP_IN),
            (B, Y_EXIT + t)]


def inner_screen(gap, t=3.0, r0=15.0):
    """Tube from the back plate down to `gap` mm short of the pole cap's back
    face (y = -28.32), r r0..r0+t (coil OD 14.29, keep-out from 20.62)."""
    return [(r0, Y_COIL_TOP), (r0 + t, Y_COIL_TOP), (r0 + t, Y_FP_IN + gap), (r0, Y_FP_IN + gap)]


def outer_screen(gap, t=3.0, r0=40.0):
    """Tube from the back plate down to `gap` mm short of the front plate's
    inner face, r r0..r0+t (outer winding sheet starts at 49.3)."""
    return [(r0, Y_COIL_TOP), (r0 + t, Y_COIL_TOP), (r0 + t, Y_FP_IN + gap), (r0, Y_FP_IN + gap)]


def sweep(study, cases, title, map_names=(), q_mode="best"):
    """Solve cases {name: params}; score each at the anode-|B|-optimal ratio
    (q_mode="best") or the spec ratio; profile + field-line figures."""
    rows, curves, fl = [], [], []
    for name, params in cases.items():
        res = solve_case(name, params, map_=name in map_names)
        q = best_q(res) if q_mode == "best" else Q_SPEC
        m, prof = evaluate(res, q)
        m["case"] = name
        rows.append(m)
        curves.append((name, prof, m))
        if name in map_names:
            fl.append((f"{name} (q={q:.2f})", res, m["NI_inner"], m["NI_outer"], params))
    table(rows, study)
    profile_fig(curves, study, title)
    if fl:
        fieldline_fig(fl, study)
    return rows


def study_poles():
    cases = {"baseline": {}}
    for t in (6.0, 4.0, 2.5):
        cases[f"chamfer t{t:g}"] = dict(pole_cap_poly=pole_cap(t_tip=t), front_plate_poly=front_plate(t_tip=t))
    for e in (2.0, 4.0):
        cases[f"chamfer t4 ext{e:g}"] = dict(pole_cap_poly=pole_cap(t_tip=4, ext=e),
                                             front_plate_poly=front_plate(t_tip=4, ext=e))
    cases["narrow gap"] = dict(pole_cap_poly=pole_cap(r=KEEP_R[0]), front_plate_poly=front_plate(bore=KEEP_R[1]))
    cases["narrow gap chamfer t4"] = dict(pole_cap_poly=pole_cap(t_tip=4, r=KEEP_R[0], w=4),
                                          front_plate_poly=front_plate(t_tip=4, bore=KEEP_R[1]))
    sweep("poles", cases, "Pole-tip shaping", map_names=("baseline", "chamfer t2.5", "narrow gap chamfer t4"))


def study_screens():
    cases = {"baseline": {}}
    for g in (3.0, 6.0, 10.0, 15.0):
        cases[f"screens gap{g:g}"] = dict(extra_iron=[inner_screen(g), outer_screen(g)])
    cases["inner screen only gap6"] = dict(extra_iron=[inner_screen(6.0)])
    cases["outer screen only gap6"] = dict(extra_iron=[outer_screen(6.0)])
    sweep("screens", cases, "Magnetic screens (tubes on the back plate)",
          map_names=("screens gap6", "screens gap15"))


def study_noplates():
    """Front/back plates non-magnetic (user, 2026-10-15: they may be 316L,
    which conflicts with MATERIALS.md's A36). Compares the field at the spec
    excitation and the A-turns needed for 300 G, with and without plates."""
    spec = (352.0, 223.4)
    cases = {"baseline": {}, "no plates": dict(plates=False),
             "no plates (thin ring)": dict(plates=False, ring_mode="thin"),
             "front plate only": dict(back_plate=False),
             "back plate only": dict(front_plate=False)}
    rows, curves, fl = [], [], []
    from metrics import combine, shape_metrics
    for name, params in cases.items():
        res = solve_case(name, params, map_=name != "no plates (thin ring)")
        y, Br, By = combine(res, *spec)
        ye = Y_EXIT
        sg = np.sign(np.interp(ye, y, Br)) or 1
        pk = np.flatnonzero((y >= ye - 15) & (y <= CHANNEL["y_anode"]))
        print(f"{name:24s} at spec 352/223.4 A-turns: B_r exit {sg * np.interp(ye, y, Br):6.1f} G, "
              f"peak {np.max(sg * Br[pk]):6.1f} G")
        for q in (Q_SPEC, best_q(res)):
            m, prof = evaluate(res, q)
            m["case"] = f"{name} q={q:.2f}"
            rows.append(m)
        curves.append((name, prof, m))
        if name in ("baseline", "no plates", "front plate only", "back plate only"):
            fl.append((f"{name}, spec ratio", res, rows[-2]["NI_inner"], rows[-2]["NI_outer"], params))
    table(rows, "noplates")
    profile_fig(curves, "noplates", "Front/back plates iron vs non-magnetic")
    fieldline_fig(fl, "noplates")


def coil_power(ni_in, ni_out):
    """Hot coil power (W) on the spec windings (EMAG_SPEC.md): inner 230 t,
    0.643 ohm; outer 4 x 184 t, 0.494 ohm each."""
    return (ni_in / 230) ** 2 * 0.643 + 4 * (ni_out / 184) ** 2 * 0.494


def study_ratio_plates():
    """Inner/outer ratio sweep with and without the front plate (back plate
    absent in both, as if it is 316L), iron-plate baseline for reference."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cases = {"baseline": {}, "front plate only": dict(back_plate=False), "no plates": dict(plates=False)}
    qs = np.round(np.arange(0.0, 3.001, 0.05), 3)
    show = (0.0, 0.3, Q_SPEC, 1.0, 1.5, 2.0, 3.0)
    rows = []
    fig, axs = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    for j, (name, params) in enumerate(cases.items()):
        res = solve_case(name, params)
        ms = []
        for q in qs:
            m, _ = evaluate(res, q)
            m["P_W"] = coil_power(m["NI_inner"], m["NI_outer"])
            ms.append(m)
        for q in show:
            m, prof = evaluate(res, q)
            m["P_W"] = coil_power(m["NI_inner"], m["NI_outer"])
            m["case"] = f"{name} q={q:.2f}"
            rows.append(m)
            y, Br, By = prof
            axs[0, j].plot(y - Y_EXIT, Br, label=f"q={q:.2f}: anode {m['anode_Br_pct']:.0f}%, "
                                                  f"pk {m['peak_offset_mm']:+.1f} mm, {m['P_W']:.0f} W")
        L = CHANNEL["y_anode"] - Y_EXIT
        a = axs[0, j]
        a.axvspan(0, L, color="0.93", zorder=0)
        a.axhline(0, color="k", lw=0.6)
        a.set_title(f"{name}: B_r, peak scaled to 300 G")
        a.set_xlabel("distance from exit (mm), >0 into channel")
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
        b = axs[1, j]
        b.plot(qs, [m["anode_Br_pct"] for m in ms], label="anode B_r % of peak")
        b.plot(qs, [m["anode_B_pct"] for m in ms], label="anode |B| % of peak")
        b.plot(qs, [m["exit_Br_pct"] for m in ms], label="exit B_r % of peak")
        b.axhline(0, color="k", lw=0.6)
        b.set_xlabel("outer/inner A-turn ratio q (outer per coil)")
        b.set_ylabel("%")
        b2 = b.twinx()
        b2.plot(qs, [m["P_W"] for m in ms], "k--", label="coil power for 300 G (W)")
        b2.set_ylabel("coil power for 300 G peak (W), spec windings")
        b.grid(alpha=0.3)
        b.legend(fontsize=7, loc="upper left")
        b2.legend(fontsize=7, loc="upper right")
        b.set_title(f"{name}: ratio sweep")
    fig.savefig(FIG / "ratio_plates.png", dpi=120)
    plt.close(fig)
    COLS.append("P_W")
    table(rows, "ratio_plates")
    for name in cases:
        sub = [r for r in rows if r["case"].startswith(name + " q")]
        print(f"{name}: min power {min(r['P_W'] for r in sub):.1f} W over the shown ratios")


def study_front_best():
    """Front plate only (back plate non-magnetic): field maps and channel
    profiles at the minimum-power ratio and the minimum-anode-|B| ratio."""
    res = solve_case("front plate only", dict(back_plate=False), map_=True)
    ms = []
    for q in Q_SCAN:
        m, prof = evaluate(res, q)
        m["P_W"] = coil_power(m["NI_inner"], m["NI_outer"])
        ms.append((m, prof))
    picks = [("min power", min(ms, key=lambda t: t[0]["P_W"])),
             ("min anode |B|", min(ms, key=lambda t: t[0]["anode_B_pct"]))]
    rows, curves, fl = [], [], []
    for lab, (m, prof) in picks:
        name = f"front only, {lab}, q={m['q']:.2f}, {m['P_W']:.0f} W"
        m["case"] = f"front only {lab}"
        rows.append(m)
        curves.append((name, prof, m))
        fl.append((f"{name}\n{m['NI_inner']:.0f} / {m['NI_outer']:.0f} A-turns for 300 G peak",
                   res, m["NI_inner"], m["NI_outer"], {}))
    base = solve_case("baseline", {}, map_=True)
    mb, pb = evaluate(base, 0.2)
    mb["P_W"] = coil_power(mb["NI_inner"], mb["NI_outer"])
    mb["case"] = "both plates iron, q=0.20"
    rows.append(mb)
    curves.append((f"both plates iron (reference), q=0.20, {mb['P_W']:.1f} W", pb, mb))
    fl.append((f"both plates iron (reference), q=0.20\n{mb['NI_inner']:.0f} / {mb['NI_outer']:.0f} A-turns",
               base, mb["NI_inner"], mb["NI_outer"], {}))
    if "P_W" not in COLS:
        COLS.append("P_W")
    table(rows, "front_best")
    profile_fig(curves, "front_best", "Front plate only")
    fieldline_fig(fl, "front_best")


def study_taper():
    """Tapered front plate (user request 2026-10-15): back face of the bore
    edge tapered from tip thickness t at the bore to the full 9.525 mm at
    bore + w; the front face stays flat on the exit plane. Both plates iron,
    and front plate only. Scored at each case's minimum-power ratio."""
    if "P_W" not in COLS:
        COLS.append("P_W")
    configs = {"both plates": {}, "front only": dict(back_plate=False)}
    geoms = {"no taper": None}
    for w in (6.0, 15.0):
        for t in (6.0, 3.0, 1.5):
            geoms[f"taper t{t:g} w{w:g}"] = (t, w)
    rows, fl = [], []
    for cname, cp in configs.items():
        curves = []
        for gname, tw in geoms.items():
            name = {"both plates": "baseline", "front only": "front plate only"}[cname] if tw is None \
                else f"{cname} {gname}"
            params = dict(cp) if tw is None else dict(cp, front_plate_poly=front_plate(t_tip=tw[0], w=tw[1]))
            show = gname in ("no taper", "taper t1.5 w15")
            res = solve_case(name, params, map_=show)
            ms = []
            for q in Q_SCAN:
                m, prof = evaluate(res, q)
                m["P_W"] = coil_power(m["NI_inner"], m["NI_outer"])
                ms.append((m, prof))
            m, prof = min(ms, key=lambda t_: t_[0]["P_W"])
            m["case"] = f"{cname}, {gname}"
            rows.append(m)
            curves.append((f"{gname}, q={m['q']:.2f}, {m['P_W']:.1f} W", prof, m))
            if show:
                fl.append((f"{cname}, {gname}\nq={m['q']:.2f}, {m['P_W']:.1f} W", res, m["NI_inner"],
                           m["NI_outer"], params))
        profile_fig(curves, f"taper_{cname.replace(' ', '_')}", f"Tapered front plate, {cname}")
    table(rows, "taper")
    fieldline_fig(fl, "taper")


STUDIES = {"taper": study_taper, "front_best": study_front_best, "ratio_plates": study_ratio_plates, "noplates": study_noplates, "ratio": study_ratio, "poles": study_poles, "screens": study_screens}

if __name__ == "__main__":
    STUDIES[sys.argv[1]]()
