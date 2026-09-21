"""
make_report.py

Builds the BPL-700 magnet-circuit report: figures (report/figures/),
generated LaTeX tables and numbers (report/generated/), then compiles
report/bpl700_magnet_report.tex with Tectonic (tools/tectonic/, gitignored).

Inputs, all produced elsewhere in the repo:
  gmsh_getdp/assembly/emag_sweep_results.csv    per-geometry magnetics (emag_sweep.py)
  gmsh_getdp/assembly/sweep_runs/*/extract.npz  channel profiles + r-y slices
  sizing/winding_design.py                      A-turns -> winding, current, voltage, power, temperature
  report/nonlinear_checks.json                  nonlinear confirmation runs (optional)

Usage (from repo root):  python report/make_report.py [--no-pdf]
"""
import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle, Polygon
from scipy.interpolate import griddata

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"
FIG = REPORT / "figures"
GEN = REPORT / "generated"
ASSEMBLY = ROOT / "gmsh_getdp" / "assembly"
sys.path.insert(0, str(ROOT / "sizing"))
sys.path.insert(0, str(ASSEMBLY))

from winding_design import WIRE, design, windows_for  # noqa: E402
from generate_regions import AISI1008_B, AISI1008_H, FNAL1020_B, FNAL1020_H  # noqa: E402

MU0 = 4e-7 * math.pi
IN = 25.4  # mm

# --- style (dataviz reference palette; slots 1-3 validate all-pairs) ------
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8983", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
MARKERS = ["o", "s", "^"]
SEQ = ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
SEQ_ORD = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue", SEQ)
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.6, "lines.markersize": 5,
    "savefig.bbox": "tight", "savefig.dpi": 300,
})

CORE_LABEL = {1.0: "core 1.000 in", 0.875: "core 0.875 in", 0.75: "core 0.750 in"}
STOCK_SPACERS = {0.0, 0.1875, 0.25, 0.3125, 0.375, 0.5}   # A36 sheet (SendCutSend), and 0 = no part
KNEE_LO, KNEE_HI = 1.3, 1.5

# Fixed geometry (mm), measured from the Onshape STEP exports.
Y_TOP_OUT, Y_TOP_IN = -32.51, -22.99
Y_ANODE_FACE = 5.59
R_POLE, R_BORE_TOP, R_PLATE = 17.46, 34.92, 76.2
R_CH_IN, R_CH_OUT = 22.098, 31.75
R_CHAMBER_OUT, R_CHAMBER_EXIT_OUT = 41.28, 34.93
R_OUTER_AXIS = 62.86
PLATE_T = 9.525


def load_rows():
    rows = []
    with open(ASSEMBLY / "emag_sweep_results.csv") as f:
        for r in csv.DictReader(f):
            for k, v in r.items():
                try:
                    r[k] = float(v)
                except ValueError:
                    r[k] = {"True": True, "False": False}.get(v, v)
            rows.append(r)
    return rows


def main_grid(rows):
    return [r for r in rows if abs(r["outer_coil_od_in"] - 1.5) < 1e-6]


def windings(r, awg):
    inner, outer = windows_for(r["emag_height_in"], r["inner_core_dia_in"], r["outer_coil_od_in"])
    return design(r["req_inner_At"], r["req_outer_At"], inner, outer, awg,
                  envelope_area_m2=r["envelope_area_m2"])


def extract(tag):
    return np.load(ASSEMBLY / "sweep_runs" / tag / "extract.npz")


# --- geometry outlines (r-y half section, mm) ------------------------------
def parts(h_in, core_dia_in, outer_od_in=1.5, outer_core_in=1.0):
    y_top = Y_TOP_IN + h_in * IN
    rc = core_dia_in * IN / 2
    ro_core, ro_coil = outer_core_in * IN / 2, outer_od_in * IN / 2
    P = {"iron": [], "coil": [], "chamber": [], "spacer": []}
    P["iron"] += [(R_BORE_TOP, Y_TOP_OUT, R_PLATE, Y_TOP_IN),            # top plate
                  (3.57, y_top, R_PLATE, y_top + PLATE_T),               # bottom plate
                  (0, Y_TOP_OUT, R_POLE, Y_TOP_IN),                      # pole cap
                  (0, Y_TOP_IN, rc, y_top),                              # inner shaft
                  (R_OUTER_AXIS - ro_core, Y_TOP_IN, R_OUTER_AXIS + ro_core, y_top)]  # outer core
    if y_top > Y_ANODE_FACE + 1e-6:
        P["spacer"] += [(R_POLE, Y_ANODE_FACE, R_CHAMBER_OUT, y_top)]
    P["coil"] += [(rc, Y_TOP_IN, R_POLE, y_top),
                  (R_OUTER_AXIS - ro_coil, Y_TOP_IN, R_OUTER_AXIS - ro_core, y_top),
                  (R_OUTER_AXIS + ro_core, Y_TOP_IN, R_OUTER_AXIS + ro_coil, y_top)]
    P["chamber"] += [(R_POLE, Y_TOP_OUT, R_CH_IN, Y_ANODE_FACE),
                     (R_CH_OUT, Y_TOP_OUT, R_CHAMBER_EXIT_OUT, Y_TOP_IN),
                     (R_CH_OUT, Y_TOP_IN, R_CHAMBER_OUT, Y_ANODE_FACE),
                     (R_CH_IN, 3.175, R_CH_OUT, Y_ANODE_FACE)]
    return P


STYLE = {"iron": dict(fc="#b9b8b2", ec="#52514e"), "spacer": dict(fc="#d9d8d3", ec="#52514e", hatch="////"),
         "coil": dict(fc="#f3c8a8", ec="#b0582a"), "chamber": dict(fc="#f4f3ef", ec="#8a8983")}


def draw_parts(ax, P, mirror=True, filled=True, lw=0.7):
    for kind, rects in P.items():
        st = STYLE[kind]
        for (r0, y0, r1, y1) in rects:
            for sgn in ((1, -1) if mirror else (1,)):
                xs = sorted([sgn * r0, sgn * r1])
                ax.add_patch(Rectangle((xs[0], y0), xs[1] - xs[0], y1 - y0,
                                       fc=st["fc"] if filled else "none", ec=st["ec"], lw=lw,
                                       hatch=st.get("hatch") if filled else None, zorder=3 if filled else 4))


def fig_geometry():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.4), sharey=True)
    for ax, (h, title) in zip(axes, [(2.0, "emag_height 2.000 in (spacer 0.875 in)"),
                                     (1.125, "emag_height 1.125 in (spacer 0)")]):
        draw_parts(ax, parts(h, 1.0), mirror=False)
        ax.add_patch(Rectangle((R_CH_IN, -31.75), R_CH_OUT - R_CH_IN, 31.75 + 3.175, fc="none",
                               ec=SERIES[0], lw=1.2, ls="--", zorder=5))
        y_top = Y_TOP_IN + h * IN
        ax.annotate("", xy=(87, Y_TOP_IN), xytext=(87, y_top),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
        ax.text(88.5, (Y_TOP_IN + y_top) / 2, "emag\nheight", fontsize=7, color=INK2, va="center")
        if y_top > Y_ANODE_FACE:
            ax.plot([R_CHAMBER_OUT, 104], [Y_ANODE_FACE] * 2, color=MUTED, lw=0.5, ls=":")
            ax.annotate("", xy=(102, Y_ANODE_FACE), xytext=(102, y_top),
                        arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8))
            ax.text(103.5, (Y_ANODE_FACE + y_top) / 2, "spacer", fontsize=7, color=INK, va="center")
        ax.text(26.9, -14, "channel", fontsize=7, color=SERIES[0], ha="center", rotation=90, zorder=6)
        ax.set_xlim(0, 118)
        ax.set_ylim(40, -40)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("radius r (mm)")
        ax.grid(False)
    axes[0].set_ylabel("axial position y (mm), exit at top")
    lbl = [("top plate", 55, -27.7), ("pole cap", 1.5, -27.7), ("inner\ncore", 1.5, -8),
           ("inner coil", 14.0, -8), ("outer core + coils", 56, -12), ("bottom plate", 40, 32.8)]
    for t, x, y in lbl:
        axes[0].text(x, y, t, fontsize=6.5, color=INK, zorder=7,
                     rotation=90 if t == "inner coil" else 0, va="center")
    axes[0].text(29.4, 17, "spacer\n(A36)", fontsize=6.5, color=INK, zorder=7, va="center", ha="center",
                 bbox=dict(fc="white", ec="none", pad=0.8))
    fig.tight_layout()
    fig.savefig(FIG / "geometry.pdf")
    plt.close(fig)


FM_Y = (-46, 42)   # common axial window so the two geometries compare directly


def field_map(ax, tag, h, core, vmin=1e-3, vmax=1.2):
    """Half-section (r >= 0) through the axis and one outer core."""
    d = extract(tag)
    scale = float(d["scale"])
    sl = d["slice_core"]
    s, y = sl[:, 0] * 1e3, sl[:, 1] * 1e3
    Bs, By, Bm = sl[:, 2] * scale, sl[:, 3] * scale, sl[:, 4] * scale
    keep = s > -3
    s, y, Bs, By, Bm = s[keep], y[keep], Bs[keep], By[keep], Bm[keep]
    gs = np.linspace(0, 90, 361)
    gy = np.linspace(FM_Y[0], FM_Y[1], int((FM_Y[1] - FM_Y[0]) / 0.25) + 1)
    GS, GY = np.meshgrid(gs, gy)
    pts = np.column_stack([s, y])
    M = griddata(pts, Bm, (GS, GY), method="linear")
    U = griddata(pts, Bs, (GS, GY), method="linear")
    V = griddata(pts, By, (GS, GY), method="linear")
    im = ax.pcolormesh(GS, GY, np.clip(M, vmin, None), cmap=CMAP, norm=LogNorm(vmin, vmax),
                       shading="auto", rasterized=True, zorder=1)
    ax.streamplot(gs, gy, np.nan_to_num(U), np.nan_to_num(V), density=1.1, color="#3a3a38",
                  linewidth=0.4, arrowsize=0.5, zorder=2)
    draw_parts(ax, parts(h, core), mirror=False, filled=False, lw=0.7)
    ax.add_patch(Rectangle((R_CH_IN, -31.75), R_CH_OUT - R_CH_IN, 34.925, fc="none",
                           ec=SERIES[1], lw=1.1, ls="--", zorder=5))
    ax.set_xlim(0, 90)
    ax.set_ylim(FM_Y[1], FM_Y[0])
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("radius r (mm)")
    return im


def fig_field_maps(rows):
    cases = [(2.0, 1.0), (1.125, 1.0)]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.9), sharey=True)
    for ax, (h, core) in zip(axes, cases):
        r = next(r for r in rows if r["emag_height_in"] == h and r["inner_core_dia_in"] == core
                 and r["outer_coil_od_in"] == 1.5)
        im = field_map(ax, r["tag"], h, core)
        ax.set_title(f"spacer {h - 1.125:.3f} in\n{r['req_inner_At']:.0f} / {r['req_outer_At']:.0f} A-turns (inner / outer)")
    axes[0].set_ylabel("axial position y (mm), exit at top")
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("|B| (T), log scale")
    cb.outline.set_edgecolor(MUTED)
    fig.savefig(FIG / "field_maps.pdf")
    plt.close(fig)


def fig_channel_zoom(rows):
    """Channel close-up, baseline geometry: B vectors over the plasma cavity."""
    r = next(r for r in rows if r["emag_height_in"] == 2.0 and r["inner_core_dia_in"] == 1.0
             and r["outer_coil_od_in"] == 1.5)
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.6))
    for ax, (h, core) in zip(axes, [(2.0, 1.0), (1.125, 1.0)]):
        rr = next(x for x in rows if x["emag_height_in"] == h and x["inner_core_dia_in"] == core
                  and x["outer_coil_od_in"] == 1.5)
        d = extract(rr["tag"])
        sc = float(d["scale"])
        sl = d["slice_core"]
        s, y = sl[:, 0] * 1e3, sl[:, 1] * 1e3
        gs = np.linspace(12, 45, 265)
        gy = np.linspace(-40, 12, 417)
        GS, GY = np.meshgrid(gs, gy)
        pts = np.column_stack([s, y])
        M = griddata(pts, sl[:, 4] * sc * 1e4, (GS, GY), method="linear")
        U = griddata(pts, sl[:, 2] * sc, (GS, GY), method="linear")
        V = griddata(pts, sl[:, 3] * sc, (GS, GY), method="linear")
        im = ax.pcolormesh(GS, GY, M, cmap=CMAP, vmin=0, vmax=450, shading="auto", rasterized=True)
        ax.streamplot(gs, gy, np.nan_to_num(U), np.nan_to_num(V), density=1.3, color=INK,
                      linewidth=0.5, arrowsize=0.6)
        draw_parts(ax, parts(h, core), mirror=False, filled=False, lw=0.7)
        ax.add_patch(Rectangle((R_CH_IN, -31.75), R_CH_OUT - R_CH_IN, 34.925, fc="none",
                               ec=SERIES[1], lw=1.1, ls="--", zorder=5))
        ax.set_xlim(12, 45)
        ax.set_ylim(12, -40)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(f"spacer {h - 1.125:.3f} in")
        ax.set_xlabel("radius r (mm)")
    axes[0].set_ylabel("axial position y (mm), exit at top")
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label("|B| (G)")
    cb.outline.set_edgecolor(MUTED)
    fig.savefig(FIG / "channel_zoom.pdf")
    plt.close(fig)


def fig_profiles(rows):
    g = sorted([r for r in main_grid(rows) if r["inner_core_dia_in"] == 1.0], key=lambda r: r["emag_height_in"])
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.9))
    for i, r in enumerate(g):
        d = extract(r["tag"])
        sc = float(d["scale"])
        y = d["prof_y"] * 1e3
        s = (y - 3.175) / (-31.75 - 3.175)
        col = SEQ_ORD[i]
        lbl = f"{r['spacer_in']:.3f} in"
        axes[0].plot(s, d["prof_Br_G"] * sc, color=col, label=lbl)
        axes[1].plot(s, d["prof_Bax_G"] * sc, color=col, label=lbl)
    ss = np.linspace(0, 1, 100)
    axes[0].plot(ss, 30 + 270 * np.exp(-3 * (1 - ss)), color=MUTED, ls="--", lw=1.2, label="target lens")
    axes[0].set_ylabel(r"$B_r$ (G), channel-averaged")
    axes[1].set_ylabel(r"$B_{axial}$ (G), channel-averaged")
    for ax in axes:
        ax.set_xlabel("position along channel (0 = anode, 1 = exit)")
        ax.set_xlim(0, 1)
    axes[0].legend(title="spacer", ncol=1, frameon=False, fontsize=6.5, title_fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "channel_profiles.pdf")
    plt.close(fig)


def series_plot(ax, rows, ykey, label=None, fn=None):
    for i, core in enumerate([1.0, 0.875, 0.75]):
        g = sorted([r for r in main_grid(rows) if r["inner_core_dia_in"] == core], key=lambda r: r["spacer_in"])
        if not g:
            continue
        x = [r["spacer_in"] for r in g]
        y = [fn(r) if fn else r[ykey] for r in g]
        ax.plot(x, y, color=SERIES[i], marker=MARKERS[i], label=CORE_LABEL[core],
                mec="white", mew=0.8)
    ax.set_xlabel("spacer height (in)")
    if label:
        ax.set_ylabel(label)


def fig_ampere_turns(rows):
    """y-range is +/-10% of the baseline so the ~1% point-to-point mesh noise
    reads as noise (shaded band) rather than as structure."""
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.5))
    base = next(r for r in main_grid(rows) if r["emag_height_in"] == 2.0 and r["inner_core_dia_in"] == 1.0)
    for ax, key, lbl in [(axes[0], "req_inner_At", "inner coil A-turns"),
                         (axes[1], "req_outer_At", "A-turns per outer coil")]:
        b = base[key]
        ax.axhspan(0.99 * b, 1.01 * b, color="#f0efec", zorder=0)
        series_plot(ax, rows, key, lbl)
        ax.set_ylim(0.9 * b, 1.1 * b)
    axes[1].text(0.02, 0.965 * base["req_outer_At"], "shaded: baseline $\\pm$1% (mesh noise)",
                 fontsize=6.5, color=INK2)
    axes[0].legend(frameon=False, fontsize=6.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "ampere_turns.pdf")
    plt.close(fig)


def fig_core_flux(rows):
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.7), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axes[0]
    ax.axhspan(KNEE_LO, KNEE_HI, color="#f0efec", zorder=0)
    ax.text(0.02, 1.4, "1020 knee region", fontsize=7, color=INK2, va="center")
    series_plot(ax, rows, "inner_shaft_B_max_T", "core shaft |B| (T)")
    for i, core in enumerate([1.0, 0.875, 0.75]):
        g = sorted([r for r in main_grid(rows) if r["inner_core_dia_in"] == core], key=lambda r: r["spacer_in"])
        ax.plot([r["spacer_in"] for r in g], [r["outer_shaft_B_max_T"] for r in g], color=SERIES[i],
                ls=":", lw=1.1)
    ax.text(0.45, 0.28, "outer cores (dotted)", fontsize=7, color=INK2)
    ax.set_ylim(0, 1.6)
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")

    ax = axes[1]
    ax.plot(np.array(FNAL1020_H), FNAL1020_B, color=INK, label="1020, measured (cores)")
    ax.plot(np.array(AISI1008_H), AISI1008_B, color=MUTED, ls="--", label="1008 proxy (plates)")
    ax.axhspan(KNEE_LO, KNEE_HI, color="#f0efec", zorder=0)
    worst = max(main_grid(rows), key=lambda r: r["inner_shaft_B_max_T"])
    ax.axhline(worst["inner_shaft_B_max_T"], color=SERIES[2], lw=1.0)
    ax.text(60, worst["inner_shaft_B_max_T"] + 0.04, "worst core in sweep", fontsize=7, color=INK2)
    ax.set_xscale("log")
    ax.set_xlim(50, 2e5)
    ax.set_ylim(0, 2.2)
    ax.set_xlabel("H (A/m)")
    ax.set_ylabel("B (T)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "core_flux.pdf")
    plt.close(fig)


def fig_winding(rows, awg=20):
    fig, axes = plt.subplots(1, 4, figsize=(6.8, 2.3))
    series_plot(axes[0], rows, None, "series current (A)", fn=lambda r: windings(r, awg)["current_A"])
    series_plot(axes[1], rows, None, "string voltage, hot (V)", fn=lambda r: windings(r, awg)["V_hot"])
    series_plot(axes[2], rows, None, "power, hot (W)", fn=lambda r: windings(r, awg)["P_hot_W"])
    series_plot(axes[3], rows, None, "hottest coil (C)", fn=lambda r: windings(r, awg)["T_coil_hot_C"])
    for ax in axes:
        ax.set_ylim(bottom=0)
    axes[0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(FIG / "winding_awg20.pdf")
    plt.close(fig)


def fig_psu(rows):
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    P = []
    for r in main_grid(rows):
        for awg in WIRE:
            w = windings(r, awg)
            if w:
                P.append((awg, w["current_A"], w["V_hot"], w["P_hot_W"]))
    P = np.array(P)
    for p in (4, 8, 16):
        I = np.linspace(0.2, 6, 200)
        ax.plot(I, p / I, color=GRID, lw=1.0, zorder=0)
        ax.text(5.7, p / 5.7 * 1.06, f"{p} W", fontsize=7, color=MUTED, ha="right")
    ax.scatter(P[:, 1], P[:, 2], s=14, color=SERIES[0], edgecolor="white", lw=0.5, zorder=3)
    for awg in WIRE:
        m = P[:, 0] == awg
        ax.text(P[m, 1].max() * 1.05, np.median(P[m, 2]), f"AWG {awg}", fontsize=7, color=INK, va="center")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.3, 8)
    ax.set_xticks([0.5, 1, 2, 5]); ax.set_xticklabels(["0.5", "1", "2", "5"])
    ax.set_yticks([1, 2, 5, 10, 20]); ax.set_yticklabels(["1", "2", "5", "10", "20"])
    ax.minorticks_off()
    ax.set_ylim(0.5, 40)
    ax.set_xlabel("series current (A)")
    ax.set_ylabel("string voltage at operating temp. (V)")
    ax.set_title("Every swept geometry x gauge, one point each")
    fig.tight_layout()
    fig.savefig(FIG / "psu_map.pdf")
    plt.close(fig)


# --- LaTeX generation -------------------------------------------------------
def tex_num(x, fmt):
    return format(x, fmt)


def spacer_cell(sp):
    stock = any(abs(sp - s) < 1e-4 for s in STOCK_SPACERS)
    return f"{sp:.3f}" + (r"$^\dagger$" if stock and sp > 0 else "")


def with_breaks(L, g):
    """Blank space between inner-core-diameter groups (rows are in g's order)."""
    out = []
    for i, (line, r) in enumerate(zip(L, g)):
        if i and r["inner_core_dia_in"] != g[i - 1]["inner_core_dia_in"]:
            out.append(r"\addlinespace[4pt]")
        out.append(line)
    return out


def write_tables(rows, extras):
    GEN.mkdir(exist_ok=True)
    g = sorted(main_grid(rows), key=lambda r: (-r["inner_core_dia_in"], r["spacer_in"]))

    # geometry
    L = []
    for r in g:
        inner, outer = windows_for(r["emag_height_in"], r["inner_core_dia_in"], r["outer_coil_od_in"])
        L.append(" & ".join([
            spacer_cell(r["spacer_in"]), f"{r['emag_height_in']:.3f}",
            f"{r['inner_core_dia_in']:.3f}", "1.375", f"{inner.build * 1e3:.2f}",
            f"{inner.build * inner.height * 1e6:.0f}",
            "1.000", f"{r['outer_coil_od_in']:.3f}", f"{outer.build * outer.height * 1e6:.0f}",
        ]) + r" \\")
    (GEN / "tab_geometry.tex").write_text("\n".join(with_breaks(L, g)) + "\n")

    # magnetics
    L = []
    for r in g:
        L.append(" & ".join([
            spacer_cell(r["spacer_in"]), f"{r['inner_core_dia_in']:.3f}",
            f"{r['req_inner_At']:.0f}", f"{r['req_outer_At']:.0f}",
            f"{r['B_peak_G']:.0f}", f"{r['B_anode_G']:.0f}", f"{r['radial_purity_pct']:.1f}",
            f"{r['inner_shaft_B_max_T']:.2f}", f"{r['outer_shaft_B_max_T']:.2f}",
            f"{r['iron_B_p999_T']:.2f}",
        ]) + r" \\")
    (GEN / "tab_magnetics.tex").write_text("\n".join(with_breaks(L, g)) + "\n")

    # AWG 20 windings
    L = []
    for r in g:
        w = windings(r, 20)
        L.append(" & ".join([
            spacer_cell(r["spacer_in"]), f"{r['inner_core_dia_in']:.3f}",
            f"{w['inner_turns']}", f"{w['inner_layers']}", f"{w['outer_turns']}", f"{w['outer_layers']}",
            f"{w['current_A']:.2f}", f"{w['R20_ohm']:.2f}", f"{w['V_hot']:.2f}", f"{w['P_hot_W']:.1f}",
            f"{w['T_coil_hot_C']:.0f}", f"{w['total_wire_m']:.0f}",
        ]) + r" \\")
    (GEN / "tab_awg20.tex").write_text("\n".join(with_breaks(L, g)) + "\n")

    # all gauges, full appendix
    L = []
    for r in g:
        first = True
        for awg in WIRE:
            w = windings(r, awg)
            if not w:
                continue
            L.append(" & ".join([
                spacer_cell(r["spacer_in"]) if first else "", f"{r['inner_core_dia_in']:.3f}" if first else "",
                f"{awg}", f"{w['inner_turns']}/{w['outer_turns']}", f"{w['current_A']:.2f}",
                f"{w['R_hot_ohm']:.2f}", f"{w['V20']:.2f}", f"{w['V_hot']:.2f}", f"{w['P_hot_W']:.1f}",
                f"{w['T_coil_hot_C']:.0f}", f"{w['J_A_per_mm2']:.1f}",
            ]) + r" \\")
            first = False
        L.append(r"\addlinespace[2pt]")
    (GEN / "tab_all_gauges.tex").write_text("\n".join(L) + "\n")

    # power-supply requirements per gauge, over the whole sweep
    L = []
    W_J = extras["stored_energy_J"]
    for awg in WIRE:
        ws = [windings(r, awg) for r in g]
        ws = [w for w in ws if w]
        Imin, Imax = min(w["current_A"] for w in ws), max(w["current_A"] for w in ws)
        Vmax = max(w["V_hot"] for w in ws)
        Pmax = max(w["P_hot_W"] for w in ws)
        L_mH = 2 * W_J / Imax ** 2 * 1e3, 2 * W_J / Imin ** 2 * 1e3
        I_rate, V_rate = 1.25 * Imax, 2.0 * Vmax
        L.append(" & ".join([
            f"{awg}", f"{Imin:.2f}--{Imax:.2f}", f"{Vmax:.1f}", f"{Pmax:.1f}",
            f"{L_mH[0]:.0f}--{L_mH[1]:.0f}", f"$\\geq${I_rate:.1f} A, $\\geq${V_rate:.0f} V",
        ]) + r" \\")
    (GEN / "tab_psu.tex").write_text("\n".join(L) + "\n")

    # headline numbers as macros
    base = next(r for r in g if r["emag_height_in"] == 2.0 and r["inner_core_dia_in"] == 1.0)
    zero = next(r for r in g if r["emag_height_in"] == 1.125 and r["inner_core_dia_in"] == 1.0)
    worst = max(g, key=lambda r: r["inner_shaft_B_max_T"])
    wb, wz = windings(base, 20), windings(zero, 20)
    macros = {
        "BaseInnerAt": f"{base['req_inner_At']:.0f}", "BaseOuterAt": f"{base['req_outer_At']:.0f}",
        "ZeroInnerAt": f"{zero['req_inner_At']:.0f}", "ZeroOuterAt": f"{zero['req_outer_At']:.0f}",
        "AtChangePct": f"{(zero['req_inner_At'] / base['req_inner_At'] - 1) * 100:+.1f}",
        "BaseI": f"{wb['current_A']:.2f}", "ZeroI": f"{wz['current_A']:.2f}",
        "BaseP": f"{wb['P_hot_W']:.1f}", "ZeroP": f"{wz['P_hot_W']:.1f}",
        "BaseV": f"{wb['V_hot']:.2f}", "ZeroV": f"{wz['V_hot']:.2f}",
        "BaseT": f"{wb['T_coil_hot_C']:.0f}", "ZeroT": f"{wz['T_coil_hot_C']:.0f}",
        "WorstCoreB": f"{worst['inner_shaft_B_max_T']:.2f}",
        "WorstCoreDesc": f"{worst['inner_core_dia_in']:.3f}~in core, {worst['spacer_in']:.3f}~in spacer",
        "BaseShaftB": f"{base['inner_shaft_B_max_T']:.2f}",
        "StoredEnergyMJ": f"{W_J * 1e3:.0f}",
        "LinVsNLPct": f"{extras['linear_vs_nonlinear_pct']:+.1f}",
        "NRuns": f"{len(rows)}",
    }
    od = [r for r in rows if abs(r["outer_coil_od_in"] - 1.625) < 1e-6]
    if od:
        ref = next(r for r in g if r["emag_height_in"] == od[0]["emag_height_in"] and r["inner_core_dia_in"] == 1.0)
        macros["OuterODChangePct"] = f"{(od[0]['req_outer_At'] / ref['req_outer_At'] - 1) * 100:+.1f}"
        macros["OuterODInnerChangePct"] = f"{(od[0]['req_inner_At'] / ref['req_inner_At'] - 1) * 100:+.1f}"
    for k, v in extras.get("macros", {}).items():
        macros[k] = v
    (GEN / "numbers.tex").write_text("\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    extras = json.loads((REPORT / "report_inputs.json").read_text())

    fig_geometry()
    fig_field_maps(rows)
    fig_channel_zoom(rows)
    fig_profiles(rows)
    fig_ampere_turns(rows)
    fig_core_flux(rows)
    fig_winding(rows)
    fig_psu(rows)
    write_tables(rows, extras)
    print("figures + tables written")

    if not args.no_pdf:
        tectonic = ROOT / "tools" / "tectonic" / "tectonic.exe"
        subprocess.run([str(tectonic), "bpl700_magnet_report.tex"], cwd=REPORT, check=True)
        # The committed copy lives at the repo root; report/*.pdf is a build artifact.
        import shutil
        shutil.copy(REPORT / "bpl700_magnet_report.pdf", ROOT / "EMAG_REPORT.pdf")
        print(f"wrote {ROOT / 'EMAG_REPORT.pdf'}")


if __name__ == "__main__":
    main()
