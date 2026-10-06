"""
field_cross_sections.py

B-field cross-section figures for the design write-up: one direct linear
solve per geometry at the excitation it actually runs at, fields kept, then
r-y slices through the thrust axis and the channel B_r profile against the
field_quality.py target.

Geometries:
  baseline -- Onshape defaults (1.0 in cores), committed 338/214.5 ratio
              scaled to 300 G at the exit
  final    -- sizing/EMAG_SPEC.md: 0.75 in cores, inner 230 t x 1.531 A,
              outer 184 t x 1.214 A

Usage (from gmsh_getdp/assembly/):  python field_cross_sections.py
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from bpl1000_geometry import CHANNEL, IN  # noqa: E402
from build_assembly import build, EXCITATION  # noqa: E402
from em_sizing_sweep import rewrite_channel, IRON_MUR, TARGET_G  # noqa: E402
from emag_sweep import run_getdp, slice_plane  # noqa: E402
from field_quality import axial_profile, evaluate_field_quality, target_profile  # noqa: E402
from generate_regions import generate_regions_pro  # noqa: E402
from onshape.cache import DEFAULTS_IN  # noqa: E402
from params_utils import parse_params  # noqa: E402
from pos_utils import parse_pos_elements  # noqa: E402

OUT = HERE / "em_sizing_figures"
RUNS = HERE / "sweep_runs_bpl1000" / "cross_sections"

CASES = {
    "baseline": dict(cfg=dict(DEFAULTS_IN), ni=(344.15, 218.41),
                     label="Baseline: 1.0 in cores (Onshape defaults)"),
    "final": dict(cfg=dict(DEFAULTS_IN, iem_coreod=0.75, iem_coilt=0.375,
                           oem_coreod=0.75, oem_coilt=0.3125),
                  ni=(230 * 1.531, 184 * 1.214),
                  label="Final spec: 0.75 in cores, 230 t x 1.53 A / 184 t x 1.21 A"),
}


def solve(name, case):
    run_dir = RUNS / name
    run_dir.mkdir(parents=True, exist_ok=True)
    if not (run_dir / "b_assembly.pos").exists():
        build(case["cfg"], str(run_dir), None, include_chamber_spacer_as_iron=True)
        rewrite_channel(run_dir / "assembly_params.txt")
        params = parse_params(run_dir / "assembly_params.txt")
        ni_in, ni_out = case["ni"]
        exc = {"inner_coil": dict(turns=1, current=ni_in, polarity=EXCITATION["inner_coil"]["polarity"]),
               "outer_coil": dict(turns=1, current=ni_out, polarity=EXCITATION["outer_coil"]["polarity"])}
        generate_regions_pro(params, exc, out_path=str(run_dir / "assembly_regions_generated.pro"),
                             iron_mur=IRON_MUR)
        run_getdp(run_dir)
    params = parse_params(run_dir / "assembly_params.txt")
    nodes, B = parse_pos_elements(run_dir / "b_assembly.pos")
    return params, nodes.mean(axis=1), B


def grid_slice(sl, s_lim, y_lim, n=320):
    from scipy.interpolate import griddata
    s, y = sl[:, 0], sl[:, 1]
    S, Y = np.meshgrid(np.linspace(*s_lim, n), np.linspace(*y_lim, n))
    Bs = griddata((s, y), sl[:, 2], (S, Y), method="linear")
    By = griddata((s, y), sl[:, 3], (S, Y), method="linear")
    return S, Y, Bs, By


def draw(name, case, params, cent, B):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    ax_x, ax_z = params["axis_x"], params["axis_z"]
    outer = [w["loc"] for w in params["windings"] if w["pole"] == "outer_coil"]
    q = max(outer, key=lambda l: l[0] + l[1])
    nrm = math.hypot(q[0] - ax_x, q[1] - ax_z)
    diag = ((q[0] - ax_x) / nrm, (q[1] - ax_z) / nrm)
    between = ((diag[0] - diag[1]) / math.sqrt(2), (diag[0] + diag[1]) / math.sqrt(2))

    fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
    for ax, (direction, title) in zip(axs, ((diag, "slice through an outer core"),
                                            (between, "slice between outer cores"))):
        sl = slice_plane(cent, B, ax_x, ax_z, direction, half_width=0.0012, r_max=0.1)
        S, Y, Bs, By = grid_slice(sl, (-0.085, 0.085), (-0.06, 0.035))
        Bm = np.hypot(Bs, By) * 1e4
        im = ax.pcolormesh(S * 1e3, Y * 1e3, Bm, shading="auto", cmap="magma",
                           norm=LogNorm(vmin=20, vmax=8000))
        ok = np.isfinite(Bs) & np.isfinite(By)
        ax.streamplot(S[0] * 1e3, Y[:, 0] * 1e3, np.where(ok, Bs, 0), np.where(ok, By, 0),
                      color="w", linewidth=0.5, density=1.4, arrowsize=0.6)
        for sgn in (-1, 1):
            r0, r1 = CHANNEL["channel_inner_r"] * 1e3, CHANNEL["channel_outer_r"] * 1e3
            y0, y1 = CHANNEL["channel_y_exit"] * 1e3, CHANNEL["channel_y_anode"] * 1e3
            xs = sorted((sgn * r0, sgn * r1))
            ax.add_patch(plt.Rectangle((xs[0], y0), xs[1] - xs[0], y1 - y0, fill=False,
                                       ec="cyan", lw=1.4, ls="--"))
        ax.set_aspect("equal")
        ax.set_xlabel("distance from thrust axis (mm)")
        ax.set_ylabel("axial position y (mm), exit at -37.8")
        ax.set_title(title)
    fig.colorbar(im, ax=axs, label="|B| (gauss, log scale)", shrink=0.9, pad=0.01)
    fig.suptitle(f"{case['label']}  (cyan dashed = discharge channel)")
    fig.savefig(OUT / f"cross_section_{name}.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    prof = axial_profile(cent, B, ax_x, ax_z, CHANNEL["channel_inner_r"], CHANNEL["channel_outer_r"],
                         CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"], n_stations=20)
    p = dict(prof, B_r=prof["B_r"] * 1e4, Bmag=prof["Bmag"] * 1e4)
    fq = evaluate_field_quality(p, B_peak=TARGET_G)
    return p, fq


def profile_figure(profiles):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5))
    s = np.linspace(0, 1, 200)
    ax.plot(s, target_profile(s, TARGET_G, 30.0, 3.0), "k--", lw=1.2,
            label="target shape (field_quality.py)")
    for (name, (p, fq)), c in zip(profiles.items(), ("0.6", "C0")):
        ax.plot(p["s"], p["B_r"], "o-", color=c, ms=4,
                label=f"{CASES[name]['label'].split(':')[0]}: score {fq['composite_score']:.1f}, "
                      f"radial purity {fq['radial_purity_pct']:.1f}%")
    ax.axhline(TARGET_G, color="0.8", lw=0.8)
    ax.set_xlabel("position along channel, anode (0) to exit (1)")
    ax.set_ylabel("radial field B_r, channel-averaged (G)")
    ax.set_title("Channel radial-field profile vs target")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "channel_profile.png", dpi=140)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    profiles = {}
    for name, case in CASES.items():
        params, cent, B = solve(name, case)
        profiles[name] = draw(name, case, params, cent, B)
        p, fq = profiles[name]
        print(f"{name}: B_r exit {p['B_r'][-1]:.1f} G, peak {np.nanmax(p['B_r']):.1f} G, "
              f"score {fq['composite_score']:.2f}, purity {fq['radial_purity_pct']:.2f}%")
    profile_figure(profiles)
    print(json.dumps({k: {"score": v[1]["composite_score"]} for k, v in profiles.items()}))


if __name__ == "__main__":
    main()
