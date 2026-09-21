"""
emag_sweep.py

Geometry sweep over emag_height (equivalently the chamber_spacer height)
and the inner core diameter, for sizing the spacer.

Geometry facts this sweep is built on, measured from the Onshape STEP
exports (2026-09-21), not assumed:

- The chamber and top plate never move. The chamber's anode end sits at
  y = +5.59mm; emag_height moves only the bottom plate. So
      spacer height = emag_height - 1.125in
  and spacer >= 0 means emag_height >= 1.125in. At exactly 1.125in Onshape
  drops the spacer body and the bottom plate sits on the chamber.
- The chamber bore and the spacer bore are both r = 17.46mm, which is the
  inner coil's OD (1.375in). With the chamber fixed, the inner winding can
  only grow INWARD, by thinning the core (inner_coil_id = core diameter).
  The pole cap stays at r = 17.46mm whatever inner_coil_id is -- only the
  shaft under the winding thins.
- Outer coils can grow outward to ~1.70in OD before touching the chamber's
  outer wall (r = 41.28mm; outer core axes are 62.86mm off-axis).

Every point is solved with LINEAR iron (constant mu_r, IRON_MUR). Justified
because the circuit is gap-dominated and the iron runs below the B-H knee;
`--validate` checks it against the committed nonlinear solve instead of
taking that on faith. The field is then linear in the ampere-turns, so one
solve per geometry at the baseline excitation gives the A-turns needed for
TARGET_G at the exit by scaling. The core flux densities it reports are
scaled the same way; any point whose scaled core field approaches the knee
(KNEE_T) is flagged, because there the linear answer stops being valid.

Each run lives in sweep_runs/<tag>/. The ~200MB .pos dumps and the mesh
are deleted after extraction; what is kept is metrics.json (small, the
source of emag_sweep_results.csv) and extract.npz (channel elements plus
two r-y slices, for the report figures). Re-running skips finished points.

Usage (from gmsh_getdp/assembly/):
  python emag_sweep.py                 # full grid
  python emag_sweep.py --only h2.000_ci1.000_oo1.500
  python emag_sweep.py --validate      # linear vs committed nonlinear, baseline geometry
"""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "thermal"))

from build_assembly import build, EXCITATION  # noqa: E402
from field_quality import axial_profile, evaluate_field_quality  # noqa: E402
from generate_regions import generate_regions_pro  # noqa: E402
from params_utils import parse_params  # noqa: E402
from pos_utils import parse_pos_elements  # noqa: E402

GETDP_EXE = HERE.parents[1] / "tools" / "getdp-3.5.0-Windows64" / "getdp.exe"
PRO_FILE = HERE / "magnetostatics_assembly.pro"
RUNS_DIR = HERE / "sweep_runs"
RESULTS_CSV = HERE / "emag_sweep_results.csv"

IN = 0.0254
SPACER_ZERO_IN = 1.125      # emag_height at which the spacer vanishes
IRON_MUR = 1500.0           # see --validate
TARGET_G = 300.0            # radial field at the channel exit
KNEE_T = 1.5                # past this the proxy B-H curve stops being near-linear
N_STATIONS = 20

# Fixed geometry, CAD-confirmed (gmsh_getdp/README.md "Channel ROI"). The
# chamber doesn't move with any sweep parameter, so these override
# detect_channel_cavity()'s ~0.76mm-off estimate rather than trusting it.
CHANNEL = dict(channel_inner_r=0.022098, channel_outer_r=0.03175,
               channel_y_anode=0.003175, channel_y_exit=-0.03175)
Y_TOP_PLATE_INNER = -0.02299     # coils/cores start here...
INNER_POLE_R = 0.0174625         # ...and the inner pole cap has this radius
OUTER_CORE_AXIS_R = 0.0628618    # |(44.45, 44.45)mm|

HEIGHTS_IN = [1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0]
CORE_DIAS_IN = [1.0, 0.875, 0.75]
OUTER_OD_IN = 1.5
EXTRA_POINTS = [(1.5, 1.0, 1.625)]   # outer-OD sensitivity check


def tag_for(h, ci, oo):
    return f"h{h:.3f}_ci{ci:.3f}_oo{oo:.3f}"


def geometry(h, ci, oo):
    return dict(emag_height=h, inner_coil_id=ci, inner_coil_od=1.375,
                outer_coil_id=1.0, outer_coil_od=oo)


def rewrite_channel(params_path):
    """Replace build_assembly.py's detected channel bounds with the fixed
    CAD-confirmed ones."""
    lines = Path(params_path).read_text().splitlines()
    vals = dict(CHANNEL,
                channel_y_min=min(CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"]),
                channel_y_max=max(CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"]))
    out = []
    for ln in lines:
        key = ln.split(" ", 1)[0]
        out.append(f"{key} {vals[key]!r}" if key in vals else ln)
    Path(params_path).write_text("\n".join(out) + "\n")


def run_getdp(run_dir):
    shutil.copy(PRO_FILE, run_dir / PRO_FILE.name)
    t0 = time.time()
    with open(run_dir / "getdp.log", "w") as log:
        rc = subprocess.run([str(GETDP_EXE), PRO_FILE.name, "-msh", "assembly.msh",
                             "-solve", "Res_a", "-pos", "Map_b"],
                            cwd=run_dir, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"getdp failed (exit {rc}), see {run_dir / 'getdp.log'}")
    return time.time() - t0


def tet_volumes(nodes):
    a, b, c = (nodes[:, i] - nodes[:, 0] for i in (1, 2, 3))
    return np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0


def shaft_profile(cent, B, cx, cz, r_core, y_lo, y_hi, n_bins=12):
    """Mean axial flux density across a core's shaft cross-section, per
    axial bin. Flux accumulates along a core as leakage joins it, so the
    max over bins is what sets saturation, not the mid-height value."""
    r = np.hypot(cent[:, 0] - cx, cent[:, 2] - cz)
    y = cent[:, 1]
    sel = (r < 0.9 * r_core) & (y > y_lo) & (y < y_hi)
    edges = np.linspace(y_lo, y_hi, n_bins + 1)
    means = []
    for i in range(n_bins):
        m = sel & (y >= edges[i]) & (y < edges[i + 1])
        if m.sum() >= 5:
            means.append(abs(B[m, 1].mean()))
    return np.array(means)


def slice_plane(cent, B, axis_x, axis_z, direction, half_width=0.0015, r_max=0.1):
    """Elements within half_width of the plane containing the thrust axis and
    the in-plane unit vector `direction` (x,z). Returns in-plane coordinate
    s (signed distance from axis), y, and B projected onto (s, y)."""
    ux, uz = direction
    dx, dz = cent[:, 0] - axis_x, cent[:, 2] - axis_z
    s = dx * ux + dz * uz
    off = -dx * uz + dz * ux
    m = (np.abs(off) < half_width) & (np.abs(s) < r_max) & (np.abs(cent[:, 1]) < r_max)
    Bs = B[m, 0] * ux + B[m, 2] * uz
    return np.stack([s[m], cent[m, 1], Bs, B[m, 1], np.linalg.norm(B[m], axis=1)], axis=1)


def extract(run_dir, geo, excitation):
    params = parse_params(run_dir / "assembly_params.txt")
    ax, az = params["axis_x"], params["axis_z"]

    nodes, B = parse_pos_elements(run_dir / "b_assembly.pos")
    cent = nodes.mean(axis=1)
    inodes, Bi = parse_pos_elements(run_dir / "b_iron.pos")
    icent = inodes.mean(axis=1)
    del nodes, inodes

    # --- channel -----------------------------------------------------------
    prof = axial_profile(cent, B, ax, az, CHANNEL["channel_inner_r"], CHANNEL["channel_outer_r"],
                         CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"], n_stations=N_STATIONS)
    B_r_G = prof["B_r"] * 1e4
    B_exit_G = float(B_r_G[-1])
    scale = TARGET_G / B_exit_G       # linear solve: everything scales with A-turns

    scaled_prof = dict(prof, B_r=prof["B_r"] * 1e4 * scale, Bmag=prof["Bmag"] * 1e4 * scale)
    fq = evaluate_field_quality(scaled_prof, B_peak=TARGET_G)

    # --- iron --------------------------------------------------------------
    h = geo["emag_height"] * IN
    y_coil_top = Y_TOP_PLATE_INNER + h
    pad = 0.002
    inner_r = geo["inner_coil_id"] * IN / 2
    outer_r = geo["outer_coil_id"] * IN / 2
    inner_shaft = shaft_profile(icent, Bi, ax, az, inner_r, Y_TOP_PLATE_INNER + pad, y_coil_top - pad)
    outer_shafts = []
    for w in params["windings"]:
        if w["pole"] == "outer_coil":
            cx, cz = w["loc"]
            outer_shafts.append(shaft_profile(icent, Bi, cx, cz, outer_r,
                                              Y_TOP_PLATE_INNER + pad, y_coil_top - pad).max())
    Bi_mag = np.linalg.norm(Bi, axis=1)

    # --- slices for figures ------------------------------------------------
    wloc =[w["loc"] for w in params["windings"] if w["pole"] == "outer_coil"]
    # direction toward whichever outer core is in the +x+z quadrant
    q = max(wloc, key=lambda l: l[0] + l[1])
    n = math.hypot(q[0] - ax, q[1] - az)
    diag = ((q[0] - ax) / n, (q[1] - az) / n)
    between = (diag[0] - diag[1], diag[0] + diag[1])
    nb = math.hypot(*between)
    between = (between[0] / nb, between[1] / nb)

    ch = ((np.hypot(cent[:, 0] - ax, cent[:, 2] - az) >= CHANNEL["channel_inner_r"]) &
          (np.hypot(cent[:, 0] - ax, cent[:, 2] - az) <= CHANNEL["channel_outer_r"]) &
          (cent[:, 1] >= CHANNEL["channel_y_exit"]) & (cent[:, 1] <= CHANNEL["channel_y_anode"]))
    np.savez_compressed(
        run_dir / "extract.npz",
        scale=scale,
        prof_y=prof["y"], prof_Br_G=B_r_G, prof_Bax_G=prof["B_axial"] * 1e4,
        prof_Bmag_G=prof["Bmag"] * 1e4, prof_counts=prof["counts"],
        channel_cent=cent[ch], channel_B=B[ch],
        slice_core=slice_plane(cent, B, ax, az, diag),
        slice_between=slice_plane(cent, B, ax, az, between),
        axis=np.array([ax, az]),
    )

    # envelope area for the thermal model, from this geometry's own mesh
    from radiation_balance import envelope_area
    area, _ = envelope_area(str(run_dir / "assembly.msh"))

    inner_at = excitation["inner_coil"]["turns"] * excitation["inner_coil"]["current"]
    outer_at = excitation["outer_coil"]["turns"] * excitation["outer_coil"]["current"]
    m = dict(
        tag=run_dir.name,
        emag_height_in=geo["emag_height"],
        spacer_in=round(geo["emag_height"] - SPACER_ZERO_IN, 4),
        inner_core_dia_in=geo["inner_coil_id"],
        outer_coil_od_in=geo["outer_coil_od"],
        base_inner_At=inner_at, base_outer_At=outer_at,
        B_exit_G_at_base=B_exit_G,
        scale_to_target=scale,
        req_inner_At=inner_at * scale,
        req_outer_At=outer_at * scale,
        B_peak_G=float(np.nanmax(B_r_G) * scale),
        y_peak_mm=float(prof["y"][np.nanargmax(B_r_G)] * 1e3),
        B_anode_G=float(B_r_G[0] * scale),
        radial_purity_pct=fq["radial_purity_pct"],
        shape_rmse_pct=fq["shape_rmse_pct"],
        monotonic_violation_pct=fq["monotonicity_violation_pct"],
        quality_score=fq["composite_score"],
        inner_shaft_B_max_T=float(inner_shaft.max() * scale),
        inner_shaft_B_mid_T=float(inner_shaft[len(inner_shaft) // 2] * scale),
        outer_shaft_B_max_T=float(max(outer_shafts) * scale),
        iron_B_p99_T=float(np.percentile(Bi_mag, 99) * scale),
        iron_B_p999_T=float(np.percentile(Bi_mag, 99.9) * scale),
        iron_B_max_T=float(Bi_mag.max() * scale),
        envelope_area_m2=area,
        n_channel_elems=int(ch.sum()),
        min_station_count=int(prof["counts"].min()),
    )
    m["linear_valid"] = bool(m["inner_shaft_B_max_T"] < KNEE_T and m["outer_shaft_B_max_T"] < KNEE_T)
    return m


def run_point(h, ci, oo, force=False, iron_mur=IRON_MUR, keep_fields=False, suffix=""):
    geo = geometry(h, ci, oo)
    run_dir = RUNS_DIR / (tag_for(h, ci, oo) + suffix)
    mpath = run_dir / "metrics.json"
    if mpath.exists() and not force:
        return json.loads(mpath.read_text())
    run_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    build(geo, str(run_dir), None, include_chamber_spacer_as_iron=True)
    rewrite_channel(run_dir / "assembly_params.txt")
    params = parse_params(run_dir / "assembly_params.txt")
    generate_regions_pro(params, EXCITATION, out_path=str(run_dir / "assembly_regions_generated.pro"),
                         iron_mur=iron_mur)
    t_build = time.time() - t0
    t_solve = run_getdp(run_dir)
    m = extract(run_dir, geo, EXCITATION)
    m.update(build_s=round(t_build, 1), solve_s=round(t_solve, 1), iron_mur=iron_mur)
    mpath.write_text(json.dumps(m, indent=2))
    if not keep_fields:
        for f in ("b_assembly.pos", "b_iron.pos", "assembly.msh", "magnetostatics_assembly.res",
                  "magnetostatics_assembly.pre"):
            (run_dir / f).unlink(missing_ok=True)
    print(f"[{run_dir.name}] build {t_build:.0f}s solve {t_solve:.0f}s  "
          f"B_exit(base)={m['B_exit_G_at_base']:.1f} G -> x{m['scale_to_target']:.3f}  "
          f"inner shaft {m['inner_shaft_B_max_T']:.2f} T", flush=True)
    return m


def write_csv(rows):
    keys = list(rows[0].keys())
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["outer_coil_od_in"], r["inner_core_dia_in"], r["emag_height_in"])):
            w.writerow(r)


def validate():
    """Baseline geometry, linear iron, compared on the same metrics to the
    committed nonlinear solve (b_assembly.pos/b_iron.pos in this directory,
    1.30 A, A36-proxy curve). Linear with mu_r=IRON_MUR is only fit for
    the sweep if the two agree closely."""
    m = run_point(2.0, 1.0, 1.5, suffix="_validate", keep_fields=False, force=True)
    params = parse_params(HERE / "assembly_params.txt")
    nodes, B = parse_pos_elements(HERE / "b_assembly.pos")
    cent = nodes.mean(axis=1)
    prof = axial_profile(cent, B, params["axis_x"], params["axis_z"],
                         CHANNEL["channel_inner_r"], CHANNEL["channel_outer_r"],
                         CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"], n_stations=N_STATIONS)
    nl_exit = prof["B_r"][-1] * 1e4
    print(f"\nnonlinear (committed) B_r exit = {nl_exit:.1f} G at base excitation")
    print(f"linear mu_r={IRON_MUR:.0f}     B_r exit = {m['B_exit_G_at_base']:.1f} G "
          f"({(m['B_exit_G_at_base'] / nl_exit - 1) * 100:+.2f}%)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="run a single tag, e.g. h2.000_ci1.000_oo1.500")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--nonlinear", metavar="TAG",
                    help="re-solve one swept geometry with the nonlinear B-H iron and compare to its "
                         "linear result, e.g. h1.125_ci0.750_oo1.500 (the highest-core-flux point)")
    args = ap.parse_args()

    if args.validate:
        validate()
        return
    if args.nonlinear:
        h, ci, oo = (float(p[1:] if p[0] == "h" else p[2:]) for p in args.nonlinear.split("_"))
        lin = run_point(h, ci, oo)
        nl = run_point(h, ci, oo, iron_mur=None, suffix="_nonlinear", force=args.force)
        for k in ("B_exit_G_at_base", "inner_shaft_B_max_T", "outer_shaft_B_max_T", "iron_B_p999_T",
                  "radial_purity_pct"):
            # both extracted at the SAME (baseline) excitation before scaling, so
            # compare the unscaled exit field and the scaled-by-own-factor iron fields
            print(f"{k:24s} linear {lin[k]:10.4f}  nonlinear {nl[k]:10.4f}  "
                  f"({(nl[k] / lin[k] - 1) * 100:+.2f}%)")
        return

    points = [(h, ci, OUTER_OD_IN) for ci in CORE_DIAS_IN for h in HEIGHTS_IN] + EXTRA_POINTS
    if args.only:
        points = [p for p in points if tag_for(*p) == args.only]
    rows = []
    for p in points:
        try:
            rows.append(run_point(*p, force=args.force))
        except Exception as e:  # keep going; one bad geometry shouldn't sink the grid
            print(f"[{tag_for(*p)}] FAILED: {e}", flush=True)
            # A failure mid-build leaves gmsh initialised with a half-built
            # model; the next build then silently writes an empty mesh.
            import gmsh
            if gmsh.isInitialized():
                gmsh.finalize()
            shutil.rmtree(RUNS_DIR / tag_for(*p), ignore_errors=True)
        done = [json.loads(f.read_text()) for f in RUNS_DIR.glob("h*_ci*_oo*/metrics.json")
                if not f.parent.name.endswith(("_validate", "_nonlinear"))]
        if done:
            write_csv(done)
    print(f"\nWrote {RESULTS_CSV}")


if __name__ == "__main__":
    main()
