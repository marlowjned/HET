"""
em_sizing_sweep.py

BPL-1000 electromagnet sizing: core diameter x winding build (coilt) for the
inner and outer electromagnets separately, plus the inner pole-plate
thickness trade.

Why two basis solves per geometry
---------------------------------
Iron is linear here (mu_r = IRON_MUR; emag_sweep.py --validate justified it
on BPL-700: the circuit is gap-dominated and the iron runs below the knee),
so the field is linear in the two ampere-turn totals:

    B(x) = NI_inner * b_inner(x) + NI_outer * b_outer(x)

where b_inner / b_outer are the fields for 1 A-turn in the inner coil alone /
in each outer coil alone. Two solves on one mesh give the whole family, so
for every geometry this:
  1. scans the outer/inner A-turn ratio q,
  2. at each q scales to TARGET_G radial at the channel exit,
  3. scores the shape with field_quality.evaluate_field_quality,
and reports the best-scoring ratio AND the fixed reference ratio 214.5/338,
so geometries can be compared both ways. That removes the "estimate turns,
then simulate" loop entirely: the required A-turns come out of the solve.

Wire gauge never enters the FEM (it sees A-turns spread over the coil
window), so gauge is post-processing -- see sizing/em_winding_trade.py.

Runs live in sweep_runs_bpl1000/<tag>/; .pos dumps and the mesh are deleted
after extraction, metrics.json + basis.npz are kept. Finished points are
skipped on re-run.

Usage (from gmsh_getdp/assembly/):
  python em_sizing_sweep.py --grid inner|outer|plate|all
  python em_sizing_sweep.py --check-superposition   # baseline: basis sum vs direct solve
  python em_sizing_sweep.py --list --grid all        # print the grid, run nothing
"""
import argparse
import csv
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "thermal"))

from bpl1000_geometry import (CHANNEL, IN, INNER_OD_MAX_IN, OUTER_OD_MAX_IN,  # noqa: E402
                              coil_top_y, inner_coil_bottom_y, Y_FRONT_PLATE_INNER,
                              CORE_BOLT_HOLE_R)
from build_assembly import build, EXCITATION  # noqa: E402
from emag_sweep import run_getdp, rewrite_channel as _rewrite_channel_bpl700  # noqa: E402,F401
from field_quality import axial_profile, evaluate_field_quality  # noqa: E402
from generate_regions import generate_regions_pro  # noqa: E402
from onshape.cache import DEFAULTS_IN  # noqa: E402
from params_utils import parse_params  # noqa: E402
from pos_utils import parse_pos_elements  # noqa: E402

RUNS_DIR = HERE / "sweep_runs_bpl1000"
RESULTS_CSV = HERE / "em_sizing_results.csv"

IRON_MUR = 1500.0
TARGET_G = 300.0
KNEE_T = 1.5
N_STATIONS = 20
REF_RATIO = EXCITATION["outer_coil"]["turns"] * EXCITATION["outer_coil"]["current"] / (
    EXCITATION["inner_coil"]["turns"] * EXCITATION["inner_coil"]["current"])   # 214.5/338
RATIOS = np.round(np.arange(0.20, 1.801, 0.01), 3)

# 1 A-turn on one pole, 0 on the other. Polarity as in EXCITATION.
BASIS = {
    "inner": {"inner_coil": dict(turns=1, current=1.0, polarity=EXCITATION["inner_coil"]["polarity"]),
              "outer_coil": dict(turns=1, current=0.0, polarity=EXCITATION["outer_coil"]["polarity"])},
    "outer": {"inner_coil": dict(turns=1, current=0.0, polarity=EXCITATION["inner_coil"]["polarity"]),
              "outer_coil": dict(turns=1, current=1.0, polarity=EXCITATION["outer_coil"]["polarity"])},
}

# Fit margin left inside the CAD limits (radial clearance for ground wrap,
# tape and assembly), in inches of diameter.
FIT_MARGIN_IN = 0.010


def _fits(cfg):
    return (cfg["iem_coreod"] + cfg["iem_coilt"] <= INNER_OD_MAX_IN - FIT_MARGIN_IN and
            cfg["oem_coreod"] + cfg["oem_coilt"] <= OUTER_OD_MAX_IN - FIT_MARGIN_IN)


def grid(which):
    pts = []
    if which in ("inner", "all"):
        for ic in (0.5, 0.625, 0.75, 0.875, 1.0, 1.125):
            for it in (0.125, 0.25, 0.375):
                pts.append(dict(DEFAULTS_IN, iem_coreod=ic, iem_coilt=it))
    if which in ("outer", "all"):
        for oc in (0.5, 0.625, 0.75, 0.875, 1.0, 1.125):
            for ot in (0.25, 0.5, 0.75):
                pts.append(dict(DEFAULTS_IN, oem_coreod=oc, oem_coilt=ot))
    if which in ("plate", "all"):
        for p in (0.125, 0.25, 0.5):      # 0.375 is the default = baseline
            pts.append(dict(DEFAULTS_IN, iem_platet=p))
    out, seen = [], set()
    for c in pts:
        t = tag_for(c)
        if t not in seen and _fits(c):
            seen.add(t)
            out.append(c)
    return out


def tag_for(c):
    return (f"ic{c['iem_coreod']:.3f}_it{c['iem_coilt']:.3f}_ip{c['iem_platet']:.3f}"
            f"_oc{c['oem_coreod']:.3f}_ot{c['oem_coilt']:.3f}_h{c['emag_height']:.3f}")


def rewrite_channel(params_path):
    """build_assembly's detect_channel_cavity estimate -> the measured CAD
    channel (bpl1000_geometry.CHANNEL)."""
    lines = Path(params_path).read_text().splitlines()
    vals = dict(CHANNEL,
                channel_y_min=min(CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"]),
                channel_y_max=max(CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"]))
    out = []
    for ln in lines:
        key = ln.split(" ", 1)[0]
        out.append(f"{key} {vals[key]!r}" if key in vals else ln)
    Path(params_path).write_text("\n".join(out) + "\n")


def shaft_signed(cent, B, cx, cz, r_core, y_lo, y_hi, n_bins=12):
    """Signed mean B_y over a core's shaft cross-section per axial bin.
    Signed (unlike emag_sweep.shaft_profile) so two basis profiles can be
    superposed before taking the magnitude. Bolt-hole air has no iron
    elements, so it is excluded automatically; the r < 0.9 r_core cut keeps
    the coil-side skin out."""
    r = np.hypot(cent[:, 0] - cx, cent[:, 2] - cz)
    y = cent[:, 1]
    sel = (r < 0.9 * r_core) & (y > y_lo) & (y < y_hi)
    edges = np.linspace(y_lo, y_hi, n_bins + 1)
    out = np.full(n_bins, np.nan)
    for i in range(n_bins):
        m = sel & (y >= edges[i]) & (y < edges[i + 1])
        if m.sum() >= 5:
            out[i] = B[m, 1].mean()
    return out


def solve_basis(run_dir, params, name):
    generate_regions_pro(params, BASIS[name], out_path=str(run_dir / "assembly_regions_generated.pro"),
                         iron_mur=IRON_MUR)
    t = run_getdp(run_dir)
    nodes, B = parse_pos_elements(run_dir / "b_assembly.pos")
    cent = nodes.mean(axis=1)
    inodes, Bi = parse_pos_elements(run_dir / "b_iron.pos")
    icent = inodes.mean(axis=1)
    return cent, B, icent, Bi, t


def evaluate(ch_cent, ch_B, ax, az, q):
    """Per unit inner A-turn: field for ratio q, its exit value, and the
    scale that brings it to TARGET_G."""
    B = ch_B[0] + q * ch_B[1]
    prof = axial_profile(ch_cent, B, ax, az, CHANNEL["channel_inner_r"], CHANNEL["channel_outer_r"],
                         CHANNEL["channel_y_anode"], CHANNEL["channel_y_exit"], n_stations=N_STATIONS)
    B_exit_G = float(prof["B_r"][-1] * 1e4)
    s = TARGET_G / B_exit_G
    sp = dict(prof, B_r=prof["B_r"] * 1e4 * s, Bmag=prof["Bmag"] * 1e4 * s)
    fq = evaluate_field_quality(sp, B_peak=TARGET_G)
    return dict(prof=prof, scale=s, fq=fq, Br_G=sp["B_r"])


def run_point(cfg, force=False):
    run_dir = RUNS_DIR / tag_for(cfg)
    mpath = run_dir / "metrics.json"
    if mpath.exists() and not force:
        return json.loads(mpath.read_text())
    run_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    build(cfg, str(run_dir), None, include_chamber_spacer_as_iron=True)
    rewrite_channel(run_dir / "assembly_params.txt")
    params = parse_params(run_dir / "assembly_params.txt")
    t_build = time.time() - t0
    ax, az = params["axis_x"], params["axis_z"]

    basis, t_solve = {}, 0.0
    for name in ("inner", "outer"):
        cent, B, icent, Bi, t = solve_basis(run_dir, params, name)
        basis[name] = (cent, B, icent, Bi)
        t_solve += t
    cent, icent = basis["inner"][0], basis["inner"][2]
    assert np.allclose(cent, basis["outer"][0]) and np.allclose(icent, basis["outer"][2]), \
        "basis element ordering differs -- the mesh changed between solves"

    r = np.hypot(cent[:, 0] - ax, cent[:, 2] - az)
    ch = ((r >= CHANNEL["channel_inner_r"]) & (r <= CHANNEL["channel_outer_r"]) &
          (cent[:, 1] >= CHANNEL["channel_y_exit"]) & (cent[:, 1] <= CHANNEL["channel_y_anode"]))
    ch_cent = cent[ch]
    ch_B = (basis["inner"][1][ch], basis["outer"][1][ch])

    # --- shafts, per basis (signed, superposable) ---------------------------
    y_top = coil_top_y(cfg["emag_height"])
    pad = 0.002
    ri, ro = cfg["iem_coreod"] * IN / 2, cfg["oem_coreod"] * IN / 2
    yi_lo = inner_coil_bottom_y(cfg["iem_platet"])
    sh_in = [shaft_signed(icent, basis[n][3], ax, az, ri, yi_lo + pad, y_top - pad) for n in ("inner", "outer")]
    outer_locs = [w["loc"] for w in params["windings"] if w["pole"] == "outer_coil"]
    sh_out = [[shaft_signed(icent, basis[n][3], cx, cz, ro, Y_FRONT_PLATE_INNER + pad, y_top - pad)
               for n in ("inner", "outer")] for cx, cz in outer_locs]
    Bi_pair = (basis["inner"][3], basis["outer"][3])

    def at_ratio(q):
        e = evaluate(ch_cent, ch_B, ax, az, q)
        s = e["scale"]
        NI_in, NI_out = s, q * s
        inner_sh = np.nanmax(np.abs(NI_in * sh_in[0] + NI_out * sh_in[1]))
        outer_sh = max(np.nanmax(np.abs(NI_in * a + NI_out * b)) for a, b in sh_out)
        Bim = np.linalg.norm(NI_in * Bi_pair[0] + NI_out * Bi_pair[1], axis=1)
        fq = e["fq"]
        return dict(
            ratio=float(q), req_inner_At=float(NI_in), req_outer_At=float(NI_out),
            B_peak_G=float(np.nanmax(e["Br_G"])),
            y_peak_mm=float(e["prof"]["y"][np.nanargmax(e["Br_G"])] * 1e3),
            B_anode_G=float(e["Br_G"][0]),
            radial_purity_pct=fq["radial_purity_pct"], shape_rmse_pct=fq["shape_rmse_pct"],
            monotonic_violation_pct=fq["monotonicity_violation_pct"],
            quality_score=fq["composite_score"],
            inner_shaft_B_max_T=float(inner_sh), outer_shaft_B_max_T=float(outer_sh),
            iron_B_p999_T=float(np.percentile(Bim, 99.9)), iron_B_max_T=float(Bim.max()),
        ), e

    scores = []
    for q in RATIOS:
        e = evaluate(ch_cent, ch_B, ax, az, q)
        scores.append(e["fq"]["composite_score"])
    q_best = float(RATIOS[int(np.nanargmax(scores))])
    best, e_best = at_ratio(q_best)
    ref, e_ref = at_ratio(REF_RATIO)

    np.savez_compressed(run_dir / "basis.npz", ch_cent=ch_cent, ch_B_inner=ch_B[0], ch_B_outer=ch_B[1],
                        axis=np.array([ax, az]), ratios=RATIOS, scores=np.array(scores),
                        prof_y=e_best["prof"]["y"], Br_best_G=e_best["Br_G"], Br_ref_G=e_ref["Br_G"])

    from radiation_balance import envelope_area
    area, _ = envelope_area(str(run_dir / "assembly.msh"))

    m = dict(tag=run_dir.name, **{k: cfg[k] for k in cfg},
             inner_em_od_in=cfg["iem_coreod"] + cfg["iem_coilt"],
             outer_em_od_in=cfg["oem_coreod"] + cfg["oem_coilt"],
             envelope_area_m2=area, n_channel_elems=int(ch.sum()),
             build_s=round(t_build, 1), solve_s=round(t_solve, 1), iron_mur=IRON_MUR)
    m.update({f"best_{k}": v for k, v in best.items()})
    m.update({f"ref_{k}": v for k, v in ref.items()})
    m["linear_valid"] = bool(max(best["inner_shaft_B_max_T"], best["outer_shaft_B_max_T"],
                                 ref["inner_shaft_B_max_T"], ref["outer_shaft_B_max_T"]) < KNEE_T)
    mpath.write_text(json.dumps(m, indent=2))
    for f in ("b_assembly.pos", "b_iron.pos", "assembly.msh", "magnetostatics_assembly.res",
              "magnetostatics_assembly.pre"):
        (run_dir / f).unlink(missing_ok=True)
    print(f"[{run_dir.name}] build {t_build:.0f}s solve {t_solve:.0f}s  best q={q_best:.2f} "
          f"-> {best['req_inner_At']:.0f}/{best['req_outer_At']:.0f} At, score {best['quality_score']:.2f}; "
          f"ref -> {ref['req_inner_At']:.0f}/{ref['req_outer_At']:.0f} At; shafts "
          f"{best['inner_shaft_B_max_T']:.2f}/{best['outer_shaft_B_max_T']:.2f} T", flush=True)
    return m


def write_csv(rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["tag"]):
            w.writerow(r)


def check_superposition():
    """Baseline geometry: basis sum at the committed EXCITATION vs a direct
    solve at that excitation. Linear iron, so these must agree to solver
    tolerance; anything more means the basis machinery is wrong."""
    cfg = dict(DEFAULTS_IN)
    run_dir = RUNS_DIR / "superposition_check"
    run_dir.mkdir(parents=True, exist_ok=True)
    build(cfg, str(run_dir), None, include_chamber_spacer_as_iron=True)
    rewrite_channel(run_dir / "assembly_params.txt")
    params = parse_params(run_dir / "assembly_params.txt")
    NI_in = EXCITATION["inner_coil"]["turns"] * EXCITATION["inner_coil"]["current"]
    NI_out = EXCITATION["outer_coil"]["turns"] * EXCITATION["outer_coil"]["current"]
    bi = solve_basis(run_dir, params, "inner")
    bo = solve_basis(run_dir, params, "outer")
    generate_regions_pro(params, EXCITATION, out_path=str(run_dir / "assembly_regions_generated.pro"),
                         iron_mur=IRON_MUR)
    run_getdp(run_dir)
    _, Bd = parse_pos_elements(run_dir / "b_assembly.pos")
    Bs = NI_in * bi[1] + NI_out * bo[1]
    err = np.linalg.norm(Bs - Bd, axis=1).max() / np.linalg.norm(Bd, axis=1).max()
    print(f"superposition: max |B_sum - B_direct| / max|B_direct| = {err:.2e}")
    for f in ("b_assembly.pos", "b_iron.pos", "magnetostatics_assembly.res", "magnetostatics_assembly.pre"):
        (run_dir / f).unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="all", choices=["inner", "outer", "plate", "all", "baseline"])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--check-superposition", action="store_true")
    args = ap.parse_args()

    if args.check_superposition:
        check_superposition()
        return
    pts = [dict(DEFAULTS_IN)] + ([] if args.grid == "baseline" else grid(args.grid))
    pts = list({tag_for(c): c for c in pts}.values())   # baseline is also in the inner grid
    if args.list:
        for c in pts:
            print(tag_for(c), f"inner OD {c['iem_coreod'] + c['iem_coilt']:.3f}",
                  f"outer OD {c['oem_coreod'] + c['oem_coilt']:.3f}")
        print(f"{len(pts)} points")
        return
    for c in pts:
        try:
            run_point(c, force=args.force)
        except Exception as e:
            print(f"[{tag_for(c)}] FAILED: {e}", flush=True)
            import gmsh
            if gmsh.isInitialized():
                gmsh.finalize()
        done = [json.loads(f.read_text()) for f in RUNS_DIR.glob("ic*/metrics.json")]
        if done:
            write_csv(done)
    print(f"\nWrote {RESULTS_CSV}")


if __name__ == "__main__":
    main()
