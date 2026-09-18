"""
current_sweep.py

Wrapper that sweeps the center-solenoid / outer-solenoid current ratio and
reports, for each ratio, how "horizontal" (radial, i.e. perpendicular to
the thruster centerline) the B field is inside the discharge channel --
the standard Hall-thruster requirement for E x B operation.

Key physics shortcut: mur_iron=1000 is a LINEAR material (no B-H
saturation curve -- see HOW_IT_WORKS.md sec. 9), so the magnetostatics
problem is linear in the winding currents and solutions superpose. That
means the sweep does NOT need one GetDP solve per ratio: it only needs two
"basis" solves --

  B_center_unit(x) = field with 1 A in the center winding, 0 A in outer
  B_outer_unit(x)  = field with 1 A in each outer winding, 0 A in center

-- and then, for ANY (I_center, I_outer_each):

  B(x) = I_center * B_center_unit(x) + I_outer_each * B_outer_unit(x)

computed instantly in Python. Nothing about the mesh, formulation, or
.pro files changes -- this only regenerates assembly_regions_generated.pro
(via generate_regions.py) between the two basis solves and reuses the
existing assembly.msh both times.

Usage (from gmsh_getdp/assembly/):
  python current_sweep.py                          # full run: 2 solves + sweep
  python current_sweep.py --skip-solve              # reuse channel_basis.npz, just re-sweep/replot
  python current_sweep.py --ratios 0.5,1,2,4        # explicit ratio list
  python current_sweep.py --outer-current 5.0       # reference outer current (A) used for magnitude reporting
"""
import argparse
import os
import subprocess
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from params_utils import parse_params
from generate_regions import generate_regions_pro
from pos_utils import parse_pos_elements

GETDP_EXE = os.path.join("..", "..", "tools", "getdp-3.5.0-Windows64", "getdp.exe")
PRO_FILE = "magnetostatics_assembly.pro"
BASIS_NPZ = "channel_basis.npz"

BASIS_EXCITATION = {
    "center": {
        "inner_coil": dict(turns=300, current=1.0, polarity=-1),
        "outer_coil": dict(turns=200, current=0.0, polarity=+1),
    },
    "outer": {
        "inner_coil": dict(turns=300, current=0.0, polarity=-1),
        "outer_coil": dict(turns=200, current=1.0, polarity=+1),
    },
}


def run_basis_solve(name, params):
    """Regenerate assembly_regions_generated.pro for one basis excitation,
    run GetDP against the existing mesh, and stash the resulting b_assembly.pos
    under a basis-specific name (the .pro file always writes the same fixed
    filename, so this has to happen between runs to avoid overwriting)."""
    print(f"\n=== Basis solve: {name} ===")
    generate_regions_pro(params, BASIS_EXCITATION[name])

    t0 = time.time()
    result = subprocess.run(
        [GETDP_EXE, PRO_FILE, "-msh", "assembly.msh", "-solve", "Res_a", "-pos", "Map_b"],
        capture_output=True, text=True,
    )
    dt = time.time() - t0
    print(f"getdp.exe exit={result.returncode} ({dt:.1f}s)")
    if result.returncode != 0:
        print(result.stdout[-3000:])
        print(result.stderr[-3000:])
        raise RuntimeError(f"GetDP solve failed for basis '{name}'")

    dest = f"b_assembly_{name}.pos"
    os.replace("b_assembly.pos", dest)
    if os.path.exists("b_iron.pos"):
        os.remove("b_iron.pos")  # not needed for the channel sweep
    return dest


def channel_mask(centroids, params):
    x, y, z = centroids[:, 0], centroids[:, 1], centroids[:, 2]
    r = np.hypot(x - params["axis_x"], z - params["axis_z"])
    return (
        (r >= params["channel_inner_r"]) & (r <= params["channel_outer_r"]) &
        (y >= params["channel_y_min"]) & (y <= params["channel_y_max"])
    )


def tet_volumes(nodes):
    """nodes: (N,4,3) corner coordinates -> (N,) tetrahedron volumes."""
    a = nodes[:, 1] - nodes[:, 0]
    b = nodes[:, 2] - nodes[:, 0]
    c = nodes[:, 3] - nodes[:, 0]
    return np.abs(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0


def extract_channel_basis(params):
    """Parse both basis .pos files, keep only channel elements + their
    per-amp B, save a small .npz, and delete the (huge) full-domain .pos
    files afterward. Returns the loaded arrays."""
    data = {}
    for name in ("center", "outer"):
        path = f"b_assembly_{name}.pos"
        print(f"Parsing {path} ...")
        nodes, B = parse_pos_elements(path)
        centroids = nodes.mean(axis=1)
        mask = channel_mask(centroids, params)
        print(f"  {len(centroids)} elements total, {mask.sum()} in channel")
        if name == "center":
            data["centroids"] = centroids[mask]
            data["vol"] = tet_volumes(nodes[mask])
            data["B_center_unit"] = B[mask]
        else:
            # sanity check: same mesh -> same element ordering/centroids
            assert np.allclose(centroids[mask], data["centroids"], atol=1e-9), \
                "center/outer basis element ordering mismatch -- did the mesh change between solves?"
            data["B_outer_unit"] = B[mask]
        os.remove(path)  # ~200MB full-domain dump; channel subset is all we need

    np.savez(BASIS_NPZ, **data)
    print(f"Wrote {BASIS_NPZ} ({len(data['centroids'])} channel elements)")
    return data


def sweep(data, ratios, outer_current_ref):
    centroids = data["centroids"]
    vol = data["vol"]
    Bc = data["B_center_unit"]
    Bo = data["B_outer_unit"]
    total_vol = vol.sum()

    rows = []
    for r in ratios:
        I_center = r * outer_current_ref
        B = I_center * Bc + outer_current_ref * Bo
        Bmag = np.linalg.norm(B, axis=1)
        axial_frac = np.abs(B[:, 1]) / np.clip(Bmag, 1e-30, None)
        axial_frac_mean = np.sum(axial_frac * vol) / total_vol
        Bmag_mean = np.sum(Bmag * vol) / total_vol
        rows.append(dict(
            ratio=r, I_center=I_center, I_outer=outer_current_ref,
            axial_frac_mean=axial_frac_mean,
            Bmag_mean=Bmag_mean, Bmag_median=float(np.median(Bmag)),
            Bmag_p90=float(np.percentile(Bmag, 90)),
        ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-solve", action="store_true",
                     help="Reuse existing channel_basis.npz instead of re-running GetDP.")
    ap.add_argument("--ratios", type=str, default=None,
                     help="Comma-separated I_center/I_outer ratios. Default: 15 log-spaced points in [0.2, 5].")
    ap.add_argument("--outer-current", type=float, default=5.0,
                     help="Reference per-coil outer-winding current (A) used for magnitude reporting. "
                          "Doesn't affect the horizontality metric (that's scale-invariant), only Bmag columns.")
    args = ap.parse_args()

    params = parse_params()

    if args.skip_solve:
        if not os.path.exists(BASIS_NPZ):
            raise SystemExit(f"{BASIS_NPZ} not found -- run without --skip-solve first.")
        print(f"Loading cached {BASIS_NPZ}")
        data = dict(np.load(BASIS_NPZ))
    else:
        if not os.path.exists("assembly.msh"):
            raise SystemExit("assembly.msh not found -- run build_assembly.py first (mesh generation, one-time).")
        run_basis_solve("center", params)
        run_basis_solve("outer", params)
        data = extract_channel_basis(params)

    if args.ratios:
        ratios = [float(x) for x in args.ratios.split(",")]
    else:
        ratios = np.logspace(np.log10(0.2), np.log10(5.0), 15)

    rows = sweep(data, ratios, args.outer_current)

    with open("sweep_results.csv", "w") as f:
        cols = list(rows[0].keys())
        f.write(",".join(cols) + "\n")
        for row in rows:
            f.write(",".join(f"{row[c]:.6g}" for c in cols) + "\n")
    print("Wrote sweep_results.csv")

    best = min(rows, key=lambda row: row["axial_frac_mean"])
    print(f"\nMost horizontal (radial-dominant) field at ratio I_center/I_outer = {best['ratio']:.4g} "
          f"(I_center={best['I_center']:.3g} A, I_outer={best['I_outer']:.3g} A each)")
    print(f"  channel axial-fraction (mean |B_y|/|B|, volume-weighted) = {best['axial_frac_mean']:.4f}")
    print(f"  channel |B| mean={best['Bmag_mean']*1e3:.3f} mT  median={best['Bmag_median']*1e3:.3f} mT  "
          f"p90={best['Bmag_p90']*1e3:.3f} mT")
    print("  -> to hit a target channel |B|, scale BOTH currents by (target / channel |B| mean above), "
          "keeping this ratio fixed.")

    ratios_arr = np.array([row["ratio"] for row in rows])
    axial_arr = np.array([row["axial_frac_mean"] for row in rows])
    bmag_arr = np.array([row["Bmag_mean"] for row in rows]) * 1e3

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    ax1.plot(ratios_arr, axial_arr, "o-")
    ax1.axvline(best["ratio"], color="r", ls="--", alpha=0.6, label=f"best r={best['ratio']:.3g}")
    ax1.set_ylabel("channel axial fraction  mean(|B_y|/|B|)\n(lower = more radial/horizontal)")
    ax1.set_xscale("log")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.plot(ratios_arr, bmag_arr, "o-", color="darkorange")
    ax2.axvline(best["ratio"], color="r", ls="--", alpha=0.6)
    ax2.set_xlabel("I_center / I_outer")
    ax2.set_ylabel(f"channel |B| mean (mT)\n(at I_outer={args.outer_current:g} A each)")
    ax2.set_xscale("log")
    ax2.grid(alpha=0.3)

    fig.suptitle("BPL-700 channel field vs. center/outer current ratio")
    fig.tight_layout()
    fig.savefig("sweep_ratio.png", dpi=150)
    print("Wrote sweep_ratio.png")


if __name__ == "__main__":
    main()
