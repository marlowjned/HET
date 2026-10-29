"""Vacuum electrostatic field of Design B at the 300 V discharge, 2-D axisymmetric (gmsh + GetDP; no bplsim).

Anode at V_d; the chamber wall (far boundary) at 0 V stands in for the cathode / ground reference. Iron and coils
are one conductor, floating (net charge 0) or grounded. Two per-volt bases give every case by superposition:
  phi_A: anode 1 V, iron 0 V            phi_I: iron 1 V, anode 0 V            (chamber wall 0 V in both)
  grounded iron: phi = V_d phi_A        floating iron: phi = V_d (phi_A + k phi_I),  k = -C_IA / C_II
A direct floating solve checks the superposition. Writes figures, out/fields_rz.npz (the per-volt bases on a regular
r-z grid; fields.py reads it) and out/summary.json.

  python run.py            (about a minute; --quick skips the chamber-size and mesh-convergence studies)
"""
import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from geometry import design_b as g  # noqa: E402
import mesh_rz  # noqa: E402

GETDP = shutil.which("getdp") or str(Path.home() / "tools/getdp/bin/getdp")
OUT = HERE / "out"
V_D = 300.0                  # discharge voltage [V]
EPS_R_BN = 4.0               # PLACEHOLDER
REGION = {"vacuum": 0, "bn": 1, "anode": 2, "iron": 3}     # region codes in the exported grid ('coil' -> iron)
GRID = dict(r=(0.0, 120.0), z=(-60.0, 150.0), h=0.25)       # export grid [mm]
COLORS = {"bn": "#cfe3f7", "anode": "#d95f02", "iron": "#9a9a9a", "coil": "#e6b35c"}


# ---------------------------------------------------------------------------------------------- solving
def solve(work, msh, **numbers):
    """Run estatic.pro in `work`; returns the global quantities (charges in C, potentials in V, energy in J)."""
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(HERE / "estatic.pro", work / "estatic.pro")         # GetDP writes its outputs next to the .pro
    numbers.setdefault("eps_r_bn", EPS_R_BN)
    args = [GETDP, "estatic.pro", "-msh", str(Path(msh).resolve()), "-solve", "ES", "-pos", "Post", "-v", "2"]
    for k, v in numbers.items():
        args += ["-setnumber", k, repr(float(v))]
    r = subprocess.run(args, cwd=work, capture_output=True, text=True)
    (work / "getdp.log").write_text(r.stdout + r.stderr)
    if r.returncode or not (work / "v_elem.txt").exists():
        raise RuntimeError(f"getdp failed in {work}:\n{(r.stdout + r.stderr)[-2000:]}")
    last = lambda f: float((work / f).read_text().split()[-1])        # noqa: E731
    return {k: last(f"{k}.txt") for k in ("Q_anode", "V_anode", "Q_iron", "V_iron", "Q_chamber", "energy")}


def bases(work, msh):
    """Both per-volt bases; returns (capacitance matrix [[C_AA, C_AI], [C_IA, C_II]] in F, raw results)."""
    a = solve(work / "basis_anode", msh, V_anode=1, Iron_floating=0, V_iron=0)
    i = solve(work / "basis_iron", msh, V_anode=0, Iron_floating=0, V_iron=1)
    c = np.array([[a["Q_anode"], i["Q_anode"]], [a["Q_iron"], i["Q_iron"]]])
    return c, {"anode": a, "iron": i}


def load_nodal(work):
    """GetDP's element-nodal potential (SimpleTable rows: x0 y0 z0 x1 y1 z1 x2 y2 z2 v0 v1 v2, metres) ->
    unique nodes (n, 2), triangles (m, 3), nodal values (n,)."""
    a = pd.read_csv(work / "v_elem.txt", sep=r"\s+", header=None).to_numpy(float)
    xy = a[:, [0, 1, 3, 4, 6, 7]].reshape(-1, 2)
    v = a[:, 9:12].reshape(-1)
    _, idx, inv = np.unique(np.round(xy * 1e9).astype(np.int64), axis=0, return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    spread = np.max(np.abs(v - v[idx][inv]))           # one value per node (the potential is continuous)
    assert spread < 1e-9 * max(1.0, np.abs(v).max()), spread
    return xy[idx], inv.reshape(-1, 3), v[idx]


# ---------------------------------------------------------------------------------------------- grid export
def region_map(R, Z):
    """Region code of every grid point (mm); conductors drawn over BN over vacuum."""
    reg = np.full(R.shape, REGION["vacuum"], np.int8)
    pts = np.column_stack([R.ravel(), Z.ravel()])
    for kind in ("bn", "anode", "iron", "coil"):
        code = REGION["iron" if kind == "coil" else kind]
        for k2, poly in g.parts().values():
            if k2 == kind:
                inside = MPath(poly).contains_points(pts, radius=1e-9).reshape(R.shape)
                reg[inside] = code
    return reg


def sample(tri, vals, R_m, Z_m):
    """Potential and E = -grad(phi) of the first-order solution at grid points (metres); masked outside the mesh."""
    f = mtri.LinearTriInterpolator(tri, vals)
    phi = f(R_m, Z_m)
    gx, gy = f.gradient(R_m, Z_m)
    return phi, -gx, -gy


def export_grid(tri, basis_vals, reg, R, Z):
    """Per-volt basis fields on the grid, conductor interiors filled with their potential and E = 0."""
    out = {}
    fill = {"A": {REGION["anode"]: 1.0, REGION["iron"]: 0.0}, "I": {REGION["anode"]: 0.0, REGION["iron"]: 1.0}}
    for b, vals in basis_vals.items():
        phi, er, ez = sample(tri, vals, R * 1e-3, Z * 1e-3)
        phi, er, ez = (np.ma.filled(x, np.nan) for x in (phi, er, ez))
        for code, v in fill[b].items():
            m = reg == code
            phi[m], er[m], ez[m] = v, 0.0, 0.0
        holes = np.isnan(phi)                                  # points on a conductor boundary outside every triangle
        if holes.any():
            from scipy.spatial import cKDTree
            nodes = tri.x * 1e3, tri.y * 1e3
            _, j = cKDTree(np.column_stack(nodes)).query(np.column_stack([R[holes], Z[holes]]))
            phi[holes] = vals[j]
            er[holes] = ez[holes] = 0.0
        out[f"phi_{b}"], out[f"Er_{b}"], out[f"Ez_{b}"] = (x.astype(np.float32) for x in (phi, er, ez))
        out[f"holes_{b}"] = int(holes.sum())
    return out


# ---------------------------------------------------------------------------------------------- figures
LABELS = {"bn_inner": "BN inner wall", "bn_outer": "BN body", "anode": "anode", "front_plate": "front plate",
          "back_plate": "back plate", "pole_cap": "pole cap", "inner_core": "inner core", "inner_coil": "inner coil",
          "outer_posts": "outer posts (ring)", "outer_coils_in": "outer coils", "outer_coils_out": "outer coils"}


def fig_geometry(path, tri):
    """Left: the parts as modelled (z horizontal). Right: the mesh at the narrowest gap."""
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axs[0]
    for name, (kind, poly) in g.parts().items():
        p = np.array(poly)[:, ::-1]
        ax.add_patch(Polygon(p, closed=True, fc=COLORS[kind], ec="k", lw=0.5))
        c = p.mean(axis=0)
        ax.text(c[0], c[1], LABELS[name], fontsize=7, ha="center", va="center",
                rotation=90 if name in ("front_plate", "back_plate", "anode", "bn_inner") else 0)
    ax.set_xlim(-65, 15)
    ax.set_ylim(0, 90)
    ax.set_aspect("equal")
    ax.set_xlabel("z [mm] (exit plane at 0, plume to the right)")
    ax.set_ylabel("r [mm]")
    ax.set_title("geometry (iron grey, coils yellow, BN blue, anode orange)", fontsize=10)
    ax = axs[1]
    ax.triplot(tri.y * 1e3, tri.x * 1e3, lw=0.15, color="k")
    for kind, poly in g.parts().values():
        ax.add_patch(Polygon(np.array(poly)[:, ::-1], closed=True, fc=COLORS[kind], ec="k", lw=0.5, alpha=0.6))
    ax.set_xlim(-12, 4)
    ax.set_ylim(14, 26)
    ax.set_aspect("equal")
    ax.set_xlabel("z [mm]")
    ax.set_title("mesh near the pole cap edge face (0.12 mm gap to the BN bore)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _map(ax, R, Z, field, reg, title, levels=None, cmap="viridis", norm=None, label=""):
    f = np.where(reg >= REGION["anode"], np.nan, field)
    im = ax.pcolormesh(Z, R, f, shading="auto", cmap=cmap, norm=norm)
    if levels is not None:
        cs = ax.contour(Z, R, f, levels=levels, colors="w", linewidths=0.6)
        ax.clabel(cs, fontsize=6, fmt="%g")
    for kind, poly in g.parts().values():
        ax.add_patch(Polygon(np.array(poly)[:, ::-1], closed=True, fc=COLORS[kind] if kind != "bn" else "none",
                             ec="k", lw=0.6, zorder=3))
    ax.set_aspect("equal")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8, label=label)


def figures(R, Z, f, reg, k, outdir):
    phi_g = V_D * f["phi_A"]
    phi_f = V_D * (f["phi_A"] + k * f["phi_I"])
    emag_g = V_D * np.hypot(f["Er_A"], f["Ez_A"])
    emag_f = V_D * np.hypot(f["Er_A"] + k * f["Er_I"], f["Ez_A"] + k * f["Ez_I"])
    lv = [1, 2, 5, 10, 20, 50, 100, 150, 200, 250, 290]
    for name, (lo_z, hi_z, hi_r) in {"potential": (-55, 60, 100), "potential_channel": (-50, 15, 50)}.items():
        fig, axs = plt.subplots(1, 2, figsize=(15, 6.2 if hi_r > 60 else 5.0))
        _map(axs[0], R, Z, phi_g, reg, f"potential, iron grounded (anode {V_D:g} V, chamber wall 0 V)", lv,
             label="V")
        _map(axs[1], R, Z, phi_f, reg, f"potential, iron floating: iron at {V_D * k:.1f} V", lv, label="V")
        for ax in axs:
            ax.set_xlim(lo_z, hi_z)
            ax.set_ylim(0, hi_r)
        fig.tight_layout()
        fig.savefig(outdir / f"{name}.png", dpi=150)
        plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.2))
    norm = LogNorm(1e2, 2e5)
    _map(axs[0], R, Z, np.maximum(emag_g, 1e2), reg, "|E|, iron grounded", cmap="inferno", norm=norm, label="V/m")
    _map(axs[1], R, Z, np.maximum(emag_f, 1e2), reg, "|E|, iron floating", cmap="inferno", norm=norm, label="V/m")
    for ax in axs:
        ax.set_xlim(-55, 60)
        ax.set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(outdir / "efield.png", dpi=150)
    plt.close(fig)

    # lines along the channel mid-radius and across the exit plane
    rm = g.channel_mid_radius()
    zc = Z[0]
    j = np.argmin(np.abs(R[:, 0] - rm))
    keep = zc >= g.ANODE_Z1
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.5))
    for lab, phi in (("grounded", phi_g), ("floating", phi_f)):
        axs[0, 0].plot(zc[keep], phi[j, keep], label=f"iron {lab}")
    axs[0, 0].set_title(f"potential along the channel mid-radius r = {R[j, 0]:.2f} mm", fontsize=10)
    axs[0, 0].set_ylabel("V")
    ez_g = V_D * f["Ez_A"][j]
    ez_f = V_D * (f["Ez_A"][j] + k * f["Ez_I"][j])
    axs[1, 0].plot(zc[keep], ez_g[keep] / 1e3, label="iron grounded")
    axs[1, 0].plot(zc[keep], ez_f[keep] / 1e3, label="iron floating")
    axs[1, 0].set_title("axial field E_z along the mid-radius", fontsize=10)
    axs[1, 0].set_ylabel("kV/m")
    for ax in axs[:, 0]:
        ax.axvline(0, color="k", lw=0.6, ls="--")
        ax.set_xlabel("z [mm] (anode face at -41.43, exit plane at 0)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    for zcut, ls in ((-20.0, "-"), (0.0, "--"), (10.0, ":")):
        i = np.argmin(np.abs(zc - zcut))
        rr = R[:, i]
        for lab, phi in (("grounded", phi_g), ("floating", phi_f)):
            line = np.where(reg[:, i] <= REGION["bn"], phi[:, i], np.nan)      # gaps where the line crosses a conductor
            axs[0, 1].plot(rr, line, ls, label=f"z = {zc[i]:+.1f} mm, iron {lab}",
                           color="C0" if lab == "grounded" else "C1")
    axs[0, 1].set_title("potential across the channel / plume", fontsize=10)
    axs[0, 1].set_xlabel("r [mm]")
    axs[0, 1].set_ylabel("V")
    axs[0, 1].set_xlim(0, 100)
    axs[0, 1].legend(fontsize=7)
    axs[0, 1].grid(alpha=0.3)
    dphi = phi_f - phi_g
    _map(axs[1, 1], R, Z, dphi, reg, "floating minus grounded potential [V]",
         levels=[5, 10, 20, 50, 100], cmap="magma", label="V")
    axs[1, 1].set_xlim(-55, 60)
    axs[1, 1].set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(outdir / "lines.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------- main
def probes(tri, vals):
    """A few fixed points (mm) for convergence checks: potential [V/V] and E_z [V/m per V]."""
    pts = {"anode_face+1": (g.channel_mid_radius(), g.ANODE_Z1 + 1.0), "mid_channel": (g.channel_mid_radius(), -20.0),
           "exit_mid": (g.channel_mid_radius(), 0.0), "plume_10": (g.channel_mid_radius(), 10.0)}
    r = np.array([p[0] for p in pts.values()]) * 1e-3
    z = np.array([p[1] for p in pts.values()]) * 1e-3
    phi, er, ez = sample(tri, vals, r, z)
    return {k: {"phi": float(phi[i]), "Ez": float(ez[i])} for i, k in enumerate(pts)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="skip the chamber-size and mesh-convergence studies")
    ap.add_argument("--chamber", type=float, default=1000.0, help="chamber wall radius [mm] (default 1000)")
    args = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    work = OUT / "work"
    msh = work / "design_b_rz.msh"
    info = mesh_rz.build(msh, chamber_r=args.chamber)
    print(f"mesh: {info['triangles']} triangles, {info['nodes']} nodes ({time.time() - t0:.1f} s)")

    C, raw = bases(work, msh)
    k = -C[1, 0] / C[1, 1]
    recip = abs(C[0, 1] - C[1, 0]) / abs(C[1, 0])
    direct = solve(work / "floating_direct", msh, V_anode=V_D, Iron_floating=1, Q_iron=0)
    print(f"capacitance matrix [pF]: C_AA {C[0, 0] * 1e12:.4f}  C_AI {C[0, 1] * 1e12:.4f}  C_IA {C[1, 0] * 1e12:.4f}  "
          f"C_II {C[1, 1] * 1e12:.4f}  (reciprocity {recip:.1e})")
    print(f"floating iron: superposition {V_D * k:.4f} V, direct solve {direct['V_iron']:.4f} V")

    pts, tris, vA = load_nodal(work / "basis_anode")
    pts2, tris2, vI = load_nodal(work / "basis_iron")
    pts3, tris3, vF = load_nodal(work / "floating_direct")
    assert np.array_equal(pts, pts2) and np.array_equal(pts, pts3)
    sup_err = np.max(np.abs(V_D * (vA + k * vI) - vF))
    print(f"superposition vs direct floating solve: max nodal difference {sup_err:.2e} V")
    tri = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)

    r = np.arange(GRID["r"][0], GRID["r"][1] + 1e-9, GRID["h"])
    z = np.arange(GRID["z"][0], GRID["z"][1] + 1e-9, GRID["h"])
    Z, R = np.meshgrid(z, r)                                 # arrays indexed [r, z]
    reg = region_map(R, Z)
    f = export_grid(tri, {"A": vA, "I": vI}, reg, R, Z)
    np.savez_compressed(OUT / "fields_rz.npz", r=(r * 1e-3).astype(np.float32), z=(z * 1e-3).astype(np.float32),
                        region=reg, **{k2: v for k2, v in f.items() if not k2.startswith("holes")},
                        C=C, k_float=k, eps_r_bn=EPS_R_BN, chamber_r=args.chamber * 1e-3)
    figures(R, Z, f, reg, k, OUT)
    fig_geometry(OUT / "geometry.png", tri)

    summary = {"V_d": V_D, "eps_r_bn": EPS_R_BN, "mesh": info, "capacitance_pF": (C * 1e12).round(5).tolist(),
               "reciprocity_rel": recip, "k_float": k, "V_iron_floating": V_D * k,
               "V_iron_floating_direct": direct["V_iron"], "superposition_max_err_V": sup_err,
               "grid_holes_filled": {b: f[f"holes_{b}"] for b in "AI"},
               "probes_grounded": {p: {"phi_V": V_D * d["phi"], "Ez_V_per_m": V_D * d["Ez"]}
                                   for p, d in probes(tri, vA).items()},
               "probes_floating": {p: {"phi_V": V_D * d["phi"], "Ez_V_per_m": V_D * d["Ez"]}
                                   for p, d in probes(tri, vA + k * vI).items()}}

    if not args.quick:
        cham = {}
        for rc in (500.0, 2000.0, 4000.0):
            m2 = work / f"chamber_{rc:g}" / "m.msh"
            mesh_rz.build(m2, chamber_r=rc)
            c2, _ = bases(work / f"chamber_{rc:g}", m2)
            cham[f"{rc:g}"] = {"V_iron_floating": V_D * -c2[1, 0] / c2[1, 1], "C_II_pF": c2[1, 1] * 1e12,
                               "C_IA_pF": c2[1, 0] * 1e12, "C_AA_pF": c2[0, 0] * 1e12}
        cham[f"{args.chamber:g}"] = {"V_iron_floating": V_D * k, "C_II_pF": C[1, 1] * 1e12, "C_IA_pF": C[1, 0] * 1e12,
                                     "C_AA_pF": C[0, 0] * 1e12}
        summary["chamber_size"] = dict(sorted(cham.items(), key=lambda kv: float(kv[0])))
        m3 = work / "fine" / "m.msh"
        info3 = mesh_rz.build(m3, chamber_r=args.chamber, h_channel=0.2, h_thruster=0.75, h_gap=0.015)
        c3, _ = bases(work / "fine", m3)
        p3, t3, vA3 = load_nodal(work / "fine" / "basis_anode")
        p4, t4, vI3 = load_nodal(work / "fine" / "basis_iron")
        tri3 = mtri.Triangulation(p3[:, 0], p3[:, 1], t3)
        k3 = -c3[1, 0] / c3[1, 1]
        summary["mesh_convergence"] = {
            "fine_mesh": info3, "V_iron_floating": V_D * k3, "capacitance_pF": (c3 * 1e12).round(5).tolist(),
            "probes_floating": {p: {"phi_V": V_D * d["phi"], "Ez_V_per_m": V_D * d["Ez"]}
                                for p, d in probes(tri3, vA3 + k3 * vI3).items()}}
        print(f"chamber size: " + ", ".join(f"R {kk} mm -> {vv['V_iron_floating']:.2f} V"
                                            for kk, vv in summary["chamber_size"].items()))
        print(f"fine mesh ({info3['triangles']} triangles): floating iron {V_D * k3:.3f} V")
    summary["runtime_s"] = round(time.time() - t0, 1)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({kk: summary[kk] for kk in ("probes_grounded", "probes_floating")}, indent=1))
    print(f"done in {time.time() - t0:.1f} s; figures and fields in {OUT}")


if __name__ == "__main__":
    main()
