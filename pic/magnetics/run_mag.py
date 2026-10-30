"""Magnetic field of Design B at the 165 G design point, for the PIC model (2-D axisymmetric, gmsh + GetDP; no bplsim).

The same model as the het library's 2-D magnetostatics (sim/library/pipelines/magnetostatics_2d.yaml):
- linear iron, mu_r 1500;
- the outer posts smeared into a ring whose mu_r is scaled by their area fraction (226);
- the outer coils as two opposite current sheets;
- outer current = inner current x q x 230/184, with q = 0.635.
The currents are scaled so the channel profile (B_r averaged across the channel, r-weighted) peaks at 165 G.
Check: for Design B, the het parity case (sim/PARITY.md) needs 0.927 A inner and 0.736 A outer.

Writes, in out/:
- bfield_rz.npz: B_r, B_z on a grid
- bfield_rz.h5: openPMD thetaMode map read by WarpX RZ, and bfield_xz.h5, the same map as a Cartesian (x = r, z)
  map for the planar x-z PIC model
- bfield.png
- summary.json, which includes the |B|(z) fit used for the PIC's anomalous collision frequency

  python run_mag.py           (about 15 s)
"""
import json
import shutil
import subprocess
import sys
from math import pi
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from geometry import design_b as g  # noqa: E402

GETDP = shutil.which("getdp") or str(Path.home() / "tools/getdp/bin/getdp")
OUT = HERE / "out"
B_PEAK = 0.0165              # [T] design point: channel-profile B_r peak
Q = 0.635                    # outer A-turns per coil / inner A-turns (het assembly legacy_bn)
MU_IRON = 1500.0             # het:materials/iron_legacy
FILL = g.POST_COUNT * pi * g.POST_R ** 2 / (2 * pi * g.POST_RING_R * 2 * g.POST_R)
MU_POSTS = 1 + FILL * (MU_IRON - 1)
TAG = {"air": 1, "iron": 2, "posts": 3, "coil_in": 4, "sheet_in": 5, "sheet_out": 6, "bnd": 10}
GRID = dict(r=(0.0, 120.0), z=(-60.0, 150.0), h=0.25)      # [mm] full map
PIC_GRID = dict(r=(15.0, 40.0), z=(-50.0, 60.0), h=0.25)   # [mm] map handed to WarpX (covers the PIC domain)


def mesh(path):
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        occ = gmsh.model.occ
        R, zc = 400.0, -27.846                               # air region of the het pipeline
        dom, _ = occ.intersect([(2, occ.addDisk(0, zc, 0, R, R))],
                               [(2, occ.addRectangle(0, zc - R - 1e-3, 0, R + 1e-3, 2 * R + 2e-3))])
        kinds, tools = [], []
        for name, (kind, poly) in g.parts().items():
            if kind in ("bn", "anode"):
                continue                                     # mu_r 1: part of the air region
            k = {"iron": "posts" if name == "outer_posts" else "iron",
                 "coil": {"inner_coil": "coil_in", "outer_coils_in": "sheet_in", "outer_coils_out": "sheet_out"}[name]
                 if kind == "coil" else None}[kind]
            pts = [occ.addPoint(r, z, 0) for r, z in poly]
            lines = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
            tools.append((2, occ.addPlaneSurface([occ.addCurveLoop(lines)])))
            kinds.append(k)
        _, omap = occ.fragment(dom, tools)
        occ.synchronize()
        groups = {k: [] for k in TAG}
        used = set()
        for k, ents in zip(kinds, omap[1:]):
            groups[k] += [t for _, t in ents]
            used |= {t for _, t in ents}
        groups["air"] = sorted({t for _, t in omap[0]} - used)
        areas = {k: sum(occ.getMass(2, t) for t in groups[k]) * 1e-6 for k in ("coil_in", "sheet_in")}
        m = gmsh.model
        for k in ("air", "iron", "posts", "coil_in", "sheet_in", "sheet_out"):
            m.addPhysicalGroup(2, groups[k], TAG[k], k)
        everything = [(2, t) for k in ("air", "iron", "posts", "coil_in", "sheet_in", "sheet_out") for t in groups[k]]
        bnd = [abs(t) for _, t in m.getBoundary(everything, combined=True, oriented=False)
               if m.getBoundingBox(1, abs(t))[3] > 1e-9]
        m.addPhysicalGroup(1, bnd, TAG["bnd"], "bnd")
        f = m.mesh.field
        fine = f.add("Box")
        for k, v in dict(VIn=0.35, VOut=25.0, XMin=0, XMax=45, YMin=-47.4274, YMax=30, Thickness=5).items():
            f.setNumber(fine, k, v)
        iron = f.add("Box")
        for k, v in dict(VIn=1.2, VOut=25.0, XMin=0, XMax=86, YMin=-65, YMax=40, Thickness=20).items():
            f.setNumber(iron, k, v)
        mn = f.add("Min")
        f.setNumbers(mn, "FieldsList", [fine, iron])
        f.setAsBackgroundMesh(mn)
        for k, v in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
                     "Mesh.MeshSizeFromCurvature": 0, "Mesh.Algorithm": 6, "Mesh.ScalingFactor": 1e-3}.items():
            gmsh.option.setNumber(k, v)
        m.mesh.generate(2)
        gmsh.write(str(path))
        return areas
    finally:
        gmsh.finalize()


def solve(work, msh, areas, i_inner, i_outer):
    shutil.copy(HERE / "mag_rz.pro", work / "mag_rz.pro")    # GetDP writes its outputs next to the .pro
    nums = dict(I_inner=i_inner, I_outer=i_outer, A_inner=areas["coil_in"], A_sheet=areas["sheet_in"],
                mu_iron=MU_IRON, mu_posts=MU_POSTS, N_inner=g.INNER_TURNS, N_outer=g.OUTER_TURNS)
    args = [GETDP, "mag_rz.pro", "-msh", str(msh), "-solve", "MS", "-pos", "Post", "-v", "2"]
    for k, v in nums.items():
        args += ["-setnumber", k, repr(float(v))]
    r = subprocess.run(args, cwd=work, capture_output=True, text=True)
    (work / "getdp.log").write_text(r.stdout + r.stderr)
    if r.returncode or not (work / "b_elem.txt").exists():
        sys.exit(f"getdp failed:\n{(r.stdout + r.stderr)[-2000:]}")
    a = pd.read_csv(work / "b_elem.txt", sep=r"\s+", header=None).to_numpy(float)
    return a[:, [0, 1, 3, 4, 6, 7]].reshape(-1, 3, 2), a[:, [9, 10, 12, 13, 15, 16]].reshape(-1, 3, 2)


class Sampler:
    """B at arbitrary points: locate the element, interpolate its nodal values (GetDP element-nodal export)."""

    def __init__(self, P, B):
        xy = P.reshape(-1, 2)
        _, idx, inv = np.unique(np.round(xy * 1e9).astype(np.int64), axis=0, return_index=True, return_inverse=True)
        self.tri = mtri.Triangulation(xy[idx, 0], xy[idx, 1], inv.reshape(-1, 3))
        self.finder = self.tri.get_trifinder()
        self.P, self.B = P, B

    def __call__(self, r, z):
        r, z = np.broadcast_arrays(np.asarray(r, float), np.asarray(z, float))
        e = self.finder(r.ravel(), z.ravel())
        if np.any(e < 0):
            raise ValueError(f"{np.sum(e < 0)} points outside the mesh")
        p = self.P[e]
        v0, v1 = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
        v2 = np.column_stack([r.ravel(), z.ravel()]) - p[:, 0]
        d00, d01, d11 = (v0 * v0).sum(1), (v0 * v1).sum(1), (v1 * v1).sum(1)
        d20, d21 = (v2 * v0).sum(1), (v2 * v1).sum(1)
        den = d00 * d11 - d01 * d01
        l1 = (d11 * d20 - d01 * d21) / den
        l2 = (d00 * d21 - d01 * d20) / den
        lam = np.column_stack([1 - l1 - l2, l1, l2])
        b = (lam[:, :, None] * self.B[e]).sum(1)
        return b[:, 0].reshape(r.shape), b[:, 1].reshape(r.shape)


def channel_profile(s):
    """The het channel profile: B_r averaged (r-weighted) over 11 radii across the channel, every 0.25 mm."""
    rr = np.linspace(g.R_IN, g.R_OUT, 11) * 1e-3
    n_z = int(np.floor((30 + g.DEPTH + 3) / 0.25 + 1e-9))
    z = (30 - 0.25 * np.arange(n_z, -1, -1)) * 1e-3
    Z, Rr = np.meshgrid(z, rr)
    br, _ = s(Rr, Z)
    return z, (br * rr[:, None]).sum(0) / rr.sum()


def bmag_model(z, b0, zp, sl, sr, fl, fr):
    """Asymmetric Gaussian with its own floor on each side (continuous at the peak). The PIC script writes the same
    expression into WarpX's parser."""
    up = z < zp
    s, fl_ = np.where(up, sl, sr), np.where(up, fl, fr)
    return fl_ + (b0 - fl_) * np.exp(-0.5 * ((z - zp) / s) ** 2)


def fit_bmag(z, bmag):
    """|B|(z) at the channel mid-radius (for the anomalous transport, which WarpX can only take as a formula)."""
    from scipy.optimize import curve_fit
    p0 = [bmag.max(), z[np.argmax(bmag)], 0.012, 0.012, 0.05 * bmag.max(), 0.1 * bmag.max()]
    p, _ = curve_fit(bmag_model, z, bmag, p0=p0, maxfev=20000)
    rms = float(np.sqrt(np.mean((bmag_model(z, *p) - bmag) ** 2)) / bmag.max())
    keys = ["B_peak_T", "z_peak_m", "sigma_up_m", "sigma_down_m", "B_floor_up_T", "B_floor_down_T"]
    return dict(zip(keys, map(float, p))), rms, bmag_model


def write_openpmd(path, r, z, br, bz, planar=False):
    """thetaMode (r, z) map for WarpX RZ, or planar=True: the same values as a Cartesian (x, z) map with x = r, for
    the planar x-z PIC model."""
    import openpmd_api as io
    if path.exists():
        path.unlink()
    s = io.Series(str(path), io.Access.create)
    it = s.iterations[0]
    B = it.meshes["B"]
    if planar:
        B.geometry = io.Geometry.cartesian
        B.axis_labels = ["x", "z"]
    else:
        B.geometry = io.Geometry.thetaMode
        B.geometry_parameters = "m=0"
        B.axis_labels = ["r", "z"]
    B.grid_spacing = [float(r[1] - r[0]), float(z[1] - z[0])]
    B.grid_global_offset = [float(r[0]), float(z[0])]
    B.grid_unit_SI = 1.0
    B.data_order = "C"
    B.unit_dimension = {io.Unit_Dimension.M: 1, io.Unit_Dimension.I: -1, io.Unit_Dimension.T: -2}
    names = ("x", "y", "z") if planar else ("r", "t", "z")
    for comp, arr in zip(names, (br, np.zeros_like(br), bz)):
        rc = B[comp]
        data = np.ascontiguousarray(arr if planar else arr[None, :, :], dtype=np.float64)
        rc.reset_dataset(io.Dataset(data.dtype, data.shape))
        rc.store_chunk(data)
        rc.position = [0.0, 0.0]
        rc.unit_SI = 1.0
    s.flush()
    s.close()


def main():
    OUT.mkdir(exist_ok=True)
    work = OUT / "work"
    work.mkdir(exist_ok=True)
    msh = work / "mag.msh"
    areas = mesh(msh)
    i_out_per_in = Q * g.INNER_TURNS / g.OUTER_TURNS
    P, B = solve(work, msh, areas, 1.0, i_out_per_in)        # linear: solve at 1 A, then scale
    s = Sampler(P, B)
    zp, prof = channel_profile(s)
    k = B_PEAK / np.max(np.abs(prof))
    i_inner, i_outer = k, k * i_out_per_in
    B *= k
    sign = np.sign(prof[np.argmax(np.abs(prof))])
    print(f"currents at {B_PEAK * 1e4:.0f} G: inner {i_inner:.5f} A, outer {i_outer:.5f} A "
          f"(het parity: 0.92733 / 0.73607 A); posts mu_r {MU_POSTS:.1f}")

    def grid(spec):
        r = np.arange(spec["r"][0], spec["r"][1] + 1e-9, spec["h"]) * 1e-3
        z = np.arange(spec["z"][0], spec["z"][1] + 1e-9, spec["h"]) * 1e-3
        Z, R = np.meshgrid(z, r)
        br, bz = s(R, Z)
        return r, z, br, bz

    r, z, br, bz = grid(GRID)
    rp, zp2, brp, bzp = grid(PIC_GRID)
    write_openpmd(OUT / "bfield_rz.h5", rp, zp2, brp, bzp)
    write_openpmd(OUT / "bfield_xz.h5", rp, zp2, brp, bzp, planar=True)
    rm = g.channel_mid_radius() * 1e-3
    zl = np.arange(-g.DEPTH, 30.0 + 1e-9, 0.25) * 1e-3
    brl, bzl = s(np.full_like(zl, rm), zl)
    bmag = np.hypot(brl, bzl)
    fit, fit_rms, ffun = fit_bmag(zl, bmag)
    _, prof = channel_profile(s)
    i_pk = np.argmax(np.abs(prof))
    np.savez_compressed(OUT / "bfield_rz.npz", r=r, z=z, Br=br.astype(np.float32), Bz=bz.astype(np.float32),
                        I_inner=i_inner, I_outer=i_outer)
    # 1-D field profile for HallThruster.jl: channel-averaged B_r (r-weighted, 11 radii) and |B| at mid-radius,
    # from the anode face to 100 mm downstream of the exit; z measured from the anode face, as HallThruster.jl does
    z1 = np.arange(-g.DEPTH, 100.0 + 1e-9, 0.5) * 1e-3
    rr = np.linspace(g.R_IN, g.R_OUT, 11) * 1e-3
    Z1, R1 = np.meshgrid(z1, rr)
    br1, _ = s(R1, Z1)
    prof1 = (br1 * rr[:, None]).sum(0) / rr.sum()
    bm1 = np.hypot(*s(np.full_like(z1, rm), z1))
    np.savetxt(OUT / "bfield_profile_1d.csv",
               np.column_stack([z1 + g.DEPTH * 1e-3, z1, np.abs(prof1), bm1]), delimiter=",", fmt="%.6e",
               header="z_from_anode_m,z_m,Br_channel_avg_T,Bmag_mid_radius_T", comments="")
    summary = {"B_peak_T": B_PEAK, "q": Q, "I_inner_A": i_inner, "I_outer_A": i_outer, "mu_r_iron": MU_IRON,
               "mu_r_posts": MU_POSTS, "profile_peak_z_mm": float(zp[i_pk] * 1e3), "profile_sign": float(sign),
               "B_r_at_exit_pct": float(100 * np.interp(0.0, zp, prof) / prof[i_pk]),
               "B_r_at_anode_pct": float(100 * np.interp(-g.DEPTH * 1e-3, zp, prof) / prof[i_pk]),
               "bmag_fit_mid_radius": fit, "bmag_fit_rms_rel": fit_rms,
               "openpmd_map": {"file": "bfield_rz.h5", "r_mm": PIC_GRID["r"], "z_mm": PIC_GRID["z"],
                               "h_mm": PIC_GRID["h"]}}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=1))

    fig, axs = plt.subplots(1, 2, figsize=(15, 5.8))
    ax = axs[0]
    Z, R = np.meshgrid(z * 1e3, r * 1e3)
    bm = np.hypot(br, bz) * 1e4
    im = ax.pcolormesh(Z, R, bm, shading="auto", cmap="viridis", vmax=400)
    ax.streamplot(z * 1e3, r * 1e3, bz, br, color="w", density=1.6, linewidth=0.5, arrowsize=0.6)
    for kind, poly in g.parts().values():
        ax.add_patch(Polygon(np.array(poly)[:, ::-1], closed=True, fc="none", ec="k", lw=0.6))
    ax.set_xlim(-65, 40)
    ax.set_ylim(0, 95)
    ax.set_aspect("equal")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title(f"|B| [G] and field lines; inner {i_inner:.3f} A, outer {i_outer:.3f} A", fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8, label="G")
    ax = axs[1]
    ax.plot(zp * 1e3, prof * 1e4, label="channel profile: B_r averaged across the channel")
    ax.plot(zl * 1e3, bmag * 1e4, label=f"|B| at mid-radius r = {rm * 1e3:.2f} mm")
    ax.plot(zl * 1e3, ffun(zl, *fit.values()) * 1e4, "--", label=f"fit used for anomalous transport (rms {fit_rms:.1%})")
    ax.axvline(0, color="k", lw=0.6, ls="--")
    ax.axvline(-g.DEPTH, color="C3", lw=0.6, ls=":", label="anode face")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("G")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("Design B at the 165 G design point", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "bfield.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
