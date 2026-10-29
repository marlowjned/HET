"""Vacuum electrostatics of Design B in 3-D (gmsh + GetDP; no bplsim).

Same formulation as the 2-D model (estatic.pro, Flag_Axi = 0), with the 3-D geometry of mesh_3d.py:
- square plates, four posts with their coils;
- a placeholder cathode;
- a spherical chamber wall, 0.5 m.

Three per-volt bases give the 3 x 3 capacitance matrix (anode, iron, cathode; chamber wall the reference), and
every case by superposition. A direct solve checks the floating case (anode 300 V, iron floating, cathode at the
chamber wall's 0 V) and supplies the potential on two planes:
- y = 0, through a post and the cathode;
- 45 degrees, between two posts.

    python run_3d.py [--threads 4] [--h-channel 1.5]      (about 10 min; out/3d/)
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from math import cos, pi, sin
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from geometry import design_b as g  # noqa: E402
import mesh_3d  # noqa: E402

GETDP = shutil.which("getdp") or str(Path.home() / "tools/getdp/bin/getdp")
OUT = HERE / "out" / "3d"
V_D = 300.0


def getdp(work, msh, post, threads, res=None, **numbers):
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(HERE / "estatic.pro", work / "estatic.pro")
    args = [GETDP, "estatic.pro", "-msh", str(msh)]
    args += ["-res", res] if res else ["-solve", "ES"]
    args += ["-pos", post, "-v", "2"]
    numbers = {"Flag_Axi": 0, "Has_cathode": 1, "eps_r_bn": 4.0, **numbers}
    for k, v in numbers.items():
        args += ["-setnumber", k, repr(float(v))]
    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OPENBLAS_NUM_THREADS": str(threads)}
    t = time.time()
    r = subprocess.run(args, cwd=work, capture_output=True, text=True, env=env)
    (work / f"getdp_{post}.log").write_text(r.stdout + r.stderr)
    if r.returncode:
        sys.exit(f"getdp failed in {work}:\n{(r.stdout + r.stderr)[-2000:]}")
    return time.time() - t


def globals_(work):
    last = lambda f: float((work / f).read_text().split()[-1])        # noqa: E731
    return {k: last(f"{k}.txt") for k in ("Q_anode", "V_anode", "Q_iron", "V_iron", "Q_cathode", "V_cathode",
                                          "Q_chamber")}


def plane(work, msh, threads, angle_deg, rmax=150.0, z0=-80.0, z1=120.0, h=1.0):
    """Potential on the half plane at azimuth angle_deg (r from -rmax to rmax through the axis), z0..z1 [mm]."""
    a = np.radians(angle_deg)
    ux, uy = cos(a), sin(a)
    n1, n2 = int(round(2 * rmax / h)), int(round((z1 - z0) / h))
    mm = 1e-3
    p = dict(P0x=-rmax * ux * mm, P0y=-rmax * uy * mm, P0z=z0 * mm, P1x=rmax * ux * mm, P1y=rmax * uy * mm,
             P1z=z0 * mm, P2x=-rmax * ux * mm, P2y=-rmax * uy * mm, P2z=z1 * mm, N1=n1, N2=n2)
    dt = getdp(work, msh, "Plane", threads, res="estatic.res", V_anode=V_D, Iron_floating=1, Q_iron=0,
               V_cathode=0, **p)
    d = np.loadtxt(work / "plane.txt")
    s = (d[:, 0] * ux + d[:, 1] * uy) * 1e3                                  # signed radius along the plane [mm]
    return s, d[:, 2] * 1e3, d[:, -1], dt


def planes_figure(out, msh, threads, v_iron):
    """Sample the floating case on the 0 and 45 degree planes, save them (planes.npz) and draw potential_3d.png.
    GetDP returns exactly 0 at points outside the dielectric mesh (conductor interiors): those are masked."""
    times, data = {}, {}
    for ang in (0.0, 45.0):
        s, z, v, dt = plane(out / "floating", msh, threads, ang)
        data[f"s_{ang:g}"], data[f"z_{ang:g}"], data[f"v_{ang:g}"] = s, z, v
        times[f"plane_{ang:g}"] = dt
    np.savez_compressed(out / "planes.npz", **data)
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    for ax, ang in zip(axs, (0.0, 45.0)):
        s, z, v = data[f"s_{ang:g}"], data[f"z_{ang:g}"], data[f"v_{ang:g}"]
        su, zu = np.unique(np.round(s, 6)), np.unique(np.round(z, 6))
        V = np.full((zu.size, su.size), np.nan)
        V[np.searchsorted(zu, np.round(z, 6)), np.searchsorted(su, np.round(s, 6))] = v
        V[V == 0.0] = np.nan                               # conductor interiors (outside the mesh)
        ax.set_facecolor("0.6")
        im = ax.pcolormesh(zu, su, V.T, shading="auto", cmap="viridis", vmin=0, vmax=V_D)
        cs = ax.contour(zu, su, V.T, levels=[20, 50, 100, 150, 180, 200, 250, 290], colors="w", linewidths=0.5)
        ax.clabel(cs, fontsize=6, fmt="%g")
        ax.set_aspect("equal")
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("signed radius along the plane [mm]")
        ax.set_title(f"iron floating ({v_iron:.1f} V), cathode 0 V; plane at {ang:g} deg"
                     + (" (through a post and the cathode)" if ang == 0 else " (between posts)")
                     + "; conductors grey", fontsize=9)
        plt.colorbar(im, ax=ax, shrink=0.8, label="V")
    fig.tight_layout()
    fig.savefig(out / "potential_3d.png", dpi=140)
    plt.close(fig)
    return times


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--h-channel", type=float, default=1.5)
    ap.add_argument("--remesh", action="store_true")
    ap.add_argument("--planes-only", action="store_true", help="redo only the plane samples and the figure")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    msh = OUT / "design_b_3d.msh"
    t0 = time.time()
    if a.remesh or not msh.exists():
        info = mesh_3d.build(msh, h_channel=a.h_channel, threads=a.threads)
        (OUT / "mesh_info.json").write_text(json.dumps(info, indent=2))
    info = json.loads((OUT / "mesh_info.json").read_text())
    print(f"mesh: {info['tetrahedra']} tetrahedra ({time.time() - t0:.0f} s)", flush=True)

    if a.planes_only:
        summary = json.loads((OUT / "summary.json").read_text())
        planes_figure(OUT, msh, a.threads, summary["V_iron_floating_direct"])
        return
    bases = {"anode": dict(V_anode=1, Iron_floating=0, V_iron=0, V_cathode=0),
             "iron": dict(V_anode=0, Iron_floating=0, V_iron=1, V_cathode=0),
             "cathode": dict(V_anode=0, Iron_floating=0, V_iron=0, V_cathode=1)}
    res, times = {}, {}
    for name, nums in bases.items():
        times[name] = getdp(OUT / name, msh, "Globals", a.threads, **nums)
        res[name] = globals_(OUT / name)
        print(f"basis {name}: {times[name]:.0f} s", flush=True)
    order = ("anode", "iron", "cathode")
    C = np.array([[res[b][f"Q_{c}"] for b in order] for c in order])     # C[i, j] = Q_i for conductor j at 1 V
    k = -C[1, 0] / C[1, 1]
    times["floating"] = getdp(OUT / "floating", msh, "Globals", a.threads, V_anode=V_D, Iron_floating=1, Q_iron=0,
                              V_cathode=0)
    direct = globals_(OUT / "floating")
    print(f"floating iron: superposition {V_D * k:.3f} V, direct {direct['V_iron']:.3f} V", flush=True)
    # cathode 20 V below the chamber wall, anode 300 V above the cathode
    v_a, v_k = V_D - 20.0, -20.0
    v_i_bias = -(C[1, 0] * v_a + C[1, 2] * v_k) / C[1, 1]

    times.update(planes_figure(OUT, msh, a.threads, direct["V_iron"]))

    summary = {"mesh": info, "capacitance_pF": (C * 1e12).round(5).tolist(), "order": order,
               "reciprocity_rel": float(np.max(np.abs(C - C.T)) / np.max(np.abs(C))),
               "k_float": k, "V_iron_floating": V_D * k, "V_iron_floating_direct": direct["V_iron"],
               "case_cathode_minus20": {"V_anode": v_a, "V_cathode": v_k, "V_iron_floating": v_i_bias},
               "solve_times_s": times, "runtime_s": round(time.time() - t0, 1)}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
