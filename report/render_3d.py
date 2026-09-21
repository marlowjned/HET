"""
render_3d.py -- cut-away 3D render of the baseline magnet circuit for the
report's opening figure.

Field: the committed nonlinear baseline solve (gmsh_getdp/assembly/
iron_field.vtu, produced by export_vtk.py; solved at
1.25 A) scaled by 1.30/1.25 to the 300 G operating point. The iron is
near-linear there, so the scaling is good to ~1% -- this is a picture,
not a measurement. Coils are drawn as plain annuli from the CAD
dimensions (the .vtu carries no region tags to pull them out by).

Usage (from repo root): python report/render_3d.py
"""
from pathlib import Path

import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "gmsh_getdp" / "assembly"
OUT = ROOT / "report" / "figures" / "render_3d.png"
SCALE = 1.30 / 1.25
CMAP = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
AIR_CMAP = LinearSegmentedColormap.from_list("seq_blue_air", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

pv.OFF_SCREEN = True


N_CUT = (1 / np.sqrt(2), 0.0, -1 / np.sqrt(2))   # plane through the axis and the (+x,+z) outer core


def cut(mesh):
    """Half-section: keep the side of the diagonal plane facing away from the
    camera, so the cut face -- through the inner pole and one outer core --
    faces the viewer."""
    return mesh.clip(normal=N_CUT, origin=(0, 0, 0), invert=True)


def annulus(cx, cz, r_in, r_out, y0, y1):
    return pv.CylinderStructured(center=(cx, (y0 + y1) / 2, cz), direction=(0, 1, 0),
                                 radius=np.linspace(r_in, r_out, 2), height=y1 - y0,
                                 theta_resolution=90, z_resolution=2).extract_surface()


def main():
    iron = pv.read(A / "iron_field.vtu")
    iron["|B| (T)"] = iron["B_magnitude"] * SCALE

    p = pv.Plotter(window_size=(1800, 1300), off_screen=True)
    p.set_background("white")
    common = dict(scalars="|B| (T)", cmap=AIR_CMAP, log_scale=True, clim=(0.02, 1.2))
    p.add_mesh(cut(iron), ambient=0.35, diffuse=0.7,
               scalar_bar_args=dict(title="|B| in iron (T), log scale", color="black",
                                    title_font_size=26, label_font_size=22, n_labels=5, fmt="%.2g",
                                    position_x=0.25, position_y=0.05, width=0.5, height=0.05),
               **common)

    y0, y1 = -0.02299, 0.02781
    coils = [annulus(0, 0, 0.0127, 0.017463, y0, y1)]
    for sx, sz in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        coils.append(annulus(sx * 0.04445, sz * 0.04445, 0.0127, 0.01905, y0, y1))
    for c in coils:
        cc = cut(c.triangulate())
        if cc.n_points:
            p.add_mesh(cc, color="#d98a4e", ambient=0.35, diffuse=0.7)

    # Camera on the removed side, looking back at the cut face, exit (top plate, -y) up.
    p.camera_position = [(0.21, -0.10, -0.17), (0.0, 0.004, 0.0), (0, -1, 0)]
    p.camera.zoom(1.15)
    p.enable_anti_aliasing("ssaa")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    p.screenshot(str(OUT))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
