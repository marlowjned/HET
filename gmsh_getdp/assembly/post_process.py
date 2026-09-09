"""
post_process.py

Summarizes the GetDP magnetostatics solve (b_assembly.pos, b_iron.pos) and
renders a cross-section |B| plot analogous to het_solenoid_bfield.m's
cross-section viewer, as a quick visual sanity check of the real-CAD
result (field concentrated at the poles, arcing across the gap between
center and outer poles since they're driven at opposite polarity).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pos_utils import parse_pos


print("Parsing b_iron.pos ...")
c_iron, B_iron = parse_pos("b_iron.pos")
mag_iron = np.linalg.norm(B_iron, axis=1)
print(f"  {len(mag_iron)} iron elements")
print(f"  |B| iron: max={mag_iron.max():.3f} T  mean={mag_iron.mean():.3f} T  "
      f"median={np.median(mag_iron):.3f} T  p90={np.percentile(mag_iron,90):.3f} T")

print("Parsing b_assembly.pos (iron+air+windings) ...")
c_all, B_all = parse_pos("b_assembly.pos")
mag_all = np.linalg.norm(B_all, axis=1)
print(f"  {len(mag_all)} total elements")
print(f"  |B| overall: max={mag_all.max():.3f} T")

with open("solve_summary.txt", "w") as f:
    f.write("BPL-700 real-CAD magnetostatics solve summary\n")
    f.write("(linear mur_iron=1000 placeholder, no B-H saturation curve)\n\n")
    f.write(f"Iron elements: {len(mag_iron)}\n")
    f.write(f"  max |B|    = {mag_iron.max():.4f} T\n")
    f.write(f"  mean |B|   = {mag_iron.mean():.4f} T\n")
    f.write(f"  median |B| = {np.median(mag_iron):.4f} T\n")
    f.write(f"  p90 |B|    = {np.percentile(mag_iron,90):.4f} T\n")
    f.write(f"  p99 |B|    = {np.percentile(mag_iron,99):.4f} T\n\n")
    f.write(f"Overall (iron+air+windings) elements: {len(mag_all)}\n")
    f.write(f"  max |B| = {mag_all.max():.4f} T\n")
print("Wrote solve_summary.txt")

# --- Cross-section plot: slab through the assembly center, plotted in the
# X-Z plane (the plane the poles/windings are arranged in -- winding axes
# are along global Y). Color = |B|, quiver = in-plane (Bx,Bz).
cy = c_all[:, 1].mean()
half_extent = (c_all[:, 1].max() - c_all[:, 1].min()) / 2
slab = 0.01  # 10mm slab half-thickness
mask = np.abs(c_all[:, 1] - cy) < slab
print(f"Cross-section slab |y-{cy:.4f}|<{slab}: {mask.sum()} elements")

x = c_all[mask, 0]
z = c_all[mask, 2]
bx = B_all[mask, 0]
bz = B_all[mask, 2]
bmag = mag_all[mask]

fig, ax = plt.subplots(figsize=(8, 8))
sc = ax.scatter(x, z, c=np.clip(bmag, 0, 3.0), cmap="inferno", s=3, vmin=0, vmax=3.0)
skip = max(1, len(x) // 400)
ax.quiver(x[::skip], z[::skip], bx[::skip], bz[::skip], color="cyan", scale=40, width=0.002, alpha=0.7)
cb = fig.colorbar(sc, ax=ax)
cb.set_label("|B| (T), clipped at 3 T")
ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")
ax.set_title(f"BPL-700 real-CAD |B| cross-section (y={cy:.4f}, slab={2*slab*1e3:.0f}mm)")
ax.set_aspect("equal")
zoom = 0.16  # assembly half-extent is ~0.15m -- crop out the coarse far-field air mesh
ax.set_xlim(-zoom, zoom)
ax.set_ylim(-zoom, zoom)
fig.tight_layout()
fig.savefig("cross_section.png", dpi=150)
print("Wrote cross_section.png")
