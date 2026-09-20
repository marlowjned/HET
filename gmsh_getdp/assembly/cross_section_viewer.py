"""
cross_section_viewer.py

Interactive cross-section |B| viewer for the real BPL-700 assembly solve,
styled after ../../matlab/het_solenoid_bfield.m's launchCrossSectionViewer:
a slider sweeps the cut angle theta (0-180 deg), showing both half-planes
(r<0 and r>0) at once, with a |B| plot, an in-plane (B_r, B_axial) quiver,
and the real discharge channel position drawn as a fixed dashed outline.

Reads assembly_field.vtu (run export_vtk.py first) and uses pyvista's
mesh.slice() to cut it -- an exact geometric plane-mesh intersection, not
an approximation. This replaces an earlier version that inverse-distance-
weight-averaged nearby element centroids on a regular query grid, which
was only ever a stand-in for real interpolation (see git history / prior
README notes) since GetDP's B field is piecewise-constant per element,
not continuously defined, so genuine "interpolation" isn't really the
right idea for it anyway -- an exact cut showing each element's true flat
value where the plane actually crosses it is more correct, not just more
convenient.

A "region of interest" checkbox switches to a renormalized view cropped
to the channel: the iron carries up to ~0.86 T at the pole faces while
the channel itself -- the physically interesting region, since that's
where the field actually acts on the plasma -- runs ~5-32 mT (50-320 G),
one to two orders of magnitude weaker. Sharing one color scale between
them makes the channel look empty; it isn't, it's just invisible next to
iron on that scale. ROI mode masks out triangles outside the channel box
and rescales color (log) to the channel's own local |B| range.

The channel box itself comes from assembly_params.txt's channel_* entries
(the real plasma cavity between the chamber's ceramic walls, NOT the
chamber part's bounding annulus), so correcting the ROI there updates
both the dashed outline and the ROI crop here without touching this file.

Run: python cross_section_viewer.py  (opens an interactive matplotlib window)
"""
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm


def load_params(path="assembly_params.txt"):
    p = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            key = parts[0]
            try:
                p[key] = float(parts[1])
            except ValueError:
                pass  # skip non-numeric rows (region-tag table, winding detail lines)
    return p


print("Loading assembly_field.vtu ...")
pv_mesh = pv.read("assembly_field.vtu")
print(f"  {pv_mesh.n_cells} cells")

params = load_params()
axis_x, axis_z = params["axis_x"], params["axis_z"]
ch_r_in, ch_r_out = params["channel_inner_r"], params["channel_outer_r"]
ch_y_min, ch_y_max = params["channel_y_min"], params["channel_y_max"]
print(f"Axis=({axis_x:.4f},{axis_z:.4f})  channel r=[{ch_r_in:.4f},{ch_r_out:.4f}] "
      f"y=[{ch_y_min:.4f},{ch_y_max:.4f}]")

# Display crop for the full (non-ROI) view -- outer poles sit at ~0.09 m
# radius (see HOW_IT_WORKS.md/plan) -- not the full coarse air domain.
R_MAX = 0.16
Y_MIN, Y_MAX = -0.06, 0.11
# Full-view color ceiling, taken from the loaded solve rather than hardcoded:
# it used to be a fixed 3.0 T, sized for the old linear no-saturation iron
# placeholder that peaked near 9.7 T. Now that the excitation is the real
# ~300 G design point the whole assembly peaks under 1 T, and a fixed 3 T
# scale renders the entire cross-section black.
B_CLIP = float(np.percentile(pv_mesh.cell_data["B_magnitude"], 99.5))
print(f"Full-view color ceiling (99.5th pctl of |B|): {B_CLIP:.3f} T")
ROI_MARGIN = 1.3  # crop margin around the channel box in ROI mode
ROI_VMAX_PCTL = 98  # clip ROI color scale at this percentile of in-box |B|,
                     # so one stray high-field element near a pole/winding
                     # edge doesn't wash out the channel's own pattern


def slice_cut(theta_deg):
    """Exact plane cut of the solved mesh at angle theta_deg / theta_deg+180
    (a full plane through the assembly's central axis, at azimuth theta --
    naturally includes both signed-r half-planes in one cut, same
    convention as MATLAB's viewer). Returns a matplotlib Triangulation in
    (r, y) plot coordinates plus per-triangle |B|, in-plane (B_r, B_y),
    and centroids (for quiver placement)."""
    theta = np.radians(theta_deg)
    normal = (np.sin(theta), 0.0, -np.cos(theta))  # perpendicular to the cut plane
    origin = (axis_x, 0.0, axis_z)  # any point on the axis line is on the plane

    sl = pv_mesh.slice(normal=normal, origin=origin).triangulate()
    pts = sl.points
    faces = sl.faces.reshape(-1, 4)[:, 1:4]  # triangulate() guarantees all-triangle faces

    r = (pts[:, 0] - axis_x) * np.cos(theta) + (pts[:, 2] - axis_z) * np.sin(theta)
    y = pts[:, 1]
    triang = mtri.Triangulation(r, y, triangles=faces)

    B = sl.cell_data["B"]
    Bmag = sl.cell_data["B_magnitude"]
    B_r = B[:, 0] * np.cos(theta) + B[:, 2] * np.sin(theta)
    B_y = B[:, 1]
    cr = r[faces].mean(axis=1)
    cy = y[faces].mean(axis=1)
    return triang, Bmag, B_r, B_y, cr, cy


fig = plt.figure(figsize=(9, 9))
ax = fig.add_axes([0.10, 0.18, 0.68, 0.75])
cax = fig.add_axes([0.82, 0.30, 0.04, 0.55])  # fixed colorbar axis -- see redraw() note

skip = 6
mesh_artist = None
quiv = None


def redraw(theta_deg):
    global mesh_artist, quiv
    triang, Bmag, B_r, B_y, cr, cy = slice_cut(theta_deg)
    roi_on = check.get_status()[0]
    in_roi = (np.abs(cr) >= ch_r_in) & (np.abs(cr) <= ch_r_out) & (cy >= ch_y_min) & (cy <= ch_y_max)

    if mesh_artist is not None:
        mesh_artist.remove()
    if quiv is not None:
        quiv.remove()
    # Clear into a fixed cax rather than repeatedly cb.remove()/fig.colorbar(ax=ax)
    # -- the latter has fig.colorbar borrow space from `ax` each call and resize
    # it back on remove(), which breaks (AttributeError on the 2nd toggle) in
    # this matplotlib version. A dedicated cax sidesteps that entirely.
    cax.cla()

    if roi_on:
        in_box = Bmag[in_roi]
        vmax = float(np.percentile(in_box, ROI_VMAX_PCTL))
        # Log scale, not linear: the channel's own field spans ~2-3 orders
        # of magnitude on its own (sub-mT mid-channel vs. several hundred
        # mT right at the pole-face edges of the box) -- a linear scale
        # over that range makes everything but the very edges look flat
        # black. Floor set from the in-box data itself (2nd percentile,
        # never below 1e-4 T) so the log range tracks what's actually there.
        pos = in_box[in_box > 0]
        vmin = max(1e-4, float(np.percentile(pos, 2)) if len(pos) else 1e-4)
        triang.set_mask(~in_roi)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        cbar_label = f"|B| (T), region-of-interest LOG scale ({vmin:.1e}-{vmax:.3f} T, {ROI_VMAX_PCTL}th pctl)"
    else:
        vmax = B_CLIP
        norm = None
        cbar_label = f"|B| (T), clipped at {B_CLIP:.3f} T (99.5th pctl)"

    mesh_artist = ax.tripcolor(triang, facecolors=Bmag, cmap="inferno", norm=norm,
                                vmin=(None if norm else 0), vmax=(None if norm else vmax), zorder=1)
    cb = fig.colorbar(mesh_artist, cax=cax)
    cb.set_label(cbar_label, fontsize=8)
    if roi_on:
        # Default LogNorm tick formatting collapses to one dense label here
        # (narrow ~2-order-of-magnitude range) -- force a handful of evenly
        # log-spaced, plainly-formatted ticks instead.
        ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
        cb.set_ticks(ticks)
        cb.set_ticklabels([f"{t:.3f}" for t in ticks])

    if roi_on:
        quiv_keep = in_roi
    else:
        # Restrict quiver candidates to the visible crop box before
        # decimating -- otherwise most decimated picks land off-screen
        # (the underlying mesh spans the full ~1.4m air sphere) and a few
        # land just outside the crop edge with a large arrow, showing up
        # as stray diagonal lines cutting across the plot.
        quiv_keep = (np.abs(cr) <= R_MAX) & (cy >= Y_MIN) & (cy <= Y_MAX)
    cr, cy, B_r, B_y = cr[quiv_keep], cy[quiv_keep], B_r[quiv_keep], B_y[quiv_keep]
    quiv = ax.quiver(
        cr[::skip], cy[::skip], B_r[::skip], B_y[::skip],
        color="cyan", scale=vmax * 10, width=0.0025, alpha=0.8, zorder=5,
    )

    if roi_on:
        ax.set_xlim(-ch_r_out * ROI_MARGIN, ch_r_out * ROI_MARGIN)
        ax.set_ylim(ch_y_min - (ch_y_max - ch_y_min) * (ROI_MARGIN - 1),
                    ch_y_max + (ch_y_max - ch_y_min) * (ROI_MARGIN - 1))
    else:
        ax.set_xlim(-R_MAX, R_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)

    title.set_text(f"BPL-700 real-CAD |B| cross-section: theta={theta_deg:.0f} / {theta_deg+180:.0f} deg"
                    + ("  [region of interest]" if roi_on else ""))
    fig.canvas.draw_idle()


# Channel outline: real chamber ceramic-wall position (constant-radius
# annular tube -- see HOW_IT_WORKS.md/plan), same rectangle drawn on both
# signed-r sides, exactly like MATLAB's two rectangle() calls. zorder=6
# keeps these on top of the plot, which is fully re-created (and
# re-added on top of existing artists) on every redraw.
for r0 in (ch_r_in, -ch_r_out):
    rect = Rectangle((r0, ch_y_min), ch_r_out - ch_r_in,
                      ch_y_max - ch_y_min, fill=False, edgecolor="white",
                      linestyle="--", linewidth=1.5, zorder=6)
    ax.add_patch(rect)

ax.set_xlabel("r (m)  [negative r = other half-plane, theta+180]")
ax.set_ylabel("y (m)  [assembly / center-pole axis]")
ax.set_aspect("equal")
title = ax.set_title("BPL-700 real-CAD |B| cross-section")

ax_slider = plt.axes([0.15, 0.05, 0.6, 0.03])
slider = Slider(ax_slider, "theta (deg)", 0, 180, valinit=0, valstep=1)

ax_check = plt.axes([0.8, 0.03, 0.18, 0.07])
check = CheckButtons(ax_check, ["region of interest"], [False])

redraw(0.0)


def on_release(event):
    redraw(slider.val)


def on_check(label):
    redraw(slider.val)


slider.on_changed(lambda v: None)  # no live redraw during drag -- see plan (matches MATLAB's ValueChangedFcn)
fig.canvas.mpl_connect("button_release_event", on_release)
fig.canvas.mpl_connect("key_release_event", on_release)  # arrow-key nudges on a focused slider
check.on_clicked(on_check)

print("Opening interactive viewer -- drag the slider and release to redraw.")
print("Check 'region of interest' to renormalize the color scale to the channel's own field range.")
plt.show()
