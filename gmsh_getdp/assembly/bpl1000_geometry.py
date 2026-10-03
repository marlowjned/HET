"""
BPL-1000 fixed geometry, measured 2026-10-03 from the Onshape STEP at the
declared defaults -- the single home for every CAD constant the sweep
drivers need, replacing emag_sweep.py's BPL-700 literals.

How it was measured: every face of every body was sampled on the B-rep
itself (gmsh.model.getValue over the face parametrisation, filtered by
isInside), and each lateral face got a least-squares circle fit. All the
cylinders came back with sub-micron residuals and land on round inch values,
so the numbers below are exact, not mesh estimates. A body-centroid axis is
NOT good enough on this STEP: occ.getCenterOfMass on its B-spline solids is
~0.17 mm off, which shows up as a spurious coil/core overlap.

Parametrisation, confirmed by single-parameter perturbations of the STEP:
- coil ID = coreod exactly; coil OD = coreod + coilt, i.e. `coilt` is a
  DIAMETRAL increment and the radial build is coilt/2.
- The inner pole cap is a fixed 1.375 in diameter disc, `iem_platet` thick,
  whose front face is pinned to the exit plane. Thinning it lengthens the
  inner coil and shaft downward by the same amount; nothing else moves.
- `emag_height` moves the back plate; the front plate, chamber and anode
  are fixed. Coil height = emag_height for the outer coils.

Coordinates: thruster axis is global Y, exit toward -Y. y values in metres.
"""
IN = 0.0254

Y_EXIT = -0.037846            # front-plate front face == pole-cap front face == chamber exit
Y_FRONT_PLATE_INNER = -0.028321   # outer coils/cores start here
PLATE_T_IN = 0.375            # front and back plates

# Discharge channel: the chamber's annular cavity, from the anode's front
# face (anode_axial, 316L, inert) to the exit plane.
CHANNEL = dict(channel_inner_r=0.022120, channel_outer_r=0.031600,
               channel_y_anode=-0.004826, channel_y_exit=Y_EXIT)

CHAMBER_BORE_R = 0.017590     # 1.385 in dia: bounds the inner EM's OD
CHAMBER_OUTER_R = 0.041021    # 3.230 in dia: bounds the outer coils' inner reach
INNER_POLE_R = 0.0174625      # 1.375 in dia, fixed
OUTER_CORE_AXIS_R = 0.0628618 # |(44.45, 44.45)| mm, 4 outer cores
CORE_BOLT_HOLE_R = 0.002553   # tapped holes in the core ends (not swept)

# Fit limits on (coreod + coilt), i.e. the coil OD, in inches.
INNER_OD_MAX_IN = 2 * CHAMBER_BORE_R / IN                          # 1.385
OUTER_OD_MAX_IN = 2 * (OUTER_CORE_AXIS_R - CHAMBER_OUTER_R) / IN   # 1.720

SPACER_ZERO_IN = 1.625        # emag_height at which the spacer body vanishes


def coil_top_y(emag_height_in):
    return Y_FRONT_PLATE_INNER + emag_height_in * IN


def inner_coil_bottom_y(iem_platet_in):
    return Y_EXIT + iem_platet_in * IN


def coil_heights_in(emag_height_in, iem_platet_in):
    """(inner, outer) winding heights in inches."""
    inner = emag_height_in + PLATE_T_IN - iem_platet_in
    return inner, emag_height_in
