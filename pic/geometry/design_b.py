"""Design B (shieldless, coils wound forward) in the r-z half-plane, shared by the electrostatic and PIC models.

Units mm. z = 0 at the exit plane, +z downstream (into the plume), r >= 0. Names follow
gmsh_getdp/axisym/figures/glossary.png.

The magnetic-circuit values are those of the het library (sim/library: assembly design_b, which extends
legacy_bn; parts bn_channel_v2, plate_6in_a36, solenoid_230t_awg20, post_coils_184t_awg20), so the magnet model
and these models share one geometry. Items marked PLACEHOLDER are not on any drawing; REPORT.md lists them.
"""
from math import ceil, pi, sqrt

# BN channel part (bn_channel_v2), channel radii at the 2026-10-29 set point ("Changes To: Thruster" slide:
# r_in 0.0205 m, r_out 0.0335 m, width 13 mm; krypton, 3.0 mg/s, 300 V, 3.53 A). Was 22.12 / 31.4452 mm.
# PLACEHOLDER: the rest of the BN part (bore, lip, body, depth) is kept from bn_channel_v2 until the new drawing.
R_IN = 20.5              # channel inner radius
R_OUT = 33.5             # channel outer radius
DEPTH = 41.4274          # exit plane to anode face
BORE_R = 17.7165         # bore radius (the pole cap sits in the bore)
LIP_R, LIP_LEN = 38.1, 9.525
BODY_R = 44.45
BN_LEN = 47.7774         # exit plane to the part's back face

# iron (all parts het:materials/iron_legacy in the magnet model)
PLATE_T = 9.525
PLATE_R = 152.4 / sqrt(pi)        # 6 in square plates -> discs of equal area, as in the magnet model
FRONT_BORE_R = 38.1
BACK_Z = -50.8                    # back plate front face
BACK_BORE_R = 3.4
CAP_R, CAP_T, CAP_EDGE, CAP_CHAMFER = 17.6, 9.525, 1.5, 6.0   # pole cap: radius, thickness, edge face length, chamfer width
CORE_R = 9.525
POST_RING_R, POST_R, POST_COUNT = 63.5, 9.525, 4

# coils: AWG 20 MW 16-C, 0.89 mm insulated; layered windings as in the het part types
WIRE_D = 0.89
INNER_TURNS, INNER_PER_LAYER = 230, 29
OUTER_TURNS, OUTER_PER_LAYER = 184, 21
INNER_BUILD = ceil(INNER_TURNS / INNER_PER_LAYER) * WIRE_D     # 8 layers: 7.12
INNER_LEN = INNER_PER_LAYER * WIRE_D                           # 25.81
OUTER_BUILD = ceil(OUTER_TURNS / OUTER_PER_LAYER) * WIRE_D     # 9 layers: 8.01
OUTER_LEN = OUTER_PER_LAYER * WIRE_D                           # 18.69

# PLACEHOLDER anode / gas distributor: a stainless ring filling the bottom of the channel slot. BN_REPORT: "the
# bottom 6.35 mm of the channel is filled, so the anode face is 6.35 mm above the part's back face".
ANODE_R_IN, ANODE_R_OUT = R_IN, R_OUT
ANODE_Z0, ANODE_Z1 = -BN_LEN, -DEPTH                           # 6.35 mm thick

# PLACEHOLDER hollow cathode (outside the thruster, off axis; only the 3-D model can place it)
CATHODE_R, CATHODE_Z = 110.0, 20.0     # orifice position (r, z) of a cathode beside the front plate


def rect(r0, r1, z0, z1):
    return [(r0, z0), (r1, z0), (r1, z1), (r0, z1)]


def parts(posts="ring"):
    """name -> (kind, polygon [(r, z), ...] in mm). kind: 'bn', 'iron', 'coil', 'anode'.

    posts='ring': the four outer posts smeared into a ring of the same radial extent (as the 2-D magnet model).
    The BN part is split into its inner ring and outer body: with the anode filling the bottom of the slot, the
    slot runs through the part's back face (BN_REPORT: 0.375 + 1.506 in deep = the part length)."""
    cap = [(0, 0), (CAP_R, 0), (CAP_R, -CAP_EDGE), (CAP_R - CAP_CHAMFER, -CAP_T), (0, -CAP_T)]
    p = {
        "bn_inner": ("bn", rect(BORE_R, R_IN, -BN_LEN, 0)),
        "bn_outer": ("bn", [(R_OUT, 0), (LIP_R, 0), (LIP_R, -LIP_LEN), (BODY_R, -LIP_LEN), (BODY_R, -BN_LEN),
                            (R_OUT, -BN_LEN)]),
        "anode": ("anode", rect(ANODE_R_IN, ANODE_R_OUT, ANODE_Z0, ANODE_Z1)),
        "front_plate": ("iron", rect(FRONT_BORE_R, PLATE_R, -PLATE_T, 0)),
        "back_plate": ("iron", rect(BACK_BORE_R, PLATE_R, BACK_Z - PLATE_T, BACK_Z)),
        "pole_cap": ("iron", cap),
        "inner_core": ("iron", rect(0, CORE_R, BACK_Z, -CAP_T)),
        "inner_coil": ("coil", rect(CORE_R, CORE_R + INNER_BUILD, -CAP_T - INNER_LEN, -CAP_T)),
    }
    if posts == "ring":
        p["outer_posts"] = ("iron", rect(POST_RING_R - POST_R, POST_RING_R + POST_R, BACK_Z, -PLATE_T))
        a, b = POST_RING_R - POST_R, POST_RING_R + POST_R
        p["outer_coils_in"] = ("coil", rect(a - OUTER_BUILD, a, -PLATE_T - OUTER_LEN, -PLATE_T))
        p["outer_coils_out"] = ("coil", rect(b, b + OUTER_BUILD, -PLATE_T - OUTER_LEN, -PLATE_T))
    return p


def channel_mid_radius():
    return 0.5 * (R_IN + R_OUT)


def summary():
    lines = [f"{k:16s} {kind:6s} r {min(x for x, _ in poly):8.3f} .. {max(x for x, _ in poly):8.3f}   "
             f"z {min(z for _, z in poly):9.4f} .. {max(z for _, z in poly):8.4f}"
             for k, (kind, poly) in parts().items()]
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
    print(f"inner coil build {INNER_BUILD:.3f} length {INNER_LEN:.3f}; outer coil build {OUTER_BUILD:.3f} "
          f"length {OUTER_LEN:.3f}; plate radius {PLATE_R:.3f}; anode thickness {ANODE_Z1 - ANODE_Z0:.4f}")
