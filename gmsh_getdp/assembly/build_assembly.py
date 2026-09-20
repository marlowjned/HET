"""
build_assembly.py

Imports the BPL-700 magnet-circuit assembly directly from a whole-assembly
Onshape STEP export (via onshape.cache.get_step_for_config), for a given
5-parameter geometry configuration (emag_height, inner/outer_coil_id/od --
all real Onshape Configuration parameters, see onshape/cache.py's
docstring for why all 5 now go through the API instead of 4 being resized
locally). This replaces the old per-part-STEP + hardcoded-TRANSFORMS
pipeline; promotes the incremental-fuse approach validated in
scratch_whole_assembly_import.py.

Body naming, confirmed via a real Onshape pull (2026-09-18, inspecting the
imported STEP's entity names directly -- not assumed):
  iron:    top_plate (x1), bottom_plate (x1), inner_emag_core (x1),
           outer_emag_core (x4). Also chamber_spacer (x1) -- its magnetic
           relevance is genuinely unconfirmed (not in the old per-part
           pipeline at all); treated as iron for now per instruction "we
           can try both" -- flip INCLUDE_CHAMBER_SPACER_AS_IRON to test
           the inert alternative.
  coil:    inner_coil (x1), outer_coil (x4) -- real winding solids exported
           by Onshape at exactly the configured ID/OD. No more
           occ.addCylinder/occ.cut primitive construction.
  other:   chamber (x1) -- magnetically inert (ceramic), kept only for the
           channel-outline metadata block at the end of this script, never
           added to the FEM domain (same reasoning as the old pipeline,
           see the "Known geometry exclusions" note below).
  injector is NOT present in this whole-assembly export at all (confirmed
  missing). Fine -- it was already excluded from the FEM solid model for
  being magnetically inert (mu_r ~= 1), so its absence changes nothing.

Known geometry exclusions (same physics reasoning as every prior version
of this pipeline): chamber is ceramic, magnetically inert -- Ampere's law
doesn't distinguish mu_r=1 from empty space, so leaving it unmeshed and
letting that volume fall into Air is physically correct, not a shortcut.

Winding axis convention: every winding/core occurrence's own rotational
axis is global Y -- confirmed from real geometry (each coil's X/Z bounding-
box span matches its own OD; sanity-asserted below, not just assumed).
Current-density SIGN per winding is a wiring-direction design choice, not
something recoverable from a solid's geometry (a hollow cylinder of
revolution has no handedness -- flipping it 180 degrees about any diameter
looks identical). Unlike the old TRANSFORMS-driven pipeline -- which had
outer poles at two different heights needing a per-occurrence sign flip --
the 4 outer_coil occurrences here are true 90-degree-rotated copies at THE
SAME y-span (a ring around the center pole), so sign is now one constant
per pole TYPE (all outer poles same polarity, center opposite), see
EXCITATION below.

A_cross (winding cross-section area, for the Js[] current-density
formula) uses the plain annulus formula from the INPUT id/od parameters
directly (pi*((od/2)^2-(id/2)^2)), not a measurement of the imported
solid's actual volume/height -- those parameters are the exact values sent
to Onshape, so they're more trustworthy than re-deriving area from a solid
that might carry small fillets/rounds Onshape adds at real edges.
"""
import argparse
import math
import os
import sys
from pathlib import Path

import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root, for `import onshape`
from onshape.cache import get_step_for_config  # noqa: E402
from generate_regions import generate_regions_pro  # noqa: E402
from params_utils import parse_params  # noqa: E402

STEP_SCALE = 0.001  # mm-labeled-as-m bug, see HOW_IT_WORKS.md sec 3 -- still present in this export
INCH_TO_M = 0.0254

# Onshape's own live defaults (GET .../configuration), used as this script's
# CLI defaults so `python build_assembly.py` with no args reproduces "the"
# baseline configuration.
ONSHAPE_DEFAULTS_IN = dict(
    emag_height=2.0, inner_coil_id=1.0, inner_coil_od=1.375,
    outer_coil_id=1.0, outer_coil_od=1.5,
)

IRON_CORE_NAMES = {"top_plate", "bottom_plate", "inner_emag_core", "outer_emag_core"}
CHAMBER_SPACER_NAME = "chamber_spacer"
CHAMBER_NAME = "chamber"
COIL_NAMES = {"inner_coil", "outer_coil"}

# Current sized for the ~300 G channel-exit target (field_quality.py's
# target definition), replacing the inherited 5 A placeholder: a lumped
# reluctance check (23.9 kA/m across the 37.4mm inner-core-to-outer-core
# air path, iron drop negligible at ~0.1 T core flux) needs only ~450-900
# A-turns, vs. the 2500 A-turns 5 A was supplying. The linear FEM sweep
# agrees independently (channel |B| mean 80.6 mT at 5 A -> 1.86 A for
# 300 G). Driving 5 A also put peak iron |B| at 2.16 T, past the B-H knee,
# which is what forced the 16-stage load-stepping the nonlinear solve
# needed; at 2 A the peak is ~0.86 T and the iron stays effectively linear.
# Turns unchanged. Opposite polarity so flux arcs across the gap
# between inner and outer poles instead of just adding axially -- see
# HOW_IT_WORKS.md sec. 5 for the physical reasoning.
EXCITATION = {
    "inner_coil": dict(turns=300, current=2.0, polarity=-1),
    "outer_coil": dict(turns=200, current=2.0, polarity=+1),
}


def detect_channel_cavity(pts, axis_x: float, axis_z: float, bin_m: float = 2.5e-4):
    """Locate the open plasma cavity inside the chamber's annular ceramic.

    `pts` must be chamber VOLUME element centroids, not surface nodes: a
    surface mesh puts no points inside solid material, so the interior of a
    ceramic wall looks identical to open space and this picks the wall
    instead of the cavity (hit empirically -- it returned r 31.98-40.94mm,
    the outer wall's interior).

    Taking the chamber's plain min/max radius (what this script did
    originally) returns the whole part's bounding annulus -- both ceramic
    walls plus, at the exit end, volume that real iron occupies -- which is
    NOT the plasma channel: measured on the 2026-09-18 geometry that mask
    was r 17.30-41.44mm and contained 5392 iron elements, inflating channel
    field stats by sampling pole material.

    The chamber is a solid of revolution, so classify each radial bin by
    what fraction of the axial range has material in it: a wall spans
    nearly the full length, while the cavity has material only where the
    ceramic floor closes the annulus at the anode end. The cavity is the
    widest contiguous run of low-axial-coverage bins.

    Returns (r_in, r_out, y_closed, y_open) -- cavity radii, and the axial
    positions of the closed (anode) and open (exit) ends.
    """
    import numpy as np

    r = np.hypot(pts[:, 0] - axis_x, pts[:, 2] - axis_z)
    y = pts[:, 1]
    n_r = max(8, int(math.ceil((r.max() - r.min()) / bin_m)))
    r_edges = np.linspace(r.min(), r.max(), n_r + 1)
    y_edges = np.linspace(y.min(), y.max(), 21)

    r_idx = np.clip(np.digitize(r, r_edges) - 1, 0, n_r - 1)
    y_idx = np.clip(np.digitize(y, y_edges) - 1, 0, len(y_edges) - 2)
    coverage = np.zeros(n_r)
    for i in range(n_r):
        sel = r_idx == i
        coverage[i] = len(np.unique(y_idx[sel])) / (len(y_edges) - 1) if sel.any() else 0.0

    is_cavity = coverage < 0.5
    best_len, best = 0, None
    i = 0
    while i < n_r:
        if is_cavity[i]:
            j = i
            while j < n_r and is_cavity[j]:
                j += 1
            if j - i > best_len:
                best_len, best = j - i, (i, j)
            i = j
        else:
            i += 1
    assert best is not None and best_len >= 4, (
        f"no plasma cavity found in the chamber's radial profile (coverage={coverage!r}) -- "
        f"chamber geometry may not be the expected annular liner"
    )
    # Trim edge bins that straddle a wall face: the 0.5 run threshold is
    # deliberately loose so the run is found at all, but a partially-filled
    # boundary bin left in would put ceramic inside the "cavity".
    lo, hi = best
    while lo < hi - 1 and coverage[lo] > 0.15:
        lo += 1
    while hi > lo + 1 and coverage[hi - 1] > 0.15:
        hi -= 1
    r_in, r_out = float(r_edges[lo]), float(r_edges[hi])

    # The only chamber material at cavity radii is the floor closing the
    # annulus; whichever axial end it sits at is the anode, the other is
    # the open exit. Probe only the middle 60% of the cavity width -- at
    # the full width a fraction of a bin of wall overlap reads as material
    # at every y and destroys the result (hit empirically: gave a 0.01mm
    # long channel).
    pad = 0.2 * (r_out - r_in)
    in_cav = (r > r_in + pad) & (r < r_out - pad)
    assert in_cav.any(), "cavity band contains no chamber nodes -- cannot locate the closed end"
    y_floor = y[in_cav]
    if abs(y_floor.mean() - y.max()) < abs(y_floor.mean() - y.min()):
        y_closed, y_open = float(y_floor.min()), float(y.min())
    else:
        y_closed, y_open = float(y_floor.max()), float(y.max())
    return r_in, r_out, y_closed, y_open


def base_name(entity_name: str) -> str:
    """Gmsh/OCC entity names come back as e.g.
    'Shapes/BPL-700 Assembly/outer_coil/outer_coil' -- last path segment
    is the real Onshape part name."""
    return entity_name.rstrip("/").split("/")[-1]


def build(params_in: dict, out_dir: str, step_path: str | None,
          include_chamber_spacer_as_iron: bool) -> None:
    os.makedirs(out_dir, exist_ok=True)

    if step_path is None:
        step_path = str(get_step_for_config(params_in))
    print(f"Using STEP: {step_path}")

    gmsh.initialize()
    gmsh.model.add("bpl700_assembly")
    occ = gmsh.model.occ

    dimtags = occ.importShapes(step_path)
    occ.synchronize()
    occ.dilate(dimtags, 0, 0, 0, STEP_SCALE, STEP_SCALE, STEP_SCALE)
    occ.synchronize()
    print(f"Imported {len(dimtags)} solid(s)")

    # --- Classify solids by name --------------------------------------------
    iron_dimtags = []
    coil_entries = []  # (base_name, dimtag)
    chamber_dimtag = None
    for d, t in dimtags:
        name = base_name(gmsh.model.getEntityName(d, t))
        if name in IRON_CORE_NAMES:
            iron_dimtags.append((d, t))
        elif name == CHAMBER_SPACER_NAME:
            if include_chamber_spacer_as_iron:
                iron_dimtags.append((d, t))
            # else: leave unclassified -> falls into Air, same as chamber
        elif name in COIL_NAMES:
            coil_entries.append((name, (d, t)))
        elif name == CHAMBER_NAME:
            assert chamber_dimtag is None, "expected exactly 1 chamber solid"
            chamber_dimtag = (d, t)
        else:
            raise AssertionError(f"unrecognized body name {name!r} (tag {t}) -- "
                                  f"update IRON_CORE_NAMES/COIL_NAMES/CHAMBER_* above")

    n_expected_iron = 7 + (1 if include_chamber_spacer_as_iron else 0)
    assert len(iron_dimtags) == n_expected_iron, \
        f"expected {n_expected_iron} iron bodies, got {len(iron_dimtags)}"
    assert len(coil_entries) == 5, f"expected 5 coil bodies, got {len(coil_entries)}"
    assert chamber_dimtag is not None, "expected 1 chamber solid for channel-outline metadata"
    print(f"Classified: {len(iron_dimtags)} iron, {len(coil_entries)} coil, 1 chamber (metadata only)")

    # --- Windings: real coil solids straight from the STEP, no primitive ---
    # construction needed. Winding0 = inner_coil (center), Winding1..4 =
    # the 4 outer_coil occurrences (order among the 4 is arbitrary --
    # symmetric by design). Measured BEFORE the iron fuse below: OCC's
    # boolean ops can renumber/retag entities not even involved in the
    # operation, so bounding-box/center-of-mass queries on dimtags held
    # across a fuse call aren't reliable (hit empirically -- see git
    # history of this script for the stale-tag bbox mismatch this caused).
    coil_entries.sort(key=lambda e: 0 if e[0] == "inner_coil" else 1)
    winding_tags = [dt for _, dt in coil_entries]
    winding_pole = [name for name, _ in coil_entries]  # "inner_coil"/"outer_coil" per winding

    # NOTE on occ.getBoundingBox(): confirmed empirically (see git history of
    # this line) that it returns a LOOSE/conservative box on these curved,
    # STEP-imported solids -- e.g. inner_coil's true X/Z extent (measured via
    # actual mesh nodes) is exactly its OD (34.9mm) on both axes, but
    # getBoundingBox() reports 34.9mm on X and 60.5mm on Z for the SAME
    # solid. This is a known OCC quirk (BRepBndLib's box bounds a curved
    # face's underlying surface, not its trimmed boundary tightly) -- it is
    # NOT a stale-tag artifact from the iron fuse below (reproduced before
    # any fuse call at all), so don't use occ.getBoundingBox() for anything
    # needing X/Z precision on these bodies. Y (the extrusion axis, bounded
    # by flat annular end-caps, not curved surface) IS reliable -- only the
    # curved lateral surface's X/Z bound is loose. getMass()/getCenterOfMass()
    # are volume integrals over the true trimmed solid, not boundary-based,
    # and are trustworthy (cross-checked against the same mesh-node data).
    winding_meta = []  # (index, pole_name, A_cross, loc_xz)
    for i, (name, dt) in enumerate(zip(winding_pole, winding_tags)):
        d, t = dt
        _, ymin, _, _, ymax, _ = occ.getBoundingBox(d, t)
        height = ymax - ymin
        id_in = params_in["inner_coil_id"] if name == "inner_coil" else params_in["outer_coil_id"]
        od_in = params_in["inner_coil_od"] if name == "inner_coil" else params_in["outer_coil_od"]
        A_cross = math.pi * (((od_in * INCH_TO_M) / 2) ** 2 - ((id_in * INCH_TO_M) / 2) ** 2)
        vol = occ.getMass(d, t)
        vol_expected = A_cross * height
        # 5% tolerance, not tighter: real CAD coil solids carry small edge
        # fillets/rounds the ideal-annulus formula doesn't capture (~1%
        # discrepancy observed on the real geometry) -- this check is a
        # sanity check on the constant-cross-section/axis-along-Y
        # assumption, not a precision validation (A_cross itself always
        # comes from the exact analytic id/od values, not this measurement).
        assert abs(vol - vol_expected) / vol_expected < 0.05, (
            f"winding{i} ({name}) volume ({vol:.6e}) doesn't match the analytic annulus "
            f"volume from id={id_in}in/od={od_in}in ({vol_expected:.6e}) -- not a simple "
            f"constant-cross-section extrusion along Y as assumed"
        )
        cx, _, cz = occ.getCenterOfMass(d, t)
        winding_meta.append(dict(index=i, pole=name, A_cross=A_cross, loc=(cx, cz)))
    print(f"Built {len(winding_tags)} winding entries (from real CAD, no Gmsh primitives)")

    # --- Fuse iron bodies INCREMENTALLY -------------------------------------
    # One-shot occ.fuse(iron[0], iron[1:]) gives 3 disjoint solids on this
    # geometry (OCC's fuzzy boolean union is order/grouping-dependent) --
    # see scratch_whole_assembly_import.py's docstring. Incremental fuse
    # (fold one part in at a time) gives the correct single solid.
    gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-6)
    iron_current = [iron_dimtags[0]]
    for dt in iron_dimtags[1:]:
        iron_current, _ = occ.fuse(iron_current, [dt])
        occ.synchronize()
    assert len(iron_current) == 1, f"iron fuse did not produce a single solid: {iron_current}"
    iron_tag = iron_current[0]
    print(f"Iron fuse: {len(iron_dimtags)} parts -> 1 solid ({iron_tag})")

    # --- Outer air domain ----------------------------------------------------
    all_solid_tags = [iron_tag] + winding_tags
    bbox_pts = [occ.getBoundingBox(d, t) for (d, t) in all_solid_tags]
    xmin = min(b[0] for b in bbox_pts); ymin = min(b[1] for b in bbox_pts); zmin = min(b[2] for b in bbox_pts)
    xmax = max(b[3] for b in bbox_pts); ymax = max(b[4] for b in bbox_pts); zmax = max(b[5] for b in bbox_pts)
    cx, cy, cz = (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2
    extent = math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2) / 2
    air_radius = 8 * extent
    print(f"Assembly bbox center=({cx:.4f},{cy:.4f},{cz:.4f}) half-extent={extent:.4f} air_radius={air_radius:.4f}")

    air_sphere = occ.addSphere(cx, cy, cz, air_radius)
    occ.synchronize()

    # --- Fragment everything together (iron blob, 5 windings, air sphere) --
    frag_inputs = [(3, air_sphere), iron_tag] + winding_tags
    frag_names = ["AirSphere", "Iron"] + [f"Winding{i}" for i in range(len(winding_tags))]

    out, out_map = occ.fragment([frag_inputs[0]], frag_inputs[1:])
    occ.synchronize()

    print("\n--- Fragment result ---")
    for nm, inp, res in zip(frag_names, frag_inputs, out_map):
        print(f"{nm:12s} {inp} -> {res}")
        assert len(res) == 1 or nm == "AirSphere", f"{nm} split unexpectedly: {res}"

    vol_before = sum(occ.getMass(d, t) for d, t in frag_inputs)
    vol_after = sum(occ.getMass(d, t) for d, t in gmsh.model.getEntities(3))
    print(f"\nTotal volume before fragment: {vol_before:.6e}  after: {vol_after:.6e}")
    # 1e-4 relative tolerance, not the old pipeline's 1e-9: that tolerance
    # was tuned for simple-primitive geometry (Gmsh-built cylinders + a
    # sphere). Real curved CAD fused/fragmented with a fuzzy boolean
    # tolerance (Geometry.ToleranceBoolean=1e-6 above) shows a small but
    # real ~1e-5 relative volume residual here -- still a meaningful
    # conservation check, just not exact-to-machine-precision on this input.
    assert abs(vol_before - vol_after) / vol_before < 1e-4

    # --- Physical groups ------------------------------------------------------
    solid_tags_final = {"Iron": iron_tag[1]}
    for i, w in enumerate(winding_tags):
        solid_tags_final[f"Winding{i}"] = w[1]

    TAG = {"Air": 1, "Iron": 2,
           "Winding0": 3, "Winding1": 4, "Winding2": 5, "Winding3": 6, "Winding4": 7,
           "OuterBnd": 8}

    air_only_tag = [t for (d, t) in out_map[0] if t not in solid_tags_final.values()]
    assert len(air_only_tag) == 1, f"expected 1 leftover air volume, got {air_only_tag}"

    for name, vtag in solid_tags_final.items():
        pg = gmsh.model.addPhysicalGroup(3, [vtag], tag=TAG[name])
        gmsh.model.setPhysicalName(3, pg, name)
    pg_air = gmsh.model.addPhysicalGroup(3, air_only_tag, tag=TAG["Air"])
    gmsh.model.setPhysicalName(3, pg_air, "Air")

    air_boundary = gmsh.model.getBoundary([(3, air_only_tag[0])], oriented=False)
    embedded_boundary_tags = set()
    for vtag in list(solid_tags_final.values()):
        b = gmsh.model.getBoundary([(3, vtag)], oriented=False)
        embedded_boundary_tags.update(t for (d, t) in b)
    outer_surf_tags = [t for (d, t) in air_boundary if t not in embedded_boundary_tags]
    print(f"\nOuter boundary surfaces: {outer_surf_tags}")
    pg_outer = gmsh.model.addPhysicalGroup(2, outer_surf_tags, tag=TAG["OuterBnd"])
    gmsh.model.setPhysicalName(2, pg_outer, "OuterBnd")

    print("\nPhysical group tags (for magnetostatics_assembly.pro):")
    for name, t in sorted(TAG.items(), key=lambda kv: kv[1]):
        print(f"  {t:2d}  {name}")

    # --- Mesh ------------------------------------------------------------
    min_winding_wall = min(
        (params_in["inner_coil_od"] - params_in["inner_coil_id"]) / 2,
        (params_in["outer_coil_od"] - params_in["outer_coil_id"]) / 2,
    ) * INCH_TO_M
    gmsh.option.setNumber("Mesh.MeshSizeMin", min_winding_wall / 2)
    gmsh.option.setNumber("Mesh.MeshSizeMax", air_radius / 4)

    refine_surfs = list(embedded_boundary_tags)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", refine_surfs)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", min_winding_wall / 2)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", air_radius / 4)
    gmsh.model.mesh.field.setNumber(2, "DistMin", min_winding_wall)
    gmsh.model.mesh.field.setNumber(2, "DistMax", extent * 1.5)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    msh_path = os.path.join(out_dir, "assembly.msh")
    gmsh.write(msh_path)

    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    print(f"\nFinal mesh: {n_nodes} nodes -> {msh_path}")

    gmsh.finalize()

    # --- Channel-outline metadata for cross_section_viewer.py ----------------
    # chamber is magnetically inert and not part of the FEM domain, but its
    # real position is still useful as a plotted reference -- the discharge
    # channel is physically where chamber's ceramic wall is. Re-import the
    # SAME whole-assembly STEP in a completely fresh gmsh session (not
    # reusing the one above -- calling mesh.generate() again after a
    # background mesh-size field was set and a full 3D mesh already written
    # crashed the native library outright, same finding as the old script),
    # pull out just the chamber body, and profile its radius vs. axial
    # position about the center pole's own axis (global Y, per the winding-
    # axis convention established above).
    axis_x, axis_z = winding_meta[0]["loc"]  # inner_coil's axis = the assembly's central axis

    gmsh.initialize()
    gmsh.model.add("chamber_profile")
    occ2 = gmsh.model.occ
    dimtags2 = occ2.importShapes(step_path)
    occ2.synchronize()
    occ2.dilate(dimtags2, 0, 0, 0, STEP_SCALE, STEP_SCALE, STEP_SCALE)
    occ2.synchronize()
    chamber2 = next(dt for dt in dimtags2 if base_name(gmsh.model.getEntityName(*dt)) == CHAMBER_NAME)
    # Remove every other imported solid before meshing -- gmsh.model.mesh.
    # getNodes() with no args returns EVERY meshed entity's nodes, not just
    # chamber's (hit empirically: without this, channel_inner_r came out
    # ~0, picking up points from inner_emag_core's on-axis centerline
    # instead of chamber's real ~17mm inner wall).
    occ2.remove([dt for dt in dimtags2 if dt != chamber2], recursive=True)
    occ2.synchronize()
    # Volume mesh, not just the surface: detect_channel_cavity() needs
    # points inside solid ceramic to tell a wall's interior from open
    # space. Size-capped so the ~4mm walls get several elements across.
    import numpy as np
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0015)
    gmsh.model.mesh.generate(3)
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coord_of = dict(zip(node_tags, np.array(node_coords).reshape(-1, 3)))
    _, elem_tags, elem_nodes = gmsh.model.mesh.getElements(3)
    tets = np.array(elem_nodes[0]).reshape(-1, 4)
    chamber_pts = np.array([[coord_of[n] for n in tet] for tet in tets]).mean(axis=1)
    gmsh.finalize()

    # The plasma cavity, NOT the chamber part's bounding annulus -- see
    # detect_channel_cavity()'s docstring for why the difference matters.
    channel_inner_r, channel_outer_r, channel_y_anode, channel_y_exit = \
        detect_channel_cavity(chamber_pts, axis_x, axis_z)
    channel_y_min = min(channel_y_anode, channel_y_exit)
    channel_y_max = max(channel_y_anode, channel_y_exit)
    print(f"\nChannel cavity (from chamber geometry): inner_r={channel_inner_r:.4f} "
          f"outer_r={channel_outer_r:.4f} anode_y={channel_y_anode:.4f} exit_y={channel_y_exit:.4f} "
          f"about axis (x={axis_x:.4f}, z={axis_z:.4f})")

    params_path = os.path.join(out_dir, "assembly_params.txt")
    with open(params_path, "w") as f:
        for name, t in sorted(TAG.items(), key=lambda kv: kv[1]):
            f.write(f"{t} {name}\n")
        f.write(f"n_windings {len(winding_tags)}\n")
        for w in winding_meta:
            f.write(f"winding{w['index']} pole={w['pole']} A_cross={w['A_cross']!r} "
                     f"loc=({w['loc'][0]!r}, {w['loc'][1]!r})\n")
        f.write(f"axis_x {axis_x!r}\n")
        f.write(f"axis_z {axis_z!r}\n")
        f.write(f"channel_inner_r {channel_inner_r!r}\n")
        f.write(f"channel_outer_r {channel_outer_r!r}\n")
        f.write(f"channel_y_min {channel_y_min!r}\n")
        f.write(f"channel_y_max {channel_y_max!r}\n")
        # Orientation, so downstream analysis doesn't have to infer which
        # end is the exit (field_quality.py's anode->exit profile depends
        # on it, and getting it backwards silently inverts the metric).
        f.write(f"channel_y_anode {channel_y_anode!r}\n")
        f.write(f"channel_y_exit {channel_y_exit!r}\n")
    print(f"Wrote {params_path}")

    # --- Generate the GetDP region/current-source include file ---------------
    # Delegated to generate_regions.py -- it used to be duplicated inline
    # here, but that stopped being reasonable once Iron's material
    # definition grew to a ~100-number B-H curve (see that module's
    # docstring); re-parsing the params file we just wrote is cheap and
    # keeps one source of truth instead of two copies that must stay in
    # sync by hand.
    pro_path = os.path.join(out_dir, "assembly_regions_generated.pro")
    generate_regions_pro(parse_params(params_path), EXCITATION, out_path=pro_path)
    print(f"Wrote {pro_path} ({len(winding_tags)} windings)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emag-height", type=float, default=ONSHAPE_DEFAULTS_IN["emag_height"])
    ap.add_argument("--inner-coil-id", type=float, default=ONSHAPE_DEFAULTS_IN["inner_coil_id"])
    ap.add_argument("--inner-coil-od", type=float, default=ONSHAPE_DEFAULTS_IN["inner_coil_od"])
    ap.add_argument("--outer-coil-id", type=float, default=ONSHAPE_DEFAULTS_IN["outer_coil_id"])
    ap.add_argument("--outer-coil-od", type=float, default=ONSHAPE_DEFAULTS_IN["outer_coil_od"])
    ap.add_argument("--out-dir", type=str, default=".")
    ap.add_argument("--step-path", type=str, default=None,
                     help="Bypass Onshape entirely and import this STEP file directly "
                          "(must still match the ID/OD args above for A_cross/mesh sizing).")
    ap.add_argument("--exclude-chamber-spacer", action="store_true",
                     help="Treat chamber_spacer as magnetically inert (excluded from the iron "
                          "fuse) instead of the default iron treatment -- its real material "
                          "isn't confirmed yet, see module docstring.")
    args = ap.parse_args()

    params_in = dict(
        emag_height=args.emag_height,
        inner_coil_id=args.inner_coil_id,
        inner_coil_od=args.inner_coil_od,
        outer_coil_id=args.outer_coil_id,
        outer_coil_od=args.outer_coil_od,
    )
    build(params_in, args.out_dir, args.step_path,
          include_chamber_spacer_as_iron=not args.exclude_chamber_spacer)


if __name__ == "__main__":
    main()
