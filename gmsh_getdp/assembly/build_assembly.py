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
           outer_emag_core (x4). Also chamber_spacer (x1) -- confirmed
           2026-09-21 to be the same A36 as the plates, so iron is the
           default; --exclude-chamber-spacer still tests the inert case.
           Its height is emag_height - 1.125in; at 1.125in Onshape omits
           the body entirely.
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

A_cross (winding cross-section area, for the Js[] current-density formula)
is the coil's r-y SECTION -- radial thickness x height -- taken via
Pappus's theorem from the solid's own volume. It is NOT the annulus
pi*((od/2)^2-(id/2)^2), which is the area for AXIAL current flow and was
used here in error until 2026-09-19; the current is azimuthal. See
detect_channel_cavity's sibling note and gmsh_getdp/README.md
"The A_cross bug".
"""
import argparse
import math
import os
import re
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
    emag_height=1.625, iem_coreod=1.0, iem_coilt=0.25, iem_platet=0.375,
    oem_coreod=1.0, oem_coilt=0.5,
)

IRON_CORE_NAMES = {"front_plate", "back_plate", "iem_core", "oem_core"}
CHAMBER_SPACER_NAME = "spacer"
CHAMBER_NAME = "chamber"
COIL_NAMES = {"coil"}

# Magnetically inert bodies that are present in the STEP but kept out of the
# FEM domain entirely, the same way `chamber` is: Ampere's law cannot tell
# mu_r=1 from empty space, so letting their volume fall into Air is exact,
# not an approximation.
#   anode_axial -- AISI 316L austenitic stainless (confirmed 2026-10-02).
#     mu_r ~ 1.00 annealed; cold work can raise it slightly via strain-
#     induced martensite, still negligible against the iron circuit.
#   Revolve1 -- 4 small stepped tubes (bore 8.1 mm, OD 10.0/14.2 mm), two
#     stacked at each of the two gas-feed holes through the back plate
#     (axes +/-26.86 mm from the thrust axis, measured 2026-10-03). Treated
#     as inert PROVISIONALLY: their material is not confirmed. They sit
#     outside the magnetic circuit, behind the back plate.
#   Revolve2, Fillet1, Imported1 -- bodies that appeared in the document
#     transiently while it was being edited mid-sweep (2026-10-03, ~15:30).
#     Treated as inert on the user's instruction; not present at the
#     microversion current at 15:52.
INERT_NAMES = {"anode_axial", "Revolve1", "Revolve2", "Fillet1", "Imported1"}

# Body names build() treated as inert because it did not recognise them.
LAST_UNKNOWN_BODIES: list[str] = []

# BPL-1000 renamed every body relative to BPL-700. Kept so a stale STEP or an
# old run directory is diagnosed by name instead of failing as "unrecognized".
BPL700_RENAMES = {
    "top_plate": "back_plate", "bottom_plate": "front_plate",
    "inner_emag_core": "iem_core", "outer_emag_core": "oem_core",
    "chamber_spacer": "spacer", "inner_coil": "coil", "outer_coil": "coil",
}

# Sized for the ~300 G channel-exit target (field_quality.py's target
# definition), replacing the inherited 5 A / 300t / 200t placeholder.
# 1.30 A measured against the AISI 1008 iron curve; peak iron |B| is
# 1.023 T, which is past 1008's peak-permeability point (0.865 T) and so
# the solve is noticeably stiffer than it was on the old generic curve.
# Opposite polarity so flux arcs across the gap between inner and outer
# poles rather than adding axially -- see HOW_IT_WORKS.md sec. 5.
#
# Turns are a real winding design, not the inherited 300/200 placeholder:
# AWG 20 heavy-build polyimide (NEMA MW 16-C), both coils wound FULL on two
# independent constant-current channels: 5 layers of 56 inner, 6 of 56
# outer, square-packed, over a 0.15mm ground-insulation allowance. From
# sizing/winding_design.py; wire dimensions from sizing/wire_data.py.
#
# The ampere-turn products are UNCHANGED from the previous entry here
# (260 x 1.30 = 338.0 and 165 x 1.30 = 214.5), so every solved field in this
# repo remains valid -- only the A-turn product reaches the solve. What
# changed is that these turns and currents describe a wire that can be
# bought, at a temperature class that is real, and the outer coils are now
# wound full instead of part-wound to 45% of their window, which is what
# made the uniform-current-density assumption below honest (see
# emag_sweep.py --fill-check).
EXCITATION = {
    "inner_coil": dict(turns=280, current=1.2071, polarity=-1),
    "outer_coil": dict(turns=336, current=0.6384, polarity=+1),
}


def detect_channel_cavity(centroids, nodes, axis_x: float, axis_z: float,
                           bin_m: float = 2.5e-4):
    """Locate the open plasma cavity inside the chamber's annular ceramic.

    Takes both halves of the same chamber volume mesh, because the two
    steps need different things:

    - `centroids` (volume element centroids) CLASSIFY which radial band is
      the cavity. This step needs points inside solid material, so surface
      nodes can't do it: a surface mesh puts nothing inside a wall, making
      the wall's interior look identical to open space -- with surface
      nodes this confidently returned r 31.98-40.94mm, the outer wall's
      guts, as the "cavity".
    - `nodes` (the same mesh's nodes) then MEASURE the bounds, because
      centroids sit half an element inside the true faces and bias every
      bound inward by ~0.2-0.3mm at a 1.5mm element size.

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

    r = np.hypot(centroids[:, 0] - axis_x, centroids[:, 2] - axis_z)
    y = centroids[:, 1]
    rn = np.hypot(nodes[:, 0] - axis_x, nodes[:, 2] - axis_z)
    yn = nodes[:, 1]
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
    r_in_approx, r_out_approx = float(r_edges[lo]), float(r_edges[hi])

    # Axial bounds first. The only chamber material at cavity radii is the
    # floor closing the annulus; whichever axial end it sits at is the
    # anode, the other is the open exit. Probe only the middle 60% of the
    # cavity width -- at the full width a fraction of a bin of wall overlap
    # reads as material at every y and destroys the result (hit
    # empirically: gave a 0.01mm long channel).
    pad = 0.2 * (r_out_approx - r_in_approx)
    in_cav = (rn > r_in_approx + pad) & (rn < r_out_approx - pad)
    assert in_cav.any(), "cavity band contains no chamber nodes -- cannot locate the closed end"
    y_floor = yn[in_cav]
    if abs(y_floor.mean() - yn.max()) < abs(y_floor.mean() - yn.min()):
        y_closed, y_open = float(y_floor.min()), float(yn.min())
    else:
        y_closed, y_open = float(y_floor.max()), float(yn.max())

    # Then the cavity radii, exactly. Nodes lie ON the ceramic faces, and in
    # the half of the channel clear of the floor the cavity is strictly
    # empty -- so the inner wall's outer face is simply the largest node
    # radius below the cavity, and the outer wall's inner face the smallest
    # above it. (Mode-of-histogram snapping was tried here first and lands
    # ~0.15mm off: face nodes aren't dense enough at this element size to
    # dominate a histogram bin.)
    r_mid = 0.5 * (r_in_approx + r_out_approx)
    span = y_closed - y_open
    open_half = ((yn - y_open) / span) < 0.5
    below, above = open_half & (rn < r_mid), open_half & (rn > r_mid)
    assert below.any() and above.any(), "cannot bracket the cavity radii on the open side"
    r_in, r_out = float(rn[below].max()), float(rn[above].min())
    return r_in, r_out, y_closed, y_open


def coil_radii(occ, d: int, t: int, height: float) -> tuple[float, float]:
    """(inner radius, outer radius) of a coil solid of revolution about Y,
    measured from the solid's own boundary faces.

    A full cylindrical lateral face of radius r and height h has area exactly
    2*pi*r*h, so r = area / (2*pi*h). occ.getMass() on a face is a true
    surface integral over the trimmed face, so this is exact -- unlike
    occ.getBoundingBox(), which is documented in build() below to return a
    LOOSE bound on exactly these curved STEP-imported solids.

    Deliberately NOT derived from the configuration parameters. On BPL-1000
    the coil OD is coreod + coilt and the ID tracks the core, but the as-built
    solids do not match that arithmetic to better than ~0.5 mm (measured
    2026-10-02: inner coil build 3.77 mm against a nominal 3.175 mm). Since
    r_mean sets A_cross and A_cross sets the ampere-turns the FEM actually
    applies -- the exact quantity the 2026-09-19 A_cross bug got wrong by
    ~1.9x -- this takes the geometry's word over the parameter arithmetic.

    Faces spanning the full height are the lateral ones; end-cap annuli have
    zero Y extent, and edge fillets span less than the full height.
    """
    lateral = []
    for sd, st in gmsh.model.getBoundary([(d, t)], oriented=False):
        sd, st = abs(sd), abs(st)
        _, sy0, _, _, sy1, _ = occ.getBoundingBox(sd, st)
        if abs((sy1 - sy0) - height) < 1e-9:
            lateral.append(occ.getMass(sd, st) / (2 * math.pi * height))
    assert len(lateral) == 2, (
        f"expected 2 full-height lateral faces on coil solid {t}, found "
        f"{len(lateral)} (radii {[round(r*1e3, 3) for r in lateral]} mm) -- the "
        f"solid is probably not the plain ring about Y this assumes")
    return min(lateral), max(lateral)


SMALL_HOLE_R_MAX = 0.004   # m. Fills 2mm lead holes, 5.1mm tapped holes, 7.1mm
                           # bolt clearances; keeps the 14.6mm gas-feed holes.


def _face_samples(tag: int, n: int = 16):
    """Points ON a (trimmed) face, from the B-rep parametrisation."""
    import numpy as np
    lo, hi = gmsh.model.getParametrizationBounds(2, tag)
    U, V = np.meshgrid(np.linspace(lo[0], hi[0], n), np.linspace(lo[1], hi[1], n))
    P = np.array(gmsh.model.getValue(2, tag, np.stack([U.ravel(), V.ravel()], 1).ravel())).reshape(-1, 3)
    keep = [gmsh.model.isInside(2, tag, list(p)) for p in P]
    return P[np.array(keep, bool)]


def fill_small_holes(occ, dimtag, r_max: float = SMALL_HOLE_R_MAX):
    """Remove holes of radius < r_max (axis along Y) from one iron body.

    BPL-1000's iron carries bolt clearances, tapped core ends and 2mm lead
    holes. Resolving them needs sub-mm elements, which took the mesh to
    ~1.3M DOFs and the direct LU solve past 31 GB of RAM (2026-10-03). They
    are magnetically negligible -- and the bolt holes are filled with steel
    bolts in the real hardware, so solid is closer to it than empty -- so
    they are filled with OCC's defeaturing (the user chose this over an
    iterative solver or a coarser mesh).

    A hole is found from its lateral face: an exact circle (least-squares
    residual < 1um) of radius < r_max. Every face of the body lying wholly
    within that radius of the same axis -- lateral wall, bottom, drill-point
    cone -- is part of the hole and goes to occ.defeature together.
    """
    import numpy as np
    occ.synchronize()
    faces = [abs(t) for _, t in gmsh.model.getBoundary([dimtag], oriented=False)]
    samples = {f: _face_samples(f) for f in faces}
    holes = []
    for f, P in samples.items():
        if len(P) < 10 or np.ptp(P[:, 1]) < 1e-6:
            continue
        x, z = P[:, 0], P[:, 2]
        A = np.c_[2 * x, 2 * z, np.ones_like(x)]
        (cx, cz, c), *_ = np.linalg.lstsq(A, x * x + z * z, rcond=None)
        r = math.sqrt(max(c + cx * cx + cz * cz, 0.0))
        if r < r_max and np.abs(np.hypot(x - cx, z - cz) - r).max() < 1e-6:
            holes.append((cx, cz, r))
    if not holes:
        return dimtag
    hole_faces = set()
    for cx, cz, r in holes:
        for f, P in samples.items():
            if len(P) and np.hypot(P[:, 0] - cx, P[:, 2] - cz).max() < r + 1e-5:
                hole_faces.add(f)
    out = occ.defeature([dimtag[1]], sorted(hole_faces))
    occ.synchronize()
    vols = [dt for dt in out if dt[0] == 3]
    assert len(vols) == 1, f"defeaturing body {dimtag} gave {vols}"
    print(f"  filled {len(holes)} hole(s) (r = "
          f"{sorted({round(h[2] * 1e3, 2) for h in holes})} mm, {len(hole_faces)} faces) in body {dimtag[1]}")
    return vols[0]


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
    inert_names_seen = []
    LAST_UNKNOWN_BODIES.clear()
    # emag_height = 1.625in is spacer height 0 on BPL-1000 (it was 1.125in on
    # BPL-700): Onshape then drops the body entirely rather than exporting a
    # zero-thickness solid. 1.625in is also the parameter's declared minimum,
    # so the spacer is absent at the default configuration.
    has_spacer = False
    for d, t in dimtags:
        name = base_name(gmsh.model.getEntityName(d, t))
        if name in IRON_CORE_NAMES:
            iron_dimtags.append((d, t))
        elif name == CHAMBER_SPACER_NAME:
            has_spacer = True
            if include_chamber_spacer_as_iron:
                iron_dimtags.append((d, t))
            # else: leave unclassified -> falls into Air, same as chamber
        elif name in COIL_NAMES:
            coil_entries.append((name, (d, t)))
        elif name == CHAMBER_NAME:
            assert chamber_dimtag is None, "expected exactly 1 chamber solid"
            chamber_dimtag = (d, t)
        elif name in INERT_NAMES:
            inert_names_seen.append(name)   # not meshed -- falls into Air
        elif name in BPL700_RENAMES:
            raise AssertionError(
                f"body {name!r} is a BPL-700 name; this pipeline targets BPL-1000 "
                f"(it would be {BPL700_RENAMES[name]!r} here). The STEP is from the "
                f"old document -- check onshape/cache.py's DID/EID.")
        else:
            # Unknown bodies are inert by the user's instruction (2026-10-03):
            # the document was being edited during sweeps and kept gaining
            # new non-magnetic bodies (Revolve2, Fillet1, Imported1,
            # Extrude1...). Recorded in LAST_UNKNOWN_BODIES so each run can
            # report what it ignored; the iron and coil counts are still
            # asserted below, so a renamed MAGNETIC body still fails loudly.
            print(f"WARNING: unrecognized body {name!r} (tag {t}) treated as inert")
            inert_names_seen.append(name)
            LAST_UNKNOWN_BODIES.append(name)

    n_expected_iron = 7 + (1 if include_chamber_spacer_as_iron and has_spacer else 0)
    assert len(iron_dimtags) == n_expected_iron, \
        f"expected {n_expected_iron} iron bodies, got {len(iron_dimtags)}"
    assert len(coil_entries) == 5, f"expected 5 coil bodies, got {len(coil_entries)}"
    assert chamber_dimtag is not None, "expected 1 chamber solid for channel-outline metadata"
    print(f"Classified: {len(iron_dimtags)} iron, {len(coil_entries)} coil, "
          f"1 chamber (metadata only), {len(inert_names_seen)} inert "
          f"({', '.join(sorted(inert_names_seen)) or 'none'})")

    # --- Windings: real coil solids straight from the STEP, no primitive ---
    # construction needed. Winding0 = inner_coil (center), Winding1..4 =
    # the 4 outer_coil occurrences (order among the 4 is arbitrary --
    # symmetric by design). Measured BEFORE the iron fuse below: OCC's
    # boolean ops can renumber/retag entities not even involved in the
    # operation, so bounding-box/center-of-mass queries on dimtags held
    # across a fuse call aren't reliable (hit empirically -- see git
    # history of this script for the stale-tag bbox mismatch this caused).
    # BPL-1000 names all five coil solids "coil" -- unlike BPL-700, which had
    # distinct inner_coil/outer_coil names -- so the pole has to come from
    # geometry. The 4 outer coils sit on a symmetric ring about the thrust
    # axis and the inner one sits on it, so the mean of all 5 centres IS the
    # axis to within the ring's symmetry, and the inner coil is the one at
    # (essentially) zero radius from it. Asserted, not assumed: the split has
    # to come out 1 + 4 with a clear radial gap, or the geometry is not what
    # this pipeline expects.
    #
    # The POLE LABELS stay "inner_coil"/"outer_coil": they are what EXCITATION
    # and generate_regions.py key on, and they describe the magnetic role,
    # which did not change with the CAD rename.
    coil_coms = [occ.getCenterOfMass(d, t) for _, (d, t) in coil_entries]
    mean_x = sum(c[0] for c in coil_coms) / len(coil_coms)
    mean_z = sum(c[2] for c in coil_coms) / len(coil_coms)
    radii = [math.hypot(c[0] - mean_x, c[2] - mean_z) for c in coil_coms]
    order = sorted(range(len(coil_entries)), key=lambda i: radii[i])
    assert radii[order[0]] < 1e-3, (
        f"no coil sits on the thrust axis (closest is {radii[order[0]]*1e3:.2f} mm off) "
        f"-- cannot identify the inner winding")
    assert radii[order[1]] > 10 * max(radii[order[0]], 1e-6), (
        f"inner/outer coil split is not clean: radii (mm) = "
        f"{[round(radii[i]*1e3, 2) for i in order]}")
    winding_tags = [coil_entries[i][1] for i in order]
    winding_pole = ["inner_coil"] + ["outer_coil"] * (len(order) - 1)
    print(f"Coil poles from geometry: 1 inner on-axis, {len(order)-1} outer at "
          f"r = {radii[order[1]]*1e3:.2f} mm")

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
        r_in_m, r_out_m = coil_radii(occ, d, t, height)
        r_mean = (r_in_m + r_out_m) / 2
        vol = occ.getMass(d, t)

        # A_cross is the cross-section the winding current actually crosses:
        # the coil's r-y section (radial thickness x height), NOT the annulus
        # pi*(ro^2-ri^2). Current here is AZIMUTHAL, so it passes through a
        # plane CONTAINING the axis; the annulus is the area for axial flow
        # and is ~1.8-1.9x too large on these coils, which silently applied
        # about half the intended ampere-turns. (The toy validation could not
        # catch this: build_toy.py fed the same wrong area into both the FEM
        # and its "analytic" reference, so the two agreed at 10.19 mT when
        # the true answer was 28.02 mT.)
        #
        # Taken via Pappus's theorem (V = 2*pi*r_centroid*A) rather than
        # dr*height: it needs only getMass(), a true volume integral, and it
        # stays correct for the small edge fillets the real CAD coils carry,
        # which a nominal dr*height misses by ~1%.
        A_cross = vol / (2 * math.pi * r_mean)
        A_nominal = (r_out_m - r_in_m) * height
        assert abs(A_cross - A_nominal) / A_nominal < 0.05, (
            f"winding{i} ({name}) section area from Pappus ({A_cross:.6e}) disagrees with "
            f"the nominal (od-id)/2 * height ({A_nominal:.6e}) by more than 5% -- the solid "
            f"is probably not the constant-section ring about Y that this assumes"
        )
        cx, _, cz = occ.getCenterOfMass(d, t)
        winding_meta.append(dict(index=i, pole=name, A_cross=A_cross, loc=(cx, cz),
                                 build=r_out_m - r_in_m))
    print(f"Built {len(winding_tags)} winding entries (from real CAD, no Gmsh primitives)")

    # --- Fuse iron bodies INCREMENTALLY -------------------------------------
    # One-shot occ.fuse(iron[0], iron[1:]) gives 3 disjoint solids on this
    # geometry (OCC's fuzzy boolean union is order/grouping-dependent) --
    # see scratch_whole_assembly_import.py's docstring. Incremental fuse
    # (fold one part in at a time) gives the correct single solid.
    iron_dimtags = [fill_small_holes(occ, dt) for dt in iron_dimtags]

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
    # Capped at the baseline inner build (0.1875in). Uncapped, a thinner core
    # (thicker inner winding) coarsened the whole mesh, which (a) left the
    # bottom plate's 2mm lead holes unresolvable -- "Invalid boundary mesh
    # (overlapping facets)" at inner_coil_id=0.875in -- and (b) made mesh
    # resolution vary across a geometry sweep, i.e. mesh noise in the
    # comparison. The cap keeps every sweep point at the baseline resolution.
    # Builds come from coil_radii()'s measurement of the real solids. (On
    # BPL-1000 they match coreod/coilt arithmetic exactly -- a 2026-10-02
    # note here claiming a ~0.5mm departure was a centroid-axis artefact,
    # see bpl1000_geometry.py -- but measuring stays the safer source.)
    # FIXED at the BPL-700 baseline build (0.1875in), not min()-ed with the
    # measured builds as it was: BPL-1000's default inner build is 0.125in,
    # which drove the near-surface size to 1.6mm, 884k nodes and 1.54M DOFs,
    # and the direct LU solve then needed ~58 GB on a 31 GB machine
    # (2026-10-03). A fixed size also keeps resolution identical across a
    # sweep, which is what the cap above was for; the thinnest swept coils
    # (0.0625in build) are then ~1 element thick, which a uniform-J source
    # region tolerates -- em_sizing_sweep.py's mesh-noise check bounds it.
    # The measured builds are still recorded in winding_meta.
    min_winding_wall = 0.1875 * INCH_TO_M
    # BPL-1000's iron carries small holes the uniform floor cannot resolve:
    # 2mm holes through the back plate, 5.1mm tapped holes in every core end
    # and 7.1mm bolt clearances. At a ~1.6mm floor their surface meshes come
    # out invalid and tetgen fails with "PLC Error: A segment and a facet
    # intersect". Curvature sizing refines exactly those holes (and nothing
    # flat), so the floor is lowered to let it work; the Threshold field
    # below still sets the resolution everywhere else.
    gmsh.option.setNumber("Mesh.MeshSizeMin", 2.5e-4)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
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

    try:
        gmsh.model.mesh.generate(3)
    except Exception as e:
        # Report WHICH surface failed and what it belongs to. "Invalid
        # boundary mesh (overlapping facets)" means the 2D mesh on that face
        # self-intersects, which on this geometry means the size field is too
        # coarse for a small feature -- but the message alone does not say
        # which feature, and guessing wastes a multi-minute rebuild per guess.
        print(f"\n3D mesh FAILED: {e}")
        for stag in {int(m) for m in re.findall(r"surface (\d+)", str(e))}:
            try:
                bb = gmsh.model.occ.getBoundingBox(2, stag)
                area = gmsh.model.occ.getMass(2, stag)
                dx, dy, dz = bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2]
                owners = [v for (_d, v) in gmsh.model.getEntities(3)
                          if (2, stag) in [(abs(a), abs(b)) for a, b in
                                           gmsh.model.getBoundary([(3, v)], oriented=False)]]
                print(f"  surface {stag}: area={area:.3e} m^2  extent="
                      f"({dx*1e3:.2f}, {dy*1e3:.2f}, {dz*1e3:.2f}) mm  "
                      f"bounds volumes {owners}")
                print(f"    bbox x[{bb[0]*1e3:.3f}, {bb[3]*1e3:.3f}] "
                      f"y[{bb[1]*1e3:.3f}, {bb[4]*1e3:.3f}] "
                      f"z[{bb[2]*1e3:.3f}, {bb[5]*1e3:.3f}] mm")
                curves = gmsh.model.getBoundary([(2, stag)], oriented=False)
                lens = sorted(gmsh.model.occ.getMass(1, abs(ct)) * 1e3 for _cd, ct in curves)
                print(f"    {len(curves)} bounding curve(s), lengths (mm): "
                      f"{[round(x, 3) for x in lens[:8]]}")
                print(f"    size field here = "
                      f"{gmsh.model.mesh.field.eval(2, (bb[0]+bb[3])/2, (bb[1]+bb[4])/2, (bb[2]+bb[5])/2)*1e3:.3f} mm"
                      if hasattr(gmsh.model.mesh.field, "eval") else "")
            except Exception as inner:
                print(f"  surface {stag}: could not characterise ({inner})")
        gmsh.finalize()  # leaving gmsh initialised makes the NEXT build write an empty mesh
        raise

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    msh_path = os.path.join(out_dir, "assembly.msh")
    gmsh.write(msh_path)

    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    assert n_nodes > 0, "empty mesh"
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
    chamber_nodes = np.array(node_coords).reshape(-1, 3)
    coord_of = dict(zip(node_tags, chamber_nodes))
    _, elem_tags, elem_nodes = gmsh.model.mesh.getElements(3)
    tets = np.array(elem_nodes[0]).reshape(-1, 4)
    chamber_pts = np.array([[coord_of[n] for n in tet] for tet in tets]).mean(axis=1)
    gmsh.finalize()

    # The plasma cavity, NOT the chamber part's bounding annulus -- see
    # detect_channel_cavity()'s docstring for why the difference matters.
    channel_inner_r, channel_outer_r, channel_y_anode, channel_y_exit = \
        detect_channel_cavity(chamber_pts, chamber_nodes, axis_x, axis_z)
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
    # One flag per real Onshape Configuration parameter, generated from the
    # defaults dict so the two cannot drift apart (they did across the
    # BPL-700 -> BPL-1000 retarget, when the parameter set changed wholesale).
    for pid, default in ONSHAPE_DEFAULTS_IN.items():
        ap.add_argument(f"--{pid.replace('_', '-')}", type=float, default=default)
    ap.add_argument("--out-dir", type=str, default=".")
    ap.add_argument("--step-path", type=str, default=None,
                     help="Bypass Onshape entirely and import this STEP file directly "
                          "(must still match the ID/OD args above for A_cross/mesh sizing).")
    ap.add_argument("--exclude-chamber-spacer", action="store_true",
                     help="Treat the spacer as magnetically inert (excluded from the iron "
                          "fuse) instead of the default iron treatment.")
    args = ap.parse_args()

    params_in = {pid: getattr(args, pid) for pid in ONSHAPE_DEFAULTS_IN}
    build(params_in, args.out_dir, args.step_path,
          include_chamber_spacer_as_iron=not args.exclude_chamber_spacer)


if __name__ == "__main__":
    main()
