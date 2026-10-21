"""
geometry.py

Parametric 2-D axisymmetric model of the BPL-1000 magnetic circuit, for fast
exploration of the channel B-field SHAPE (seconds per solve, vs ~6 min for the
3-D STEP pipeline in ../assembly/). Gmsh plane: x = r, y = thrust axis, same y
origin as the 3-D model (exit plane at y = -37.846 mm, downstream is -y).

Dimensions are the final spec (sizing/EMAG_SPEC.md), measured from the frozen
STEP on 2026-10-15 (r about the inner-core axis):
  inner core shaft    r 0..9.525,          y -28.321..12.954
  inner pole cap      r 0..17.4625, 9.525 thick, front face on the exit plane
  inner coil          r 9.525..14.2875,     y -28.321..12.954
  front plate         6 in square, bore r 34.925 (2.75 in), y -37.846..-28.321
  back plate          6 in square, centre hole r 3.40,      y 12.954..22.479
  outer cores (x4)    r 9.525 at R 62.862,  y -28.321..12.954
  outer coils (x4)    r 9.525..13.494 about each core
  chamber, anode      non-magnetic -> air

Two approximations make it axisymmetric:
- The square plates become discs of equal area (r = 152.4 / sqrt(pi) mm).
- The 4 outer posts become a ring at R with the same total iron area. Each
  post's coil becomes a pair of concentric winding sheets at the posts' radial
  extremes carrying -NI_out (inside) and +NI_out (outside). An Ampere loop
  through the ring and back outside it then encloses NI_out, exactly as a loop
  through one post encloses that post's coil. ring_mode picks how the ring is
  smeared:
    "thin"    full-mu ring of equal area (2.9 mm wide)
    "smeared" ring as wide as a post (19.05 mm), mu_eff = 1 + f (mu_r - 1),
              f = fill fraction -- right axial reluctance, but radially it is
              far more permeable than 4 separate posts
  calibrate.py compares both against the 3-D solve.

Extra iron (screens, pole extensions) is passed as polygons in mm via
params["extra_iron"]; the pole cap and front-plate bore can be replaced by
polygons via params["pole_cap_poly"] / params["front_plate_poly"]. Plate
thicknesses, rear air gaps (core or posts to back plate) and the axial extent
of each coil are parameters too (see DEFAULTS); layout() returns every shape
as a polygon, which the figures use to draw the geometry.
"""
import math
from pathlib import Path

import gmsh

MM = 1e-3
IN_MM = 25.4

Y_EXIT = -37.846
Y_FP_IN = -28.321          # front plate inner face == coil/core bottom
Y_COIL_TOP = 12.954        # back plate inner face
PLATE_T = 9.525

CHANNEL = dict(r_in=22.120, r_out=31.600, y_anode=-4.826, y_exit=Y_EXIT)

DEFAULTS = dict(
    mu_r=1500.0,
    plate_r=152.4 / math.sqrt(math.pi),
    fp_bore_r=34.925,
    bp_hole_r=3.40,
    core_r=9.525,             # inner core shaft
    pole_cap_r=17.4625,
    pole_cap_t=9.525,
    icoil_ro=14.2875,
    ring_R=62.8618,
    ocore_r=9.525,
    ocoil_build=13.494 - 9.525,
    ring_mode="smeared",
    plates=True,              # False: front/back plates non-magnetic (316L) -> air
    front_plate=True,         # per-plate switches (both also need plates=True)
    back_plate=True,
    trim=None,                # trim coil window (r0, r1, y0, y1) mm, own excitation NI_TRIM
    extra_iron=[],            # list of [(r, y), ...] polygons, mm
    pole_cap_poly=None,       # replaces the pole cap rectangle if given
    front_plate_poly=None,    # replaces the front plate rectangle if given
    fp_t=PLATE_T,             # front plate thickness (exit face fixed, inner face moves)
    bp_t=PLATE_T,             # back plate thickness (inner face fixed)
    core_gap=0.0,             # air gap between the inner core's rear end and the back plate
    ring_gap=0.0,             # air gap between the outer posts' rear ends and the back plate
    icoil_y=None,             # inner coil axial window (y0, y1) mm; default fills the core
    ocoil_y=None,             # outer coil axial window (y0, y1) mm; default fills the posts
    air_R=400.0,
    h_fine=0.35,              # mesh size in the channel/pole-gap box (mm)
    h_iron=1.2,
    h_far=25.0,
)

# Region ids shared with axisym.pro
R_AIR, R_IRON, R_WIN, R_WOI, R_WOO, R_RING, R_TRIM = 1, 2, 3, 4, 5, 6, 7
R_BND = 10


def ring_geometry(p):
    """(r0, r1, mu_r) of the smeared outer-core ring."""
    R, rc = p["ring_R"], p["ocore_r"]
    area = 4 * math.pi * rc ** 2
    if p["ring_mode"] == "thin":
        w = area / (2 * math.pi * R)
        return R - w / 2, R + w / 2, p["mu_r"]
    f = area / (2 * math.pi * R * 2 * rc)
    return R - rc, R + rc, 1 + f * (p["mu_r"] - 1)


def _poly(occ, pts):
    pts = [q for i, q in enumerate(pts) if q != pts[i - 1]]   # drop repeated vertices
    tags = [occ.addPoint(r * MM, y * MM, 0) for r, y in pts]
    lines = [occ.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
    return occ.addPlaneSurface([occ.addCurveLoop(lines)])


def layout(p):
    """Every material shape of the model as (region, [(r, y), ...] mm), in
    priority order: earlier wins where shapes overlap. Also used to draw the
    geometry on field maps."""
    p = dict(DEFAULTS, **p)
    y_fp = Y_EXIT + p["fp_t"]                       # front plate inner face
    y_bp = Y_COIL_TOP                                # back plate inner face
    ring_r0, ring_r1, _ = ring_geometry(p)
    oc0, oc1 = p["ring_R"] - p["ocore_r"], p["ring_R"] + p["ocore_r"]
    b = p["ocoil_build"]

    def rect(r0, r1, y0, y1):
        y0, y1 = min(y0, y1), max(y0, y1)
        return [(r0, y0), (r1, y0), (r1, y1), (r0, y1)]

    out = [(R_IRON, rect(0, p["core_r"], Y_EXIT, y_bp - p["core_gap"])),
           (R_IRON, p["pole_cap_poly"] or rect(0, p["pole_cap_r"], Y_EXIT, Y_EXIT + p["pole_cap_t"]))]
    if p["plates"] and p["front_plate"]:
        out.append((R_IRON, p["front_plate_poly"] or rect(p["fp_bore_r"], p["plate_r"], Y_EXIT, y_fp)))
    if p["plates"] and p["back_plate"]:
        out.append((R_IRON, rect(p["bp_hole_r"], p["plate_r"], y_bp, y_bp + p["bp_t"])))
    out += [(R_IRON, poly) for poly in p["extra_iron"]]
    out.append((R_RING, rect(ring_r0, ring_r1, y_fp, y_bp - p["ring_gap"])))
    iy0, iy1 = p["icoil_y"] or (Y_EXIT + p["pole_cap_t"], y_bp)
    out.append((R_WIN, rect(p["core_r"], p["icoil_ro"], iy0, iy1)))
    oy0, oy1 = p["ocoil_y"] or (y_fp, y_bp)
    out.append((R_WOI, rect(oc0 - b, oc0, oy0, oy1)))
    out.append((R_WOO, rect(oc1, oc1 + b, oy0, oy1)))
    if p["trim"]:
        r0, r1, y0, y1 = p["trim"]
        out.append((R_TRIM, rect(r0, r1, y0, y1)))
    return out


def build(p, out_dir):
    """Mesh the model into out_dir/axi.msh and write out_dir/regions.pro.
    Returns a dict of derived quantities (winding areas, ring mu)."""
    p = dict(DEFAULTS, **p)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("axi")
    occ = gmsh.model.occ

    _, _, ring_mu = ring_geometry(p)
    shapes = [(reg, _poly(occ, poly)) for reg, poly in layout(p)]

    disk = occ.addDisk(0, -10 * MM, 0, p["air_R"] * MM, p["air_R"] * MM)
    half = occ.addRectangle(0, (-10 - p["air_R"] - 1) * MM, 0, (p["air_R"] + 1) * MM,
                            (2 * p["air_R"] + 2) * MM)
    air, _ = occ.intersect([(2, disk)], [(2, half)])
    occ.synchronize()

    inputs = [(2, t) for _, t in shapes] + air
    _, omap = occ.fragment(inputs, [])
    occ.synchronize()

    owner = {}
    for i, (reg, _) in enumerate(shapes):
        for d, t in omap[i]:
            owner.setdefault(t, reg)
    for d, t in omap[len(shapes)]:
        owner.setdefault(t, R_AIR)

    by_reg = {}
    for t, reg in owner.items():
        by_reg.setdefault(reg, []).append(t)
    for reg, tags in by_reg.items():
        gmsh.model.addPhysicalGroup(2, tags, reg)

    # outer boundary: curves not on the axis (x = 0) that bound only the air
    bnd = []
    for d, t in gmsh.model.getBoundary([(2, t) for t in owner], combined=True, oriented=False):
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(1, t)
        if xmax > 1e-9:
            bnd.append(t)
    gmsh.model.addPhysicalGroup(1, bnd, R_BND)

    # mesh size: fine box over channel + pole gap + downstream, medium in the
    # magnet, coarse far field
    f = gmsh.model.mesh.field
    box = f.add("Box")
    f.setNumber(box, "VIn", p["h_fine"] * MM)
    f.setNumber(box, "VOut", p["h_far"] * MM)
    f.setNumber(box, "XMin", 0.0)
    f.setNumber(box, "XMax", 45 * MM)
    f.setNumber(box, "YMin", (Y_EXIT - 30) * MM)
    f.setNumber(box, "YMax", 0.0)
    f.setNumber(box, "Thickness", 60 * MM)
    fields = [box]
    gaps = []
    if p["core_gap"] > 0:
        gaps.append((0.0, p["core_r"], Y_COIL_TOP - p["core_gap"], p["core_gap"]))
    if p["ring_gap"] > 0:
        r0, r1, _ = ring_geometry(p)
        gaps.append((r0, r1, Y_COIL_TOP - p["ring_gap"], p["ring_gap"]))
    for r0, r1, y0, g in gaps:            # resolve thin rear air gaps
        bg = f.add("Box")
        f.setNumber(bg, "VIn", min(0.3, g / 2) * MM)
        f.setNumber(bg, "VOut", p["h_far"] * MM)
        f.setNumber(bg, "XMin", r0 * MM)
        f.setNumber(bg, "XMax", r1 * MM)
        f.setNumber(bg, "YMin", (y0 - 1) * MM)
        f.setNumber(bg, "YMax", (y0 + g + 1) * MM)
        f.setNumber(bg, "Thickness", 10 * MM)
        fields.append(bg)
    box2 = f.add("Box")
    f.setNumber(box2, "VIn", p["h_iron"] * MM)
    f.setNumber(box2, "VOut", p["h_far"] * MM)
    f.setNumber(box2, "XMin", 0.0)
    f.setNumber(box2, "XMax", (p["plate_r"] + 5) * MM)
    f.setNumber(box2, "YMin", (Y_EXIT - 40) * MM)
    f.setNumber(box2, "YMax", (Y_COIL_TOP + p["bp_t"] + 5) * MM)
    f.setNumber(box2, "Thickness", 80 * MM)
    mn = f.add("Min")
    f.setNumbers(mn, "FieldsList", fields + [box2])
    f.setAsBackgroundMesh(mn)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(2)
    gmsh.write(str(out / "axi.msh"))

    def area(reg):
        return sum(gmsh.model.occ.getMass(2, t) for t in by_reg.get(reg, []))

    derived = dict(A_win=area(R_WIN), A_woi=area(R_WOI), A_woo=area(R_WOO), A_trim=area(R_TRIM) or 1.0,
                   ring_mu=ring_mu, n_elem=int(sum(len(e) for e in gmsh.model.mesh.getElements(2)[1])))
    gmsh.finalize()

    (out / "regions.pro").write_text(
        "// generated by geometry.py\n"
        "Group {\n"
        f"  Air = Region[{R_AIR}]; Iron = Region[{R_IRON}]; Ring = Region[{R_RING}];\n"
        f"  WInner = Region[{R_WIN}]; WOutIn = Region[{R_WOI}]; WOutOut = Region[{R_WOO}];\n"
        f"  WTrim = Region[{R_TRIM}]; Bnd = Region[{R_BND}];\n"
        "}\n"
        f"A_WIN = {derived['A_win']!r};\nA_WOI = {derived['A_woi']!r};\nA_WOO = {derived['A_woo']!r};\n"
        f"A_TRM = {derived['A_trim']!r};\n"
        f"MUR_IRON = {p['mu_r']!r};\nMUR_RING = {ring_mu!r};\n")
    return derived
