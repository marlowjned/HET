"""
build_assembly.py

Imports the real BPL-700 iron parts (per-part STEP exports in
../../matlab/new_export/), applies the same per-occurrence placement
transforms extracted from the whole-assembly STEP export (see
../../matlab/new_export/apply_transforms.m for the original MATLAB version
and derivation), fuses them into one solid, builds the 5 winding sleeves,
and fragments everything (iron + windings + outer air sphere) into one
conformal mesh.

chamber and injector are deliberately NOT included as solid regions: both
are magnetically inert (mu_r ~= 1, same as air/vacuum -- ceramic and
non-magnetic stainless respectively), so for a magnetostatic solve their
exact shape is physically irrelevant. Iron already carries its own correct
real shape (with its true bore/cavity where chamber/injector sit)
regardless of whether those parts are separately meshed -- omitting them
just leaves that space classified as Air, which is physically correct.
This also sidesteps two real geometry problems hit when they WERE included
(see matlab/README.md and this script's git history):
  - chamber and iron genuinely overlap in a thin 3D sliver at their mating
    face (not just touch), which needs a real inflate/clearance fix to
    mesh cleanly, not just a tolerant boolean;
  - injector's own STEP tessellation is independently dirty (matlab's
    scratch_check_tessellation.m found 53 degenerate triangles on it,
    the worst of any part) and produced self-intersecting facets when
    meshed no matter how it was placed.
If per-region material fidelity for chamber/injector is ever needed (e.g.
for a thermal or structural model, not magnetostatics), reintroduce them
with an explicit inflate-and-subtract clearance gap sized to the mesh's
min element size -- see the SHRINK logic removed from this script's
history for a starting point.

The win this script DOES demonstrate: fusing all 7 real iron parts (top
plate, bottom plate, center pole, 4 outer poles) into one solid succeeds
on the first try via Gmsh/OpenCASCADE's fuzzy-tolerance boolean kernel --
this is the exact interface that broke MATLAB's fegeometry.union() on 5
of 9 real touching interfaces (matlab/README.md, "Blocker hit at step 1").
"""
import gmsh
import numpy as np
import math
import os

STEP_DIR = os.path.join("..", "..", "matlab", "new_export")

# Same transform table as matlab/new_export/apply_transforms.m: name,
# translation, local Z-axis direction (global), local X-refdir (global).
# Extracted from the NAUO->IDT chain in BPL-700-assem-SIMS-V1.step. Only the
# 7 iron occurrences are placed below (see module docstring for why chamber/
# injector are excluded); chamber's own transform is kept separately, below,
# since it's still needed to derive the channel-outline metadata used by
# cross_section_viewer.py (chamber is meshed there but never added to the
# FEM domain -- see the "Channel-outline metadata" block near the end).
TRANSFORMS = [
    ("outer_solenoid",  (-0.125033506959677, -0.0217433018267154,  0.0770308774188161), (0, 1, 0), (1, 0, 0)),
    ("outer_solenoid",  ( 0.00196649304032328, 0.0798566981732845,  0.0770308774188161), (0, -1, 0), (1, 0, 0)),
    ("center_solenoid", (-0.0615335069596768,  0.0798566981732845,  0.0135308774188161), (0, -1, 0), (1, 0, 0)),
    ("bottom_plate",    (-0.0615335069596767,  0.0893816981732846,  0.0135308774188161), (0, 0, 1),  (1, 0, 0)),
    ("outer_solenoid",  ( 0.00196649304032332, 0.0798566981732845, -0.0499691225811839), (0, -1, 0), (1, 0, 0)),
    ("top_plate",       (-0.0615335069596767, -0.0217433018267154,  0.0135308774188161), (0, 0, 1),  (1, 0, 0)),
    ("outer_solenoid",  (-0.125033506959677,  -0.0217433018267154, -0.0499691225811839), (0, 1, 0),  (1, 0, 0)),
]

CHAMBER_TRANSFORM = ("chamber", (-0.0615926647537434, 0.0798566981732845, 0.0135553813794386), (0, -1, 0), (1, 0, 0))

STEP_FILE = {
    "top_plate":       "BPL-700 Assembly - top_plate-1__Body1.step",
    "bottom_plate":    "BPL-700 Assembly - bottom_plate-1__Body1.step",
    "center_solenoid": "BPL-700 Assembly - center_solenoid-1__Body1.step",
    "outer_solenoid":  "BPL-700 Assembly - outer_solenoid-1__Body1.step",
}

# Winding sleeve dims, same as matlab/new_export/scratch_windings.m -- local
# frame id/od/z0/z1 per pole TYPE (center pole is stepped; see
# scratch_center_gap.m), not a generic formula.
clearance = 0.002
outer_dims = dict(
    id=2 * (0.01588 + clearance),
    od=2 * (0.01588 + clearance) + 2 * 0.008,
    z0=0.00762,
    z1=0.00762 + 0.10160 * 0.85,
)
center_dims = dict(
    id=2 * (0.01588 + clearance),
    od=2 * (0.0306 - clearance),
    z0=0.005,
    z1=0.096,
)


def axis_angle_from_R(R):
    """Port of apply_transforms.m's rotmat2axisangle -- returns (axis, angle_rad)."""
    c = max(-1.0, min(1.0, (np.trace(R) - 1) / 2))
    angle = math.acos(c)
    if angle < 1e-9:
        return np.array([1.0, 0.0, 0.0]), 0.0
    if abs(math.pi - angle) < 1e-6:
        M = (R + np.eye(3)) / 2
        idx = int(np.argmax(np.diag(M)))
        v = np.sqrt(np.maximum(M[:, idx], 0))
        others = [k for k in range(3) if k != idx]
        for k in others:
            if M[idx, k] < 0:
                v[k] = -v[k]
        return v / np.linalg.norm(v), math.pi
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / (2 * math.sin(angle))
    return v / np.linalg.norm(v), angle


def rotation_matrix(zdir, xref):
    Zax = np.array(zdir, dtype=float); Zax /= np.linalg.norm(Zax)
    Xax = np.array(xref, dtype=float); Xax /= np.linalg.norm(Xax)
    Yax = np.cross(Zax, Xax)
    return np.column_stack([Xax, Yax, Zax])  # columns = global dirs of local X,Y,Z


def place(dimtags, loc, zdir, xref):
    """Rotate dimtags about origin per (zdir,xref), then translate by loc -- same order as apply_transforms.m."""
    R = rotation_matrix(zdir, xref)
    axis, angle = axis_angle_from_R(R)
    if angle > 1e-6:
        gmsh.model.occ.rotate(dimtags, 0, 0, 0, axis[0], axis[1], axis[2], angle)
    gmsh.model.occ.translate(dimtags, loc[0], loc[1], loc[2])


gmsh.initialize()
gmsh.model.add("bpl700_assembly")
occ = gmsh.model.occ

# Import each unique iron part once as a template, then copy per occurrence.
# NOTE: these STEP files' SI_UNIT header declares METRE, but the raw
# coordinate values are physically millimeters (e.g. top_plate imports
# with a ~190-unit bounding box -- a 190 mm plate is a sane size for this
# thruster, a 190 m one is not; this also matches the ~0.1 m-scale
# translations in TRANSFORMS above, which were extracted assuming mm->m).
# MATLAB's fegeometry STEP importer evidently assumes mm unconditionally
# and auto-converts; Gmsh's OCC-based reader respects the (here, wrong)
# declared unit and doesn't rescale -- so rescale explicitly here.
STEP_SCALE = 0.001
templates = {}
for name, fname in STEP_FILE.items():
    path = os.path.join(STEP_DIR, fname)
    dimtags = occ.importShapes(path)
    assert len(dimtags) == 1 and dimtags[0][0] == 3, f"{name}: expected 1 solid, got {dimtags}"
    occ.dilate(dimtags, 0, 0, 0, STEP_SCALE, STEP_SCALE, STEP_SCALE)
    templates[name] = dimtags[0]
occ.synchronize()

iron_tags = []
for name, loc, zdir, xref in TRANSFORMS:
    src = templates[name]
    copy = occ.copy([src])
    place(copy, loc, zdir, xref)
    iron_tags.append(copy[0])

# Remove the unplaced templates (still sitting at local-frame origin)
occ.remove(list(templates.values()), recursive=True)
occ.synchronize()

print(f"Placed {len(iron_tags)} iron parts")

# --- Fuse the 7 iron parts into one solid ----------------------------------
# This is the step that broke MATLAB's fegeometry.union() on 5 of 9 real
# touching interfaces (sliver faces, unmeshable). Set a fuzzy boolean
# tolerance to give OCC's kernel slack on near-but-not-exactly coincident
# mating faces (the outer_solenoid_2/_7 <-> top_plate placement/mate-chain
# issue documented in matlab/README.md).
gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-6)

iron_result, iron_map = occ.fuse([iron_tags[0]], iron_tags[1:])
occ.synchronize()
print(f"Iron fuse: {len(iron_tags)} parts -> {len(iron_result)} solid(s): {iron_result}")
assert len(iron_result) == 1, f"iron fuse did not produce a single solid: {iron_map}"
iron_tag = iron_result[0]

# --- Build the 5 winding sleeves -------------------------------------------
# Same geometry as matlab/new_export/scratch_windings.m: a hollow cylinder
# built in the pole's LOCAL frame (id/od/z0/z1 along local Z), then rotated
# + translated by the same per-occurrence transform as its host pole.
pole_occurrences = [t for t in TRANSFORMS if t[0] in ("outer_solenoid", "center_solenoid")]
winding_tags = []
for name, loc, zdir, xref in pole_occurrences:
    d = center_dims if name == "center_solenoid" else outer_dims
    h = d["z1"] - d["z0"]
    outer_cyl = occ.addCylinder(0, 0, d["z0"], 0, 0, h, d["od"] / 2)
    inner_cyl = occ.addCylinder(0, 0, d["z0"], 0, 0, h, d["id"] / 2)
    cut = occ.cut([(3, outer_cyl)], [(3, inner_cyl)])
    sleeve = cut[0][0]
    place([sleeve], loc, zdir, xref)
    winding_tags.append(sleeve)
occ.synchronize()
print(f"Built {len(winding_tags)} winding sleeves: {winding_tags}")

# --- Outer air domain -------------------------------------------------------
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

# --- Fragment everything together: iron blob, 5 windings, and the air
# sphere -- into one conformal complex. All solids are non-overlapping
# (touching at most), so each input should map to exactly one output
# fragment; the air sphere is the only one expected to split (into "the
# leftover air after removing every embedded solid").
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
print(f"\nTotal volume before fragment: {vol_before:.6e}  after: {vol_after:.6e}  "
      f"(should match -- fragment only re-partitions, doesn't add/remove material)")
assert abs(vol_before - vol_after) / vol_before < 1e-9

# --- Physical groups ---------------------------------------------------
# Tags pinned explicitly (matches magnetostatics_assembly.pro). GetDP's
# Region[] only understands integer tags, not physical-name strings.
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

# Outer boundary surface = boundary of the air volume not shared with any
# embedded solid
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

# --- Mesh --------------------------------------------------------------
# Global size from the air domain; refine near the iron/winding surfaces
# so thin winding walls and pole gaps are resolved.
min_winding_wall = min(
    (outer_dims["od"] - outer_dims["id"]) / 2,
    (center_dims["od"] - center_dims["id"]) / 2,
)
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
gmsh.write("assembly.msh")

n_nodes = len(gmsh.model.mesh.getNodes()[0])
print(f"\nFinal mesh: {n_nodes} nodes")

gmsh.finalize()

# --- Channel-outline metadata for cross_section_viewer.py ------------------
# chamber is magnetically inert and not part of the FEM domain (see module
# docstring), but its real position is still useful as a plotted reference
# -- the discharge channel is physically where chamber's ceramic wall is.
# Import+place it exactly like the iron parts above, profile its radius vs.
# axial position about the center pole's own axis (which is the assembly's
# natural central axis -- every pole occurrence's local axis is global Y,
# per TRANSFORMS), and keep only the resulting numbers, not the geometry.
# This mirrors matlab/new_export/scratch_center_gap.m's profiling logic,
# redone here since chamber isn't part of this pipeline's live geometry.
# Run in a completely fresh gmsh session (not reusing the one above) --
# calling mesh.generate() again after a background mesh-size field was set
# and a full 3D mesh already written crashed the native library outright.
gmsh.initialize()
gmsh.model.add("chamber_profile")
occ2 = gmsh.model.occ

chamber_name, chamber_loc, chamber_zdir, chamber_xref = CHAMBER_TRANSFORM
center_entry = next(t for t in TRANSFORMS if t[0] == "center_solenoid")
axis_loc = center_entry[1]  # (x, y, z) -- axis_x, axis_z below; axis runs along global Y

chamber_dt = occ2.importShapes(os.path.join(STEP_DIR, "BPL-700 Assembly - chamber-1__Body1.step"))
occ2.dilate(chamber_dt, 0, 0, 0, STEP_SCALE, STEP_SCALE, STEP_SCALE)
place(chamber_dt, chamber_loc, chamber_zdir, chamber_xref)
occ2.synchronize()
gmsh.model.mesh.generate(2)  # surface only -- fast, sufficient for a radius profile
_, chamber_node_coords, _ = gmsh.model.mesh.getNodes()
chamber_pts = np.array(chamber_node_coords).reshape(-1, 3)
gmsh.finalize()

axis_R = rotation_matrix(center_entry[2], center_entry[3])
local = (axis_R.T @ (chamber_pts - np.array(axis_loc)).T).T
r_local = np.sqrt(local[:, 0] ** 2 + local[:, 1] ** 2)
y_local_from_axis = -local[:, 2]  # local Z is global -Y here (see rotation_matrix); flip for a plain global-Y span

channel_inner_r = float(r_local.min())
channel_outer_r = float(r_local.max())
channel_y_min = float((axis_loc[1] + y_local_from_axis).min())
channel_y_max = float((axis_loc[1] + y_local_from_axis).max())
print(f"\nChannel outline (from chamber geometry): inner_r={channel_inner_r:.4f} "
      f"outer_r={channel_outer_r:.4f} y=[{channel_y_min:.4f},{channel_y_max:.4f}] "
      f"about axis (x={axis_loc[0]:.4f}, z={axis_loc[2]:.4f})")

with open("assembly_params.txt", "w") as f:
    for name, t in sorted(TAG.items(), key=lambda kv: kv[1]):
        f.write(f"{t} {name}\n")
    f.write(f"n_windings {len(winding_tags)}\n")
    for i, (name, loc, zdir, xref) in enumerate(pole_occurrences):
        d = center_dims if name == "center_solenoid" else outer_dims
        A_cross = math.pi * ((d["od"] / 2) ** 2 - (d["id"] / 2) ** 2)
        f.write(f"winding{i} pole={name} A_cross={A_cross!r} loc={loc} zdir={zdir} xref={xref}\n")
    f.write(f"axis_x {axis_loc[0]!r}\n")
    f.write(f"axis_z {axis_loc[2]!r}\n")
    f.write(f"channel_inner_r {channel_inner_r!r}\n")
    f.write(f"channel_outer_r {channel_outer_r!r}\n")
    f.write(f"channel_y_min {channel_y_min!r}\n")
    f.write(f"channel_y_max {channel_y_max!r}\n")

# --- Generate the GetDP region/current-source include file -----------------
# Reasonable placeholder excitation, same turns/current as the parametric
# check in matlab/het_solenoid_bfield.m: center 300 turns, outer 200 turns,
# both at 5 A. Center and outer driven at opposite magnetic polarity (field
# lines arc across the annular gap instead of just adding axially) -- same
# physical choice as het_solenoid_bfield.m's centerCoil.sign/outerCoil.sign.
#
# Every pole occurrence here happens to have its local winding axis (zdir)
# along the GLOBAL Y axis (either +Y or -Y -- see TRANSFORMS), so the
# current density is built directly as a scaled cross product zdir x (P -
# loc): this is the azimuthal direction around that axis (GetDP has no
# built-in cross-product operator, so it's hand-expanded into components
# below), with magnitude Jmag = turns*current/A_cross. Its SIGN is chosen
# per winding so that the resulting B direction is consistent across all 4
# outer poles in GLOBAL coordinates (not just relative to each winding's
# own, possibly mirrored, local frame) and opposite for the center pole:
# sign = polarity_of(pole_type) * zdir_y (since zdir_y in {+1,-1} flips
# which geometric rotation direction corresponds to "B along global +Y").
mu0 = 4 * math.pi * 1e-7
EXCITATION = {
    "center_solenoid": dict(turns=300, current=5.0, polarity=-1),
    "outer_solenoid":  dict(turns=200, current=5.0, polarity=+1),
}

lines = []
lines.append("// AUTO-GENERATED by build_assembly.py -- do not hand-edit.")
lines.append("// Region tags and winding current-density sources for the real")
lines.append("// BPL-700 assembly. See build_assembly.py for derivation.")
lines.append("")
lines.append("Group {")
for name, t in sorted(TAG.items(), key=lambda kv: kv[1]):
    lines.append(f"  {name} = Region[{t}];")
winding_names = [f"Winding{i}" for i in range(len(winding_tags))]
lines.append(f"  Windings = Region[{{{', '.join(winding_names)}}}];")
lines.append("  DomainC  = Region[{Windings}];")
lines.append("  DomainCC = Region[{Air, Iron}];")
lines.append("  Domain   = Region[{DomainC, DomainCC}];")
lines.append("}")
lines.append("")
lines.append("Function {")
lines.append(f"  mu0 = {mu0!r};")
lines.append("  mur_iron = 1000;  // reasonable placeholder (linear) -- see README caveat: soft iron")
lines.append("                    // needs a real nonlinear B-H curve for accuracy near saturation")
lines.append("  nu[Air]  = 1 / mu0;")
lines.append("  nu[Iron] = 1 / (mur_iron * mu0);")
lines.append("  nu[Windings] = 1 / mu0;  // copper, non-magnetic")
lines.append("")
for i, (name, loc, zdir, xref) in enumerate(pole_occurrences):
    d = center_dims if name == "center_solenoid" else outer_dims
    A_cross = math.pi * ((d["od"] / 2) ** 2 - (d["id"] / 2) ** 2)
    exc = EXCITATION[name]
    zy = zdir[1]
    assert zdir[0] == 0 and zdir[2] == 0 and abs(zy) == 1, \
        f"winding {i} ({name}) local axis isn't along global Y ({zdir}) -- Js[] formula below assumes it is"
    Jmag = exc["polarity"] * zy * exc["turns"] * exc["current"] / A_cross
    lines.append(f"  // Winding{i}: pole={name} turns={exc['turns']} I={exc['current']}A "
                  f"A_cross={A_cross:.6e} zdir_y={zy:+d} polarity={exc['polarity']:+d} -> Jmag={Jmag:.6e} A/m^2")
    lines.append(f"  Js[Winding{i}] = ({Jmag!r} / Sqrt[(X[]-({loc[0]!r}))^2 + (Z[]-({loc[2]!r}))^2]) * "
                  f"Vector[(Z[]-({loc[2]!r})), 0, -(X[]-({loc[0]!r}))];")
lines.append("}")
lines.append("")

with open("assembly_regions_generated.pro", "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"\nWrote assembly_regions_generated.pro ({len(winding_tags)} windings)")
