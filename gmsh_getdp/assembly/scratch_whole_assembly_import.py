"""
One-off validation: can Gmsh/OCC import the WHOLE-ASSEMBLY STEP directly
(matlab/BPL-700-assem-SIMS-V1.step) and skip the per-part-STEP-export +
apply_transforms.m hardcoded-placement-table workflow that build_assembly.py
currently uses?

Answer: yes, with one caveat (see "INCREMENTAL FUSE" below). Kept here as a
record of the finding, not part of the build pipeline -- build_assembly.py
still uses the per-part + hardcoded-transform approach as of this writing.
This becomes the basis for the parametric (per-sweep-point) CAD pipeline:
each sweep point can export ONE whole-assembly STEP from Onshape and skip
maintaining any transform table, since Onshape's own mates already encode
correct placement and Gmsh recovers it exactly.

Checks performed:
1. Import the whole assembly, list the resulting solids and confirm their
   real-world placement matches apply_transforms.m's hand-extracted table
   (matched by part name; X/Z compared directly, Y compared after
   accounting for the known local-origin-vs-centroid offset -- see
   inline comments). All 9 solids' X/Z match to numerical precision.
2. Fuse the 7 iron parts (everything except chamber/injector) into one
   solid and mesh it, confirming the pipeline's core geometry step still
   works when sourced this way.

FINDING -- INCREMENTAL FUSE REQUIRED: fusing all 7 iron parts in one batch
call, `occ.fuse([iron_tags[0]], iron_tags[1:])` (as build_assembly.py does
today), produced 3 disjoint solids instead of 1 on this file -- NOT because
placement is wrong, but because OCC's fuzzy boolean union is order/grouping
-dependent (matlab/README.md already documents this same non-symmetry for
MATLAB's kernel; Gmsh/OCC has a milder version of it). Fusing
INCREMENTALLY -- folding one part at a time into a growing result --
produced a clean single solid (73 boundary faces, meshes with no errors).
Any future parametric pipeline sourcing iron from a whole-assembly STEP
should fuse incrementally, not in one batch call.

GOTCHA -- POST-IMPORT dilate() UNIT FIX IS UNRELIABLE ON ASSEMBLY SOLIDS:
this file's raw coordinates are in mm despite the SI_UNIT(METRE) header
(same bug as the per-part exports, see HOW_IT_WORKS.md sec 3), so a
0.001 unit fix is still needed. Applying it via occ.dilate() on solids
pulled out of an ASSEMBLY-structured STEP (as opposed to a lone single-
solid STEP) does not behave as a pure linear scale here -- comparing raw
vs. post-dilate centroids shows each part shifting by a small (~8 mm),
part-specific amount along one axis, not the 0 shift a true scale-about-
origin should give. Likely cause: OCC often stores an assembly solid's
placement as a separate "Location" transform layered on the shape's raw
geometry, and a blanket dilate() may not compose with that the same way
it does for an already-flattened single-solid import. This did NOT
corrupt the fuse result above (73 faces, 1 solid, clean mesh -- the
RELATIVE placement between the 7 iron parts came out correct), so this
script's checks still pass, but it means dilate()-based unit correction
should not be trusted in general on assembly-sourced solids. The clean
fix for the real Phase 2 pipeline: don't reintroduce this mm/m mismatch
at all -- configure the new parametrized Onshape export to emit correct
real-world units directly, rather than relying on a post-import scale
correction.
"""
import gmsh

STEP_FILE = r"../../matlab/BPL-700-assem-SIMS-V1.step"
STEP_SCALE = 0.001  # same mm-labeled-as-m bug as the per-part exports (see HOW_IT_WORKS.md sec 3)

PART_NAMES = ("outer_solenoid", "center_solenoid", "injector", "chamber",
              "bottom_plate", "top_plate")

# From matlab/new_export/apply_transforms.m -- translation applied to each
# part's LOCAL ORIGIN after rotation (not its centroid).
EXPECTED = [
    ("outer_solenoid",  [-0.125033506959677, -0.0217433018267154,  0.0770308774188161]),
    ("outer_solenoid",  [ 0.00196649304032328, 0.0798566981732845,  0.0770308774188161]),
    ("center_solenoid", [-0.0615335069596768,  0.0798566981732845,  0.0135308774188161]),
    ("injector",        [-0.0615926647537434,  0.0766816981732846,  0.0135553813794385]),
    ("chamber",         [-0.0615926647537434,  0.0798566981732845,  0.0135553813794386]),
    ("bottom_plate",    [-0.0615335069596767,  0.0893816981732846,  0.0135308774188161]),
    ("outer_solenoid",  [ 0.00196649304032332, 0.0798566981732845, -0.0499691225811839]),
    ("top_plate",       [-0.0615335069596767, -0.0217433018267154,  0.0135308774188161]),
    ("outer_solenoid",  [-0.125033506959677,  -0.0217433018267154, -0.0499691225811839]),
]

gmsh.initialize()
gmsh.model.add("scratch_whole_assembly_import")

gmsh.model.occ.importShapes(STEP_FILE)
gmsh.model.occ.synchronize()
volumes = gmsh.model.getEntities(dim=3)
gmsh.model.occ.dilate(volumes, 0, 0, 0, STEP_SCALE, STEP_SCALE, STEP_SCALE)
gmsh.model.occ.synchronize()

print(f"Imported {len(volumes)} volume(s) from the whole-assembly STEP\n")

parts = []  # (tag, short_name, centroid, extent)
for dim, tag in volumes:
    label = gmsh.model.getEntityName(dim, tag)
    short = next(n for n in PART_NAMES if n in label)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    centroid = ((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
    extent = (xmax - xmin, ymax - ymin, zmax - zmin)
    parts.append((tag, short, centroid, extent))

print("--- Placement check vs. apply_transforms.m (matched by part name) ---")
used = set()
for tag, short, centroid, extent in parts:
    best_i, best_d = None, float("inf")
    for i, (exp_name, exp_loc) in enumerate(EXPECTED):
        if i in used or exp_name != short:
            continue
        d = ((centroid[0]-exp_loc[0])**2 + (centroid[2]-exp_loc[2])**2) ** 0.5
        if d < best_d:
            best_d, best_i = d, i
    used.add(best_i)
    _, exp_loc = EXPECTED[best_i]
    dy = centroid[1] - exp_loc[1]
    half_extent_y = extent[1] / 2
    y_ok = abs(abs(dy) - half_extent_y) < 0.002
    print(f"  tag={tag:2d} {short:16s} XZ_residual={best_d:.6f} m  "
          f"dY={dy:+.4f} m vs half-extent-Y={half_extent_y:.4f} m  "
          f"{'OK (local-origin offset)' if y_ok else 'unexplained (informational only -- excluded from iron solve if injector)'}")

print("\n--- Iron fuse check (7 parts: both plates, center + 4x outer pole) ---")
iron_tags = [tag for tag, short, _, _ in parts if short not in ("chamber", "injector")]
assert len(iron_tags) == 7

gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-6)

current = [(3, iron_tags[0])]
for t in iron_tags[1:]:
    current, _ = gmsh.model.occ.fuse(current, [(3, t)])  # INCREMENTAL -- see module docstring
    gmsh.model.occ.synchronize()

print(f"Fuse result: {len(current)} solid(s) (want 1)")
assert len(current) == 1, "iron fuse did not produce a single solid"

faces = gmsh.model.getBoundary(current, oriented=False)
print(f"Fused solid has {len(faces)} boundary faces")

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.02)
gmsh.model.mesh.generate(3)
nodes = gmsh.model.mesh.getNodes()
elems = gmsh.model.mesh.getElements(dim=3)
n_tets = sum(len(e) for e in elems[1])
print(f"Meshed OK: {len(nodes[0])} nodes, {n_tets} volume elements")

gmsh.finalize()
