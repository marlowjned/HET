# HET B-field simulation

Goal: an accurate 3-D magnetostatic simulation of the BPL-700 Hall Effect
Thruster's magnetic circuit, using the real Onshape CAD geometry and
per-part material properties (relative permeability), solved in MATLAB's
PDE Toolbox.

## Status

- **Parametric solenoid model**: working (`het_solenoid_bfield.m`).
- **Real CAD geometry import**: working — see `new_export/`.
- **Material assignment**: not yet wired into a solve. Materials are
  known (see below); the multi-domain classification pipeline that
  assigns them is prototyped in `new_export/scratch_classify.m`, but its
  first step (a clean `union` of iron + chamber + injector) hit a real
  geometry-kernel limitation and the plan is being revised — see
  "Next step".

## Files

### `het_solenoid_bfield.m`
Preliminary magnetostatic model using idealized parametric geometry (no
CAD): one center-pole solenoid + four outer-pole solenoids spaced around
the discharge channel opening, driven at opposite polarity. No thruster
body/core material (relative permeability = 1 everywhere), so `|B|` here
is a topology/sanity check, not a final result. Includes an interactive
cross-section viewer. Confirms the field-generation pipeline (geometry,
`cellLoad` current density, boundary conditions, solve, visualization)
works before the real CAD geometry and materials are added.

### `new_export/` — real CAD geometry, correctly positioned
The current, working import of the actual BPL-700 assembly from Onshape:

- `BPL-700 Assembly - *.step` — one flat STEP file per named part
  (`chamber`, `injector`, `top_plate`, `bottom_plate`, `center_solenoid`,
  `outer_solenoid` x4). Each imports cleanly into `fegeometry` with no
  cleanup needed. These carry correct *shape* but not correct assembled
  *position* — Onshape's per-part STEP export flattens patterned
  instances to one local placement (see Gotchas below).
- `apply_transforms.m` — applies the real per-occurrence rotation +
  translation (extracted once from a whole-assembly STEP export,
  `BPL-700-assem-SIMS-V1.step` in the repo root) to each part above,
  producing correctly assembled geometry. Run this first; it writes
  `placed_parts.mat`. The transform values are hardcoded in this script;
  the one-time extraction script that produced them
  (`scratch_check_transforms.py`) has served its purpose and lives in
  `junk/` — only needed again if the assembly is re-exported with
  different placements.
- `placed_parts.mat` — output of `apply_transforms.m`: a cell array of
  9 correctly-positioned `fegeometry` objects (`placed`) with matching
  identity labels (`labels`, e.g. `outer_solenoid_1`).
- `check_placed.m` — loads `placed_parts.mat`, prints the overall
  assembly bounding box, and reports every pair of parts that touch or
  overlap (via `addCell` probing). Useful after any re-export to confirm
  the geometry still looks like a sane assembly.
- `scratch_center_gap.m` — one-off profiling of the chamber's radial bore
  vs. axial position in the center-pole's local frame. Confirmed the
  center pole is stepped: it's only "flush" against the chamber over the
  last ~11mm of its length (pole radius steps from 15.9mm to 30.8mm right
  there); across the other ~100mm the chamber's own bore stays at a
  constant ~30.6mm radius, leaving a genuine ~15mm annular gap around the
  thin pole section. That gap is the center winding's natural home — the
  chamber bore itself defines the winding OD, no invented wall thickness
  needed.
- `scratch_windings.m` — builds a hollow cylindrical winding-sleeve solid
  for each of the 5 poles (center + 4 outer), sized to fit the real
  clearance gap without touching iron/chamber/injector (the center one
  sized specifically to the thin-section gap found by
  `scratch_center_gap.m`), and verifies no winding-to-structure or
  winding-to-winding collisions. Geometry stage confirmed complete: all 5
  windings float freely with zero collisions. Writes `windings.mat`.
- `windings.mat` — output of `scratch_windings.m`: the 5 winding-sleeve
  `fegeometry` solids plus their labels/rotation/translation.
- `scratch_isolate_mesh.m`, `scratch_keepbound.m` — one-off checks that
  `union`-ing the iron/chamber/injector solids together stays meshable
  (incrementally, and with/without `KeepBoundaries`) before committing to
  the union-based multi-domain approach below.
- `scratch_classify.m` — first working prototype of the full multi-domain
  plan described in "Next step" below: unions iron + chamber + injector
  into one blob, adds the 5 windings as separate cells, meshes the whole
  9-cell geometry once, then classifies every element of the iron/chamber/
  injector blob into its correct region via `findCell` (with a
  nearest-solid fallback for elements it can't place directly). Prototype
  only — confirms the approach works on this geometry, but its output
  (`classified.mat`) isn't wired into a `femodel` solve yet and isn't
  checked into the repo.

### `het_view_assembly.m` (repo root)
Loads `placed_parts.mat` and renders all 9 parts together, color-coded by
part identity, with centroid labels, so the assembly can be rotated and
visually verified (correct part count, no stray duplicates, plausible
contacts, sane proportions) before trusting it for the solve. Confirmed
against the real geometry: two pole plates connected by 4 outer poles and
1 center pole, chamber nested in the middle.

### `BPL-700-assem-SIMS-V1.step` (repo root)
The whole-assembly STEP export used as the source of truth for real
per-instance placement transforms (see `apply_transforms.m`). Kept in
case the transform extraction ever needs to be redone or double-checked.

### `junk/`
Superseded exports, one-time extraction/diagnostic scripts, and abandoned
approaches from getting to the current pipeline (see history below).
Nothing load-bearing for the current workflow lives here.

## Materials (known, not yet wired into the solve)

| Part                          | Material         |
|--------------------------------|------------------|
| `injector`                     | Stainless steel  |
| `chamber`                      | Ceramic          |
| everything else (both plates, center + outer solenoids) | Soft iron |

For a magnetostatic solve, the only material property that matters is
`RelativePermeability`. Soft iron needs a nonlinear B-H curve for
accuracy near saturation (MATLAB supports this via a function handle,
see `pde/NonlinearRelativePermeabilityExample`) — the specific soft iron
alloy is still needed to get a real B-H curve. **Gotcha**: MATLAB's
built-in material catalog (`materialProperties(Material="steel")`)
silently returns `RelativePermeability = 1` — it has no real magnetic
data — so it must not be used for the iron parts.

## Next step: multi-domain material assignment

`addCell` (used by `check_placed.m` to probe contacts) only supports
nesting one solid strictly inside another — it cannot build a geometry
where two solids share a real conformal boundary, which is what a
correctly-mated assembly needs. The original approach, prototyped
end-to-end in `new_export/scratch_classify.m` (see above):

1. `union` the iron + chamber + injector parts into one solid "structure
   blob" (union tolerates touching/overlapping solids; validated
   incrementally in `scratch_isolate_mesh.m`/`scratch_keepbound.m`).
2. Build the 5 winding-sleeve solids separately (`scratch_windings.m`,
   `windings.mat`) — these float in the air gap and stay their own cells,
   not unioned into the blob.
3. `addCell` the structure blob and all 5 windings into the air domain
   sphere (7 cells total) and mesh once, conformally.
4. For each element of the structure-blob cell, classify which original
   part (iron / chamber / injector) it falls inside via `findCell`
   (vectorized point-in-solid test), with a nearest-solid distance
   fallback for elements `findCell` can't place directly.
5. Rebuild via `fegeometry(mesh, ElementIDToRegionID)` — MATLAB's
   sanctioned mechanism for a multi-material geometry with real shared
   (conformal) boundaries (`pde/DMultidomainFegeometryObjectExample`) —
   and wire the result into a `femodel` magnetostatic solve with real
   per-region `RelativePermeability` and winding `cellLoad` currents.

### Blocker hit at step 1: `union` isn't clean across every interface

Isolating the failure (per-part-count meshability checks in
`scratch_isolate_mesh.m`) found that `union`-ing the 8 iron parts
together works perfectly: 73 faces, meshes fine, ~103k elements. But the
moment `chamber` (16 faces on its own) gets unioned into that iron blob,
the result jumps to 234 faces and fails to mesh across every
`Hmax`/`Hmin` combination tried — not a resolution problem, a geometry-
quality one. The likely cause: chamber's mating surface is near-but-not-
exactly coincident with the iron's, so the boolean kernel produces a pile
of degenerate sliver faces instead of one clean shared boundary, rather
than failing outright. (Also worth noting: `union` for `fegeometry` is a
brand-new MATLAB feature, introduced in R2025a.)

This matters because the whole mesh-once-then-`findCell`-classify plan
above depended on that union being clean. It isn't, at least for the
iron/chamber interface.

**Proposed fix (decision pending as of this writing):** since chamber and
injector are magnetically inert (relative permeability ≈ 1, same as air)
and their touching interface with the iron only matters geometrically —
not electromagnetically — for this solve, give chamber and injector a
small clearance gap from the iron (and from each other) via an
inflate-and-subtract, rather than a zero-gap union. The iron itself stays
one true zero-gap union (needed for flux continuity across it). This also
simplifies the pipeline: every region becomes its own natively separate
`addCell`'d cell (like the windings already are), making the
mesh-once-then-`findCell`-reclassify step (steps 3-5 above) unnecessary
entirely. The alternative not yet tried: re-tessellating/simplifying
chamber's mating face before the union, to see if a clean conformal union
is achievable after all.

## Gotchas encountered (for future re-exports)

- Onshape's raw "Export as STEP" wraps every part in an assembly
  structure that MATLAB's STEP reader rejects
  (`STEPGeometryMustBeSolid`). Exporting each part individually as its
  own flat STEP file (as in `new_export/`) avoids this entirely.
- A STEP file can only contain **one** solid shape for `fegeometry` /
  `importGeometry` — a multi-body assembly STEP cannot be imported
  directly, flattened or not.
- Per-part STEP export does **not** preserve patterned-instance assembly
  placement — all 4 `outer_solenoid` copies exported to the identical
  local position. The fix was to separately export the whole assembly to
  STEP once, extract each occurrence's real placement transform from the
  `NEXT_ASSEMBLY_USAGE_OCCURRENCE` / `ITEM_DEFINED_TRANSFORMATION` chain,
  and apply it in MATLAB (`apply_transforms.m`) to the reliable per-part
  geometry.
- `pdegplot` silently resets `hold` to `off` after every call — re-assert
  `hold on` before each subsequent plot/legend call in a multi-part
  figure or later parts silently wipe out earlier ones.
- `union` for `fegeometry` (new in R2025a) can silently produce
  degenerate sliver faces instead of a clean shared boundary when two
  solids' mating surfaces are near-but-not-exactly coincident — it
  doesn't error, but the result then fails to mesh at any `Hmax`/`Hmin`.
  Symptom to watch for: a disproportionate face-count jump after adding
  one part to a union (see "Blocker hit at step 1" below).

## Future features / things to consider

- Full thruster material included (in progress, see above)
- Force acting on solenoids
- E field sim
- Trace electron trajectory, quantify stability
- Coil thermals?

Requires MATLAB R2023b+ with the Partial Differential Equation Toolbox.
