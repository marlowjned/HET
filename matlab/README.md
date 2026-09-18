# HET B-field simulation

Goal: an accurate 3-D magnetostatic simulation of the BPL-700 Hall Effect
Thruster's magnetic circuit, using the real Onshape CAD geometry and
per-part material properties (relative permeability), solved in MATLAB's
PDE Toolbox.

## Status

- **Parametric solenoid model**: working (`het_solenoid_bfield.m`).
- **Real CAD geometry import**: working — see `new_export/`.
- **Material assignment**: not yet wired into a MATLAB solve. Materials
  are known (see below); the multi-domain classification pipeline that
  assigns them is prototyped in `new_export/scratch_classify.m`, but its
  first step (a clean `union` of iron + chamber + injector) hit a real
  geometry-kernel limitation. Re-diagnosed in detail since this was first
  written: it's not one bad interface, it's 5 of the assembly's 9 real
  touching interfaces (see "Blocker hit at step 1" under "Next step"),
  and the plan is being revised.
- **A second, independent pipeline outside MATLAB now exists and is
  further along**: `../gmsh_getdp/` (Gmsh for geometry/meshing, GetDP for
  the FEM solve — both compiled, not Python) fuses all 7 real iron parts
  cleanly on the first try (vs. MATLAB's 4/9 clean interfaces) and has a
  working, validated (1.4% vs. an exact analytic check) magnetostatic
  solve on the real assembly, including reasonable placeholder materials
  and winding currents. See `../gmsh_getdp/README.md` for status and
  caveats (still linear iron, no B-H saturation curve — same gap as here).

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
- `scratch_pairwise_diagnose.m` — isolates every part pair that actually
  touches (rediscovers the contact graph itself via the same `addCell`
  probing `check_placed.m` uses) and unions just that pair alone, checking
  whether the resulting face count is exactly additive (clean shared
  boundary) or inflated (sliver faces), and whether the pair meshes at
  all. Localizes exactly which mating interfaces are broken instead of
  only knowing "the accumulated blob fails somewhere." See "Blocker hit at
  step 1" below for results.
- `scratch_check_tessellation.m` — checks each part's own STEP-derived
  boundary triangulation, independent of any union, for degenerate
  near-zero-area triangles. Distinct failure mode from the mating-pair
  check above: a part can be individually clean and still produce a bad
  union with a neighbor (a placement/tolerance problem), or be
  individually dirty regardless of who it's unioned with (an export-
  tessellation problem). See "Blocker hit at step 1" below for results.

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

**Re-diagnosed since the paragraph above was first written — the
original "iron works, only chamber is broken" story does not hold up.**
Re-running `scratch_isolate_mesh.m` on this machine (MATLAB R2025b — the
only release installed here; `union` for `fegeometry` was introduced in
R2025a, so it's possible the original "73 faces, meshes fine" result was
produced on R2025a and the boolean kernel's tolerance behavior changed —
unconfirmed, but consistent with everything below), the iron-only chain
now fails outright at the 6th of 8 parts with `Unable to combine the
geometries`, and even the 5-part blob just before that already fails to
mesh at every `Hmax`/`Hmin` tried.

Isolating every pair of parts that actually touches (the same `addCell`
contact graph `check_placed.m` computes; automated in
`scratch_pairwise_diagnose.m`) shows the problem isn't "iron blob vs.
chamber" — it's specific individual mating interfaces, and it's worse
than one blocker:

| interface | result |
|---|---|
| `outer_solenoid_1` ↔ `top_plate` | clean (exactly additive faces), meshes |
| `outer_solenoid_1` ↔ `bottom_plate` | clean, meshes |
| `outer_solenoid_9` ↔ `top_plate` | clean, meshes |
| `outer_solenoid_9` ↔ `bottom_plate` | clean, meshes |
| `outer_solenoid_7` ↔ `top_plate` | mildly inflated (+3–10 spurious faces depending on union order), unmeshable |
| `outer_solenoid_2` ↔ `top_plate` | badly inflated (13–17 spurious faces), unmeshable at any resolution |
| `center_solenoid` ↔ `chamber` | union "succeeds" (23→49 faces) but unmeshable at any resolution |
| `injector` ↔ `chamber` | **hard failure**: `Unable to create a valid geometry. The triangulation contains degenerate triangles.` |
| `chamber` ↔ `top_plate` | same hard failure as injector↔chamber |

Of the 9 real touching interfaces in the whole assembly, 4 are clean and
5 are broken — including all 3 of chamber's contacts (only 1 of which,
iron, was previously suspected) and 2 of the 4 outer-pole-to-top-plate
mates. Two more findings narrow this further:

- **The broken pole mates correlate with placement, not part shape.**
  All 4 `outer_solenoid` copies are the same STEP file placed by
  `apply_transforms.m`'s hardcoded per-occurrence transform table.
  `outer_solenoid_1`/`_9` (clean) share `top_plate`'s exact translation
  coordinate on one axis, consistent with having been mated directly to
  `top_plate` in Onshape. `outer_solenoid_2`/`_7` (broken) instead share
  that coordinate with `chamber`/`center_solenoid`, consistent with their
  placement being derived through a *different*, longer mate chain that
  only reaches `top_plate` indirectly — any tolerance accumulated along
  that chain would surface exactly as this kind of near-but-not-quite
  coincident mating face. This is a placement/mate-chain problem, not
  something inherent to the outer-pole shape.
- **`injector`'s own STEP tessellation is unusually dirty independent of
  any union.** `scratch_check_tessellation.m` checks each part's own
  boundary triangulation for degenerate (near-zero-area) triangles:
  `injector` has 53 triangles under 1e-8 area on a 39-face part, far more
  than any other part (the 4 outer-solenoid copies and center_solenoid
  all have a small, *identical* handful — 20 or fewer, plausibly benign
  cylinder-cap artifacts common to revolved STEP tessellation — and
  `chamber`/`bottom_plate`/`top_plate` have zero). `chamber` and
  `top_plate` are each individually clean, so their mutual union failure
  is purely a placement/tolerance issue introduced by the boolean itself
  — but `injector`'s failure against `chamber` may be compounded by
  `injector`'s own dirty export, not placement alone.
- **`union(A,B)` is not symmetric.** The exact face-count inflation for a
  broken pair varies depending on which operand is unioned into which
  (e.g. `outer_solenoid_2`+`top_plate` produced 35–39 faces across
  otherwise-identical reruns depending on argument order) — the pair
  stays broken either way, but the specific symptom isn't deterministic
  across orderings. Worth keeping in mind for any future fix attempt or
  checker script: test both orderings, don't trust one.

This matters because the whole mesh-once-then-`findCell`-classify plan
above depended on every one of these unions being clean. Most aren't.

**Proposed fix for the inert parts (chamber, injector) — decision pending
as of this writing:** since chamber and injector are magnetically inert
(relative permeability ≈ 1, same as air) and their touching interfaces
only matter geometrically — not electromagnetically — for this solve,
give them a small clearance gap from whatever they touch (and from each
other) via an inflate-and-subtract, rather than a zero-gap union. This
also simplifies the pipeline: every such region becomes its own natively
separate `addCell`'d cell (like the windings already are), making the
mesh-once-then-`findCell`-reclassify step (steps 3-5 above) unnecessary
for those regions entirely.

**This fix does not cover the `outer_solenoid_2`/`_7` ↔ `top_plate`
mates.** Those are iron-to-iron interfaces that need a true zero-gap
union for flux continuity (per the reasoning above), so they can't be
gapped away the same way. Since the likely cause is a placement/mate-
chain issue in the source assembly, the real fix is on the Onshape side:
correct the pattern/mates so occurrences 2 and 7 sit exactly flush
against `top_plate` the way 1 and 9 already do, then re-extract
transforms. Not yet done — needs a decision on whether to fix this in
Onshape now or work around it in MATLAB first (e.g. re-tessellating the
mating faces before union, unconfirmed whether that would even work).

## CAD software choice: Onshape vs. Fusion 360 vs. SolidWorks

Came up while scoping the Onshape-side fix above (a parametric gap and/or
corrected mates, re-exported programmatically): does the CAD package
matter for where this pipeline is headed? Yes, mainly on how automatable
headless parameter-change-and-export is — not just "does it have an API"
generically.

- **Onshape** (current choice) — best fit, and already in use. A real
  cloud REST API: API-key auth, read/write FeatureScript parameters,
  trigger STEP export, all callable from any language/OS (including this
  Mac) with no local install or license seat consumed just to automate.
  Exactly the shape this pipeline already assumes.
- **Fusion 360** — has a genuine headless path via Autodesk Platform
  Services (Design Automation + Data Management APIs), architecturally
  similar to Onshape's, but noticeably more complex to stand up (OAuth
  app registration, "AppBundles"/"Activities") and less mature/documented
  for this use case. Its simpler API is in-app Python scripting, which
  needs a running Fusion session — at least that works natively on Mac,
  since Fusion (unlike SolidWorks) isn't Windows-only.
- **SolidWorks** — weakest fit for this need specifically. Automation is
  COM-based and desktop-bound: requires a licensed SolidWorks instance
  actually running, on Windows only. No lightweight remote-callable REST
  layer like Onshape's; batch export exists (Task Scheduler, Excel-linked
  design tables) but "change a parameter and pull a fresh STEP from a
  script on this Mac" isn't how it's built to work.

If "best for the future of the simulation architecture" means what this
pipeline already does today — source-parametric gaps, scripted
re-export, eventually sweeping a dimension programmatically — Onshape is
the clear architectural fit and SolidWorks the biggest step backward.
The real tradeoff is team familiarity/existing licenses, which is a
legitimate cost but a separate axis from API capability.

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
- **A standalone CAD-checker script**, likely necessary at some point:
  automate what `scratch_pairwise_diagnose.m` and
  `scratch_check_tessellation.m` currently do by hand into one script
  that runs after every re-export/re-transform and reports, per touching
  interface, whether the union is clean and meshable, plus per-part
  tessellation quality — rather than re-discovering broken mates one at a
  time. Not urgent today, but as the assembly gets re-exported from
  Onshape (mate fixes, added parts, etc.) this turns a manual multi-probe
  diagnosis into a single pass/fail gate.
- Force acting on solenoids
- E field sim
- Trace electron trajectory, quantify stability
- Coil thermals?

Requires MATLAB R2023b+ with the Partial Differential Equation Toolbox.
