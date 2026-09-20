# HET B-field simulation: Gmsh + GetDP pipeline

A second, independent magnetostatics pipeline for the same BPL-700
geometry as `../matlab/`, built after MATLAB's `fegeometry.union()` proved
unable to cleanly union 5 of the assembly's 9 real touching interfaces
(sliver faces, unmeshable at any resolution -- see `../matlab/README.md`,
"Blocker hit at step 1"). Both mesh generation (Gmsh) and the FEM
assembly/solve (GetDP) are compiled native code -- Python here is only
orchestration (build geometry, write the mesh, shell out to `getdp.exe`,
read results back), not the numerical engine.

## Status

- **Toy validation**: done and passing. `toy_test/` builds a single
  parametric solenoid coil (same dimensions as
  `../matlab/het_solenoid_bfield.m`'s center coil) and solves it with
  GetDP's magnetic-vector-potential formulation. Result: 10.05 mT on-axis
  at coil center, vs. 10.19 mT from the exact analytic on-axis integral
  for a finite thick solenoid (1.4% agreement) -- confirms the mesh,
  formulation, and gauge are all correct before trusting them on the real
  CAD.
- **Real assembly**: done and solving, now sourced from a **whole-assembly
  Onshape STEP export** (via `../onshape/cache.py`) rather than the
  original per-part-STEP + hardcoded-transform approach. `assembly/
  build_assembly.py` imports all 14 named bodies directly (7 iron + 5 real
  coil solids + chamber + chamber_spacer), fuses the iron bodies
  **incrementally** into one solid, keeps the 5 coil solids as real
  windings (no more `occ.addCylinder`/`occ.cut` primitive construction),
  meshes, and solves. **All real touching interfaces mesh cleanly** (vs.
  MATLAB's 4/9) -- see "The geometry win" below. All 5 real Onshape
  Configuration parameters (`emag_height`, `inner/outer_coil_id/od`) now
  drive this pull -- see `../onshape/cache.py`'s docstring for why (the
  soft-iron core's diameter IS the coil ID it wraps, so there's no
  separate core-diameter parameter). `chamber_spacer` is treated as iron
  for now (`build_assembly.py --exclude-chamber-spacer` to test the inert
  alternative) -- its real material isn't confirmed.
- **Center/outer current-ratio sweep**: done at the old geometry and the
  old linear-iron material (pole names `center_solenoid`/`outer_solenoid`,
  since renamed to `inner_coil`/`outer_coil` -- `current_sweep.py` itself
  is updated for the rename, but hasn't been re-run against the new
  whole-assembly geometry, and its channel statistics predate the ROI
  correction below). Found the most radial (E x B - correct) channel
  field at `I_center = I_outer` -- see "Current ratio sweep" below. **Its
  2-basis-solve superposition trick assumes a linear material**, which
  looked fatal once Iron became nonlinear -- but at the real ~300 G design
  point the iron peaks at 0.86 T, well below the B-H knee, so it behaves
  linearly and the trick is valid again there (confirmed: the nonlinear
  peak matches a linear-scaled prediction to 0.4%). Re-check that if a
  future design point pushes the iron toward saturation.
- **Materials**: nonlinear iron (generic soft-steel B-H curve) replaces
  the `mur_iron = 1000` linear placeholder, and **converges** -- 11
  Picard iterations, 19.6 min wall, 805 MB peak. See "Nonlinear iron
  (B-H curve)" below for how the excitation, not the method, was what
  made that work.
- **Excitation**: 2 A per coil (inner 300 t, outer 200 t), sized to the
  ~300 G channel-exit target rather than inherited from the old MATLAB
  check. Measured 292 G at the exit plane. See "Operating point" below.
- **Channel ROI**: the channel bounds in `assembly_params.txt` are the
  real plasma cavity between the chamber's ceramic walls, not the chamber
  part's bounding annulus (which contained iron). See "Channel ROI" below.
- **chamber**: deliberately excluded from the solid model (see
  `assembly/build_assembly.py` docstring) -- magnetically inert (mu_r ~=
  1, same as air), so omitting it doesn't affect the B-field solve.
  **injector** isn't present in the whole-assembly STEP export at all
  (confirmed) -- harmless, it was already excluded for the same reason.

## Setup

```
pip install gmsh meshio pyvista        # geometry/meshing API, VTK export, VTK-based slicing (all C++ core underneath)
# GetDP: download the Windows binary from https://getdp.info/ (bin/Windows/getdp-3.5.0-Windows64c.zip)
# and unzip it -- this repo expects it at ../tools/getdp-3.5.0-Windows64/getdp.exe
# (../tools/ is gitignored; not checked in)
```

## Running

```
cd toy_test
python build_toy.py
"../../tools/getdp-3.5.0-Windows64/getdp.exe" magnetostatics_toy.pro -msh toy.msh -solve Res_a -pos Map_b
# check b_center.txt against the printed B_ideal

cd ../assembly
python build_assembly.py                       # writes assembly.msh + assembly_regions_generated.pro + assembly_params.txt
"../../tools/getdp-3.5.0-Windows64/getdp.exe" magnetostatics_assembly.pro -msh assembly.msh -solve Res_a -pos Map_b
python post_process.py                          # writes solve_summary.txt + cross_section.png (quick static check)
python export_vtk.py                             # writes assembly_field.vtu + iron_field.vtu (standard VTK) -- run before the viewer below
python cross_section_viewer.py                  # interactive angle-sweep viewer (opens a window), reads assembly_field.vtu
```

## The geometry win

MATLAB's `union()` broke on 5 of 9 real touching interfaces: all 3 of
chamber's contacts (hard failure, degenerate triangles) and 2 of the 4
outer-pole-to-top-plate mates (sliver faces, unmeshable). Gmsh/OCC's
fuzzy-tolerance boolean kernel (`Geometry.ToleranceBoolean`) fused all 7
real iron parts (top plate, bottom plate, center pole, 4 outer poles) into
a single solid **on the first try**, including the exact
outer_solenoid_2/_7 <-> top_plate mates that broke MATLAB. This is the
main validation that switching toolchains actually fixes the root
problem, not just relocates it.

## Result file formats

- **`assembly.msh`** -- the mesh (nodes, elements, physical-group region
  tags). Gmsh legacy ASCII v2.2. No field data.
- **`b_assembly.pos`, `b_iron.pos`** -- GetDP's native output (Gmsh's
  `.pos` ASCII post-processing view format): one `VS(...){...}` line per
  tetrahedron, its 4 corner coordinates and constant field value. Plain
  text, large (~200MB for the full assembly).
- **`assembly_field.vtu`, `iron_field.vtu`** (`export_vtk.py`) -- the same
  field data as standard VTK (XML/binary, `UnstructuredGrid`), for
  ParaView or any other VTK-standard consumer -- the interoperable format
  if this ever gets layered against other simulation tools. Built via
  `meshio` directly from the `.pos` element data (each tet keeps its own
  4 corner points rather than being deduplicated against a shared mesh --
  see `export_vtk.py`'s docstring for why: **`gmsh.view.write(tag,
  "*.vtk")` looked like it worked (no error) but silently wrote a
  non-standard format** -- no `# vtk DataFile` header, none of the
  expected `POINTS`/`CELLS` sections -- so don't reach for it as a
  shortcut here). Verified by a full round-trip read: peak `|B|` in the
  `.vtu` matches the known `.pos` value (9.671343399912887 T) exactly.
- `post_process.py` still reads `.pos` directly via `pos_utils.py` (its
  own quick static check, unchanged). `cross_section_viewer.py` now reads
  `assembly_field.vtu` via `pyvista` instead -- see "Known limitations"
  below for why that's a real accuracy improvement, not just a format
  swap, and note it means `export_vtk.py` must be run before the viewer.

## Current ratio sweep

`assembly/current_sweep.py` sweeps the center-solenoid / outer-solenoid
current ratio and finds which one makes the discharge-channel field most
radial (perpendicular to the thruster centerline) -- the standard Hall
thruster requirement for E x B electron confinement.

**Key physics shortcut**: with the linear `mur_iron = 1000` placeholder
(no B-H saturation), the magnetostatics problem is linear in the winding
currents, so solutions superpose. The script solves GetDP only twice --
once with 1 A in the center winding alone, once with 1 A in each outer
winding alone -- caches both per-element channel fields to
`channel_basis.npz`, then evaluates `B(x) = I_center * B_center_unit(x) +
I_outer * B_outer_unit(x)` for any current pair instantly in Python (no
remeshing or re-solving). Use `--skip-solve` to re-sweep or re-plot from
the cached basis without paying for the two GetDP solves again.

**Result** (`sweep_results.csv`, `sweep_ratio.png`; 25 ratios from 0.4 to
2.5, `I_outer` fixed at the reference 5 A/coil): the channel field is most
radial -- lowest volume-weighted axial fraction, `mean(|B_y|/|B|) =
0.4848` -- at **`I_center / I_outer = 1.0`**, i.e. equal current in the
center and each outer coil (`I_center = I_outer = 5 A`). This happens to
coincide with the placeholder excitation already used for the base solve
(`center 300t/5A, outer 200t/5A`) -- not by design, since the turns
differ 300 vs. 200 per coil, but the ratio metric is fairly flat near its
minimum (axial fraction only rises from 0.485 to ~0.54 across the whole
swept range), so this isn't a sharp optimum. At that point:

- channel `|B|` mean = 80.6 mT, median = 8.5 mT, p90 = 33.3 mT (the mean
  is pulled up by a small high-field volume near the pole faces; the
  median is more representative of the bulk channel field).
- to hit a different target channel field while keeping the field
  direction optimal, scale **both** currents together by
  `target / 80.6 mT`, keeping the 1:1 ratio fixed.

**Those three magnitude numbers are stale** -- they were computed with
the old channel mask (the chamber part's bounding annulus, which
included iron; see "Channel ROI" above), so they average pole-face field
into the channel and the spread between mean/median/p90 is largely that
contamination rather than real structure. Over the corrected plasma
cavity the field is far more uniform and almost purely radial (97.7%).
The *ratio* result is unaffected -- that metric is scale-invariant and
the contamination was symmetric across ratios -- but don't use the mT
figures above to size current; use "Operating point" instead.

This is a field-direction optimum only -- it says nothing about the
field *magnitude* needed for a real operating point (that depends on
electron Larmor radius / discharge voltage requirements, not modeled
here), and the underlying solve is the linear/no-saturation placeholder
-- though see "Nonlinear iron" above for why superposition turns out to
be legitimate at the real design point anyway.

## Solenoid winding specs (current configuration)

Geometry from `assembly/build_assembly.py` (winding sleeve dims) and
`assembly/generate_regions.py`/`assembly_params.txt` (turns/current
excitation). The solve itself only ever uses a bulk current density
(`turns * current / A_cross`, ../gmsh_getdp/HOW_IT_WORKS.md sec. 9) --
it has no notion of individual wire strands -- so the per-wire numbers
below (gauge, resistance, voltage, power) are **not part of the FEM
model**; they're a separate hand-calc, assuming round AWG 12 magnet wire
(2.053 mm bare copper diameter) and copper resistivity at two reference
temperatures, added here to size the actual coil/power supply.

| | outer coil (x4) | center coil |
|---|---|---|
| ID / OD | 35.76 / 51.76 mm | 35.76 / 57.20 mm |
| winding height | 86.4 mm | 91.0 mm |
| window cross-section | 1099.8 mm^2 | 1565.3 mm^2 |
| turns | 200 | 300 |
| mean turn length | 137.5 mm | 146.0 mm |
| total wire length | 27.5 m | 43.8 m |
| assumed wire | AWG 12 (2.053 mm dia, 3.31 mm^2) | AWG 12 (same) |
| implied packing factor | 60.2% | 63.4% |
| resistance @ 20 C | 139.5 mOhm | 222.3 mOhm |
| resistance @ 100 C (est. operating) | 183.4 mOhm | 292.2 mOhm |

AWG 12 is a guess, not a spec -- but it's notable that it lands at a
plausible 60-63% packing factor (round wire in a wound coil tops out
around 78-91% theoretical, 60-75% is normal once insulation and
hand-winding slack are counted) in **both** windings despite their
different turn counts and window areas, which is at least consistent
with one wire gauge being used for the whole coil set rather than
suggesting a modeling error.

At the design operating point (`I_center = I_outer = 2 A`, see "Operating
point" below -- the 1:1 ratio is also the most radial-field ratio found
in the sweep above):

| | outer coil (each) | center coil | total (1 center + 4 outer) |
|---|---|---|---|
| voltage @ 20 C | 0.28 V | 0.44 V | -- |
| voltage @ 100 C | 0.37 V | 0.58 V | -- |
| power @ 20 C | 0.56 W | 0.89 W | 3.1 W |
| power @ 100 C | 0.73 W | 1.17 W | 4.1 W |

If all 5 coils are wired in series on one supply at this operating point
(same 2 A through all of them, true only at this 1:1 ratio), the whole
magnet circuit needs **~1.6-2.1 V at 2 A, ~3-4 W total** depending on
winding temperature -- a trivially small bench-supply load, nothing like
the discharge (anode) supply. Moving off the 1:1 ratio breaks the series
assumption, since `I_center != I_outer` then requires either two
independent supplies or a shunt/trim resistor on one leg; scaling within
that constraint, power grows with the square of current on whichever
coil's current changes.

## Operating point

Target: **~300 G radially at the channel exit**, tapering toward the
anode (the magnetic-lens shape `assembly/field_quality.py` scores
against). The excitation that meets it is **2 A through every coil**,
turns unchanged at inner 300 / outer 200.

Two independent derivations agree:

- **Lumped reluctance**: 300 G in air needs H = 23.9 kA/m; across the
  37.4mm inner-core-to-outer-core air path that's ~450-900 A-turns
  depending on how much of the path holds full field. The iron leg
  contributes ~2.5 A-turns (227mm at ~0.1 T), i.e. nothing -- this is a
  gap-dominated circuit.
- **The solve itself**: at 2 A the measured field is 292 G at the exit
  plane, peaking at 319 G about 4mm inboard of it.

Pin 300 G *at the exit plane* and you want 2.05 A; pin it as the channel
*peak* and you want 1.88 A. Either way 2 A is within a few percent, so
it's the baseline.

For contrast, the previous placeholder was 5 A -- inherited from
`../matlab/het_solenoid_bfield.m`'s parametric check, never derived from
a target. That supplies 2500 A-turns, roughly 3x the requirement, and is
what drove the iron past saturation and made the nonlinear solve
intractable (see "Nonlinear iron" above).

## Channel ROI

The channel bounds in `assembly_params.txt` (`channel_inner_r`,
`channel_outer_r`, `channel_y_anode`, `channel_y_exit`) define the **real
plasma cavity** -- the open annulus between the chamber liner's ceramic
walls. Confirmed against the CAD: ID 1.74 in, OD 2.5 in, 1.375 in axial
span, exit at y = -1.25 in, anode end at y = +0.125 in.

This used to be the chamber *part's* bounding annulus (r 17.30-41.44mm,
y -32.51..+5.59mm) -- i.e. the whole ceramic liner including both walls.
That is not the plasma channel, and it mattered: that box contained 5392
iron elements at the exit end, so channel field statistics were averaging
in pole-face material. Correcting it moved radial purity from 82.5% to
97.7% and removed a spurious 2.2x jump between adjacent axial stations.

`detect_channel_cavity()` in `build_assembly.py` finds the cavity from
the chamber geometry so this survives a geometry sweep. Two traps it
documents, both hit while building it:

- It needs chamber **volume** element centroids, not surface nodes. A
  surface mesh puts no points inside solid material, so a wall's interior
  reads identically to open space -- with surface nodes it confidently
  returned the outer wall's guts (r 31.98-40.94mm) as the "cavity".
- The detected band must be trimmed clear of the wall faces before
  locating the closed end, or a fraction of a bin of wall overlap reads
  as material at every y and collapses the channel to 0.01mm long.

**Known gap**: the detector lands ~0.76mm off the CAD-confirmed bounds
(it reads centroids, which sit half an element inside the true faces), so
re-running `build_assembly.py` will overwrite hand-confirmed values with
slightly shifted ones. Either snap its output to the nearest real chamber
face or take the channel dimensions as explicit parameters.

## Whole-assembly STEP import (now the live pipeline)

`assembly/scratch_whole_assembly_import.py` was the spike that validated
this; `build_assembly.py` now uses it for real, replacing the old
per-part-STEP-export + hardcoded-placement-table pipeline entirely
(`TRANSFORMS`, ported from `matlab/new_export/apply_transforms.m`, is
gone). Findings that carried over from the spike into the real pipeline:

- Gmsh/OCC recovers Onshape's own real per-occurrence placement exactly
  and preserves part names as entity labels -- no manual per-file
  transform bookkeeping needed. Confirmed body names (2026-09-18 pull):
  `top_plate`, `bottom_plate`, `inner_emag_core`, `outer_emag_core` (x4)
  -- iron; `inner_coil`, `outer_coil` (x4) -- real winding solids;
  `chamber`, `chamber_spacer`.
- Fusing all iron parts in **one batch call** only produced 3 disjoint
  solids on this file, not 1 -- OCC's fuzzy boolean union is
  order/grouping-dependent (same non-symmetry already noted in
  `../matlab/README.md` for MATLAB's kernel). Fusing **incrementally**
  (one part folded into the growing result at a time) gives the correct
  single solid, and is what `build_assembly.py` does now.
- The usual mm-labeled-as-m unit bug is present in this STEP too --
  `build_assembly.py` applies the `0.001` scale to all imported solids at
  once, right after import (before any boolean op), which the spike found
  more reliable than scaling pre-placed per-part templates.
- **New gotcha found wiring this in for real** (not seen in the spike,
  which only checked placement, not areas/volumes): `occ.getBoundingBox()`
  returns a **loose/conservative box on curved STEP-imported solids** --
  e.g. a coil's true X/Z extent (confirmed via actual mesh nodes) equals
  its OD on both axes, but `getBoundingBox()` reports the right value on
  one axis and up to ~1.75x too large on the other, for the *same* solid,
  reproducible with no boolean ops involved at all. Only the axis bounded
  by curved (not flat end-cap) surface is affected. `build_assembly.py`
  no longer uses `getBoundingBox()` for anything needing X/Z precision on
  these bodies -- `getMass()`/`getCenterOfMass()` (volume integrals over
  the true trimmed solid, not boundary-based) are trustworthy and used
  instead.

## Nonlinear iron (B-H curve) -- converged

Replaced the `mur_iron = 1000` linear placeholder with a real saturating
B-H curve, since the linear solve was reporting unphysical peak |B|
(9.7 T under the old geometry, 2.16 T under the new whole-assembly
geometry at the old placeholder excitation -- both well past where
real soft iron saturates, ~1.5-2 T).

**The short version of everything below**: two attempts at this failed
while the excitation was still the inherited 5 A placeholder, and the
fix turned out not to be a better nonlinear method but a correct
operating point. At 5 A the poles are driven to 2.16 T, past the B-H
knee, and no amount of load-stepping made that cheap. At the real
~300 G design current (2 A) peak iron is 0.86 T, the material is
effectively linear, and the solve converges in 11 Picard iterations /
19.6 min / 805 MB with a 2-stage ramp. The debugging history below is
kept because the failure modes are informative, not because the
16-stage ramp is still needed.

- **Material**: generic soft steel, using GetDP's own bundled
  `SteelGeneric` dataset (`tools/getdp-3.5.0-Windows64/templates/
  Lib_Materials.pro`, 49-point H/B table) -- copied inline into
  `generate_regions.py` (not `Include`d from `tools/`, which is
  gitignored/machine-local) rather than sourcing a real alloy datasheet,
  since the actual BPL-700 pole material isn't specified. Reluctivity is
  built as a piecewise-linear interpolation of `H/B` vs `B^2`
  (`InterpolationLinear[SquNorm[$1]]{...}`, the same derivation
  `Lib_Materials.pro` itself uses), assigned to `Iron` only --
  `Air`/`Windings` stay linear.
- **Formulation**: `magnetostatics_assembly.pro`'s single
  `nu[] * Dof{d a}` Galerkin term is now split into a linear part
  (`Air+Windings`) and a nonlinear (Picard) part over `Iron`
  (`nu[{d a}]`), evaluated at the previous iteration's field -- same
  `Vol_L_Mag`/`Vol_NL_Mag` split GetDP's own `Lib_Magnetostatics_a_phi.pro`
  template uses. Full Newton (with a `dhdb[]` tangent term) was the
  template's own default, but its exact tangent-derivative macro syntax
  (`SquDyadicProduct[#1]`-based) wasn't confidently reproducible without
  risking a silently-wrong-not-obviously-broken linearization, so Picard
  was chosen for v1 as the lower-risk option.
- **Plain Picard diverges from a cold start**: tried first, confirmed
  empirically -- residual dropped 613 -> 96 over 2 iterations, then
  exploded to 1.7e5 on the 3rd. Root cause: this B-H curve's slope gets
  extremely steep near saturation (H jumps from 3.1e5 to 7.6e5 A/m over
  the last few tenths of a Tesla), and the linear solve's own peak
  (2.16 T) starts the very first nonlinear iteration right in that steep
  region -- too large a jump for unrelaxed Picard to track.
- **Fix: load-stepping (continuation)**. Excitation is scaled by a
  runtime `$IFrac` (in `Js[]`, `generate_regions.py`), ramped up in
  stages from a small fraction to full current in
  `magnetostatics_assembly.pro`'s Resolution, Picard-converging fully at
  each stage before stepping up -- every stage continues from the
  previous stage's already-converged field (no `InitSolution` between
  stages), so each individual Picard step only tracks a small change
  instead of jumping from zero into saturation.
- **First attempt (8 uniform stages, `NL_iter_max=20`, `NL_tol_rel=1e-6`)
  ran to completion but did not converge.** Confirmed empirically
  (~97 min wall / ~10.3 CPU-hr): stages up to 37.5% converged cleanly (2,
  4, 9 iterations); the 50% stage was still improving *monotonically*
  every iteration (never oscillating) but hit the 20-iteration cap at rel.
  residual ~2e-3 -- genuinely converging, just under-budgeted, not
  diverging. Because that stage got cut off unconverged, the next stage
  (62.5%) inherited a not-quite-right field and started oscillating; by
  75-100% the residual was exploding 1-3 orders of magnitude per
  iteration. The saved (garbage) result had peak |B| = 6.5 T, confirming
  it wasn't a real solution. Root cause: iteration/tolerance budget, not
  an unreachable solution.
- **Second attempt: finer ramp (16 non-uniform stages, denser above 40%),
  `NL_iter_max=60`, `NL_tol_rel` relaxed to 1e-4.** This behaved much
  better -- stages 0-3 (10/20/30/40%) all converged cleanly and
  monotonically (2, 4, 6, 9 iterations respectively, well inside the new
  budget), no oscillation at all through 40% (vs. the first attempt's
  visible strain already building by 37.5%). **Killed by the harness at
  ~53 min in (partway into stage 4, 45%) due to host-level memory
  pressure** -- not a divergence or a bug in the solve itself. Both
  attempts showed the same pattern: `getdp.exe`'s resident memory
  oscillates (grows during a Generate/Solve/factorization, drops after)
  rather than climbing monotonically, but the *peaks* got large enough
  (workstation has 31 GB RAM; observed available-memory dips as low as
  ~4.5 GB on the second run) that the environment's own low-memory
  protection stepped in before stage 4 could finish. Likely cause:
  MUMPS's LU factorization data isn't being released between the many
  repeated `Generate[]`/`Solve[]` calls this manual iteration loop makes
  (16 stages x up to several iterations each, all against a 589k-DOF
  system) -- plausible but not root-caused.
- **Resolution: the excitation was wrong, not the solver.** Both failures
  above were chasing a saturation regime the design never called for. A
  lumped reluctance check against the 300 G channel target (see
  "Operating point" below) needs ~450-900 A-turns; 5 A was supplying
  2500. Dropping to 2 A puts peak iron at 0.859 T -- below the knee,
  where this B-H curve is nearly straight -- and the same formulation,
  mesh and Picard iteration then converge easily:

  | | 5 A, 16 stages | 2 A, 2 stages |
  |---|---|---|
  | Picard iterations | hit caps, oscillated | 11 total, monotonic |
  | wall time | ~97 min / killed at 53 min | 19.6 min |
  | peak memory | host dipped to ~4.5 GB free | 805 MB |
  | peak iron \|B\| | 2.157 T (past knee) | 0.859 T |
  | result | garbage (6.5 T) / unfinished | converged, rel 9.9e-5 |

  Contraction was clean throughout: stage 0 converged in 3 iterations
  (50x/28x/24x per step), stage 1 in 8, never oscillating. The host
  memory pressure that killed the second attempt simply never arises with
  ~1/8th as many large factorizations.
- **A useful corollary**: at 0.859 T the nonlinear peak matches a
  linear-scaled prediction from the old `mur_iron=1000` run to 0.4%.
  That is not a coincidence -- this is a gap-dominated circuit (iron drop
  ~2.5 of ~900 A-turns), so the answer is insensitive to iron
  permeability as long as it's large. Practical consequence:
  superposition is valid at the design point, so `current_sweep.py`'s
  cheap 2-basis-solve sweep is trustworthy again.
- **If a future design point does need saturation**, restore a fine ramp
  (git history has the 16-stage schedule) and budget accordingly, or
  switch to full Newton with the `dhdb[]` tangent -- the Picard/
  load-stepping combination is genuinely slow there.

## Known limitations / next steps
- **chamber/injector geometry not in the solve.** Fine for magnetostatics
  (inert, mu_r ~= 1) but would need reintroducing -- with a real inflate/
  clearance-gap fix, not a knife-edge boolean cut -- for any future
  thermal or structural model where their exact shape matters.
- **Excitation magnitude is now derived, the turns count isn't.** 2 A
  comes from the 300 G channel-exit target (see "Operating point"), and
  the 1:1 center/outer ratio is separately the most radial-field point
  the sweep found. What's still inherited is the 300/200 turns split
  itself, and the per-wire gauge/resistance/power numbers in "Solenoid
  winding specs" remain a downstream hand-calc the FEM solve doesn't
  model. The 300 G target also hasn't been traced back to a
  discharge-voltage/Larmor-radius requirement for this specific
  thruster -- `../sizing/basic_het_sizing.py` assumes 150 G for its
  magnetization check, and reconciling those two numbers is open.
- **`cross_section_viewer.py` now cuts the mesh exactly**, via
  `pyvista`'s `mesh.slice()` on `assembly_field.vtu` -- a real plane-mesh
  intersection, not an approximation. An earlier version inverse-distance
  -weight-averaged the 4 nearest element centroids on a regular query
  grid (a stand-in for MATLAB's `interpolateMagneticFlux`); that's gone
  now. Since `B` is piecewise-constant per element (see `HOW_IT_WORKS.md`
  sec. 9), each slice triangle just shows its parent element's true flat
  value (`ax.tripcolor(..., facecolors=...)`), which is actually the
  physically correct thing to show, not merely more convenient. One minor
  cosmetic side effect of no longer smoothing: the quiver arrows are
  noisier in the near-zero-field air region, since real per-element
  values (including small numerical near-cancellation vectors right at
  source/material boundaries) are shown unsmoothed -- harmless, but
  visually busier than the old averaged version.
- **The vacuum/air region was never missing** -- `Air` already covers
  the whole non-iron, non-winding space (including the channel) and was
  always part of the solve and `b_assembly.pos`. It just didn't *look*
  present: the channel's real field (~5-32 mT at the design point) is
  one to two orders of magnitude weaker than the iron near the poles
  (0.86 T), so sharing one color scale made it look empty.
  `cross_section_viewer.py`'s "region of interest" checkbox fixes the
  *display*, not the physics: it crops to the channel box and rescales
  color to a log range fit to the channel's own local |B|, which is what
  actually reveals the field structure there. Both that viewer and
  `post_process.py` now take their full-view color ceiling from the
  loaded solve rather than a hardcoded constant -- the old fixed 3 T
  ceiling was sized for the 9.7 T linear era and renders the current
  0.86 T solve almost entirely black.
- **Mesh resolution is a single pass**, not refinement-checked -- no
  convergence study yet (compare a coarser/finer mesh's peak/median |B|).
