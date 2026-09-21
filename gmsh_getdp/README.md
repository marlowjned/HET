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
  GetDP's magnetic-vector-potential formulation. Result: 27.63 mT on-axis
  at coil center, vs. 28.02 mT from the exact analytic on-axis integral
  for a finite thick solenoid (1.4% agreement) -- confirms the mesh,
  formulation, and gauge are all correct before trusting them on the real
  CAD. Run `python check_toy.py` after the solve: it asserts the
  comparison rather than leaving it as a number quoted here, which is what
  let a real error hide for months -- see "The A_cross bug" below.
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
  separate core-diameter parameter). `chamber_spacer` is iron (confirmed
  2026-09-21: the same A36 as the plates; `--exclude-chamber-spacer`
  still tests the inert case).
- **Center/outer current-ratio sweep**: done at the old geometry and the
  old linear-iron material (pole names `center_solenoid`/`outer_solenoid`,
  since renamed to `inner_coil`/`outer_coil` -- `current_sweep.py` itself
  is updated for the rename, but hasn't been re-run against the new
  whole-assembly geometry, and its channel statistics predate the ROI
  correction below). Found the most radial (E x B - correct) channel
  field at `I_center = I_outer` -- see "Current ratio sweep" below. **Its
  2-basis-solve superposition trick assumes a linear material**, which
  looked fatal once Iron became nonlinear -- but at the real ~300 G design
  point the iron peaks at 1.023 T, below the B-H knee, so it behaves
  near-linearly and the trick is roughly valid there. Note this got
  *weaker* when the material became A36-proxy rather than the old generic
  curve: 1.023 T is past that curve's peak-permeability point (0.865 T),
  so superposition is now an approximation rather than near-exact.
  Re-check it if a design point pushes the iron further.
- **Materials**: nonlinear iron (ASTM A36, run on a proxy B-H curve) replaces
  the `mur_iron = 1000` linear placeholder, and **converges** -- 33
  Picard iterations, 68 min wall, 806 MB peak. See "Nonlinear iron
  (B-H curve)" below for how the excitation, not the method, was what
  made that work. The iteration count roughly tripled when the material
  went from the generic curve to the A36 proxy: the design's peak iron
  field sits past that curve's permeability peak, where Picard contracts
  slowly. Worth budgeting for in the geometry sweep.
- **Excitation**: 1.30 A per coil (inner 260 t, outer 165 t) -- a real
  winding design sized to the ~300 G channel-exit target, not the old
  inherited 300/200 at 5 A. Measured 288.6 G at the exit plane at 1.25 A,
  hence 1.30 A for 300 G. See "Operating point" below.
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
python check_toy.py                             # asserts the FEM against the exact analytic field; exits non-zero on failure

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

## Spacer-height sweep (`emag_sweep.py`)

Sizes the `chamber_spacer`. Full write-up, figures and every coil geometry:
**`../EMAG_REPORT.pdf`**. Geometry facts it rests on,
measured from the STEP exports:

- **spacer = emag_height - 1.125 in.** The chamber and top plate never
  move; emag_height moves only the bottom plate. At 1.125 in Onshape drops
  the spacer body entirely (`build_assembly.py` accepts that).
- The chamber bore and spacer bore are both the inner coil's OD
  (r 17.46 mm), so the inner winding can only grow inward by thinning the
  core (`inner_coil_id`). The pole cap stays 1.375 in regardless.

Grid: emag_height 1.125-2.0 in (spacer 0-0.875 in, 0.125 steps) x inner
core 1.000/0.875/0.750 in, plus one outer-OD 1.625 in check -- 25 points.
**Linear iron** (mu_r 1500): +1.25% on exit field vs the committed nonlinear
solve at equal current (`--validate`), ~2.5 min/point instead of ~70.

Results (`assembly/emag_sweep_results.csv`), all scaled to 300 G at the exit:

- A-turns for 300 G: inner 326-339, outer 207-215 each, across the whole
  grid. Spacer height moves it ~1% (= mesh noise); a 0.750 in core adds ~2%.
- Core shaft flux 0.34-0.74 T, highest at the **longest** spacer (the steel
  spacer wraps the inner coil and pulls leakage flux through the core).
  Knee on the measured 1020 curve is 1.3-1.5 T.
- The spacer mostly changes winding room: AWG 20 at the 1.000 in core goes
  1.27 A / 5.1 W (spacer 0.875) -> 2.25 A / 10.0 W (spacer 0). A 0.750 in
  core with no spacer is 1.28 A / 5.4 W -- today's numbers.

Worst-case nonlinear check (`--nonlinear h2.000_ci0.750_oo1.500`, 0.74 T
core): the 50% stage converged in 9 Picard iterations, but the 100% stage
**stalled** at rel. residual ~0.11 for all 60 iterations (2.4 h) -- not
converged. Unconverged, it matches linear to 0.3% at the exit and puts the
inner core 3.4% lower; its p99.9 iron |B| is 13% higher (a few corner
elements, plausibly where Picard stalls). A real check needs a finer ramp
or Newton.

Two pipeline fixes came out of it: `build_assembly.py`'s mesh fine size is
now capped at the baseline value (a thinner core coarsened the whole mesh,
which failed on the bottom plate's 2 mm lead holes and made resolution vary
across the sweep), and a meshing exception now finalizes gmsh (leaving it
initialized made the *next* build silently write an empty mesh).

## Operating point

Target: **~300 G radially at the channel exit**, tapering toward the
anode (the magnetic-lens shape `assembly/field_quality.py` scores
against). The design point that meets it is **338 A-turns on the inner
pole and 215 on each outer pole** -- realized as 260 / 165 turns at
1.30 A, wired in series.

Two independent derivations agree:

- **Lumped reluctance**: 300 G in air needs H = 23.9 kA/m; across the
  37.4mm inner-core-to-outer-core air path that's ~450-900 A-turns
  depending on how much of the path holds full field. The iron leg
  contributes ~2.5 A-turns (227mm at ~0.1 T), i.e. nothing -- this is a
  gap-dominated circuit.
- **The solve itself** (A36-proxy iron): 288.6 G at the exit plane at
  1.25 A, peaking at 308.5 G about 3.5mm inboard of it. Scaling to 300 G
  gives 1.30 A. Note the iron is no longer strictly linear at this point,
  so that scaling is slightly optimistic.

For contrast, the original placeholder was 5 A at 300/200 turns,
inherited from `../matlab/het_solenoid_bfield.m`'s parametric check and
never derived from a target. That is ~3x the required A-turns, and is
what drove the iron past saturation and made the nonlinear solve
intractable (see "Nonlinear iron" above).

## Solenoid winding design

Only the **A-turn product** reaches the FEM (as a bulk current density,
see HOW_IT_WORKS.md sec. 2), so the wire choice below changes the
current/voltage split and nothing about the field. Every row delivers the
same 325/206 A-turns and therefore the same 300 G.

Windows measured from the real CAD coil solids: **4.76 x 50.8 mm**
(inner) and **6.35 x 50.8 mm** (outer); mean turn lengths 94.8 and
99.7 mm. Turn counts assume fiberglass-served wire, hexagonally nested,
with the inner coil wound full and the outer coils wound to 0.63x that so
one series current gives the 1.58 A-turn ratio the field wants.

| AWG | N inner | N outer | layers in/out | current | V @ 20 C | V @ 60 C | total P |
|---|---|---|---|---|---|---|---|
| 16 | 105 | 67 | 3 / 2 | 3.10 A | 1.46 V | 1.69 V | 5.2 W |
| 18 | 172 | 109 | 4 / 3 | 1.89 A | 2.30 V | 2.67 V | 5.0 W |
| **20** | **260** | **165** | **5 / 4** | **1.30 A** | **3.81 V** | **4.42 V** | **5.8 W** |
| 22 | 384 | 244 | 6 / 4 | 0.85 A | 5.84 V | 6.76 V | 5.7 W |
| 24 | 624 | 396 | 8 / 6 | 0.52 A | 9.26 V | 10.72 V | 5.6 W |

AWG 20 is what `build_assembly.py`'s `EXCITATION` is set to. Switching
rows needs no re-solve -- just the turns/current pair.

**Dissipation is ~5 W in every row**, which is the useful structural
result here:

```
P = (N*I)^2 * rho * L_turn / (k * A_window)
```

Wire gauge cancels: thinner wire raises resistance exactly as fast as the
lower current reduces I^2. The only levers on heat are window area,
packing factor, turn length, and the A-turn requirement itself. Choosing
a gauge is a power-supply matching decision, not a thermal one.

Conductor current density is 2.4-2.7 A/mm^2 across the table --
comfortable even with vacuum/radiation-only cooling. A standalone
radiation balance (all 5 W leaving a ~0.070 m^2 envelope by radiation
alone, no conduction path) puts the assembly at roughly 42-60 C depending
on surface emissivity, so the magnet circuit is not thermally
interesting on its own.

**Soft inputs**: the insulated-diameter figures are representative of
single glass serving rather than taken from a datasheet, and they set the
turn counts -- a real wire spec could move N by 10-20%, which moves the
current, not the field or the power. Wiring in series is what forces
equal current through coils of unequal resistance; running the coils at
different currents needs two supplies or a trim resistor on one leg.

## Iron material data (ASTM A36)

The poles and plates are ASTM A36, from SendCutSend hot-rolled
pickled-and-oiled stock -- their 1008 cold rolled tops out at 0.135 in,
too thin for anything here, and A36 is stocked at exactly the 0.375 in
plate thickness. See `../MATERIALS.md` for the full property table and
sourcing. Magnetic data is a 19-point DC
magnetization curve inline in `assembly/generate_regions.py`; thermal and
physical properties are in `../thermal/radiation_balance.py`.

| property | value | use |
|---|---|---|
| saturation (B at 1000 Oe) | 2.165 T | magnetics (proxy curve) |
| peak relative permeability | 2164, at H = 318 A/m | magnetics (proxy curve) |
| mu_r at this design's peak iron field | ~1935 at 1.023 T | magnetics (proxy curve) |
| thermal conductivity | 50 W/m-K | justifies the isothermal-iron assumption |
| specific heat | 470 J/kg-K | transient model, when one exists |
| density | 7900 kg/m^3 | as above |
| max service temperature | 400 C | thermal limit |
| Curie point | ~770 C | magnetic circuit fails well before this |

**Read the provenance before trusting the magnetic curve.** It is **AISI
1008 data used as a proxy for A36**, deliberately not relabelled -- no
numeric A36 table could be obtained. The 1008 points themselves come from
the Ansys Maxwell SV material library, transcribed via a public
engineering forum; the primary source could not be re-fetched (403).
Neither a measurement of A36 nor of the actual BPL-700 stock.

A36 carries up to 0.26% carbon against 1008's ~0.08% plus up to 1.2% Mn,
so **real A36 is less permeable than this curve and the model is
optimistic about the iron** -- plausibly a few percent of channel field,
i.e. a design current nearer 1.35 A. See `../MATERIALS.md`.

What supports them being genuine digitized data: every H value is an exact
Oersted multiple (2, 4, 6, 8, 10, 20, 40 ... 5000 Oe), which is how
pre-SI magnetization tables were published; the curve is monotonic in both
variables; differential permeability rises then falls as a real
magnetization curve must; and 2.165 T at 1000 Oe sits squarely in the
2.1-2.2 T band expected of low-carbon steel.

What should keep you cautious:

- **1008 is not an electrical steel.** ASTM imposes no magnetic
  requirement on it at all, so two suppliers' 1008 can differ
  substantially and nothing obliges either to match this curve.
- **Processing history dominates.** Peak mu_r of 2164 is characteristic of
  as-rolled rather than annealed material; annealed low-carbon steel
  reaches 3000-5000. If the real poles are annealed they will outperform
  this model slightly.

Fortunately this matters less here than it would in most magnetic designs,
because the circuit is gap-dominated: the iron carries ~2.5 of the ~900
A-turns. Switching from the old generic curve to 1008 cut permeability at
the operating point by 2.7x (mu_r 5700 -> 2150) and moved the channel
field by well under a percent. A measured curve on the actual stock is
still the only way to do better, but it would buy very little.

## The A_cross bug (fixed 2026-09-19, worth knowing about)

Both `build_assembly.py` and `build_toy.py` computed the winding's
current-density area as the **annulus**, `pi*(ro^2 - ri^2)`. The current
is azimuthal, so it crosses a plane *containing* the axis: the correct
area is the coil's r-z section (radial thickness x height). The annulus
is the area for *axial* flow, and it is 1.84x (inner) / 1.94x (outer) too
large here, so every solve silently applied about half its nameplate
A-turns. The toy case was 2.75x off.

**The toy validation could not catch it**, because `build_toy.py` fed the
same wrong area into both the FEM source term and its analytic reference.
The two agreed to 1.4% at 10.19 mT while the true answer was 28.02 mT --
a comparison that cannot fail is not a validation. That is why
`check_toy.py` now exists and asserts, and why `A_cross` is derived via
Pappus's theorem (`A = V / 2*pi*r_mean`) from the solid's own volume
rather than from nominal dimensions.

Second-order lesson from the same fix: `magnetostatics_toy.pro` held a
hand-copied `Jmag` literal, so the first corrected re-solve silently
reused the old value. It now `Include`s a file `build_toy.py` generates,
matching how the assembly pipeline already worked.

Everything *shape*-based from before the fix survived it -- the ROI work,
the lens profile, radial purity, the saturation conclusion -- because a
scalar error in J rescales the field without changing its geometry. What
changed is the label on the current axis.

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
~300 G design current (1.30 A, 338/215 A-turns) peak iron is 1.023 T,
the material is effectively linear, and the solve converges in 12 Picard
iterations / 28 min / 806 MB with a 2-stage ramp. The debugging history below is
kept because the failure modes are informative, not because the
16-stage ramp is still needed.

- **Material**: **ASTM A36**, modelled with a 19-point AISI 1008 proxy H/B table inline in
  `generate_regions.py`. This replaced GetDP's bundled generic
  `SteelGeneric` dataset, which was a stand-in from when the pole material
  was unspecified. Reluctivity is built as a piecewise-linear
  interpolation of `H/B` vs `B^2`
  (`InterpolationLinear[SquNorm[$1]]{...}`, the derivation GetDP's own
  `Lib_Materials.pro` uses), assigned to `Iron` only -- `Air`/`Windings`
  stay linear. See "Iron material data" below for provenance and the
  caveats, which are real.
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
  2500. Dropping to the real design A-turns puts peak iron at 1.023 T -- below the knee,
  where this B-H curve is nearly straight -- and the same formulation,
  mesh and Picard iteration then converge easily:

  | | 5 A placeholder, 16 stages | design point, 2 stages |
  |---|---|---|
  | Picard iterations | hit caps, oscillated | 12 total, monotonic |
  | wall time | ~97 min / killed at 53 min | 28 min |
  | peak memory | host dipped to ~4.5 GB free | 806 MB |
  | peak iron \|B\| | 2.157 T (past knee) | 1.023 T |
  | result | garbage (6.5 T) / unfinished | converged, rel 9.9e-5 |

  Contraction was clean throughout: stage 0 converged in 3 iterations
  (48x/27x/23x per step), stage 1 in 9, never oscillating. The host
  memory pressure that killed the second attempt simply never arises with
  ~1/8th as many large factorizations.
- **A useful corollary**: on the old generic curve the nonlinear peak matched a
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
- **The excitation is derived now; the 300 G target itself isn't.** Both
  the A-turns and the turns/current split come from the field requirement
  plus the real coil windows (see "Operating point" and "Solenoid winding
  design"). What has *not* been traced back to first principles is the
  300 G number -- `../sizing/basic_het_sizing.py` assumes 150 G for its
  magnetization check, and reconciling those two is open. The per-wire
  resistance/voltage/power figures also remain a hand-calc the FEM
  doesn't model, and rest on representative rather than datasheet
  insulated-wire diameters.
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
