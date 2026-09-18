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
- **Real assembly**: done and solving. `assembly/` imports the real
  BPL-700 iron parts, places them with the same per-occurrence transforms
  as `../matlab/new_export/apply_transforms.m`, fuses them into one solid,
  builds the 5 windings, meshes, and solves. **All 9 real touching
  interfaces mesh cleanly** (vs. MATLAB's 4/9) -- see "The geometry win"
  below.
- **Center/outer current-ratio sweep**: done. `assembly/current_sweep.py`
  found the most radial (E x B - correct) channel field at
  `I_center = I_outer` -- see "Current ratio sweep" below.
- **Materials**: linear placeholder only (`mur_iron = 1000`, no B-H
  saturation curve). See "Known limitations" below.
- **chamber / injector**: deliberately excluded from the solid model (see
  `assembly/build_assembly.py` docstring) -- they're magnetically inert
  (mu_r ~= 1, same as air), so omitting them doesn't affect the B-field
  solve, and it sidesteps two more geometry problems (a genuine 3D overlap
  sliver between chamber and iron, and injector's independently dirty STEP
  tessellation -- both already flagged in `../matlab/README.md`).

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

This is a field-direction optimum only -- it says nothing about the
field *magnitude* needed for a real operating point (that depends on
electron Larmor radius / discharge voltage requirements, not modeled
here), and the underlying solve is still the linear/no-saturation
placeholder (see "Known limitations" below).

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

At the swept operating point (`I_center = I_outer = 5 A`, the most
radial-field ratio found above):

| | outer coil (each) | center coil | total (1 center + 4 outer) |
|---|---|---|---|
| voltage @ 20 C | 0.70 V | 1.11 V | -- |
| voltage @ 100 C | 0.92 V | 1.46 V | -- |
| power @ 20 C | 3.49 W | 5.56 W | 19.5 W |
| power @ 100 C | 4.59 W | 7.31 W | 25.7 W |

If all 5 coils are wired in series on one supply at this operating point
(same 5 A through all of them, true only at this 1:1 ratio), the whole
magnet circuit needs **~3.9-5.1 V at 5 A, ~20-26 W total** depending on
winding temperature -- a low-voltage, low-power circuit, consistent with
running off a small bench supply rather than anything resembling the
discharge (anode) power supply. Moving off the 1:1 ratio breaks the
series assumption, since `I_center != I_outer` then requires either two
independent supplies or a shunt/trim resistor on one leg; scaling within
that constraint, power grows with the square of current on whichever
coil's current changes (e.g. `I_center` swept from 2-12.5 A above swings
center-coil power alone from ~0.9 W to ~35-46 W at fixed `I_outer = 5A`).

## Whole-assembly STEP import (spike, not yet in the build pipeline)

`assembly/scratch_whole_assembly_import.py` validates a simpler alternative
to the current per-part-STEP-export + hardcoded-placement-table pipeline
(`build_assembly.py`'s `TRANSFORMS`, ported from
`matlab/new_export/apply_transforms.m`): importing the **whole-assembly**
STEP directly. Findings, relevant to any future parametric/sweep pipeline:

- Gmsh/OCC recovers Onshape's own real per-occurrence placement exactly
  (verified against `apply_transforms.m`'s table) and preserves part names
  as entity labels -- no manual per-file transform bookkeeping needed.
- Fusing all 7 iron parts in **one batch call** (as `build_assembly.py`
  does today) only produced 3 disjoint solids on this file, not 1 --
  OCC's fuzzy boolean union is order/grouping-dependent (same non-symmetry
  already noted in `../matlab/README.md` for MATLAB's kernel). Fusing
  **incrementally** (one part folded into the growing result at a time)
  gave the correct single solid (73 boundary faces, meshes cleanly).
- The usual mm-labeled-as-m unit bug is present in this STEP too, but
  fixing it via `occ.dilate()` on solids pulled from an *assembly*-
  structured STEP is unreliable (produces small, part-specific placement
  drift -- see the script's docstring for the suspected cause). Didn't
  break this particular result, but shouldn't be trusted going forward.
  The real fix: have any future parametrized Onshape export emit correct
  real-world units directly, rather than relying on a post-import scale
  correction.

Not wired into `build_assembly.py` yet -- kept as a validated reference
for whichever CAD pipeline design gets built next.

## Known limitations / next steps

- **Linear iron, no saturation.** `mur_iron = 1000` is a reasonable
  placeholder, not a real B-H curve -- the solve shows |B| up to ~9.7 T in
  the iron near the windings (median ~2.7 T), well past where real soft
  iron saturates (~1.5-2 T). GetDP supports nonlinear `nu[]` via a B-H
  interpolation function; needs the actual alloy's B-H curve (same gap
  noted in `../matlab/README.md`).
- **chamber/injector geometry not in the solve.** Fine for magnetostatics
  (inert, mu_r ~= 1) but would need reintroducing -- with a real inflate/
  clearance-gap fix, not a knife-edge boolean cut -- for any future
  thermal or structural model where their exact shape matters.
- **Excitation values are placeholders**, not a designed operating point:
  same turns/current as `../matlab/het_solenoid_bfield.m`'s parametric
  check (center 300t/5A, outer 200t/5A each, opposite polarity). The
  current-ratio sweep (see above) confirms this 1:1 current ratio also
  happens to be the most radial-field point found, but the magnitude (5 A)
  is still arbitrary -- no discharge-voltage/Larmor-radius requirement has
  been used to pick it, and the per-wire gauge/resistance/power numbers in
  "Solenoid winding specs" above are a downstream hand-calc, not something
  the FEM solve itself models.
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
  present: the channel's real field (~1-30 mT) is ~2-3 orders of
  magnitude weaker than the saturated iron near the poles (~9.7 T), so
  sharing one 0-3T color scale made it look empty. `cross_section_viewer
  .py`'s "region of interest" checkbox fixes the *display*, not the
  physics: it crops to the channel box and rescales color to a log
  range fit to the channel's own local |B|, which is what actually
  reveals the field structure there.
- **Mesh resolution is a single pass**, not refinement-checked -- no
  convergence study yet (compare a coarser/finer mesh's peak/median |B|).
