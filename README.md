# HET — Hall Effect Thruster simulation

Simulating the magnetic circuit (and eventually more) of the BPL-700 Hall
Effect Thruster, from first-pass sizing through real-CAD magnetostatic FEM
solves. Two independent geometry/FEM toolchains exist side by side because
the first one (MATLAB) hit a real geometry-kernel limitation partway
through; rather than being redundant, the second (Gmsh+GetDP) is the one
that's currently ahead and validated. This file is a map of the repo, not a
status report — each subfolder's own README/notes carry current status.

## Layout

- **`sizing/`** — `basic_het_sizing.py`: zero-CAD, textbook top-level
  sizing (Goebel & Katz scaling laws) — discharge current, channel
  geometry, thrust/Isp/efficiency, and a Larmor-radius magnetization
  check, from power/voltage targets. Upstream of everything else: this is
  where target channel dimensions and required B-field would come from if
  a new design point is needed, rather than the current placeholder
  numbers used downstream.
- **`matlab/`** — first magnetostatics pipeline: idealized parametric
  coil model (working), plus a real-CAD (Onshape STEP) import and
  placement pipeline (working). Blocked partway through multi-material
  assignment by a `fegeometry.union()` limitation on 5 of 9 real
  part-to-part interfaces (sliver/degenerate faces) — see
  `matlab/README.md` for the full diagnosis. Still the source of the real
  per-part transforms and material list used by the other pipeline.
- **`gmsh_getdp/`** — second, independent magnetostatics pipeline (Gmsh
  for geometry/meshing, GetDP for the FEM solve, both compiled native
  code) built specifically because it doesn't hit MATLAB's union
  blocker. Currently the most complete part of the project: validated
  against an analytic solenoid check, solves the full real assembly,
  and has run a center/outer current-ratio sweep. See
  `gmsh_getdp/README.md` for status/results and `gmsh_getdp/HOW_IT_WORKS.md`
  for the physics-to-code walkthrough (read that one to understand *why*
  any given step exists, not just what it does).
- **`tools/`** — gitignored, machine-local third-party binaries (the
  Windows GetDP executable). Not checked in; see `gmsh_getdp/README.md`
  "Setup" for how to repopulate it on a new machine.
- **`onshape/`** — HMAC API-key auth layer + caching STEP pull for a
  parametrized Onshape CAD document, keyed on all 5 real Configuration
  parameters (`emag_height`, `inner/outer_coil_id/od`). Auth and the
  multi-param STEP export are both validated end-to-end and wired
  directly into `gmsh_getdp/assembly/build_assembly.py`'s geometry
  ingestion. See `onshape/README.md`.

## How the pieces fit together

`sizing/` → target channel geometry and B-field (not yet actually wired
into the CAD dimensions below — the CAD is the existing real BPL-700
hardware, not a from-scratch design derived from this script). CAD
geometry now comes directly from a whole-assembly Onshape STEP export
(`onshape/cache.py` → `gmsh_getdp/assembly/build_assembly.py`), which
recovers Onshape's own real per-occurrence placement and part names
directly — the `matlab/` pipeline's hand-extracted per-part transforms
are no longer used by `gmsh_getdp/`. Material list (which parts are iron
vs. inert) originates in `matlab/README.md`'s diagnosis but is
re-confirmed against the real STEP body names in
`gmsh_getdp/assembly/build_assembly.py`. From here, active development is
in `gmsh_getdp/` — `matlab/` is not being advanced further unless the
union blocker gets resolved on the Onshape side.

## Current overall status

- Real-geometry magnetostatic B-field solve: **working**, in
  `gmsh_getdp/assembly/`, sourced from a parametrized whole-assembly
  Onshape STEP export (all 5 geometry params) rather than a fixed
  hardcoded geometry. Nonlinear iron (generic soft-steel B-H curve,
  replacing the old linear `mur_iron=1000` placeholder) is implemented
  but not yet confirmed convergent — see `gmsh_getdp/README.md`
  "Nonlinear iron (B-H curve)". Winding currents are still placeholders.
- Top-level thruster sizing (power/thrust/Isp from design targets):
  **working**, standalone, in `sizing/`.
- Not yet started: wiring `sizing/`'s targets into the CAD dimensions,
  the geometry/thermal sweep runsheet and orchestrator across the 5
  Onshape parameters, electric-field/plasma simulation, electron
  trajectory/stability, thermal, and structural/force analysis.
