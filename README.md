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
- **`onshape/`** — HMAC API-key auth layer for scripting against a
  parametrized Onshape CAD document (electromagnet height/spacer/iron
  radius as configuration variables), in support of a geometry sweep. Auth
  validated end-to-end; STEP export per configuration not yet built. See
  `onshape/README.md`.

## How the pieces fit together

`sizing/` → target channel geometry and B-field (not yet actually wired
into the CAD dimensions below — the CAD is the existing real BPL-700
hardware, not a from-scratch design derived from this script). CAD
geometry and per-part placement transforms originate in `matlab/`
(`new_export/apply_transforms.m`) and are reused as hardcoded values in
`gmsh_getdp/assembly/build_assembly.py`, since re-deriving them from the
STEP files a second way isn't necessary. Material list (which parts are
iron vs. inert) is documented once, in `matlab/README.md`, and consumed by
both pipelines. From here, active development is in `gmsh_getdp/` —
`matlab/` is not being advanced further unless the union blocker gets
resolved on the Onshape side.

## Current overall status

- Real-geometry magnetostatic B-field solve: **working**, in
  `gmsh_getdp/assembly/`, with a placeholder linear-iron material model
  (no B-H saturation curve yet) and placeholder winding currents.
- Top-level thruster sizing (power/thrust/Isp from design targets):
  **working**, standalone, in `sizing/`.
- Not yet started: wiring `sizing/`'s targets into the CAD dimensions,
  nonlinear iron (B-H curve), electric-field/plasma simulation, electron
  trajectory/stability, thermal, and structural/force analysis.
