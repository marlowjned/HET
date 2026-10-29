# HET field and PIC models: assumptions and issues (2026-10-29)

Overnight work on the request: "GO ahead with placeholders for now… Keep working until you get a successful
example. Then work on getting a working 2D PIC example, adding placeholders as needed. Write a report of any
assumptions you make, or any issues you encounter."

Everything is standalone in `HET/pic` (no bplsim). Names follow `gmsh_getdp/axisym/figures/glossary.png`. New
terms are defined where they first appear: anode, chamber wall, cathode plane, iron.

## Summary

| piece | status | headline |
|---|---|---|
| vacuum E-field (`estatics/`) | **done, verified** | iron floating: **198 V** (0.66 of the anode voltage); iron grounded: the vacuum potential dies away inside the channel (1.9 V at the exit) |
| magnet field map (`magnetics/`) | **done, verified** | reproduces the het parity case (inner coil 0.9271 A vs 0.9273 A at 165 G) and writes the openPMD map WarpX reads |
| 2-D PIC (`pic2d/`) | **runs** | see [PIC results](#pic-results); a first-pass model with many placeholders |
| 3-D vacuum E-field (`estatics/run_3d.py`) | **done** | square plates, four posts with coils, a placeholder cathode: iron floating at **184.8 V** (cathode at 0 V), 171.5 V with the cathode 20 V below ground |

## 1. Assumptions shared by all models

| # | assumption | value used | replace with |
|---|---|---|---|
| A1 | Design B magnetic circuit and BN channel part | het library values (`sim/library`: assembly `design_b` over `legacy_bn`; parts `bn_channel_v2`, `plate_6in_a36`, `solenoid_230t_awg20`, `post_coils_184t_awg20`). Channel r 22.12-31.4452 mm (r_in assumed, as in BN_REPORT), depth 41.4274 mm | CAD |
| A2 | **Anode / gas distributor (placeholder)** | one stainless ring filling the bottom 6.35 mm of the channel slot: r 22.12-31.4452 mm, z -47.7774 to -41.4274 mm. BN_REPORT says "the bottom 6.35 mm of the channel is filled". **This differs from the "3 mm thick" anode I proposed before you went to bed.** I used the filled 6.35 mm because the report says the slot bottom is filled. The slot is 1.881 in deep, the full part length, so with the anode in it the BN is two pieces in r-z (inner wall ring, outer body). The anode's back face then sees the back plate across the 3.02 mm gap of the magnet model | the real anode/distributor drawing; how it is held and fed |
| A3 | Iron is one conductor | pole cap, inner core, back plate, outer posts and front plate are bolted together. Square plates as discs of equal area; posts as a solid ring in 2-D (as the magnet model) | 3-D model has the real square plates and four posts |
| A4 | **Coils lumped with the iron** | the windings touch the iron in the model, and the magnet supply is assumed referenced to the thruster body | if the coil supply floats separately, give the coils their own conductor. They sit inside the iron, so the channel field changes little |
| A5 | Cathode not in 2-D | it is mounted outside the thruster, off axis, which a 2-D axisymmetric model cannot hold. The chamber wall at 0 V stands in for ground and the cathode | 3-D model (`Has_cathode` already in `estatic.pro`); real cathode position and its potential to ground (typically -10 to -30 V) |
| A6 | BN relative permittivity | 4.0 (hot-pressed BN is about 4.0-4.2) | grade data sheet |
| A7 | "Large chamber" | grounded half-circle, radius 1 m | see the sensitivity in §2 |
| A8 | Argon flow | 2 mg/s, evenly injected over the anode face: 3.015e19 atoms/s, mass-flow equivalent current 4.831 A | |

## 2. Vacuum E-field (`estatics/`): results and what floating does

Method, figures and numbers: [`estatics/README.md`](estatics/README.md). The formulation (`estatic.pro`) is
checked against concentric spheres with a floating shell (exact solution). The floating potential agrees to
1.2e-4 and the charges to 1.8e-4 and 3.8e-4.

- **Floating iron sits at 198.0 V** (anode 300 V, chamber wall 0 V).
  - The anode couples almost only to the iron: of its 17.54 pF, 0.001 pF reaches the chamber.
  - The iron has 9.0 pF to the chamber, so floating it divides the 300 V: 102 V anode-to-iron and 198 V
    iron-to-chamber.
- **Sensitivity:**
  - chamber size: 192.2 / 198.0 / 200.7 / 202.0 V at 0.5 / 1 / 2 / 4 m;
  - mesh: -0.05 % when the elements are halved;
  - anode-to-back-plate gap (A2): this gap carries a large part of the anode-to-iron capacitance, so the 198 V
    depends on the real anode mount.
- **Grounded iron:** the vacuum potential falls from 300 V at the anode face to 28 V mid-channel and 1.9 V at the
  exit (e-folding about 7 mm). The iron behind the 4.4 mm BN walls screens the channel.
- **Superposition** works as you expected. Two per-volt solves (anode basis, iron basis) give every case. A direct
  floating solve matches the superposition to 5e-12 V, and the capacitance matrix is symmetric to 6e-16.
  - With plasma, the plasma's space charge adds a third term. The basis fields are useful to check the PIC field
    solver and to separate the applied field from the plasma's own field in PIC output.
  - A floating conductor's potential, though, changes with the plasma, as below.

**What difference floating makes.** The 198 V holds only without plasma: before ignition, or if the discharge
goes out. The whole thruster body then sits about 200 V above the chamber, with arcing risk to grounded
hardware.

In operation, a floating surface collects electron and ion current from the plasma. It charges until the two
balance, a few electron temperatures below the local plasma potential. For the iron, that means the pole cap
front face and the front plate front face, both in the near plume. Three consequences:

1. The floating potential is set by the near-plume plasma, not by capacitance. It is typically within some tens
   of volts of the cathode.
2. Ions reach the pole faces with that sheath drop. Pole erosion therefore depends on the choice: floating,
   grounded, or tied to the cathode. Tying the body to the cathode is the usual choice on modern magnetically
   shielded thrusters to limit pole erosion.
3. The front-plate potential shapes the near-plume field, and so the beam divergence.

Answering this needs a model with the iron surfaces in the plasma and a current-balance condition. The
2-D PIC below has neither (§4). The 3-D or an embedded-boundary PIC model could, with the "spacecraft charging"
technique: WarpX's example adjusts a conductor's potential each step from the charge it collects.

### 3-D vacuum E-field (`estatics/run_3d.py`, `mesh_3d.py`)

The same `estatic.pro` with `Flag_Axi = 0`. The geometry has the real 6 in square plates (with their bores), the
four outer posts with one coil each, and the 2-D profiles revolved for the axisymmetric parts. The chamber wall
is a 0.5 m sphere, and the mesh has 1.55 M tetrahedra (1 s to mesh, 47 s per solve on 4 threads).

The **cathode is a placeholder**: a 6.35 mm radius cylinder parallel to the thruster axis, at x = 110 mm, from
z = -40 mm to its orifice at z = +10 mm.

Capacitance matrix in pF (chamber wall the reference; symmetric to machine precision):

| | anode | iron | cathode |
|---|---|---|---|
| anode | 17.799 | -17.791 | -0.001 |
| iron | -17.791 | 28.887 | -1.323 |
| cathode | -0.001 | -1.323 | 2.009 |

- **Iron floating, cathode at the chamber wall's 0 V:** 184.8 V. Superposition (k = 0.6159) and the direct
  floating solve agree to 1e-12 V. The 2-D model with the same 0.5 m chamber gives 192.2 V. The difference comes
  from the cathode's 1.3 pF next to the iron and from the real plate and post shapes.
- **Cathode 20 V below the chamber wall, anode 300 V above the cathode:** iron at 171.5 V.
- Figure `estatics/out/3d/potential_3d.png`: the potential on the plane through a post and the cathode, and on
  the plane between two posts.
- The anode still couples almost only to the iron (17.80 pF against 17.54 pF in 2-D), so the 2-D model is a
  good stand-in for the channel field.

## 3. Magnet field for the PIC (`magnetics/`)

| # | assumption | value |
|---|---|---|
| M1 | The het 2-D magnetostatics model | linear iron, mu_r 1500; posts smeared into a ring of mu_r 226; outer coils as two opposite current sheets; outer current = inner x 0.635 x 230/184 |
| M2 | Design point | channel-profile B_r peak 165 G. Gives inner coil 0.9271 A, outer 0.7358 A (het parity: 0.9273 / 0.7361 A); peak at z = -2.25 mm, 95.7 % at the exit, 2.76 % at the anode, as the parity table |
| M3 | Static, vacuum B | the plasma's own field is neglected: the Hall current is small next to the coil ampere-turns |

## 4. 2-D PIC (`pic2d/het2d.py`): assumptions and placeholders

WarpX 26.10 (conda-forge, CPU, OpenMP), PICMI Python input. **Planar x-z geometry**, all three velocity
components kept.

The channel annulus is unrolled into a slab: x = r, z axial, y = azimuthal. The y extent is the circumference at
the channel mid-radius (2 pi x 26.78 mm = 168.3 mm), which reproduces the channel area exactly: 2 pi r_mid w =
pi (r_out^2 - r_in^2). Reported currents are totals over that circumference. The 1/r geometric factors are lost,
which is acceptable since w / r_mid = 0.35.

I started axisymmetric (r-z) but switched: WarpX 26.10 leaves the radial Poisson boundary conditions of an r-z
grid unset when the inner radius is above 0 (I5). `--geometry rz` is kept for when that is fixed.

| # | item | placeholder used | why / what it should become |
|---|---|---|---|
| P1 | **Domain** | the channel annulus (x = r from 22.12 to 31.4452 mm), from the anode face (z = -41.43 mm) to a **cathode plane** 30 mm downstream of the exit. 32 x 240 cells (0.29 x 0.30 mm) | WarpX cannot model the BN walls as dielectrics (I4), so the channel walls are the domain's radial boundaries. The plume is therefore confined to the channel radii |
| P2 | **Walls** (both radial boundaries, channel and plume) | field normal derivative zero ("Neumann"). Ions absorbed. **Electrons reflected when their wall-normal energy is below 30 eV, absorbed above**: a fixed-cutoff stand-in for the sheath of a floating BN wall. No secondary electrons | Fully absorbing walls (no charging) sent the 6 A of cathode electrons straight to the walls along the plume field lines (I10). A real floating wall's sheath is a few Te. Needs: wall charging (a dielectric solver, I4), BN secondary emission (yield about 1 near 40-50 eV) |
| P3 | Electrodes | anode face at 300 V over the full channel width; cathode plane at 0 V | real anode shape and cathode coupling |
| P4 | **Neutrals: fixed background**, no depletion | 6.81e19 m^-3 in the channel: 2 mg/s as a free-molecular half-Maxwellian at 600 K through the 1569 mm^2 channel area. Downstream of the exit, n_ch (2/pi) atan(w/2z) (view factor of a slot of the channel width) + 3e17 m^-3 facility background | kinetic neutrals (WarpX DSMC supports electron-impact ionization of a neutral species), which also gives propellant depletion. This placeholder is the main reason the discharge current is not self-limited (see results) |
| P5 | Collisions | warpx-data argon cross sections. Electrons: elastic (Zatsarinny-Bartschat), one lumped excitation at 11.5 eV (Yamabe-Buckman-Phelps), ionization at 15.76 eV (Biagi/MagBoltz). Ions: isotropic + backward elastic (Phelps). Collisions every 2 steps for electrons (nu dt <= 0.02), every 20 for ions | more excitation levels; charge exchange as its own process |
| P6 | **Anomalous electron transport** | Bohm-type collision frequency nu = alpha omega_ce: alpha = 1/160 in the channel, 1/16 in the plume, switching over 1 mm at the exit. \|B\| from a fit of the mid-radius profile (1.8 % rms), uniform in r. Implemented as a fictitious elastic collision with sigma ∝ 1/v (so nu does not depend on energy) and a 1e-15 kg target, which turns electrons without taking energy | the coefficients are the usual starting values. Calibrate against the measured discharge current (and later thrust); 2-D r-z PIC needs this because the azimuthal turbulence that carries the current is not resolved |
| P7 | **Permittivity scaling** | eps0 x 400 (Debye length and plasma period x 20). Done by dividing the deposited charge by 400 just before every Poisson solve: WarpX's eps0 is a compile-time constant. Densities stay physical, so all collision rates are real. Checked active: A/B run, 5 V vs 107 V at the exit after 15 ns | standard in HET PIC; check the sensitivity (200 / 400 / 800) once the physics is settled |
| P8 | **Cathode electrons** | 6 A one-way thermal flux at 5 eV from the cathode plane. More than needed: the excess returns to the boundary through a small potential dip (space-charge limited), so the net current sets itself | a current-controlled injection (inject what the beam takes out) is the other common choice |
| P9 | Start | uniform plasma 5e16 m^-3, electrons 10 eV, ions 600 K: the initial condition of the LANDMARK benchmark (Charoy et al. 2019). Leveling-thinning resampling keeps each species below 80 macroparticles per cell on average | starting from a 1e15 m^-3 seed instead, the cathode electrons piled up at the magnetic barrier and only about 2 mA crossed in 1.3 µs (I7) |
| P10 | Numerics | dt 5e-11 s. At 1e18 m^-3, scaled: omega_pe dt = 0.14 and Debye length 0.66 mm at 20 eV, against 0.3 mm cells. omega_ce dt = 0.14 at 160 G. Linear particle shape, 12 OpenMP threads | |

## PIC results

(filled in after the example run; see below)

## 5. Issues encountered

| # | issue | what I did |
|---|---|---|
| I1 | GetDP writes its output files next to the `.pro` file, not in the working directory | each solve runs on a copy of the `.pro` in its own folder |
| I2 | GetDP's axisymmetric Jacobian integrates per radian, so global charges came out as Q / 2 pi | `estatic.pro` scales charges and energies by 2 pi (and divides an imposed charge); the sphere check confirms totals |
| I3 | A gmsh Threshold size field without `StopAtDistMax` forced 0.4 mm elements over the whole 1 m domain: 9.5 GB and swapping | fixed (`StopAtDistMax = 1`); the mesh now takes 1 s |
| I4 | **WarpX's electrostatic solver has no dielectric.** The `MacroscopicProperty` docstring says it works with the electrostatic solver, but the source shows the solve uses a constant eps0. Embedded boundaries are always Dirichlet (fixed potential) | the BN walls are the domain boundaries with a zero-normal-field condition (P1, P2). Real BN walls need one of: a Python Poisson solver through WarpX's `poissonsolver` hook; a small WarpX patch passing eps_r as the multigrid coefficient (the effective-potential solver already uses a variable one); or a different code |
| I5 | **WarpX bug: r-z electrostatics with an inner radius above 0.** `PoissonBoundaryHandler::DefinePhiBCs` (Source/FieldSolver/ElectrostaticSolvers/PoissonBoundaryHandler.cpp) sets the radial boundary conditions only when r_lo = 0, and its loop then starts at z. With r_lo > 0, `lobc[0]`, `hibc[0]` and `dirichlet_flag[0..1]` are never initialised. **Symptoms:** (1) multigrid diverges whenever the grid can coarsen (32 x 240; residual x10^20, abort); (2) on a non-coarsenable grid (31 x 238) it converges but the radial walls act as about 0 V. That gave 78 V at 5 mm from the anode and 0.25 V at the exit, where Neumann walls should give about 280 V and 126 V | planar x-z geometry (§4), where all four boundaries are set: 284 / 138 / 63 V, the expected linear vacuum profile. Worth an upstream issue with this minimal case: Neumann radial boundaries at r_lo > 0 in RZ electrostatics |
| I6 | WarpX's `MCCCollisions` does not accept constants for its density expression | numbers are written into the expressions |
| I7 | Ignition from a weak seed is slow. From 1e15 m^-3, the cathode electrons collected in a thin layer at the B peak (z = -3 to 0 mm) and along the inner boundary of the plume. About 2 mA reached the anode after 1.3 µs. **Correction:** I first blamed a dense initial plasma for holding the 300 V in an anode sheath; that symptom was the r-z boundary bug (I5), as the identical pattern at 1e15 and 3e16 m^-3 shows | restarted from the LANDMARK initial plasma (P9). The seed run is kept in `pic2d/runs/ignition_attempt` |
| I8 | Hyper-threads hurt: 24 threads ran at 35 ms/step against 24 ms/step on 12 | 12 threads |
| I9 | bplsim open issue O1 (3-D nonlinear magnet solve) finished overnight: the final 100 % ramp step converges, but the 85 % step does not | recorded with options in `bplsim/PROGRESS.md`; not part of this work |
| I10 | With fully absorbing radial walls, 5.7 A of the 6 A of injected cathode electrons went to the radial walls (4.9 A inner, 0.9 A outer) along the plume field lines. Nothing reached the channel: such a wall never charges, so nothing turns electrons back | the electron wall cutoff (P2). It is also why electrons collect along the inner boundary in the plume, where the real field lines run on to the pole cap. A real plume region (embedded boundaries, item 2 of the next steps) removes this artefact |

## Next steps

1. **Neutral depletion** (P4): kinetic argon neutrals with WarpX's DSMC ionization. This makes the 2 mg/s a real
   boundary condition and limits the discharge current.
2. **Real geometry** with embedded boundaries: the pole cap, the front plate and an expanding plume. This needs
   dielectric walls (I4).
3. **BN secondary electron emission** at the walls.
4. **Calibrate the anomalous coefficients** (P6) against a measured discharge current.
5. **3-D:** the vacuum 3-D model runs (`run_3d.py`). It needs the real cathode position and mount.
   3-D PIC means WarpX 3-D with the CAD as an STL embedded boundary, on the GPU (the RTX 5080). The conda WarpX
   build is CPU-only, so this needs a CUDA build.
6. **Floating vs tied iron in operation:** an embedded-boundary model with the iron as a conductor whose
   potential follows its collected current.
