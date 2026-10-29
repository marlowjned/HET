# Parity of bplsim with the legacy 2-D code

bplsim's `getdp_axisym` backend with the `het` library (this folder) against the legacy
`gmsh_getdp/axisym` code, for the two BN designs, at 165 G peak with q = 0.635 (2026-10-29).

Legacy numbers were produced by the legacy code itself (`AXISYM_CONFIG=bn`, `target_sigma.score`), not
copied from reports. Design B is `results/bn_noshield/best.json` without the post shim.

| Quantity | Design A legacy | Design A bplsim | Design B legacy | Design B bplsim |
|---|---|---|---|---|
| misfit RMS (% of peak) | 1.1727 | 1.1782 | 1.9355 | 1.9435 |
| misfit max (% of peak) | 5.0450 | 5.0716 | 4.3850 | 4.3702 |
| sigma (mm) | 15.2 | 15.2 | 15.4 | 15.4 |
| B_r at anode (% of peak) | 1.7348 | 1.7051 | 2.7592 | 2.7354 |
| \|B\| at anode (% of peak) | 3.3741 | 3.3709 | 3.2436 | 3.2236 |
| B_r at exit (% of peak) | 94.955 | 94.928 | 95.615 | 95.630 |
| peak position (mm) | -2.75 | -2.75 | -2.25 | -2.75 |
| inner current (A) | 0.91300 | 0.91342 | 0.92680 | 0.92733 |
| outer current (A) | 0.72469 | 0.72502 | 0.73565 | 0.73607 |
| power (W) | 1.5737 | 1.5755 | 1.8641 (fixed turns) | 1.8668 |
| saturation ratio | 0.2881 | 0.2930 | 0.2271 | 0.2235 |

**Tolerances** (`tests/test_parity.py`) come from the legacy model's own mesh sensitivity: re-meshing
Design A at 0.8x and 1.2x the mesh size moves the misfit by +-0.3 %, B_r at the anode by +-0.3 % and the
currents by 0.05 %. The test allows 0.1 percentage points on %-quantities, 1 % on physical quantities,
one 0.5 mm grid step on the peak position, and 3 % on the saturation ratio.

**Known, explained differences**

- *Mesh.* bplsim meshes each part as its own region, so the triangulation differs slightly from the legacy
  mesh (82 k vs 82 k elements, different layout). This is the source of the 0.05-0.5 % differences above.
- *Peak position, Design B.* The peak is flat to within the mesh noise over 0.5 mm; both values are grid
  points 0.5 mm apart.
- *Power, Design B.* The legacy `P165_W` scales coil resistance with the inverse window area (it assumed
  the turns follow the window). bplsim uses the fixed-turns formula (I^2 R from the actual turns, mean turn
  length, wire area and temperature), which matches the legacy `power_fixed_turns` to 0.15 % at 50 °C.
  Design A uses the spec windows, so both legacy formulas agree there.
- *Saturation.* The legacy model treats the core and pole cap as one iron region for the 95th percentile;
  bplsim has them as separate parts (the core starts behind the cap).
