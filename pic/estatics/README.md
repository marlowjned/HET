# Vacuum electrostatic field, Design B (2-D axisymmetric)

The field of the thruster with no plasma: anode at the 300 V discharge voltage, chamber wall at 0 V, and the iron
(with the coils) either floating or grounded. Standalone gmsh and GetDP, no bplsim. Assumptions and open questions
are in [`../REPORT.md`](../REPORT.md).

```bash
source pic/activate.sh            # once per shell (pic/setup_env.sh builds the environment)
cd pic/estatics
python verify_spheres.py          # formulation check against concentric spheres (exact solution), 2 s
python run.py                     # Design B: bases, floating/grounded cases, figures, grid export, 40 s
python run.py --quick             # without the chamber-size and mesh-convergence studies, 9 s
```

| file | content |
|---|---|
| `estatic.pro` | GetDP formulation: scalar potential with one unknown per conductor (`BF_GroupOfNodes`); each conductor at a set potential or floating (set net charge). 2-D axisymmetric or 3-D |
| `mesh_rz.py` | gmsh mesh of the r-z half plane: vacuum + BN, conductor interiors cut out, grounded half-circle chamber wall |
| `../geometry/design_b.py` | the Design B dimensions (from the het library), shared with the PIC model |
| `run.py` | solves, superposition check, figures, `out/fields_rz.npz`, `out/summary.json` |
| `fields.py` | reads `fields_rz.npz` and returns any case by superposition |
| `verify_spheres.py` | exact-solution check of `estatic.pro` |

## Method

Two solves per volt give every case (the chamber wall is at 0 V in both):

- `phi_A`: anode 1 V, iron 0 V
- `phi_I`: iron 1 V, anode 0 V

The iron-grounded case is `300 phi_A`. The iron-floating case is `300 (phi_A + k phi_I)`, where
`k = -C_IA / C_II` is set by zero net charge on the iron. A direct floating solve agrees with the superposition to
5e-12 V. In operation the plasma adds its own space-charge part.

## Results (chamber wall radius 1 m, BN eps_r 4)

| quantity | iron grounded | iron floating |
|---|---|---|
| iron potential | 0 V | **198.0 V** |
| potential 1 mm from the anode face, channel mid-radius | 271.0 V | 290.1 V |
| potential mid-channel (z = -20 mm) | 27.8 V | 207.3 V |
| potential at the exit plane, mid-radius | 1.9 V | 194.3 V |
| E_z 1 mm from the anode face | 28.5 kV/m | 9.7 kV/m |

Capacitance matrix in pF (anode, iron; chamber wall the reference):

| | anode | iron |
|---|---|---|
| anode | 17.544 | -17.543 |
| iron | -17.543 | 26.574 |

The anode couples almost only to the iron. Its capacitance to the chamber is C_AA + C_AI = 0.001 pF. Floating iron
therefore sits at C_IA/C_II = 66 % of the anode voltage.

- **Grounded iron:** the vacuum potential dies away inside the channel, with an e-folding length of about 7 mm.
  The iron behind the thin BN walls screens it.
- **Floating iron:** the whole thruster body sits near 200 V. The potential then falls off slowly into the plume,
  toward the chamber wall.

**Checks:**

- **Chamber size:** the floating potential is 192.2 / 198.0 / 200.7 / 202.0 V for a chamber radius of
  0.5 / 1 / 2 / 4 m.
- **Mesh refinement** (half the element size, 352k triangles): the floating potential changes by -0.05 %, the
  capacitances by -0.2 %, and the probe potentials and fields by at most 1.5 %.

The 198 V applies only without plasma, for example at start-up. Once the discharge runs, the plasma sets the
floating potential by current balance; see `../REPORT.md`.

Figures, all in `out/`:

- `geometry.png`: the parts as modelled, and the mesh at the pole-cap gap
- `potential.png`, `potential_channel.png`: potential maps
- `efield.png`: |E| maps
- `lines.png`: profiles along the channel and the floating-minus-grounded map

`out/fields_rz.npz` holds the per-volt bases on a 0.25 mm grid, r 0-120 mm and z -60 to +150 mm, with a region
map.
