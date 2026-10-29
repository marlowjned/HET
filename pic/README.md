# HET field and PIC models (standalone, no bplsim)

Simulations of the BPL-1000 Design B (shieldless, coils wound forward) beyond the magnet: the vacuum electric
field, the magnet field map, and a first 2-D particle-in-cell discharge model. Everything here runs without bplsim
and can be moved into it later. Assumptions, placeholders and issues: [`REPORT.md`](REPORT.md).

| folder | what | status |
|---|---|---|
| `geometry/design_b.py` | Design B r-z geometry shared by all models (values from the het library), plus placeholders for the anode and cathode | done |
| `estatics/` | vacuum E-field at 300 V (gmsh + GetDP): floating or grounded iron, superposition bases, grid export | done, verified against an exact solution |
| `magnetics/` | Design B B-field at the 165 G design point (gmsh + GetDP), openPMD map for WarpX | done, matches the het parity case to 0.03 % |
| `pic2d/` | 2-D PIC-MCC discharge model (WarpX, planar x-z): argon, 2 mg/s, 300 V | runs end to end; the discharge is not yet self-sustaining (REPORT.md) |

## Environment (WSL, no sudo)

```bash
bash pic/setup_env.sh      # micromamba + conda-forge env "hetpic": WarpX 26.10, gmsh, openPMD-api, numpy, ...
source pic/activate.sh     # every shell; also puts GetDP (~/tools/getdp) on PATH
```

The PIC model also needs the argon cross sections from
[warpx-data](https://github.com/BLAST-WarpX/warpx-data), cloned to `~/tools/warpx-data`:
`git clone --depth 1 https://github.com/BLAST-WarpX/warpx-data ~/tools/warpx-data`.

## Running

```bash
cd pic/estatics && python run.py                 # 40 s
cd pic/magnetics && python run_mag.py            # 15 s; writes the B map the PIC reads
cd pic/pic2d && python het2d.py --tag example && python analyse.py runs/example     # 6 us, about 20 min
```
