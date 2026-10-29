# het: the HET project for bplsim

The BPL-1000 Hall thruster's design library and trade studies, run with
[bplsim](../../bplsim) (`Projects/bplsim`). This folder replaces the per-study scripts in
`gmsh_getdp/axisym/`; those stay until every study is ported.

## Setup (Windows, PowerShell or Git Bash)

```
cd ../../bplsim
uv venv --python 3.12 .venv
uv pip install --python .venv/Scripts/python.exe -e ".[dev]"
```

GetDP: either put `tools/getdp-3.5.0-Windows64` on PATH, set `BPLSIM_TOOL_GETDP` to `getdp.exe`, or add a
site profile `~/.config/bplsim/sites/local.yaml`:

```yaml
name: local
data_root: ~/bplsim_data
executor: local
tools: {getdp: C:/Users/<you>/Projects/HET/tools/getdp-3.5.0-Windows64/getdp.exe}
```

## Layout

| Folder | Contents |
|---|---|
| `library/materials` | iron (legacy proxy and real grades), BN, copper, AWG 20 wire |
| `library/part_types` | channel insulator, plates, core, pole cap, screens, post array, solenoid, post coils |
| `library/parts` | real parts: BN channel v2, 6 in plates, caps, the spec coils |
| `library/assemblies` | `legacy_bn` (as the legacy 2-D model), `design_a`, `design_b`, `bpl1000_bn` (illustrative) |
| `library/probes`, `metrics`, `references`, `figures` | channel profile/ring/wall, field shape, iron saturation, coil power, half-Gaussian target, profile figures |
| `library/measurements`, `pipelines` | what to measure; GetDP 2-D with the legacy mesh settings |
| `studies` | parity checks, the gap sweep, the Design A search |
| `reports` | report templates (`*.md.j2`) and their generated output |
| `tests` | parity with the legacy code (`PARITY.md`) |

## Common commands (from this folder)

```
bplsim check                                         # validate the library
bplsim show het:assemblies/design_a                  # every resolved value and where it came from
bplsim study plan studies/design_a_gap_sweep.yaml    # dry run
bplsim study run studies/design_a_gap_sweep.yaml     # run (or resume)
bplsim report build reports/design_a_gap_sweep.md.j2 # regenerate the report from the results
bplsim query "params.inner_screen_gap >= 20 mm" --metrics misfit.rms_pct,power.total
python -m pytest tests -m slow                       # parity with the legacy 2-D code
```
