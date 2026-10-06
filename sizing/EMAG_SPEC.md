# BPL-1000 electromagnet specification

Final design, 2026-10-05. Numbers come from `sizing/em_final_spec.py`, which
reads the combined linear solve of exactly this geometry
(`gmsh_getdp/assembly/sweep_runs_bpl1000/ic0.750_it0.375_ip0.375_oc0.750_ot0.312_h1.625`,
frozen Onshape copy, microversion `d75ed0ea…`). Supersedes the BPL-700 design
point in `WIRE_AND_SUPPLY.md`.

## Specification

Common to all five magnets:
- **Wire:** AWG 20 heavy-build aromatic polyimide, NEMA MW 16-C, 240 °C. Max OD is 0.892 mm. All fits use the maximum OD, not the nominal.
- **Cores:** 0.75 in diameter, ASTM A108 (1018 cold-drawn). All five come from one bar stock.
- **Ground insulation:** 0.15 mm on each core.
- **Winding:** square packing, 46 turns per layer over the 1.625 in winding height.
- **Supply:** two constant-current channels. The four outer coils are in series on one channel.

| | Inner (×1) | Outer (×4, in series) |
|---|---|---|
| Core diameter | 0.750 in | 0.750 in |
| Winding | 5 layers × 46 = **230 turns** | 4 layers × 46 = **184 turns** per coil |
| Wound copper OD | 1.113 in | 1.043 in |
| CAD coil OD (window) | 1.125 in | 1.0625 in |
| Fit limit | 1.385 in (chamber bore) | 1.720 in (chamber outer wall) |
| Required ampere-turns (300 G at exit) | 352.0 | 223.4 per coil |
| **Channel current** | **1.53 A** | **1.21 A** |
| Resistance at 20 °C / hot | 0.557 / 0.643 Ω | 0.429 / 0.494 Ω per coil |
| Channel voltage (hot, coils only) | 0.98 V | 2.40 V (4 in series) |
| Power | 1.51 W | 0.73 W per coil, 2.91 W total |
| Coil temperature (steady state) | 59.6 °C | 58.9 °C |
| Wire per coil | 17.2 m (80 g) | 13.3 m (62 g) |
| Peak core flux | 0.66 T | 0.32 T |

**Totals:**
- 4.42 W, with the iron at 58.4 °C steady state. The thermal model is radiation-only at envelope emissivity 0.30, which is its dominant assumption.
- 70.2 m of AWG 20 wire, about 330 g. A 1 lb spool (~98 m) is enough; 2 lb leaves room for practice winds.

**Field:**
- Quality score 86.9 at the committed outer/inner ratio (214.5/338 = 0.635).
- Re-optimising the ratio for this geometry gives 0.51, which needs 382/195 ampere-turns and scores 87.0. It isn't worth a separate design. With two independent channels, the ratio can be set on the bench.

**Supply:** both channels sit in the requested 1–2 A band, at under 3 V plus lead and feedthrough drop. The two-channel 30 V / 3 A bench supply in `WIRE_AND_SUPPLY.md` covers this with wide margin.

### Onshape configuration to set

| Parameter | Value (in) |
|---|---|
| emag_height | 1.625 |
| iem_coreod | 0.75 |
| iem_coilt | 0.375 |
| iem_platet | 0.375 |
| oem_coreod | 0.75 |
| oem_coilt | 0.3125 |

`coilt` is a diametral increment: coil OD = coreod + coilt. Each window leaves clearance outside the wound copper: 0.006 in radial on the inner coil and 0.010 in on the outer.

## Why the exact gauge and insulation don't matter much

The FEM only sees ampere-turns spread over the coil window. Gauge and insulation therefore don't change the field at all. They change only how many turns fit, and from that the current, voltage, power and wound diameter.

1. **Insulation thickness barely moves the diameter.** Holding the turn count at spec (230 / 184):

   | Film | Inner wound OD | Outer wound OD |
   |---|---|---|
   | Bare copper (zero film, physical lower bound) | 1.082 in | 1.018 in |
   | Heavy build, max OD (spec) | 1.113 in | 1.043 in |
   | Double heavy film | 1.221 in | 1.144 in |

   From no insulation at all to twice the specified film, the diameter moves by −0.03 / +0.11 in. Every case fits well inside the limits.

2. **Gauge barely moves the power.** In the spec windows, every gauge from AWG 16 to 28 gives 4.4–5.2 W total, with the iron at 58–63 °C. What gauge really sets is the current:

   | AWG | Inner current | Outer current | Total power |
   |---|---|---|---|
   | 16 | 4.05 A | 3.85 A | 5.19 W |
   | 18 | 2.38 A | 2.01 A | 4.50 W |
   | **20** | **1.53 A** | **1.21 A** | **4.42 W** |
   | 22 | 1.03 A | 0.78 A | 4.64 W |
   | 24 | 0.71 A | 0.52 A | 4.96 W |

   The scatter comes from whole-layer quantisation, not the gauge itself. AWG 20 is the only gauge that puts both channels in 1–2 A.

3. **Diameter vs power has diminishing returns.** Adding layers past the spec saves little:

   | Layers | Inner OD | Total power | Outer OD | Total power |
   |---|---|---|---|---|
   | 2 | 0.90 in | 6.5 W | 0.90 in | 7.3 W |
   | 3 | 0.97 in | 5.3 W | 0.97 in | 5.3 W |
   | 4 | 1.04 in | 4.8 W | **1.04 in** | **4.4 W** |
   | 5 | **1.11 in** | **4.4 W** | 1.11 in | 3.9 W |
   | 6 | 1.18 in | 4.2 W | 1.18 in | 3.5 W |

   Bold rows are the spec. Going smaller than this costs power and current fast, while going larger saves under 0.5 W per added layer.

Figure: `gmsh_getdp/assembly/em_sizing_figures/final_spec_sensitivity.png`.
Full sweep data: `gmsh_getdp/assembly/em_sizing_results.csv` and `sizing/em_winding_trade.csv`.

## Known biases and open items

- **Linear iron (μr 1500).** Peak core flux is 0.66 T, well below the knee. On BPL-700, linear iron read +1.2% against a nonlinear solve; this geometry hasn't been checked nonlinearly. That bias, the A36 B-H proxy, and the FEM spreading current over the whole window (−0.8% on BPL-700) all mean the real hardware will want slightly more current than modelled. The two-channel supply absorbs this.
- **Packing.** The design assumes square packing. Hand winding usually lands between square and hex packing and may fit more turns. Count what you actually wind, then set each channel to required ampere-turns ÷ turns.
- **Non-magnetic bodies.** `Revolve1` (the gas-feed fittings) and anything else the pipeline didn't recognise are treated as non-magnetic, at the user's instruction.
- **Arcing and isolation.** The film handles turn-to-turn insulation. Size the coil-to-core ground insulation, the over-wrap and the lead insulation to the potential difference set by how the magnet supply is referenced, not to the coils' few volts. The discharge voltage and grounding scheme are still to be decided.
