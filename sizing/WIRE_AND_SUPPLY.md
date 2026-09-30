# Magnet wire and power supply specification

What wire to buy, how many turns of it fit, what currents that implies, and
what supply drives it. Written 2026-09-30, when the previous design point's
wire turned out not to be purchasable.

Computed by `winding_design.py` (run it for the full gauge trade); wire
dimensions live in `wire_data.py`; material provenance and confidence
ratings in [`../MATERIALS.md`](../MATERIALS.md).

---

## Specified wire

**AWG 20 heavy-build aromatic polyimide magnet wire, NEMA MW 16-C, 240 °C.**

| | |
|---|---|
| bare copper | 0.813 mm (0.0320 in) nominal |
| max overall dia | **0.892 mm (0.0351 in)** — NEMA MW 1000 maximum, not nominal |
| temperature class | 240 °C |
| availability | stocked, 2 oz to 10 lb spools |
| example source | Remington Industries, 20 AWG 240 °C polyimide |

The number that matters is the **maximum** overall diameter. A fit
calculation wants the worst case a supplier is allowed to ship, and using
nominal instead is how you lose a layer on the bench.

## Why not fiberglass

The project previously specified "single-glass-served, 200 °C, AWG 20." Three
things are wrong with that, and they compound:

1. **It isn't a standard size.** NEMA MW 1000's round served-wire dimension
   table (polyester-glass, i.e. "Dacron glass") runs 1/0 → AWG 12 and no
   finer. Glass serving is built from filaments whose thickness doesn't
   scale with the conductor — it adds ~0.15–0.21 mm to the OD at any gauge —
   which is exactly why it stops being offered on fine wire.
2. **It isn't stocked.** No US retailer carries round glass-served magnet
   wire in any gauge. S&W Wire and Bridgeport Insulated Wire both make it,
   to order. Coarsening to AWG 12 gets you into the standard's table but
   still not onto anyone's shelf.
3. **The round product is Class 155, not 200 °C.** Dacron Glass 155 is
   NEMA MW 45. The 200 °C variant (Dacron Glass High Temp, MW 51) is made
   "only in shaped or heavy round wire sizes." The old 200 °C figure in the
   thermal model was therefore optimistic for the wire it named — and that
   entry used to be the *first* limit reached.

Polyimide is better on every axis that matters here: it's **thinner** than
the glass estimate (so more turns, not fewer), its dimensions come from a
standard rather than an estimate, it's **240 °C**, and being a solid film it
has no porous covering to trap gas — which matters in a vacuum chamber.

The original rationale for glass was plasma arc resistance. The coils sit
outside the discharge chamber with the PTFE wall between them and the
plasma, so arc resistance at the wire is a second-order concern; vacuum
compatibility and temperature class are not.

**Fiberglass remains specified as a fallback** at AWG 12 (the finest
standard round size) in `wire_data.py`, so it can be ordered if a
requirement ever demands it. It needs a made-to-order run, and it pushes the
inner channel to ~7.3 A.

## The coil windows, and what eats them

Baseline geometry (`emag_height` 2.000 in, inner core 1.000 in):

| | inner | outer (×4) |
|---|---|---|
| bore radius (bare core) | 12.70 mm | 12.70 mm |
| ground insulation allowance | 0.15 mm | 0.15 mm |
| coil OD / 2 | 17.46 mm | 19.05 mm |
| **usable build** | **4.612 mm** | **6.200 mm** |
| height | 50.8 mm | 50.8 mm |

The inner winding **cannot grow outward** — the chamber bore sits on it. It
grows only inward, by thinning the core, which is what the `emag_height`
sweep trades against core saturation. That makes the inner window the
binding one, and it is what a wire diameter has to be chosen against.

The 0.15 mm ground insulation is two half-lapped wraps of 0.05 mm polyimide
tape over the bare steel. It is an **allowance, not a datasheet value**, and
it is broken out rather than buried because it costs ~3% of the inner build —
the same order as the OD tolerance sitting next to it.

## Turn counts, both packing models

`winding_design.py` reports every result under two packing models, because
which one is real depends on how the coil is wound:

- **hex** — each layer settles into the valleys of the one below, radial
  pitch `d·√3/2`. Best case; needs a smooth, dimensionally consistent wire.
- **square** — layers stack directly, pitch `d`. Conservative, and what a
  hand-wound coil tends to do.

AWG 20 polyimide, baseline geometry:

| | inner | outer |
|---|---|---|
| turns, hex | 280 (5 × 56) | 392 (7 × 56) |
| turns, square | **280 (5 × 56)** | **336 (6 × 56)** |
| unused build, square | 0.152 mm | 0.848 mm |

**The inner coil holds 5 layers either way.** That is the single best reason
to prefer this wire: the turn count no longer depends on the packing
assumption at all. AWG 16 and 18 share that property; the old glass-served
AWG 20 did not — it needed hex nesting for its 5th layer and missed the
square-packed one by **0.037 mm**, a 20% swing in turns (260 → 208) resting
on an assumption about a wire type least likely to nest.

Run `python winding_design.py` for the margin column across AWG 12–24. Some
gauges are genuinely knife-edge (AWG 15 and 17 square-pack to within
0.004 mm and 0.000 mm of losing a layer); AWG 20's 0.152 mm is 17% of a wire
diameter, which is the kind of slack a build tolerance can absorb.

## Design point

Two independent constant-current channels, **both coils wound full**:

| | inner | outer (×4, series on one channel) |
|---|---|---|
| turns | 280 | 336 (square) / 392 (hex) |
| current | **1.207 A** | **0.638 A** (square) / 0.547 A (hex) |
| ampere-turns | 338.0 | 214.5 |
| resistance, 20 °C | 0.859 Ω | 4.243 Ω (all four) |
| voltage, 20 °C | 1.04 V | 2.71 V |
| dissipation | 1.25 W | 1.73 W |

Total **2.98 W at 20 °C**, 3.33 W hot. Steady-state temperature 50 °C, with
the first material limit (wire insulation, 240 °C) at 4.17 A — **3.5×**
design current. See [`../thermal/README.md`](../thermal/README.md).

Ampere-turns are unchanged from the previous design point (338.0 / 214.5),
so **every field already solved in this repo remains valid**. Only the
current/turns split changed.

> **Wind first, then set the current.** With two channels the turn count is
> a current setting, not a design error: count what you actually wound and
> set each channel to `req_At / N`. Nothing upstream has to be re-derived.

## Two channels vs one series string

| | series string | two channels |
|---|---|---|
| outer turns | 178 of 336 (53% of window) | 336 (full) |
| current | 1.207 A through everything | 1.207 A / 0.638 A |
| dissipation, 20 °C | 4.34 W | **2.98 W** |
| A-turn ratio | frozen in copper at winding time | a dial |

Two channels win on three counts:

1. **1.46× less dissipation** (1.57× under hex packing), because `P ∝ 1/N`
   and both coils get wound full instead of the outer being part-wound to
   hit a ratio.
2. **The 1.58 A-turn ratio stays adjustable.** That ratio came from a linear
   solve with a proxy B-H curve; freezing it in copper on a first article is
   the wrong thing to make permanent.
3. **It fixes an FEM assumption.** `build_assembly.py` applies a uniform
   current density over the whole coil solid. With the outer coils
   part-wound to 45% of their window that was 1.84× off. Wound full, the
   remaining error is the un-fillable remainder of the window, measured by
   `emag_sweep.py --fill-check` at **−0.80%** on the channel-exit field —
   see below.

### What the fill-factor error actually costs

`emag_sweep.py --fill-check` re-solves the baseline with the outer coil
solid shrunk to the copper's real outer radius (1.433 in instead of the CAD
window's 1.500 in): same ampere-turns, same iron, current confined to where
it really sits.

| | exit B_r |
|---|---|
| full CAD window (what the FEM does) | 306.52 G |
| real copper extent | 304.08 G |
| | **−0.80%** |

So the A-turn requirement is ~0.8% higher than the sweep reports, i.e. the
inner channel wants ~1.217 A rather than 1.207 A. Left uncorrected in the
committed design point, because it is under the 1.25% linear-vs-nonlinear
iron discrepancy already accepted — but noted as a **systematic bias in a
known direction**, and it does not cancel against the other one: the
A36-proxy B-H curve is also optimistic (see `../MATERIALS.md`). Both say the
real hardware needs slightly more current than modelled.

Only the outer coil is varied, because `outer_coil_od` touches nothing but
the outer coil solids, whereas `inner_coil_od` also drives the chamber and
spacer bore and would change the iron. The inner coil needs no test: its
copper is pushed 0.15 mm *out* at the bore by the ground insulation and sits
0.152 mm *in* at the top, so its mean radius moves by ~1 µm.

A pre-solve estimate — extrapolating from the existing
`outer_coil_od` 1.500→1.625 sensitivity pair — put this at 0.17%, about 4.6×
too small. The sensitivity is not linear in mean coil radius. Worth
remembering before trusting a similar extrapolation.

## Power supply

### Requirements, independent of gauge

- **Two independent constant-current channels.** The field follows current,
  not voltage, and copper resistance rises ~0.4%/K as the coils warm, so a
  constant-voltage supply would let the field drift down by the same amount.
- **Floating (isolated) outputs.** The coils sit next to the discharge; an
  isolated output keeps the magnet circuit independent of the discharge
  supply's reference. Verify channel-to-channel and channel-to-earth
  isolation on the datasheet — don't assume it.
- **Flyback diode across each coil string.** Stored field energy at the
  design point is 65.9 mJ (set by the field, so the same at every gauge), so
  an open circuit or a relay opening under current would otherwise generate
  a voltage spike. Energy this small makes diode selection trivial — any
  1 A / 100 V part covers it. What the inductance actually constrains is
  ramp rate, not the diode.
- **Ramp, don't step.** `L/R` is of order 100 ms, so ramping over ≳1 s makes
  the turn-on quasi-static. Ripple is heavily filtered by the same time
  constant and does not matter.

### Envelope at the specified AWG 20

Maximum over **all 25 swept geometries × both packing models**, so a supply
chosen here covers whichever spacer height is picked:

| | inner channel | outer channel |
|---|---|---|
| max current | 2.04 A | 1.10 A |
| max voltage, hot | 1.19 V | 3.23 V |
| max total power | 4.86 W | |

**A two-channel 30 V / 3 A bench supply covers every AWG 20 design** with
≥1.4× current margin and large voltage margin. Candidates meeting the
ratings: Rigol DP832 (2 × 30 V/3 A plus a 5 V rail, ~$429) or Siglent
SPD3303X-E (2 × 32 V/3.2 A, ~$349–459). Two single-channel units work
equally well and may be cheaper.

**Current requirement is gauge-dependent, and coarser wire is what strains
the supply**, since fewer turns need more current for the same A-turns:

| AWG | max inner current, all geometries | max outer | supply needed |
|---|---|---|---|
| 16 | 5.44 A | 2.65 A | ≥6.8 A channel |
| 18 | 3.26 A | 1.69 A | ≥4.1 A channel |
| **20** | **2.04 A** | **1.10 A** | **≥2.5 A channel** |
| 22 | 1.36 A | 0.66 A | ≥1.7 A channel |
| 24 | 0.95 A | 0.43 A | ≥1.2 A channel |

Total dissipation is 4.8–5.2 W in every row: gauge cancels, because thinner
wire raises resistance exactly as fast as the lower current reduces `I²`.
**Choosing a gauge is a supply-matching decision, not a thermal one.** AWG 20
sits where a commodity 3 A supply has comfortable margin and the turn count
is still hand-windable.

Low-voltage/high-current operation is the regime where lead and contact
resistance start to matter: at 5.4 A into a 0.5 V load, a couple of metres
of lead is a large fraction of the total. CC regulation handles it, but it
argues against the coarse end independently of the supply rating.

## Bill of materials — wire

AWG 20 polyimide, baseline geometry. Order to the **hex** case: if the
winding nests you need more wire, not less.

| | length | |
|---|---|---|
| inner coil | 26.1 m (86 ft) | 280 turns |
| outer coil, each | 38.5 m (126 ft) | 392 turns, ×4 |
| **total** | **179.9 m (590 ft)** | |
| +25% winding waste | 225 m (738 ft) | leads, start/finish, one botched coil |
| +40% | 252 m (826 ft) | |

**Order a 2.5 lb spool (786 ft / 240 m)** — covers the hex case with 33%
margin. Step up to 5 lb (1573 ft) if a rewind is likely; the wire is the
cheapest part of this assembly and running short mid-coil wastes a whole
winding.

Also needed: 0.05 mm polyimide tape for ground insulation and layer
interleaving, and a way to count turns.

## Soft inputs

- **Ground insulation is an allowance, not a spec** (0.15 mm). If the real
  build uses thicker tape or adds interlayer insulation, the inner window
  loses a layer before anything else does. It is the input most likely to
  bite.
- **Packing model.** Both are reported; the design point uses square. If the
  coils are machine-wound and actually nest, the outer channel drops to
  0.547 A and nothing else changes.
- **The 300 G exit target is still not derived** — `basic_het_sizing.py`
  assumes 150 G for its magnetization check, and nothing reconciles the two.
  That would move the A-turn requirement, and hence every current here, by a
  factor this document cannot bound. See `../MATERIALS.md` "Open questions".
- **Self-inductance per channel is not computed.** Only the total stored
  energy (65.9 mJ) is known; splitting it into two self-inductances plus a
  mutual term needs one solve per channel energised alone. It affects ramp
  timing only, not the supply rating or the flyback diode.

## Sources

- NEMA Standards Publication MW 1000, via the [magnet wire / winding wire engineering data handbook (12/2009 Rev.1)](https://www.casa.co.nz/cables/wire/magnet/EngData-Book_66p.pdf) — round heavy-build film dimensions (p.18) and round served polyester-glass dimensions (p.36)
- [MWS Wire — magnet wire insulation guide](https://mwswire.com/wp-content/uploads/2018/01/Insulation-Guide-Types-and-Characteristics-v2.pdf) — Dacron Glass 155 (MW 45) and Dacron Glass High Temp (MW 51) size restrictions and temperature classes
- [Remington Industries — 20 AWG 240 °C polyimide](https://www.remingtonindustries.com/magnet-wire/magnet-wire-240-c-20-awg-polyimide-7-spool-sizes-available/) — MW 16-C, 0.0351 in overall, spool sizes
- [S&W Wire](https://www.swwireco.com/) and [Bridgeport Insulated Wire](https://www.thebridgeportwireco.com/products) — glass/Daglas-served magnet wire, made to order
