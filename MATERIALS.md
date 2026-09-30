# Material properties

Every material property used anywhere in this repo, what uses it, where it
came from, and how much it should be trusted. Added because several
numbers in this project were placeholders that looked like data, and one
of them (the iron B-H curve) survived months of use before anyone asked
where it came from.

## How to read the confidence column

| tag | meaning |
|---|---|
| **measured** | measured on the actual BPL-700 hardware. Nothing here qualifies yet. |
| **verified** | the cited source was fetched and the value read directly out of it |
| **corroborated** | consistent across several independent secondary sources; no primary read |
| **secondhand** | single secondary source, primary inaccessible |
| **estimate** | representative engineering value, no specific source |
| **constant** | defined or standard physical constant |
| **derived** | computed from another row in this file |

---

## ASTM A36 structural steel

Top plate, bottom plate and `chamber_spacer` (confirmed 2026-09-21: the
spacer is the same soft iron as the plates). The emag cores are a
different grade -- see "ASTM A108 cold-finished bar" below.

### Procurement source

Parts come from **SendCutSend hot-rolled mild steel**, specified as
**A36/A1018 HRP&O** (hot rolled, pickled and oiled). They also stock 1008
cold rolled, but only to 0.135 in -- too thin for anything in this
assembly, which is what selects A36.

| property | value | confidence |
|---|---|---|
| grade | A36 / A1018 HRP&O | **verified** |
| carbon content | 0 to 0.26% | **verified** |
| manganese | up to ~1.2% | corroborated |
| density | 7.9 g/cm^3 (supplier: 7.85) | **verified** |
| tensile strength (ultimate) | 480 MPa (69 ksi) | **verified** |
| yield strength | 290 MPa (41 ksi) | **verified** |
| available thicknesses | .187 / .250 / .313 / .375 / .500 in | **verified** |
| laser cutting tolerance | +/-.005 in | **verified** |
| max part size | 30 x 44 in (30 x 56 in on custom quote) | **verified** |

Sources: [SendCutSend -- mild steel][scs], [makeitfrom.com ASTM A36][mifa36].

**Manufacturability.** The 0.375 in plates are a direct hit: A36 is
stocked at exactly .375 in, so each plate is a single laser-cut piece. No
lamination, and none of the interlayer-reluctance concern that a laminated
stack would have brought.

| part | thickness / dia | A36 sheet? |
|---|---|---|
| top plate | 9.53 mm (0.375 in) | **yes, exact stock thickness** |
| bottom plate | 9.53 mm (0.375 in) | **yes, exact stock thickness** |
| chamber_spacer | emag_height - 1.125 in (0 to 22.2 mm) | **yes at 0.187 / 0.250 / 0.313 / 0.375 / 0.500 in**, i.e. emag_height 1.312 / 1.375 / 1.438 / 1.500 / 1.625 in |
| outer emag core (x4) | 25.4 mm dia shaft | no -- A108 bar, see below |
| inner emag core | 25.4 mm dia shaft (= inner_coil_id), 34.9 mm pole cap | no -- A108 bar, see below |

(An earlier version of this table gave the inner core as 42.9 mm dia; the
CAD measures a 25.4 mm shaft under the winding and a 34.9 mm pole cap.)

The spacer is an annulus (r 17.46-41.28 mm) whose height tracks
emag_height, so a spacer height that lands on a stock thickness is a single
laser-cut part, and one that doesn't needs a stack or machining. Worth
weighing in the spacer-height trade (`gmsh_getdp/assembly/emag_sweep.py`).

Laser cutting leaves a heat-affected zone with altered microstructure at
every cut edge; the Fermilab measurement cited below deliberately sampled
"far away from the flame cutting region" for this reason. On parts this
size the HAZ is a non-trivial fraction of the whole.

### Magnetic

| property | value | used by | confidence |
|---|---|---|---|
| DC magnetization curve | 19 points, H 0-397887 A/m, B 0-2.585 T | `gmsh_getdp/assembly/generate_regions.py` | **secondhand** |
| saturation (B at 1000 Oe) | 2.165 T | -- | derived |
| peak relative permeability | 2164, at H = 318 A/m | -- | derived |
| mu_r at this design's peak iron field | ~1935 at 1.023 T | -- | derived |
| Curie temperature | ~770 C | `thermal/radiation_balance.py` limit | estimate |

**The curve in the model is AISI 1008 data, used as a proxy for A36.** No
numeric A36 B-H table could be obtained; the published A36 magnetization
work is transformer-tank literature behind paywalls. It is deliberately
not relabelled.

A36 carries up to 0.26% carbon against 1008's ~0.08%, plus up to 1.2%
manganese -- more pearlite, more domain-wall pinning, so **real A36 is
less permeable than this curve and the model is optimistic about the
iron**. Partly offsetting: this A36 is hot-rolled, a more relaxed
microstructure than cold-rolled 1008. Empirically, swapping the generic
GetDP curve for 1008 dropped mu_r 2.7x and moved the channel field 3.9%;
if A36 is another ~1.8x down, expect a few percent more, i.e. a design
current nearer 1.35 A than 1.30 A. That is an extrapolation, not a
result.

The B-H points originate in the Ansys Maxwell SV material library and were
transcribed via a public engineering forum. The primary source returned
HTTP 403 and could not be re-read, so this is the weakest-provenance
number in the project that anything actually depends on.

Reasons to believe it is genuine digitized data rather than invented:
every H value is an exact Oersted multiple (2, 4, 6, 8, 10, 20, 40, 60,
80, 100, 200, 400, 600, 800, 1000, 2000, 4000, 5000 Oe), which is how
pre-SI magnetization tables were published; the curve is monotonic in both
variables; differential permeability rises then falls as a real
magnetization curve must; and B at 1000 Oe lands in the 2.1-2.2 T band
expected of low-carbon steel.

Reasons to stay cautious:

- **AISI 1008 is not an electrical steel.** ASTM imposes no magnetic
  requirement on it at all. Two suppliers' 1008 may differ materially and
  neither is obliged to match this curve.
- **Processing dominates.** Peak mu_r of 2164 suggests as-rolled rather
  than annealed material; annealed low-carbon steel reaches 3000-5000.

**Why this matters less than it looks.** The magnetic circuit is
gap-dominated: at the design point the 37.4 mm air gap costs ~893
A-turns while the 227 mm iron path costs ~9. Halving or doubling mu_r
moves the total by well under a percent. What the iron buys is not
sensitivity to its own grade but the ~22x reduction in required A-turns
versus having no iron at all (that same path would cost ~18,800 A-turns
at mu_r = 1). Grade is a detail; presence is not.

### Thermal and physical

All from [makeitfrom.com ASTM A36 (SS400, S275)][mifa36], fetched and read
directly.

| property | value | used by | confidence |
|---|---|---|---|
| thermal conductivity | 50 W/m-K | justifies the isothermal-iron assumption | **verified** |
| specific heat capacity | 470 J/kg-K | carried for a future transient model | **verified** |
| density | 7900 kg/m^3 (supplier: 7850) | as above | **verified** |
| max service temperature | 400 C | `thermal/radiation_balance.py` limit | **verified** |
| elastic modulus | 190 GPa | unused | **verified** |
| thermal expansion | 12 um/m-K | unused | **verified** |
| electrical conductivity | 12% IACS | unused (relevant only to eddy currents) | **verified** |
| melting range | 1420-1460 C | unused | **verified** |

---

## ASTM A108 cold-finished bar (emag cores)

Inner and outer emag cores. Round turned parts, so they come from bar
stock rather than SendCutSend sheet: **McMaster "multipurpose low-carbon
steel" rod, ASTM A108**. A108 is a specification for cold-finished carbon
bar and covers several grades; for this product it is taken to be
**1018, cold drawn** (see "Open questions").

### Physical and thermal

From [makeitfrom.com -- cold drawn 1018][mif1018], fetched and read
directly.

| property | value | confidence |
|---|---|---|
| carbon | 0.15-0.20% | **verified** |
| manganese | 0.6-0.9% | **verified** |
| density | 7.9 g/cm^3 | **verified** |
| tensile strength (ultimate) | 480 MPa | **verified** |
| yield strength | 400 MPa | **verified** |
| thermal conductivity | 52 W/m-K | **verified** |
| specific heat | 470 J/kg-K | **verified** |
| max service temperature | 400 C | **verified** |
| electrical conductivity | 7.0% IACS | **verified** |

Thermally it is interchangeable with A36 (same conductivity, heat
capacity and service limit), so the thermal model's single isothermal
iron node covers both without change.

### Magnetic

| property | value | confidence |
|---|---|---|
| DC magnetization curve | measured 1020, 21 points, 0-1.9 T | **verified** (primary source read; grade is a proxy) |
| mu_r at 1.0 T | ~2170 | derived |
| mu_r at 1.2 T | ~1460 | derived |
| mu_r at 1.4 T | ~820 | derived |
| mu_r at 1.5 T | ~520 | derived |

Curve: [Fermilab TM-1197][fnal] Table 2, toroids machined from purchased
1020 stock and measured by Arnold Engineering. The report PDF was fetched
and all 21 rows read off the scanned table page, not just its OCR text.
Stored as `FNAL1020_B`/`FNAL1020_H` in
`gmsh_getdp/assembly/generate_regions.py`.

**This is a much better-provenanced curve than the plates' A36 proxy** --
a primary measurement, not a forum transcription -- but it is still 1020
standing in for 1018. The two are adjacent plain-carbon grades (1018:
0.15-0.20% C, 1020: 0.18-0.23% C), so the curve should be close. Cold
drawing adds residual strain, which lowers permeability relative to
annealed bar; FNAL describes its samples as behaving like hot-rolled
material. If the cores ever run near the knee, a stress-relief anneal
after machining is the cheap fix.

What the curve is used for: the core saturation margin in the spacer
sweep. The sweep itself runs linear iron, since the channel field is
insensitive to iron grade (gap-dominated circuit); what grade does
control is where a thinned core saturates. The knee on this curve is
around 1.3-1.5 T (mu_r falls from ~1460 at 1.2 T to ~520 at 1.5 T).

The same report's Table 1 is a measured **1008** curve -- a primary-source
candidate to replace the secondhand Maxwell 1008 data the plates' proxy
uses. Not yet transcribed.

---

## Copper (windings)

| property | value | used by | confidence |
|---|---|---|---|
| resistivity at 20 C | 1.68e-8 ohm-m | `thermal/radiation_balance.py` | constant |
| temperature coefficient | 0.00393 /K | as above | constant |
| AWG bare diameters | 16: 1.290, 18: 1.024, 20: 0.813, 22: 0.643, 24: 0.511 mm | winding trade study | **verified** |

Standard annealed-copper values; the bare diameters are now the NEMA
MW 1000 nominals ([source][nemamw]) rather than a generic AWG table, and
agree with it to within 0.001 mm. The magnetics solve does **not** use
copper properties at all -- windings enter only as a prescribed current
density, and `nu[Windings] = 1/mu0` treats them as non-magnetic.

---

## Magnet wire

**Specified: AWG 20 heavy-build aromatic polyimide, NEMA MW 16-C, 240 C.**
Full write-up, including the power supply, in
[`sizing/WIRE_AND_SUPPLY.md`](sizing/WIRE_AND_SUPPLY.md); the tables live in
`sizing/wire_data.py`.

This entry used to describe fiberglass-served wire with estimated
diameters, and was flagged here as the softest input in the project. It was
worse than soft -- it described a product that cannot be bought. Replaced
2026-09-30.

### Dimensions

All from NEMA MW 1000, and all **maximum** overall diameters rather than
nominal: a turns-that-fit calculation wants the worst case a supplier is
allowed to ship.

| property | value | used by | confidence |
|---|---|---|---|
| bare copper dia | 16: 1.290, 18: 1.024, 20: 0.813, 22: 0.643, 24: 0.511 mm | winding design | **verified** |
| max OD, heavy-build film | 16: 1.384, 18: 1.110, **20: 0.892**, 22: 0.714, 24: 0.577 mm | turns-that-fit calculation | **verified** |
| max OD, single glass serving | AWG 12 only: 2.200 mm (bare 2.052) | fiberglass fallback | **verified** |
| max OD, double glass serving | AWG 12 only: 2.276 mm | fiberglass fallback | **verified** |
| temperature class, polyimide | 240 C (MW 16-C) | `thermal/radiation_balance.py` | **verified** |
| temperature class, glass served | 155 C (MW 45) | as above, fallback | **verified** |
| ground insulation allowance | 0.15 mm radial | winding design | **estimate** |

Source: [NEMA MW 1000 via the magnet wire engineering data handbook][nemamw],
round heavy-build film dimensions p.18, round served polyester-glass p.36.
Transcribed by parsing the PDF text layer, then range-checked in
`sizing/wire_data.py`'s self-test.

Heavy-build dimensions are **insulation-chemistry-independent** -- MW 16-C
(polyimide, 240 C), MW 35-C (polyurethane/nylon, 200 C) and MW 80-C
(polyamideimide, 200 C) share this table. Only temperature class and price
change with the polymer, so the turn counts survive a change of wire
chemistry.

**One transcription error found and corrected.** The AWG 12 double-glass OD
is printed as 0.0986 in. Every other row's double serving adds 0.226-0.310 mm
over bare; 0.0986 in would add 0.452 mm. It is a digit transposition of
0.0896 in (0.0088 in over bare, matching AWG 11's 0.226 mm and AWG 10's
0.229 mm). Corrected in `wire_data.py` with the reasoning recorded there,
and the consistency check that caught it is kept as an assertion.

### Why not fiberglass

Three independent findings, any one of which would have been enough:

1. **Not a standard size.** NEMA MW 1000 tabulates round polyester-glass
   served wire from 1/0 down to **AWG 12 and no finer**. Glass serving is
   filaments, whose thickness does not scale with the conductor -- it adds
   ~0.15-0.21 mm to the OD at every gauge -- so it stops being offered on
   fine wire. The old design point was AWG 20.
2. **Not stocked, at any gauge.** No US retailer carries round glass-served
   magnet wire. [S&W Wire][swwire] and [Bridgeport Insulated Wire][bpt]
   make it to order.
3. **The round product is Class 155, not 200 C.** Round glass serving is
   Dacron Glass 155 (MW 45). The 200 C variant (Dacron Glass High Temp,
   MW 51) is made "only in shaped or heavy round wire sizes"
   ([MWS insulation guide][mwsins]). The old 200 C entry here was therefore
   optimistic for the wire it named -- and it was the *first* limit reached
   in the thermal model.

Polyimide is better on every axis that matters: thinner than the glass
estimate it replaces (so **more** turns, not fewer), dimensions from a
standard rather than an estimate, 240 C, and a solid film with no porous
covering to trap gas -- which matters in vacuum. Glass serving is
specifically a poor vacuum choice for that last reason.

The original rationale for glass was plasma arc resistance. The coils sit
outside the discharge chamber with the PTFE wall between them and the
plasma, so arc resistance at the wire is second-order; vacuum compatibility
and temperature class are not.

Fiberglass is kept as a **spec-only fallback** at AWG 12 in
`sizing/wire_data.py`, so it can be ordered if a requirement ever demands
it. It needs a made-to-order run and pushes the inner channel to ~7.3 A.

### What this changed downstream

Turn counts went **up** (the wire is thinner), and the two-channel supply
decision let both coils be wound full instead of the outer being part-wound
to hit an ampere-turn ratio:

| | before | after |
|---|---|---|
| inner / outer turns | 260 / 165 | 280 / 336 |
| current | 1.30 A series | 1.207 A / 0.638 A, two channels |
| dissipation | 5.8 W | 3.0 W |
| design-point temperature | 64 C | 50 C |
| first limit | wire at 200 C, 2.4x design current | wire at 240 C, 3.5x |

**Ampere-turns are unchanged** (338.0 / 214.5), so every field solved in
this repo remains valid -- only the current/turns split moved.

---

## PTFE (chamber)

| property | value | used by | confidence |
|---|---|---|---|
| max continuous service temperature | 260 C | `thermal/radiation_balance.py` limit | **verified** |
| thermal conductivity | 0.24-0.25 W/m-K | not used (PTFE has no source term) | corroborated |
| relative permeability | ~1 | why it is excluded from the FEM domain | constant |
| arc resistance (ASTM D495) | 700 s | not modelled; motivates the material choice | **verified** |
| transition points | 19 C and 327 C | -- | **verified** |
| decomposition | appreciable above 400 C | -- | **verified** |

Magnetically PTFE is indistinguishable from vacuum, which is why
`build_assembly.py` leaves the chamber out of the solid model entirely --
that is a recognition that the physics does not care, not a shortcut.

Thermally it matters despite having no heat source of its own: it floats
at the iron temperature and has the lowest limit of anything in the
build, so it sets the current ceiling.

Note PTFE is a poor choice once plasma is actually present (erosion,
outgassing). Presumably a test article rather than the flight part.

---

## Surface emissivity

The single most consequential thermal input: it moves the design-point
temperature from 40 C to 121 C across plausible finishes, more than any
geometric detail the model resolves. From [design1st emissivity
tables][emis], fetched and read directly.

| surface | emissivity | confidence |
|---|---|---|
| steel, polished | 0.066 | **verified** |
| iron, polished | 0.14-0.38 | **verified** |
| **mild steel** | **0.20-0.32** | **verified** |
| iron plate, rusted red | 0.61 | **verified** |
| sheet, rough oxide layer | 0.81 | **verified** |
| black anodize | 0.88 | **verified** |
| black paints (3M Black Velvet, Chemglaze Z306) | 0.91 | **verified** |

`thermal/radiation_balance.py` defaults to **0.30**, which is at the
*optimistic* end of the sourced 0.20-0.32 mild-steel band. At 0.20 the
design point runs 77 C rather than 61 C. The script sweeps emissivity for
exactly this reason -- read the sweep, not the single default.

Blackening the exterior is the cheapest thermal margin available: real
black coatings reach 0.88-0.91, better than the 0.80 used as the
"blackened" case, and would drop the design point to ~40 C.

---

## Vacuum / air

| property | value | used by | confidence |
|---|---|---|---|
| mu0 | 4*pi*1e-7 H/m | everywhere | constant |
| relative permeability of air | 1 | `nu[Air] = 1/mu0` | constant |
| Stefan-Boltzmann constant | 5.670374419e-8 W/m^2-K^4 | `thermal/radiation_balance.py` | constant |
| ambient (chamber wall) temperature | 300 K | thermal boundary condition | **estimate** |

No convection and no conduction to a mount are assumed. If the real
fixture bolts the thruster to something cold, the thermal model is
conservative.

---

## Open questions

- **Which A108 grade the core bar is.** A108 is a cold-finished-bar
  specification, not a grade; McMaster's "multipurpose low-carbon" rod is
  taken to be 1018. Worth confirming on the order/cert.
- **No numeric A36 B-H curve.** The model runs 1008 data as a proxy, which
  is optimistic about the iron. This is the single weakest input that
  anything depends on.
- **The 300 G channel-exit target is not derived.**
  `sizing/basic_het_sizing.py` assumes 150 G for its magnetization check.
  Those two numbers disagree by 2x and nothing currently reconciles them.
  Unlike most entries here, this one would move the design point, not just
  its uncertainty band.
- **Ground insulation is an allowance, not a spec** (0.15 mm radial, two
  half-lapped wraps of 0.05 mm polyimide tape). Now the softest input in
  the winding design, and the one most likely to cost the inner coil a
  layer, since that window has only 0.152 mm of unused build left at
  AWG 20. Resolved by deciding the actual build method.
- **Which packing model the wound coils achieve.** Hex nesting and square
  packing are both reported; the design point uses square. At AWG 20 the
  inner coil holds 5 layers either way, so this only sets the outer
  channel's current (0.638 A square, 0.547 A hex) -- and with two channels
  it is measured after winding rather than assumed beforehand.
- **Nothing in this file is measured on the actual hardware.** The B-H
  curve and the emissivity are the two worth measuring if anything is
  going to be. Wire dimensions are now at least from a standard, and are
  verifiable with a micrometer on arrival -- worth doing, since they are
  maxima and the real wire will run under them.

---

## Sources

- [makeitfrom.com -- ASTM A36 (SS400, S275) Structural Carbon Steel][mifa36] -- iron thermal/physical properties
- [makeitfrom.com -- SAE-AISI 1008 (G10080) Carbon Steel][mif1008] -- 1008 properties, for the proxy B-H comparison
- [design1st -- Thermal Emissivity Values][emis] -- emissivity tables
- [The Plastic Shop -- PTFE Technical Information][ptfe] -- PTFE service temperature, arc resistance
- [FEMM -- DC Magnetization Curves of Soft Magnetic Materials][femm] -- context on low-carbon steel B-H data and its origin (Metals Handbook, 8th ed., Vol. 1, ASM 1966)
- [makeitfrom.com -- cold drawn 1018 carbon steel][mif1018] -- core (A108) physical/thermal properties
- [SendCutSend -- cold rolled mild steel][scs] -- **procurement source**; grade, mechanical properties, thicknesses, tolerances
- [NEMA MW 1000, via the magnet wire engineering data handbook][nemamw] -- **magnet wire dimensions**; round heavy-build film (p.18), round served polyester-glass (p.36). Parsed from the PDF text layer into `sizing/wire_data.py`
- [MWS Wire -- magnet wire insulation guide][mwsins] -- insulation types, NEMA/IEC designations, temperature classes, and the size restrictions that rule out fine round glass-served wire
- [Remington Industries -- 20 AWG 240 C polyimide][remington] -- **procurement source** for the specified wire; MW 16-C, 0.0351 in overall, spool sizes
- [S&W Wire][swwire], [Bridgeport Insulated Wire][bpt] -- glass/Daglas-served magnet wire, made to order (no stocked retail source exists)
- [FNAL TM-1197 via INIS -- B vs H curves for 1008 and 1020 steels][fnal] -- measured toroid curves. **Table 2 (1020) is transcribed and is the core material's curve.** Table 1 (1008) is the best candidate for replacing the plates' secondhand B-H proxy; not yet transcribed.

[mif1008]: https://www.makeitfrom.com/material-properties/SAE-AISI-1008-G10080-Carbon-Steel
[mif1018]: https://www.makeitfrom.com/material-properties/Cold-Drawn-1018-Carbon-Steel
[mifa36]: https://www.makeitfrom.com/material-properties/ASTM-A36-SS400-S275-Structural-Carbon-Steel
[emis]: https://www.design1st.com/Design-Resource-Library/engineering_data/ThermalEmissivityValues.pdf
[ptfe]: https://www.theplasticshop.co.uk/ptfe-technical-information.html
[femm]: https://www.femm.info/wiki/SoftMagneticMaterials
[fnal]: https://inis.iaea.org/records/wmz02-45d12
[scs]: https://sendcutsend.com/materials/mild-steel/cold-rolled/
[nemamw]: https://www.casa.co.nz/cables/wire/magnet/EngData-Book_66p.pdf
[mwsins]: https://mwswire.com/wp-content/uploads/2018/01/Insulation-Guide-Types-and-Characteristics-v2.pdf
[remington]: https://www.remingtonindustries.com/magnet-wire/magnet-wire-240-c-20-awg-polyimide-7-spool-sizes-available/
[swwire]: https://www.swwireco.com/
[bpt]: https://www.thebridgeportwireco.com/products
