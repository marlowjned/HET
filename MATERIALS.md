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

Top plate, bottom plate, inner and outer emag cores. Possibly also
`chamber_spacer` -- see "Open questions".

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
| outer emag core (x4) | 25.4 mm dia round | no -- turned from bar |
| inner emag core | 42.9 mm dia round | no -- turned from bar |
| chamber_spacer | 22.2 mm | no |

The cores are round turned parts and must be sourced as bar stock
elsewhere, so **the cores and plates will not be the same lot and need not
be the same grade** unless that is specified on purpose. Worth deciding
rather than inheriting -- the cores carry the highest flux density in the
circuit (peak 1.023 T at the pole tips), so if either part deserves the
better magnetic material it is those.

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

## Copper (windings)

| property | value | used by | confidence |
|---|---|---|---|
| resistivity at 20 C | 1.68e-8 ohm-m | `thermal/radiation_balance.py` | constant |
| temperature coefficient | 0.00393 /K | as above | constant |
| AWG bare diameters | 16: 1.291, 18: 1.024, 20: 0.812, 22: 0.644, 24: 0.511 mm | winding trade study | constant |

Standard annealed-copper values. The magnetics solve does **not** use
copper properties at all -- windings enter only as a prescribed current
density, and `nu[Windings] = 1/mu0` treats them as non-magnetic.

---

## Fiberglass-served magnet wire

Chosen for insulation quality and plasma arc resistance rather than
temperature capability.

| property | value | used by | confidence |
|---|---|---|---|
| insulated OD, single glass serving | 16: 1.45, 18: 1.18, 20: 0.96, 22: 0.79, 24: 0.65 mm | turns-that-fit calculation | **estimate** |
| insulated OD, double glass serving | 16: 1.58, 18: 1.30, 20: 1.09, 22: 0.92, 24: 0.78 mm | as above | **estimate** |
| continuous temperature limit | 200 C | `thermal/radiation_balance.py` limit | **estimate** |

**These are the softest inputs in the winding design.** The insulated
diameters set how many turns fit, which sets the current and voltage (not
the field, and not the dissipation -- see
`gmsh_getdp/README.md` "Solenoid winding design"). A real wire datasheet
could move turn counts by 10-20%.

The 200 C limit is deliberately conservative; glass serving itself
tolerates far more and the real ceiling depends on the binder. If that is
firmed up, PTFE becomes the binding thermal constraint instead.

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

- **`chamber_spacer` material is unconfirmed.** Currently treated as iron;
  `build_assembly.py --exclude-chamber-spacer` tests the inert
  alternative. Nobody has confirmed which it is.
- **The cores cannot come from SendCutSend.** They are turned round parts,
  not sheet, so they need bar stock from another supplier -- and therefore
  a deliberate grade choice. If they end up as something other than A36,
  this file needs a second iron entry rather than an edit to the first.
- **No numeric A36 B-H curve.** The model runs 1008 data as a proxy, which
  is optimistic about the iron. This is the single weakest input that
  anything depends on.
- **The 300 G channel-exit target is not derived.**
  `sizing/basic_het_sizing.py` assumes 150 G for its magnetization check.
  Those two numbers disagree by 2x and nothing currently reconciles them.
  Unlike most entries here, this one would move the design point, not just
  its uncertainty band.
- **Nothing in this file is measured on the actual hardware.** The B-H
  curve and the emissivity are the two worth measuring if anything is
  going to be.

---

## Sources

- [makeitfrom.com -- ASTM A36 (SS400, S275) Structural Carbon Steel][mifa36] -- iron thermal/physical properties
- [makeitfrom.com -- SAE-AISI 1008 (G10080) Carbon Steel][mif1008] -- 1008 properties, for the proxy B-H comparison
- [design1st -- Thermal Emissivity Values][emis] -- emissivity tables
- [The Plastic Shop -- PTFE Technical Information][ptfe] -- PTFE service temperature, arc resistance
- [FEMM -- DC Magnetization Curves of Soft Magnetic Materials][femm] -- context on low-carbon steel B-H data and its origin (Metals Handbook, 8th ed., Vol. 1, ASM 1966)
- [SendCutSend -- cold rolled mild steel][scs] -- **procurement source**; grade, mechanical properties, thicknesses, tolerances
- [FNAL via INIS -- B vs H curves for 1008 and 1020 steels][fnal] -- a real measurement of 1008 toroids; numeric data is in the PDF and has not been extracted. **This is the best candidate for replacing the secondhand B-H curve.**

[mif1008]: https://www.makeitfrom.com/material-properties/SAE-AISI-1008-G10080-Carbon-Steel
[mifa36]: https://www.makeitfrom.com/material-properties/ASTM-A36-SS400-S275-Structural-Carbon-Steel
[emis]: https://www.design1st.com/Design-Resource-Library/engineering_data/ThermalEmissivityValues.pdf
[ptfe]: https://www.theplasticshop.co.uk/ptfe-technical-information.html
[femm]: https://www.femm.info/wiki/SoftMagneticMaterials
[fnal]: https://inis.iaea.org/records/wmz02-45d12
[scs]: https://sendcutsend.com/materials/mild-steel/cold-rolled/
