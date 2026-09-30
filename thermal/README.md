# Thermal

Steady-state thermal model of the magnet circuit **running on its own** --
coil ohmic dissipation only, no plasma load. The question it answers is
"does anything get too hot just from energising the magnets, and at what
current does that change".

## Running

```
cd thermal
python radiation_balance.py
```

No arguments. It reads the radiating envelope from
`../gmsh_getdp/assembly/assembly.msh` if that exists (convex hull of the
iron and winding nodes) and falls back to a bounding cylinder otherwise,
so it runs standalone.

## Result at the design point

AWG 20 heavy-build polyimide, 280 turns inner at 1.207 A and 336 turns
outer at 0.638 A, on two independent constant-current channels -- see
`../sizing/WIRE_AND_SUPPLY.md` and `../gmsh_getdp/README.md`
"Solenoid winding design":

| | value |
|---|---|
| dissipation | 3.3 W |
| iron / PTFE | 50 C |
| coils, bonded to core | 50 C |
| coils, radiation-coupled only | 77 C |

Nothing is close to a limit. Margin to the first one, scanning the inner
channel with the outer tracking it at a fixed ratio:

| limit | | current | vs design |
|---|---|---|---|
| wire insulation (polyimide) | 240 C | 4.17 A | 3.5x |
| PTFE chamber | 260 C | 4.42 A | 3.7x |
| ASTM A36 service limit | 400 C | 6.11 A | 5.1x |
| ASTM A36 Curie point | 770 C | 9.81 A | 8.1x |

Wire insulation is still the first limit, but with 3.5x margin instead of
the 2.4x reported before. Both numbers moved for the same reason: the
previous design point specified fiberglass-served wire at a 200 C class,
and that wire cannot be bought at this gauge -- the round glass product is
Class 155, which would have made this table *worse*, not better. The
polyimide actually specified is 240 C. See `../MATERIALS.md` "Magnet wire".

Dissipation fell from 5.8 W to 3.3 W because both coils are now wound full
on separate channels, rather than the outer being part-wound to 45% of its
window to hit an ampere-turn ratio from a single series current, and
P ~ 1/N.

Iron properties are ASTM A36 (k = 50 W/m-K, cp = 470 J/kg-K,
rho = 7900 kg/m^3, service 400 C, Curie ~770 C), from
[makeitfrom.com](https://www.makeitfrom.com/material-properties/ASTM-A36-SS400-S275-Structural-Carbon-Steel).
See `../MATERIALS.md` for sourcing and confidence.
Heat capacity and density are carried but unused by this steady-state
model -- they are what a transient warm-up version would need.

## What the model does and doesn't include

The thruster is assumed **radiatively isolated in vacuum**: no convection,
and no conduction path to a mount. Every watt therefore leaves by
radiation from the outer envelope, which makes the result a near-pure
radiation balance. If the real test fixture bolts the thruster to
something cold, this model is conservative.

Iron is treated as **one isothermal node**: conduction through it costs
~2 K at these power levels, far below the uncertainty in emissivity. The
coils are separate nodes because that interface is the only place a real
gradient forms -- in vacuum a coil that isn't bonded to the core can only
radiate across its gap, which is worth ~25 K.

**PTFE has no source term** in this configuration, so it floats at the
iron temperature. It still sets the current ceiling, because it has the
lowest limit in the build.

## The assumption that actually matters

**Envelope emissivity**, and it isn't close:

| eps | hottest part | surface |
|---|---|---|
| 0.15 | 98 C | polished steel |
| 0.30 | 65 C | bare machined steel (assumed) |
| 0.50 | 50 C | lightly oxidised |
| 0.80 | 42 C | blackened / anodised |

Sourced ranges (see `../MATERIALS.md`): mild steel is 0.20-0.32, so the
0.30 default sits at the *optimistic* end; real black coatings reach
0.88-0.91, better than the 0.80 used here.

That's a 56 C spread at the design current from one number, larger than
any geometric detail this model could resolve. Blackening the exterior is
the cheapest thermal margin available. If the real finish is unknown, read
the sweep rather than the single design-point figure.

## Not modelled

Plasma heat load to the channel walls and anode, which for a running
thruster is one to two orders of magnitude above the coil dissipation and
would dominate everything here. This model is deliberately the
magnets-only case; a plasma-loaded model needs the chamber in the domain
as a real heat-receiving body, and needs a heat flux distribution that
nothing in this repo currently produces.
