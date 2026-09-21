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

1.30 A, 260/165 turns of AWG 20 -- see `../gmsh_getdp/README.md`
"Solenoid winding design":

| | value |
|---|---|
| dissipation | 5.8 W |
| iron / PTFE | 64 C |
| coils, bonded to core | 65 C |
| coils, radiation-coupled only | 92 C |

Nothing is close to a limit. Margin to the first one:

| limit | | current | vs design |
|---|---|---|---|
| fiberglass wire insulation | 200 C | 3.09 A | 2.4x |
| PTFE chamber | 260 C | 3.76 A | 2.9x |
| ASTM A36 service limit | 400 C | 5.26 A | 4.0x |
| ASTM A36 Curie point | 770 C | 8.79 A | 6.8x |

The 200 C wire figure is deliberately conservative -- glass serving itself
tolerates far more, and the real limit depends on the binder. If that
number is firmed up, PTFE becomes the binding constraint.

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
