# How this pipeline works: code and physics side by side

This walks through `toy_test/` and `assembly/` step by step, connecting
what the code does to what it means physically. For status/setup/known
issues, see `README.md` instead — this file is about *why* each step
exists.

## 1. The physics problem

We want the magnetic field `B` produced by current-carrying coils and
shaped by iron pole pieces. The governing equations (magnetostatics — no
time variation, so no induced E-field to worry about) are:

```
curl H = J        (Ampere's law: current is the source of H)
div B  = 0         (no magnetic monopoles)
B = mu * H          (material response: B and H related by permeability)
```

Two hard parts, both handled the same way in `toy_test` and `assembly`:

- **`div B = 0` has to hold everywhere, exactly**, not just approximately.
  Directly solving for `B`'s three components with a generic FEM
  discretization doesn't guarantee that. The standard fix: introduce the
  **magnetic vector potential** `A`, defined by `B = curl A`. Since
  `div(curl(anything)) = 0` identically, `div B = 0` is now automatic —
  solved by construction, not by the mesh. Substituting into Ampere's law
  gives the actual PDE we solve:

  ```
  curl( nu * curl A ) = J        where nu = 1/mu  (reluctivity)
  ```

  This is what every `.pro` file's `Formulation` block encodes (section 6
  below).

- **Real space is infinite; a mesh isn't.** The field technically extends
  to infinity, decaying with distance. Both scripts truncate this by
  building a large sphere of "air" around the real hardware (`air_radius
  = 8 * <assembly half-extent>` in both `build_toy.py` and
  `build_assembly.py`) and forcing `A = 0` on its outer boundary — see
  section 4. Physically this says "far enough away, there's no field,"
  which is true in the limit and a good approximation at 8x the source's
  own size.

## 2. Geometry construction (`build_toy.py`, `build_assembly.py`)

### Building the conductor and iron as solids

`toy_test/build_toy.py` builds one hollow cylinder (a coil) directly with
Gmsh/OpenCASCADE primitives (`occ.addCylinder`, then `occ.cut` to hollow
it out) sitting inside an air sphere. This is the "can I trust the
formulation at all" check — no real CAD, no placement complexity, just
geometry simple enough to compare against a hand-derivable analytic
answer (see section 8).

`assembly/build_assembly.py` does the same job on the real BPL-700
hardware: imports a single whole-assembly Onshape STEP export, which
carries every part's real placement and name already (no hand-maintained
transform table), and fuses the 7 real iron parts (top plate, bottom
plate, center pole, 4 outer poles) into one solid — folding them in
**one at a time**, since fusing them in a single batch call leaves 3
disjoint solids on this geometry (see `README.md`).

**Why fuse them at all, physically?** In the real hardware these iron
pieces are bolted/mated together into one continuous magnetic circuit —
flux has to pass through the top plate into an outer pole and back out
the bottom plate without any discontinuity. For the FEM mesh to represent
that correctly, the elements on either side of each mating face need to
share exactly the same nodes there (a *conformal* mesh) — otherwise the
solver sees two disconnected iron blobs with an air gap between them,
which is wrong: a real air gap of even a fraction of a millimeter in a
magnetic circuit adds enormous reluctance (see section 6) and would
visibly change the answer. This was the exact step that broke in
MATLAB (`../matlab/README.md`, "Blocker hit at step 1") — not a physics
problem, a numerical-tolerance one: two real STEP surfaces meant to be
coincident come in a few microns apart, and MATLAB's boolean kernel
couldn't decide whether that counts as touching. Gmsh/OCC's
`Geometry.ToleranceBoolean` fuzzy tolerance is the fix — it explicitly
tells the kernel how much slack to allow before treating two surfaces as
the same one.

**Why are chamber and injector left out entirely?** Both are magnetically
inert — ceramic and non-magnetic stainless steel, `mu_r ~= 1`, the same
as air. Ampere's law doesn't distinguish a `mu_r=1` object from empty
space; leaving that volume unlabeled and letting it fall into the "Air"
region has zero effect on `B`. This isn't a shortcut around the physics,
it's recognizing the physics doesn't care — see the long comment at the
top of `build_assembly.py` for the two geometry problems this also
sidesteps.

### The winding sleeves are a *homogenized* current source, not real wire

The windings come in from the STEP as real coil solids (Onshape exports
them at the configured ID/OD), and are treated as one smeared-out
conductor rather than hundreds of individual turns of real wire. This is
standard practice:
a coil of `N` turns carrying current `I`, wound tightly enough that
individual turns are much smaller than the mesh, is magnetically
equivalent to a solid conductor carrying a smeared-out **bulk current
density**

```
J = N * I / A_cross          (A/m^2, A_cross = coil's cross-sectional area)
```

flowing azimuthally around the pole axis. Meshing every real turn would
be enormously more expensive for no accuracy gain at this scale — see
section 5 for how that current density is actually built into the
GetDP source term.

### Fragmenting everything together

Both scripts finish geometry construction with one big `occ.fragment()`
call across every solid (iron, windings, air sphere). **Physically, this
is what makes the field continuity conditions at material interfaces
exact.** At any interface between two different materials, Maxwell's
equations require:

- the *tangential* component of `H` to be continuous,
- the *normal* component of `B` to be continuous.

Between iron (`mu_r` in the thousands) and air (`mu_r=1`), that means `B` bends
sharply at the boundary — field lines refract, much like light at a
lens surface, entering the iron nearly along the surface normal
regardless of the field direction in air (this is *why* iron pole pieces
work at all: they funnel and concentrate flux). A finite-element mesh can
only represent that correctly if the mesh nodes exactly coincide at the
interface (a "conformal" mesh) — `fragment()` is what guarantees that,
by construction, for every pair of solids that touch.

## 3. Units bookkeeping (not physics)

`build_assembly.py` scales every imported STEP part by `0.001` right
after import. This has nothing to do with the physics — it's a data
bug workaround: the STEP files' header declares `SI_UNIT(METRE)` but the
actual coordinate values are physically millimeters (confirmed by
comparing the imported bounding box against the known real part size).
MATLAB's importer silently assumes mm and rescales; Gmsh's doesn't. Get
this wrong and every length in the model — and therefore every field
value, since field strength depends on geometry — is off by 1000x.

## 4. Boundary conditions: what "the edge of the model" means physically

In every `.pro` file:

```
Constraint {
  { Name MVP_0;
    Case { { Region OuterBnd; Value 0; } }
  }
  ...
}
```

This pins the vector potential `A = 0` on the outer sphere's surface. As
in section 1, this represents "no field infinitely far from an isolated
source" — a real approximation (not an identity), which is why the air
sphere is built 8x larger than the source geometry: far enough that the
field has already decayed to something the `A=0` truncation doesn't
meaningfully distort. (The toy case's 1.4% agreement with the exact
analytic answer, computed on an *infinite*-domain formula, is direct
evidence this truncation distance is adequate.)

## 5. Assigning materials and current sources (`assembly_regions_generated.pro`)

This file is generated by `build_assembly.py`, not hand-written — see
its "Generate the GetDP region/current-source include file" section.

### Reluctivity: how "eager" each material is to carry flux

```
nu[Air]  = 1 / mu0;
nu[Iron] = SteelGeneric_nu[$1];   // B-dependent, see section 9
```

`nu` (reluctivity) is the magnetic analogue of electrical resistivity —
literally the coefficient in the governing PDE `curl(nu * curl A) = J`
from section 1. High `mu_r` (several thousand for this steel below
saturation — see section 9) means **low** reluctivity: flux preferentially routes through
iron rather than air, the same way current preferentially routes through
a low-resistance wire rather than through free space. This is the entire
reason the assembly is iron pole pieces shaped a certain way, not an air
gap with coils floating in it — the iron *channels* the flux to where
it's wanted (the gap between poles), instead of letting it spread out
diffusely through the surrounding air.

### The current source: encoding "wind the coil this way, this direction"

```
Js[Winding0] = (909248.99 / Sqrt[(X[]-(-0.125...))^2 + (Z[]-(0.077...))^2])
             * Vector[(Z[]-(0.077...)), 0, -(X[]-(-0.125...))];
```

This looks dense but is one physical idea: current flowing in a circle
around each pole's own axis (an azimuthal current), with magnitude
`J = N*I/A_cross` (section 2) and direction set by the right-hand rule —
curl your right hand's fingers in the current's direction, your thumb
points along the resulting `B` inside the coil. In vector form, that
azimuthal direction is `(pole axis) x (radial direction from pole axis)`
— a cross product. GetDP has no built-in cross-product operator, so
`build_assembly.py` hand-expands it into `X[]`/`Z[]` components (every
real pole's own axis runs along the global Y direction — confirmed
directly from the imported geometry, not assumed — which is what
collapses the general 3D cross product into just these two components).

**The sign is the physically interesting part.** `build_assembly.py`
picks it so all 4 outer poles produce `B` pointing the same way in global
coordinates, and the center pole opposite — i.e. flux exits one polarity,
arcs across the gap between center and outer poles, and returns through
the other. That arcing-across-the-gap field is the entire point of a
cusped/mirror magnetic circuit like this one (same design intent as
`../matlab/het_solenoid_bfield.m`'s `centerCoil.sign`/`outerCoil.sign`,
carried over here for the real geometry). Flip that sign convention and
you'd get a physically different, much weaker configuration where center
and outer fields just cancel head-on in the gap instead of connecting
through it.

## 6. The gauge condition: bookkeeping, not physics

```
{ Name MVP_Gauge;
  Case { { Region Domain; SubRegion OuterBnd; Value 0.; } }
}
```

paired with `EntityType EdgesOfTreeIn; EntitySubType StartingOn` in the
`FunctionSpace` block. This one is subtle enough that it's worth being
explicit that **it has no physical effect on the answer** — it exists
purely to make the linear algebra solvable.

`B = curl A` doesn't pin down `A` uniquely: for *any* scalar function
`phi`, `A' = A + grad(phi)` gives `curl A' = curl A` too (`curl(grad(
anything)) = 0` identically), so `A'` produces the exact same physical
`B`. Without an extra constraint, the discrete system built from `curl-
curl` is singular — it has infinitely many solutions `A` (differing by
some `grad(phi)`), all equally valid, and a direct solver can't pick one.
The **tree-cotree gauge** fixes this by pinning `A` to zero along a
spanning tree of mesh edges anchored on the outer boundary — this
removes exactly the `grad(phi)` freedom (nothing else) from the solution
space, so what's left is unique and solvable, and `B = curl A` comes out
identical to what any other valid gauge choice would give.

This was found the hard way: the first `toy_test` solve, run without this
constraint, "succeeded" (no error) but gave physically wrong, unstable
results (spurious off-axis field components, ~3x too weak on-axis).
Adding the gauge fixed both — see the comment in
`toy_test/magnetostatics_toy.pro` for the full derivation and the GetDP
manual tutorial it's based on.

## 7. The weak form (`Formulation` block): where the PDE actually gets solved

```
Galerkin { [ nu[] * Dof{d a}, {d a} ]; In Domain; ... }
Galerkin { [ -Js[], {a} ]; In DomainC; ... }
```

This is `curl(nu * curl A) = J` (section 1) rewritten in **weak form** —
instead of solving the PDE pointwise, multiply both sides by an arbitrary
test function `a'`, integrate over the whole domain, and integrate the
left side by parts once (this is what turns a second-derivative-like
operator into two first derivatives, `curl A` and `curl a'`, which is why
`d a` — GetDP's notation for `curl` of the `Form1` quantity `a` — appears
on both sides). The two lines above are exactly the two resulting terms:

- `nu[] * Dof{d a}, {d a}` — the "energy" term: proportional to magnetic
  energy density `(1/2) H . B` integrated over space, evaluated for every
  possible test function.
- `-Js[], {a}` — the source term: work done by the driving current.

The solver's job (`Resolution` block, `Solve[Sys_a]`) is finding the `A`
that balances these for every possible test function simultaneously —
that's what "finite element solve" means concretely here. Edge elements
(`BF_Edge`, one degree of freedom per mesh edge, not per node) are the
natural discretization for a `curl`-based quantity like `A`: unlike
ordinary nodal elements, they don't force spurious continuity of every
vector component across material boundaries, which is what real physics
requires there (section 2).

## 8. Validation logic (`toy_test/`)

Before trusting any of the above on real, complex CAD, `toy_test/
build_toy.py` builds a single coil simple enough to have an exact
closed-form answer, and `post_process`-equivalent logic in
`build_toy.py`'s own header computes that answer directly:

```
B(0) = (mu0 * J / 2) * integral over the coil's r,z extent of
       r^2 / (r^2 + z^2)^(3/2)
```

— the on-axis field from a "thick" solenoid (finite length *and* finite
radial wall thickness, not the textbook infinitely-thin-shell
approximation, which would be a bad match for this coil's proportions:
its radial wall thickness is comparable to its mean radius). The FEM
result (10.05 mT) landing within 1.4% of this exact integral (10.19 mT)
is the actual evidence that sections 1, 4, 5, 6, and 7 above are all
correctly implemented — not just "the solver ran without errors."

## 9. Reading the result, and where the model is still a simplification

`post_process.py` parses GetDP's `.pos` output (`b = {d a}`, i.e. `B =
curl A`, evaluated per mesh element — piecewise-constant, since `A` is
represented by first-order edge elements and `curl` of a first-order
field is zeroth-order/constant per element) and reports summary
statistics plus a cross-section plot.

Iron is no longer the `mur_iron = 1000` linear constant this section
originally described. `generate_regions.py` now assigns a real saturating
B-H curve (GetDP's bundled `SteelGeneric`), because a constant `mu_r` has
no saturation cap and was reporting peak `|B|` near 9.7 T — a real core
would never reach that; it would saturate, its effective `mu_r` would
collapse toward 1 in the saturated region, and flux would redistribute in
response. That nonlinearity is why the `.pro`'s Galerkin term is split
into linear (Air+Windings) and Picard-iterated nonlinear (Iron) parts —
see `README.md`, "Nonlinear iron (B-H curve)".

The interesting practical result is that, at this thruster's actual
operating point, it barely matters. The magnetic circuit is
**gap-dominated**: the iron contributes ~2.5 of the ~900 A-turns needed
for a 300 G channel field, so the answer is insensitive to iron
permeability as long as it's large. At the design excitation peak iron
`|B|` is 0.859 T — below the knee, where the curve is nearly straight —
and the nonlinear result matches a linear-scaled prediction from the old
`mur_iron=1000` run to 0.4%.

The general lesson still holds, though: trust this model's field *shape*
(where flux concentrates, the pattern across the gap) more than its
peak-field numbers, since shape is far less sensitive to material
modeling than the values at the most concentrated hot spots are. And if
an excitation is ever chosen that does push the iron into saturation, the
nonlinear solve becomes dramatically more expensive — that history is
documented in `README.md`.
