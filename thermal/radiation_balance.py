"""
radiation_balance.py

Steady-state thermal model of the BPL-700 magnet circuit running on its
own -- coil ohmic dissipation only, no plasma load. Answers "does anything
get too hot just from energising the magnets", and at what current the
answer changes.

Scope and why it's this simple:

- **Vacuum.** No convection, and the thruster is assumed radiatively
  isolated (no conduction path to a mount). Every watt therefore leaves by
  radiation from the outer envelope, which makes the assembly temperature
  a near-pure radiation balance:  eps*sigma*A*(T^4 - T_amb^4) = Q.
- **Iron is effectively isothermal.** Conduction through A36 at these
  power levels costs ~2 K end to end (k = 50 W/m-K, ~5 cm path, a few
  watts), which is far below the uncertainty in emissivity, so the iron
  and everything bolted to it is one node.
- **The coils are the only heat source**, and the only place a real
  gradient can form: in vacuum a coil that is not well bonded to iron can
  only radiate across its gap. That interface is modelled explicitly
  because it is the binding constraint, and it is swept because its
  conductance is genuinely unknown.
- **PTFE has no source term** in this configuration, so it floats at the
  iron temperature. It still matters, because it has the lowest limit of
  any material in the build (~260 C continuous) and therefore sets the
  current ceiling.

The single most important input is EMISSIVITY of the outer surface, which
is why it is an explicit parameter everywhere rather than a constant
buried in a formula: it moves the answer by ~160 C at high current, more
than any geometric detail this model could resolve. Bare machined steel is
~0.2-0.3; oxidised, blackened or anodised is ~0.7-0.8. If the real
hardware's finish is unknown, read the eps sweep at the bottom, not the
single design-point number.

Usage:
    python radiation_balance.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

SIGMA = 5.670374419e-8      # Stefan-Boltzmann, W/m^2-K^4
RHO_CU_20 = 1.68e-8         # copper resistivity at 20 C, ohm-m
ALPHA_CU = 0.00393          # copper temperature coefficient, 1/K

T_AMB = 300.0               # K -- vacuum chamber wall temperature
EPS_ENVELOPE = 0.30         # outer surface emissivity; SEE MODULE DOCSTRING
EPS_GAP = 0.80              # coil <-> iron surfaces, for the radiation-only case

# Material limits (continuous service), deg C
LIMITS = {
    "PTFE chamber": 260.0,
    "fiberglass wire insulation": 200.0,   # conservative; some serve rates far higher
    "ASTM A36 service limit": 400.0,       # makeitfrom.com max operating temp
    "ASTM A36 (Curie)": 770.0,             # magnetic circuit fails well before this
}

# ASTM A36 structural steel, the pole/plate material (SendCutSend hot-rolled
# pickled-and-oiled -- the only grade they stock thick enough for the
# 0.375in plates). Source: makeitfrom.com ASTM A36 (SS400, S275).
IRON_K = 50.0          # W/m-K thermal conductivity
IRON_CP = 470.0        # J/kg-K specific heat
IRON_RHO = 7900.0      # kg/m^3 density
# Not used by this steady-state model -- CP and RHO are here for whenever a
# transient (warm-up time) version is wanted, which is the obvious next step
# if anyone cares how long a test run takes to reach these temperatures.


@dataclass
class Coil:
    """One winding. Geometry from the real CAD coil solids; turns/current
    from build_assembly.py's EXCITATION -- keep those in sync by hand, the
    magnetics pipeline needs gmsh to import and this script deliberately
    does not."""
    name: str
    count: int                 # how many of this coil in the assembly
    turns: int
    wire_dia: float            # m, BARE copper diameter of the chosen gauge
    r_in: float                # m, winding inner radius
    r_out: float               # m, winding outer radius
    height: float              # m, axial

    @property
    def window(self) -> float:
        """Winding cross-section the current crosses: the r-z rectangle.
        Same quantity as the magnetics pipeline's A_cross -- NOT the
        annulus (see gmsh_getdp/HOW_IT_WORKS.md sec. 2)."""
        return (self.r_out - self.r_in) * self.height

    @property
    def turn_length(self) -> float:
        return 2 * math.pi * (self.r_in + self.r_out) / 2

    @property
    def wire_area(self) -> float:
        return math.pi * (self.wire_dia / 2) ** 2

    @property
    def packing(self) -> float:
        """Copper fraction of the window -- DERIVED, not an input. Serves as
        a fit check: above ~0.75 the turns cannot physically be wound."""
        return self.turns * self.wire_area / self.window

    def resistance(self, T_C: float) -> float:
        """R = rho*N*L_turn/A_wire, from the actual conductor.

        Deliberately NOT the rho*N^2*L/(k*A_window) form: that one assumes
        the winding FILLS its window, and the outer coils here are wound to
        165 of a possible 364 turns so one series current gives the right
        ampere-turn ratio. Using a full-wind packing factor on a
        part-wound coil understated its resistance by 2.2x.
        """
        rho = RHO_CU_20 * (1 + ALPHA_CU * (T_C - 20.0))
        return rho * self.turns * self.turn_length / self.wire_area

    def dissipation(self, current: float, T_C: float) -> float:
        """Watts for ONE coil of this type."""
        return current ** 2 * self.resistance(T_C)

    @property
    def bore_area(self) -> float:
        """Inner cylindrical face -- the conduction path to the core."""
        return 2 * math.pi * self.r_in * self.height

    @property
    def radiating_area(self) -> float:
        """Both cylindrical faces, for the not-bonded case."""
        return 2 * math.pi * (self.r_in + self.r_out) * self.height


# The committed design point: AWG 20 fiberglass-served, inner wound full,
# outer wound to the ampere-turn ratio the field wants. See
# gmsh_getdp/README.md "Solenoid winding design" for the gauge trade.
AWG20_BARE = 0.812e-3
INNER = Coil("inner", 1, 260, AWG20_BARE, 0.0127, 0.0174625, 0.0508)
OUTER = Coil("outer", 4, 165, AWG20_BARE, 0.0127, 0.019050, 0.0508)
COILS = [INNER, OUTER]
I_DESIGN = 1.30

for _c in COILS:
    assert _c.packing < 0.75, (
        f"{_c.name} coil: {_c.turns} turns of {_c.wire_dia*1e3:.3f} mm wire implies "
        f"{_c.packing:.0%} copper fill, which will not wind"
    )


def envelope_area(msh_path: str = "../gmsh_getdp/assembly/assembly.msh") -> tuple[float, str]:
    """Radiating area of the hardware's outer envelope, taken from the real
    mesh rather than a hand-typed cylinder.

    Uses the convex hull of the iron and winding nodes. The hull is the
    right model here because internal cavities (the channel, the gaps
    around the coils) radiate to each other and net out close to zero --
    only the outward-facing surface loses heat to the chamber walls. It is
    a mild underestimate wherever the real surface is concave, which is
    conservative: less area means a hotter predicted assembly.

    Falls back to a bounding cylinder if the mesh or gmsh isn't available,
    so the script still runs standalone.
    """
    try:
        import gmsh
        from scipy.spatial import ConvexHull
    except ImportError as exc:  # pragma: no cover
        return _fallback_area(f"import failed ({exc})")

    try:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(msh_path)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        pts = np.array(coords).reshape(-1, 3)
        # Physical groups: 2 = Iron, 3..7 = Winding0..4 (assembly_params.txt)
        keep = []
        for tag in (2, 3, 4, 5, 6, 7):
            for ent in gmsh.model.getEntitiesForPhysicalGroup(3, tag):
                ntags, ncoords, _ = gmsh.model.mesh.getNodes(3, ent, includeBoundary=True)
                keep.append(np.array(ncoords).reshape(-1, 3))
        gmsh.finalize()
        solid = np.vstack(keep)
        hull = ConvexHull(solid)
        return float(hull.area), f"convex hull of {len(solid)} iron+winding nodes ({msh_path})"
    except Exception as exc:  # pragma: no cover
        try:
            gmsh.finalize()
        except Exception:
            pass
        return _fallback_area(f"mesh read failed ({exc})")


def _fallback_area(why: str) -> tuple[float, str]:
    r, h = 0.0764, 0.0691
    return 2 * math.pi * r * (r + h), f"bounding cylinder r={r*1e3:.1f} h={h*1e3:.1f} mm -- {why}"


def total_dissipation(current: float, T_coil_C: float, coils: list[Coil] | None = None) -> float:
    coils = COILS if coils is None else coils
    return sum(c.count * c.dissipation(current, T_coil_C) for c in coils)


def solve_steady(current: float, area: float, eps: float = EPS_ENVELOPE,
                 coil_interface: str = "contact", h_contact: float = 500.0,
                 t_amb: float = T_AMB, coils: list[Coil] | None = None) -> dict:
    """Steady-state temperatures.

    The iron node sheds all the heat by radiation; each coil sits above it
    by its own dissipation divided by its interface conductance. Solved by
    fixed-point iteration because both the T^4 radiation law and copper's
    R(T) are nonlinear, and they push in opposite directions (hotter iron
    radiates harder, hotter copper dissipates more).

    coil_interface:
      "contact"   -- coil bonded/clamped to the core, G = h_contact * bore area.
                     h_contact ~ 500-5000 W/m^2-K for clamped metal in vacuum.
      "radiation" -- worst case, coil thermally floating in its pocket, coupled
                     only by radiation across the gap.

    coils: defaults to the committed design (COILS); pass another set to
    evaluate a different winding, e.g. from sizing/winding_design.py.
    """
    COILS_ = COILS if coils is None else coils
    T_iron = t_amb
    T_coil = {c.name: t_amb for c in COILS_}

    for _ in range(200):
        Q = sum(c.count * c.dissipation(current, T_coil[c.name] - 273.15) for c in COILS_)
        T_iron_new = (Q / (eps * SIGMA * area) + t_amb ** 4) ** 0.25
        for c in COILS_:
            if coil_interface == "contact":
                G = h_contact * c.bore_area
            else:
                G = 4 * EPS_GAP * SIGMA * T_iron_new ** 3 * c.radiating_area
            T_coil[c.name] = T_iron_new + c.dissipation(current, T_coil[c.name] - 273.15) / G
        if abs(T_iron_new - T_iron) < 1e-9:
            T_iron = T_iron_new
            break
        T_iron = T_iron_new

    return {
        "current": current,
        "Q": sum(c.count * c.dissipation(current, T_coil[c.name] - 273.15) for c in COILS_),
        "T_iron_C": T_iron - 273.15,
        "T_ptfe_C": T_iron - 273.15,          # no source term; floats at iron
        "T_coil_C": {k: v - 273.15 for k, v in T_coil.items()},
        "T_hot_C": max(T_coil.values()) - 273.15,
    }


def _fmt_row(r: dict) -> str:
    return (f"{r['current']:6.2f} {r['Q']:8.2f} {r['T_iron_C']:9.1f} "
            f"{r['T_coil_C']['inner']:9.1f} {r['T_coil_C']['outer']:9.1f} {r['T_ptfe_C']:9.1f}")


def main() -> None:
    area, how = envelope_area()

    print("=" * 78)
    print("BPL-700 magnet circuit -- steady-state thermal, magnets only, vacuum")
    print("=" * 78)
    print(f"\nradiating envelope : {area*1e4:.1f} cm^2")
    print(f"  source           : {how}")
    print(f"ambient            : {T_AMB:.0f} K ({T_AMB-273.15:.0f} C)")
    print(f"envelope emissivity: {EPS_ENVELOPE:.2f}   <-- dominant assumption, see sweep below")
    print(f"design current     : {I_DESIGN:.2f} A")
    for c in COILS:
        print(f"  {c.name:6s} x{c.count}  {c.turns:3d} turns of {c.wire_dia*1e3:.3f} mm  "
              f"window {c.window*1e6:6.1f} mm^2  fill={c.packing:.0%}  "
              f"R(20C)={c.resistance(20):.3f} ohm")
    R_series = sum(c.count * c.resistance(20) for c in COILS)
    print(f"  series total       : {R_series:.3f} ohm -> "
          f"{I_DESIGN*R_series:.2f} V at {I_DESIGN:.2f} A")

    print("\n--- design point ---")
    for mode, label in (("contact", "coil bonded to core (h=500 W/m^2-K)"),
                        ("radiation", "coil floating, radiation-coupled only")):
        r = solve_steady(I_DESIGN, area, coil_interface=mode)
        print(f"  {label}")
        print(f"    Q = {r['Q']:.2f} W   iron/PTFE {r['T_iron_C']:.1f} C   "
              f"inner coil {r['T_coil_C']['inner']:.1f} C   outer coil {r['T_coil_C']['outer']:.1f} C")

    print("\n--- current margin (contact interface) ---")
    print(f"{'I (A)':>6} {'Q (W)':>8} {'iron C':>9} {'inner C':>9} {'outer C':>9} {'PTFE C':>9}")
    rows = [solve_steady(I, area) for I in (1.0, I_DESIGN, 2.0, 3.0, 5.0, 7.0, 10.0)]
    for r in rows:
        print(_fmt_row(r))

    print("\n--- current at which each limit is reached (contact interface) ---")
    for name, limit in LIMITS.items():
        lo, hi = 0.1, 60.0
        for _ in range(80):
            mid = (lo + hi) / 2
            if solve_steady(mid, area)["T_hot_C"] < limit:
                lo = mid
            else:
                hi = mid
        print(f"  {name:28s} {limit:6.0f} C  ->  {lo:5.2f} A  ({lo/I_DESIGN:4.1f}x design)")

    print("\n--- sensitivity to envelope emissivity (at design current) ---")
    print(f"{'eps':>6} {'iron C':>9} {'hottest C':>11}   note")
    for eps, note in ((0.15, "polished steel"), (0.30, "bare machined steel"),
                      (0.50, "lightly oxidised"), (0.80, "blackened / anodised")):
        r = solve_steady(I_DESIGN, area, eps=eps)
        print(f"{eps:6.2f} {r['T_iron_C']:9.1f} {r['T_hot_C']:11.1f}   {note}")

    print("\n--- sensitivity to the coil-to-iron interface (at design current) ---")
    print(f"{'h (W/m^2-K)':>12} {'inner coil C':>14} {'rise above iron K':>19}")
    base = solve_steady(I_DESIGN, area)["T_iron_C"]
    for h in (100.0, 500.0, 2000.0, 5000.0):
        r = solve_steady(I_DESIGN, area, h_contact=h)
        print(f"{h:12.0f} {r['T_coil_C']['inner']:14.1f} {r['T_coil_C']['inner']-base:19.1f}")
    r = solve_steady(I_DESIGN, area, coil_interface="radiation")
    print(f"{'radiation':>12} {r['T_coil_C']['inner']:14.1f} {r['T_coil_C']['inner']-base:19.1f}")


if __name__ == "__main__":
    main()
