"""
winding_design.py

Turns an ampere-turn requirement plus a coil window into a real winding:
turns that physically fit, series current, resistance, voltage, power,
wire length and steady-state temperature. The magnetics solve only ever
sees the A-turn product (gmsh_getdp/README.md "Solenoid winding design"),
so everything here is downstream of it and none of it feeds back into the
field.

Reproduces the committed design as a check: AWG 20 in the 2.0in-tall
inner window gives 52 turns/layer x 5 layers = 260, and the outer window
holds 364 -- the numbers the README and thermal/radiation_balance.py use.

Winding model: hexagonally nested round wire. Turns per layer =
floor(height / d_ins); layer pitch d_ins*sqrt(3)/2; layers =
floor((build - d_ins) / pitch) + 1. The inner coil is wound full; the outer
coils are wound to whatever turn count gives the A-turn ratio the field
wants at the SAME series current, unless the outer window is the one that
runs out first, in which case the outer coil is full and the inner is
part-wound. Mean turn length uses the radius actually wound (a part-wound
coil's copper sits at the bore), not the window's mid-radius.

Insulated diameters are single-glass-served ESTIMATES (MATERIALS.md,
"Fiberglass-served magnet wire") -- the softest input here. A real
datasheet moves turn counts by 10-20%, which moves current and voltage,
not the field and barely the power.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "thermal"))
from radiation_balance import Coil, solve_steady, RHO_CU_20, ALPHA_CU  # noqa: E402

IN = 0.0254

# AWG: (bare copper dia, single-glass-served OD), mm. MATERIALS.md.
WIRE = {
    16: (1.291, 1.45),
    18: (1.024, 1.18),
    20: (0.812, 0.96),
    22: (0.644, 0.79),
    24: (0.511, 0.65),
}


def turns_fit(build_m: float, height_m: float, d_ins_m: float) -> tuple[int, int, int]:
    """(max turns, turns per layer, layers) for a hexagonally nested winding."""
    per_layer = int(math.floor(height_m / d_ins_m + 1e-9))
    pitch = d_ins_m * math.sqrt(3) / 2
    layers = int(math.floor((build_m - d_ins_m) / pitch + 1e-9)) + 1 if build_m >= d_ins_m else 0
    return per_layer * layers, per_layer, layers


@dataclass
class CoilWindow:
    name: str
    count: int
    r_in: float      # m (= core radius)
    r_out: float     # m (= coil OD / 2)
    height: float    # m (= emag_height)

    @property
    def build(self) -> float:
        return self.r_out - self.r_in


@dataclass
class WoundCoil:
    window: CoilWindow
    awg: int
    turns: int
    layers_used: int

    @property
    def d_bare(self) -> float:
        return WIRE[self.awg][0] * 1e-3

    @property
    def d_ins(self) -> float:
        return WIRE[self.awg][1] * 1e-3

    @property
    def r_wound(self) -> float:
        """Outer radius of the copper actually wound."""
        return self.window.r_in + self.d_ins * (1 + (self.layers_used - 1) * math.sqrt(3) / 2)

    @property
    def turn_length(self) -> float:
        return 2 * math.pi * (self.window.r_in + self.r_wound) / 2

    @property
    def wire_length(self) -> float:
        return self.turns * self.turn_length

    def resistance(self, T_C: float) -> float:
        rho = RHO_CU_20 * (1 + ALPHA_CU * (T_C - 20.0))
        return rho * self.wire_length / (math.pi * (self.d_bare / 2) ** 2)

    def as_thermal_coil(self) -> Coil:
        return Coil(self.window.name, self.window.count, self.turns, self.d_bare,
                    self.window.r_in, self.r_wound, self.window.height)


def design(req_inner_At: float, req_outer_At: float, inner: CoilWindow, outer: CoilWindow,
           awg: int, envelope_area_m2: float | None = None) -> dict | None:
    """One series-wired winding set meeting both A-turn requirements at a
    common current. Returns None if the gauge doesn't fit at all."""
    d = WIRE[awg][1] * 1e-3
    n_in_max, pl_in, _ = turns_fit(inner.build, inner.height, d)
    n_out_max, pl_out, _ = turns_fit(outer.build, outer.height, d)
    if n_in_max == 0 or n_out_max == 0:
        return None

    # Lowest current at which both windows can supply their A-turns.
    I = max(req_inner_At / n_in_max, req_outer_At / n_out_max)
    n_in = min(n_in_max, math.ceil(req_inner_At / I - 1e-9))
    n_out = min(n_out_max, math.ceil(req_outer_At / I - 1e-9))
    # Integer turns overshoot slightly; set the current by the inner coil
    # (it sets the channel field most directly) and report the outer A-turns
    # actually delivered.
    I = req_inner_At / n_in

    wi = WoundCoil(inner, awg, n_in, math.ceil(n_in / pl_in))
    wo = WoundCoil(outer, awg, n_out, math.ceil(n_out / pl_out))
    R20 = wi.resistance(20) + outer.count * wo.resistance(20)

    out = dict(
        awg=awg, current_A=I,
        inner_turns=n_in, inner_turns_max=n_in_max, inner_layers=wi.layers_used,
        outer_turns=n_out, outer_turns_max=n_out_max, outer_layers=wo.layers_used,
        outer_At_delivered=n_out * I, outer_At_error_pct=(n_out * I / req_outer_At - 1) * 100,
        inner_wire_m=wi.wire_length, outer_wire_m_each=wo.wire_length,
        total_wire_m=wi.wire_length + outer.count * wo.wire_length,
        R20_ohm=R20, V20=I * R20, P20_W=I * I * R20,
        J_A_per_mm2=I / (math.pi * (wi.d_bare * 1e3 / 2) ** 2),
    )
    if envelope_area_m2 is not None:
        th = solve_steady(I, envelope_area_m2, coils=[wi.as_thermal_coil(), wo.as_thermal_coil()])
        T_hot = th["T_hot_C"]
        R_hot = wi.resistance(th["T_coil_C"]["inner"]) + outer.count * wo.resistance(th["T_coil_C"]["outer"])
        out.update(T_iron_C=th["T_iron_C"], T_coil_hot_C=T_hot,
                   R_hot_ohm=R_hot, V_hot=I * R_hot, P_hot_W=I * I * R_hot)
    return out


def windows_for(emag_height_in: float, inner_core_dia_in: float, outer_coil_od_in: float,
                inner_coil_od_in: float = 1.375, outer_core_dia_in: float = 1.0):
    h = emag_height_in * IN
    inner = CoilWindow("inner", 1, inner_core_dia_in * IN / 2, inner_coil_od_in * IN / 2, h)
    outer = CoilWindow("outer", 4, outer_core_dia_in * IN / 2, outer_coil_od_in * IN / 2, h)
    return inner, outer


if __name__ == "__main__":
    # Self-check against the committed design (README "Solenoid winding design").
    inner, outer = windows_for(2.0, 1.0, 1.5)
    assert turns_fit(inner.build, inner.height, 0.96e-3)[0] == 260
    assert turns_fit(outer.build, outer.height, 0.96e-3)[0] == 364
    r = design(338.0, 214.5, inner, outer, 20)
    print(f"AWG20 at 2.0in: {r['inner_turns']}/{r['outer_turns']} turns, {r['current_A']:.3f} A, "
          f"{r['V20']:.2f} V, {r['P20_W']:.2f} W at 20 C")
    assert r["inner_turns"] == 260 and r["outer_turns"] == 165
