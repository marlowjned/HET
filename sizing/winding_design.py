"""
winding_design.py

Turns an ampere-turn requirement plus a coil window into a real winding:
turns that physically fit, channel currents, resistance, voltage, power,
wire length and steady-state temperature. The magnetics solve only ever
sees the A-turn product (gmsh_getdp/README.md "Solenoid winding design"),
so everything here is downstream of it and none of it feeds back into the
field.

What changed on 2026-09-30, and why the old numbers are not in here
--------------------------------------------------------------------
This module used to size AWG 20 "single-glass-served" wire from estimated
insulated diameters. Two things were wrong with that:

1. The wire did not exist. NEMA MW 1000 tabulates round polyester-glass
   served wire from 1/0 down to AWG 12 and no finer, no US retailer stocks
   it in any gauge, and the round product is Class 155 rather than the
   200 C the thermal model assumed. See sizing/wire_data.py.
2. The turn count was being computed to 0.04 mm on a +/-10% estimate, with
   no allowance for ground insulation, and it depended on the wire
   hexagonally nesting -- which glass serving, being rough and
   dimensionally variable, is the least likely wire type to do. Square
   packing took the old design point from 260 turns to 208, a 20% swing on
   an assumption.

So: dimensions now come from NEMA MW 1000 maxima, ground insulation is an
explicit allowance, and every result is reported for BOTH packing models
rather than silently assuming the optimistic one.

Two-channel vs series
---------------------
The inner and outer coils need different ampere-turns (~338 and ~214.5 at
the baseline geometry, a ratio of 1.58). Two ways to deliver that:

- TWO CHANNELS (specified): independent supply currents, so both coils are
  wound FULL and each channel is set to req_At/N. Turn-count uncertainty
  stops being a design risk -- you count what you wound and dial the
  current. Both coils full also means P ~ 1/N is minimised, and the FEM's
  uniform-current-density-over-the-whole-window assumption becomes correct
  instead of being 1.8x off on the part-wound outer coils.
- SERIES (costed for comparison): one current, so the outer coils must be
  part-wound to hit the ratio. Cheaper supply, ~1.7x the dissipation, and
  the ratio is frozen in copper at winding time even though it came from a
  linear solve with a proxy B-H curve.

Winding model: turns per layer = floor(height / d_ins); layers =
floor(build / d_ins) for square packing, or floor((build - d_ins)/pitch)+1
with pitch = d_ins*sqrt(3)/2 for hexagonal nesting. Mean turn length uses
the radius actually wound, not the window's mid-radius, so a part-wound
coil is not credited with copper it does not have.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "thermal"))
from wire_data import (GROUND_INSULATION_MM, OPTIONS, POLYIMIDE,  # noqa: E402
                       Insulation)
from radiation_balance import Coil, solve_steady, RHO_CU_20, ALPHA_CU  # noqa: E402

IN = 0.0254

# The COMMITTED baseline excitation -- what gmsh_getdp/assembly/
# build_assembly.py's EXCITATION actually delivers, and what every solved
# field in this repo was solved at. It is deliberately NOT the 300 G-exact
# figure: at the baseline geometry it produces 306.5 G at the channel exit,
# about 2% above the 300 G target, and the exactly-scaled requirement is
# emag_sweep_results.csv's req_inner_At = 330.8 / req_outer_At = 209.9.
#
# These defaults are here so this module reproduces the committed design
# point standalone. Per-geometry work should pass req_inner_At/req_outer_At
# from the sweep CSV instead, which is what report/make_report.py does.
REQ_INNER_AT = 338.0
REQ_OUTER_AT = 214.5

PACKINGS = ("hex", "square")


def turns_fit(build_m: float, height_m: float, d_ins_m: float,
              packing: str = "square") -> tuple[int, int, int, float]:
    """(max turns, turns per layer, layers, radial extent of the copper).

    `packing`:
      "hex"    -- each layer settles into the valleys of the one below,
                  radial pitch d_ins*sqrt(3)/2. The best case, and it needs
                  a smooth, dimensionally consistent wire to be real.
      "square" -- layers stack directly on top of each other, pitch d_ins.
                  The conservative case, and what a hand-wound coil of
                  rough wire actually tends to do.
    """
    if packing not in PACKINGS:
        raise ValueError(f"packing must be one of {PACKINGS}, got {packing!r}")
    if build_m < d_ins_m or height_m < d_ins_m:
        return 0, 0, 0, 0.0
    per_layer = int(math.floor(height_m / d_ins_m + 1e-9))
    if packing == "hex":
        pitch = d_ins_m * math.sqrt(3) / 2
        layers = int(math.floor((build_m - d_ins_m) / pitch + 1e-9)) + 1
        extent = d_ins_m * (1 + (layers - 1) * math.sqrt(3) / 2)
    else:
        layers = int(math.floor(build_m / d_ins_m + 1e-9))
        extent = d_ins_m * layers
    return per_layer * layers, per_layer, layers, extent


def layer_margin(build_m: float, d_ins_m: float, packing: str) -> float:
    """Unused build left over after the last layer, m.

    Worth looking at: a design whose margin is a few hundredths of a
    millimetre is one datasheet revision away from losing a whole layer.
    The old AWG 20 glass-served design point missed its 5th square-packed
    layer by 0.037 mm, which is what motivated reporting this at all.
    """
    _, _, layers, extent = turns_fit(build_m, 1.0, d_ins_m, packing)
    if not layers:
        return build_m
    # Clamped at zero: turns_fit's floor() carries a 1e-9 tolerance so that a
    # layer landing exactly on the build counts, which can leave this a few
    # 1e-16 negative. A margin of zero is the honest reading -- it means the
    # last layer exactly fills the window and any real-world tolerance loses it.
    return max(0.0, build_m - extent)


@dataclass
class CoilWindow:
    """The space available to a winding, from the real CAD coil solids.

    r_core is the BARE core radius; the first turn sits a ground-insulation
    allowance outside it. The inner window cannot grow outward -- the
    chamber bore sits on it -- only inward by thinning the core, which is
    what the emag_height sweep trades.
    """
    name: str
    count: int
    r_core: float                                  # m, bare core radius
    r_out: float                                   # m, coil OD / 2
    height: float                                  # m, = emag_height
    ground_mm: float = GROUND_INSULATION_MM

    @property
    def r_in(self) -> float:
        return self.r_core + self.ground_mm * 1e-3

    @property
    def build(self) -> float:
        return self.r_out - self.r_in


@dataclass
class WoundCoil:
    window: CoilWindow
    insulation: Insulation
    awg: int
    turns: int
    layers_used: int
    packing: str

    @property
    def d_bare(self) -> float:
        return self.insulation.bare_mm(self.awg) * 1e-3

    @property
    def d_ins(self) -> float:
        return self.insulation.od_mm(self.awg) * 1e-3

    @property
    def r_wound(self) -> float:
        """Outer radius of the copper actually wound."""
        if self.packing == "hex":
            ext = self.d_ins * (1 + (self.layers_used - 1) * math.sqrt(3) / 2)
        else:
            ext = self.d_ins * self.layers_used
        return self.window.r_in + ext

    @property
    def turn_length(self) -> float:
        return 2 * math.pi * (self.window.r_in + self.r_wound) / 2

    @property
    def wire_length(self) -> float:
        return self.turns * self.turn_length

    def resistance(self, T_C: float) -> float:
        rho = RHO_CU_20 * (1 + ALPHA_CU * (T_C - 20.0))
        return rho * self.wire_length / (math.pi * (self.d_bare / 2) ** 2)

    def as_thermal_coil(self, current_ratio: float = 1.0) -> Coil:
        return Coil(self.window.name, self.window.count, self.turns, self.d_bare,
                    self.window.r_in, self.r_wound, self.window.height,
                    current_ratio=current_ratio)


def _wind_full(window: CoilWindow, ins: Insulation, awg: int, packing: str):
    n, per_layer, layers, _ = turns_fit(window.build, window.height,
                                        ins.od_mm(awg) * 1e-3, packing)
    if n == 0:
        return None
    return WoundCoil(window, ins, awg, n, layers, packing), per_layer


def design_two_channel(inner: CoilWindow, outer: CoilWindow, awg: int,
                       ins: Insulation = POLYIMIDE, packing: str = "square",
                       req_inner_At: float = REQ_INNER_AT,
                       req_outer_At: float = REQ_OUTER_AT,
                       envelope_area_m2: float | None = None) -> dict | None:
    """Both coils wound full, each on its own constant-current channel.

    Returns None if the gauge does not fit the window at all.
    """
    wi = _wind_full(inner, ins, awg, packing)
    wo = _wind_full(outer, ins, awg, packing)
    if wi is None or wo is None:
        return None
    (wi, _), (wo, _) = wi, wo

    I_in, I_out = req_inner_At / wi.turns, req_outer_At / wo.turns
    R_in20, R_out20 = wi.resistance(20), outer.count * wo.resistance(20)

    out = dict(
        topology="two-channel", awg=awg, insulation=ins.key, packing=packing,
        inner_turns=wi.turns, inner_layers=wi.layers_used,
        outer_turns=wo.turns, outer_layers=wo.layers_used,
        current_inner_A=I_in, current_outer_A=I_out,
        inner_wire_m=wi.wire_length, outer_wire_m_each=wo.wire_length,
        total_wire_m=wi.wire_length + outer.count * wo.wire_length,
        R20_inner=R_in20, R20_outer=R_out20,
        V20_inner=I_in * R_in20, V20_outer=I_out * R_out20,
        P20_W=I_in ** 2 * R_in20 + I_out ** 2 * R_out20,
        J_inner=I_in / (math.pi * (wi.d_bare * 1e3 / 2) ** 2),
        J_outer=I_out / (math.pi * (wo.d_bare * 1e3 / 2) ** 2),
        inner_margin_mm=layer_margin(inner.build, wi.d_ins, packing) * 1e3,
        outer_margin_mm=layer_margin(outer.build, wo.d_ins, packing) * 1e3,
    )
    _add_thermal(out, wi, wo, I_in, I_out / I_in, outer.count, envelope_area_m2)
    return out


def design_series(inner: CoilWindow, outer: CoilWindow, awg: int,
                  ins: Insulation = POLYIMIDE, packing: str = "square",
                  req_inner_At: float = REQ_INNER_AT,
                  req_outer_At: float = REQ_OUTER_AT,
                  envelope_area_m2: float | None = None) -> dict | None:
    """One series string: the outer coils are part-wound to hit the A-turn
    ratio at the single shared current. Kept to cost the alternative."""
    wi_ = _wind_full(inner, ins, awg, packing)
    wo_ = _wind_full(outer, ins, awg, packing)
    if wi_ is None or wo_ is None:
        return None
    (wi_full, pl_in), (wo_full, pl_out) = wi_, wo_
    n_in_max, n_out_max = wi_full.turns, wo_full.turns

    # Lowest current at which both windows can supply their A-turns, then
    # set the current by the inner coil (it drives the channel field most
    # directly) and report what the outer actually delivers.
    I = max(req_inner_At / n_in_max, req_outer_At / n_out_max)
    n_in = min(n_in_max, math.ceil(req_inner_At / I - 1e-9))
    n_out = min(n_out_max, math.ceil(req_outer_At / I - 1e-9))
    I = req_inner_At / n_in

    wi = WoundCoil(inner, ins, awg, n_in, math.ceil(n_in / pl_in), packing)
    wo = WoundCoil(outer, ins, awg, n_out, math.ceil(n_out / pl_out), packing)
    R20 = wi.resistance(20) + outer.count * wo.resistance(20)

    out = dict(
        topology="series", awg=awg, insulation=ins.key, packing=packing,
        inner_turns=n_in, inner_turns_max=n_in_max, inner_layers=wi.layers_used,
        outer_turns=n_out, outer_turns_max=n_out_max, outer_layers=wo.layers_used,
        outer_fill_frac=n_out / n_out_max,
        current_inner_A=I, current_outer_A=I,
        outer_At_delivered=n_out * I,
        outer_At_error_pct=(n_out * I / req_outer_At - 1) * 100,
        inner_wire_m=wi.wire_length, outer_wire_m_each=wo.wire_length,
        total_wire_m=wi.wire_length + outer.count * wo.wire_length,
        R20_ohm=R20, V20=I * R20, P20_W=I * I * R20,
        J_inner=I / (math.pi * (wi.d_bare * 1e3 / 2) ** 2),
    )
    _add_thermal(out, wi, wo, I, 1.0, outer.count, envelope_area_m2)
    return out


def _add_thermal(out: dict, wi: WoundCoil, wo: WoundCoil, I_ref: float,
                 ratio: float, outer_count: int,
                 envelope_area_m2: float | None) -> None:
    if envelope_area_m2 is None:
        return
    coils = [wi.as_thermal_coil(), wo.as_thermal_coil(ratio)]
    th = solve_steady(I_ref, envelope_area_m2, coils=coils)
    R_in_hot = wi.resistance(th["T_coil_C"]["inner"])
    R_out_hot = outer_count * wo.resistance(th["T_coil_C"]["outer"])
    out.update(
        T_iron_C=th["T_iron_C"], T_coil_hot_C=th["T_hot_C"], P_hot_W=th["Q"],
        T_inner_C=th["T_coil_C"]["inner"], T_outer_C=th["T_coil_C"]["outer"],
        R_hot_inner=R_in_hot, R_hot_outer=R_out_hot,
        V_hot_inner=I_ref * R_in_hot, V_hot_outer=I_ref * ratio * R_out_hot,
    )


def windows_for(emag_height_in: float, inner_core_dia_in: float,
                outer_coil_od_in: float = 1.500,
                inner_coil_od_in: float = 1.375,
                outer_core_dia_in: float = 1.0,
                ground_mm: float = GROUND_INSULATION_MM):
    """The two coil windows for one point of the emag_height sweep."""
    h = emag_height_in * IN
    inner = CoilWindow("inner", 1, inner_core_dia_in * IN / 2,
                       inner_coil_od_in * IN / 2, h, ground_mm)
    outer = CoilWindow("outer", 4, outer_core_dia_in * IN / 2,
                       outer_coil_od_in * IN / 2, h, ground_mm)
    return inner, outer


def trade_table(inner: CoilWindow, outer: CoilWindow, ins: Insulation = POLYIMIDE,
                gauges=None, **kw) -> list[dict]:
    """Every gauge x packing combination that fits, two-channel."""
    gauges = ins.gauges if gauges is None else gauges
    rows = []
    for awg in gauges:
        for packing in PACKINGS:
            r = design_two_channel(inner, outer, awg, ins, packing, **kw)
            if r:
                rows.append(r)
    return rows


def _self_check() -> None:
    """Reproduce the committed design point (thermal/radiation_balance.py
    INNER/OUTER and gmsh_getdp/assembly/build_assembly.py EXCITATION)."""
    inner, outer = windows_for(2.0, 1.0)
    assert abs(inner.build - 4.6125e-3) < 1e-6, inner.build
    assert abs(outer.build - 6.2000e-3) < 1e-6, outer.build

    # The inner coil holds 5 layers of 56 under BOTH packing models -- the
    # property that specifying a film-insulated wire bought us.
    for packing in PACKINGS:
        r = design_two_channel(inner, outer, 20, POLYIMIDE, packing)
        assert r["inner_turns"] == 280, (packing, r["inner_turns"])
        assert abs(r["current_inner_A"] - 1.2071) < 5e-4, r["current_inner_A"]

    sq = design_two_channel(inner, outer, 20, POLYIMIDE, "square")
    hx = design_two_channel(inner, outer, 20, POLYIMIDE, "hex")
    assert (sq["outer_turns"], hx["outer_turns"]) == (336, 392)
    assert abs(sq["current_outer_A"] - 0.6384) < 5e-4, sq["current_outer_A"]
    assert abs(hx["current_outer_A"] - 0.5472) < 5e-4, hx["current_outer_A"]

    # Two channels must beat the series string on dissipation, because both
    # coils get wound full and P ~ 1/N. At 20 C the advantage is 1.46x with
    # square packing and 1.57x with hex -- larger with hex because the outer
    # window holds a 7th layer that only the two-channel scheme can fill.
    ser_sq = design_series(inner, outer, 20, POLYIMIDE, "square")
    ser_hx = design_series(inner, outer, 20, POLYIMIDE, "hex")
    assert abs(ser_sq["P20_W"] / sq["P20_W"] - 1.46) < 0.02, ser_sq["P20_W"] / sq["P20_W"]
    assert abs(ser_hx["P20_W"] / hx["P20_W"] - 1.57) < 0.02, ser_hx["P20_W"] / hx["P20_W"]


_self_check()


if __name__ == "__main__":
    inner, outer = windows_for(2.0, 1.0)
    print(f"Baseline geometry: emag_height 2.000 in, inner core 1.000 in")
    print(f"  inner window  build {inner.build*1e3:.3f} mm x {inner.height*1e3:.1f} mm "
          f"(r {inner.r_in*1e3:.2f}-{inner.r_out*1e3:.2f} mm)")
    print(f"  outer window  build {outer.build*1e3:.3f} mm x {outer.height*1e3:.1f} mm "
          f"(r {outer.r_in*1e3:.2f}-{outer.r_out*1e3:.2f} mm)")
    print(f"  ground insulation allowance {GROUND_INSULATION_MM} mm on the bare core")

    for key in ("polyimide", "glass"):
        ins = OPTIONS[key]
        print(f"\n{'='*94}")
        print(f"{ins.label}  --  {ins.nema}, {ins.temp_C:.0f} C, {ins.availability}")
        print('=' * 94)
        gauges = [g for g in ins.gauges if 12 <= g <= 24]
        hdr = (f"{'AWG':>4} {'pack':>6} | {'N_in':>5} {'I_in':>6} {'N_out':>5} {'I_out':>6} "
               f"| {'P(W)':>5} {'V_in':>5} {'V_out':>6} | {'wire m':>7} "
               f"{'margin in/out (mm)':>19}")
        print(hdr, '\n' + '-' * len(hdr))
        for awg in gauges:
            for packing in PACKINGS:
                r = design_two_channel(inner, outer, awg, ins, packing)
                if not r:
                    continue
                print(f"{awg:>4} {packing:>6} | {r['inner_turns']:>5} {r['current_inner_A']:>6.3f} "
                      f"{r['outer_turns']:>5} {r['current_outer_A']:>6.3f} | {r['P20_W']:>5.2f} "
                      f"{r['V20_inner']:>5.2f} {r['V20_outer']:>6.2f} | {r['total_wire_m']:>7.1f} "
                      f"{r['inner_margin_mm']:>9.3f} {r['outer_margin_mm']:>9.3f}")

    print(f"\n{'='*94}\nSeries string, for comparison (AWG 20 polyimide, square)\n{'='*94}")
    s = design_series(inner, outer, 20, POLYIMIDE, "square")
    t = design_two_channel(inner, outer, 20, POLYIMIDE, "square")
    print(f"  series      : inner {s['inner_turns']}, outer part-wound to {s['outer_turns']}"
          f" of {s['outer_turns_max']} ({s['outer_fill_frac']:.0%} of window), "
          f"{s['current_inner_A']:.3f} A, {s['V20']:.2f} V, {s['P20_W']:.2f} W")
    print(f"  two-channel : inner {t['inner_turns']}, outer {t['outer_turns']} (full), "
          f"{t['current_inner_A']:.3f}/{t['current_outer_A']:.3f} A, "
          f"{t['V20_inner']:.2f}/{t['V20_outer']:.2f} V, {t['P20_W']:.2f} W")
    print(f"  -> series costs {s['P20_W']/t['P20_W']:.2f}x the dissipation and freezes the "
          f"A-turn ratio in copper")
