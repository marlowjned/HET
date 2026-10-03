"""
em_winding_trade.py

BPL-1000 winding trade: every swept geometry (gmsh_getdp/assembly/
em_sizing_results.csv) x wire gauge x packing -> turns, channel currents,
voltage, power, temperature, wire length and mass, and EM overall OD.

The FEM only ever sees ampere-turns over the coil window, so this is pure
post-processing of the sweep's required A-turns -- see em_sizing_sweep.py.
For a window wound full, P ~ rho (NI)^2 MLT / (fill * A_window): power is
set by the window and barely by the gauge, which only splits it between
volts and amps (plus layer quantisation). That is why gauge is a column
here and not an axis of the field sweep.

Both excitations the sweep reports are carried: `best` (outer/inner ratio
re-optimised per geometry for field shape) and `ref` (the fixed 214.5/338
ratio of the committed design).

Usage:  python sizing/em_winding_trade.py
"""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "gmsh_getdp" / "assembly"))

from winding_design import CoilWindow, PACKINGS, design_two_channel  # noqa: E402
from wire_data import POLYIMIDE  # noqa: E402
from bpl1000_geometry import IN, coil_heights_in  # noqa: E402

SWEEP_CSV = ROOT / "gmsh_getdp" / "assembly" / "em_sizing_results.csv"
OUT_CSV = ROOT / "sizing" / "em_winding_trade.csv"
GAUGES = (18, 20, 22, 24)
RHO_CU_KG_M3 = 8960.0
N_OUTER = 4


def windows(row):
    h_in, h_out = coil_heights_in(float(row["emag_height"]), float(row["iem_platet"]))
    ic, it = float(row["iem_coreod"]), float(row["iem_coilt"])
    oc, ot = float(row["oem_coreod"]), float(row["oem_coilt"])
    inner = CoilWindow("inner", 1, ic * IN / 2, (ic + it) * IN / 2, h_in * IN)
    outer = CoilWindow("outer", N_OUTER, oc * IN / 2, (oc + ot) * IN / 2, h_out * IN)
    return inner, outer


def main():
    rows = list(csv.DictReader(open(SWEEP_CSV)))
    out = []
    for row in rows:
        inner, outer = windows(row)
        for exc in ("best", "ref"):
            for awg in GAUGES:
                for packing in PACKINGS:
                    r = design_two_channel(
                        inner, outer, awg, POLYIMIDE, packing,
                        req_inner_At=float(row[f"{exc}_req_inner_At"]),
                        req_outer_At=float(row[f"{exc}_req_outer_At"]),
                        envelope_area_m2=float(row["envelope_area_m2"]))
                    if r is None:
                        continue
                    a_cu = math.pi * (POLYIMIDE.bare_mm(awg) * 1e-3 / 2) ** 2
                    out.append(dict(
                        tag=row["tag"], excitation=exc,
                        iem_coreod=row["iem_coreod"], iem_coilt=row["iem_coilt"],
                        iem_platet=row["iem_platet"],
                        oem_coreod=row["oem_coreod"], oem_coilt=row["oem_coilt"],
                        inner_em_od_in=row["inner_em_od_in"], outer_em_od_in=row["outer_em_od_in"],
                        ratio=row[f"{exc}_ratio"],
                        req_inner_At=row[f"{exc}_req_inner_At"], req_outer_At=row[f"{exc}_req_outer_At"],
                        quality_score=row[f"{exc}_quality_score"],
                        inner_shaft_B_max_T=row[f"{exc}_inner_shaft_B_max_T"],
                        outer_shaft_B_max_T=row[f"{exc}_outer_shaft_B_max_T"],
                        linear_valid=row["linear_valid"],
                        inner_build_mm=inner.build * 1e3, outer_build_mm=outer.build * 1e3,
                        **r,
                        inner_wire_kg=r["inner_wire_m"] * a_cu * RHO_CU_KG_M3,
                        outer_wire_kg_total=N_OUTER * r["outer_wire_m_each"] * a_cu * RHO_CU_KG_M3,
                    ))
    keys = []
    for r in out:
        keys += [k for k in r if k not in keys]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    print(f"{len(rows)} geometries -> {len(out)} winding rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()
