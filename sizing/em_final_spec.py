"""
em_final_spec.py

Final BPL-1000 electromagnet specification, and the evidence that it is
insensitive to the things it might be expected to be sensitive to.

Design (chosen 2026-10-05): AWG 20 heavy-build polyimide, 0.75 in cores for
all five magnets (one bar stock), two constant-current channels in the
1-2 A band the user asked for:
  inner: 5 layers x 46 = 230 turns, iem_coilt 0.375 in
  outer: 4 layers x 46 = 184 turns per coil, oem_coilt 0.3125 in
The required ampere-turns come from the combined linear solve of exactly
this geometry (em_sizing_sweep.py point FINAL_TAG), not from the
one-magnet-at-a-time sweeps.

Three sensitivity studies, all post-processing (no field solves -- the FEM
only sees ampere-turns over the window, see em_sizing_sweep.py):
  1. Insulation: overall diameter and power at the specified turn count,
     bracketed by BARE copper (zero film, the physical lower bound) and
     DOUBLE heavy-build film. Bracketing avoids importing single-build
     table values that are not in wire_data.py.
  2. Gauge: AWG 16-28 in the specified windows.
  3. Diameter: layer count (= overall diameter) vs power and temperature
     for each magnet, the other held at its specified winding.

Usage:  python sizing/em_final_spec.py      (prints tables, writes figure)
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "thermal"))
sys.path.insert(0, str(ROOT / "gmsh_getdp" / "assembly"))

from winding_design import CoilWindow, WoundCoil, turns_fit, GROUND_INSULATION_MM  # noqa: E402
from wire_data import POLYIMIDE  # noqa: E402
from radiation_balance import solve_steady  # noqa: E402

IN = 0.0254
FINAL = dict(iem_coreod=0.75, iem_coilt=0.375, oem_coreod=0.75, oem_coilt=0.3125,
             iem_platet=0.375, emag_height=1.625)
FINAL_TAG = "ic0.750_it0.375_ip0.375_oc0.750_ot0.312_h1.625"
AWG = 20
LAYERS = {"inner": 5, "outer": 4}
N_OUTER = 4
RUN_DIR = ROOT / "gmsh_getdp" / "assembly" / "sweep_runs_bpl1000" / FINAL_TAG
FIG = ROOT / "gmsh_getdp" / "assembly" / "em_sizing_figures" / "final_spec_sensitivity.png"
RHO_CU_KG_M3 = 8960.0


def final_metrics():
    return json.loads((RUN_DIR / "metrics.json").read_text())


def height_m(which):
    h = FINAL["emag_height"]
    if which == "inner":
        h += 0.375 - FINAL["iem_platet"]
    return h * IN


def coil(which, layers, awg=AWG, core_in=0.75):
    """A coil wound with `layers` full layers of heavy-build wire (max OD)."""
    d = POLYIMIDE.od_mm(awg) * 1e-3
    h = height_m(which)
    per_layer = int(math.floor(h / d + 1e-9))
    r_core = core_in * IN / 2
    r_out = r_core + GROUND_INSULATION_MM * 1e-3 + layers * d
    w = CoilWindow(which, 1 if which == "inner" else N_OUTER, r_core, r_out, h)
    c = WoundCoil(w, POLYIMIDE, awg, per_layer * layers, layers, "square")
    return c, d, 2 * r_out / IN


def od_for_turns(which, turns, film_scale, awg=AWG, core_in=0.75):
    """Wound OD (in) needed to fit `turns` with the film scaled by
    film_scale (0 = bare copper, 1 = heavy build max, 2 = double film).
    Thinner film also fits more turns per layer, so fewer layers may be
    needed -- both effects are counted."""
    bare = POLYIMIDE.bare_mm(awg) * 1e-3
    d = bare + film_scale * (POLYIMIDE.od_mm(awg) * 1e-3 - bare)
    per_layer = int(math.floor(height_m(which) / d + 1e-9))
    layers = math.ceil(turns / per_layer)
    return layers, per_layer, (core_in * IN + 2 * (GROUND_INSULATION_MM * 1e-3 + layers * d)) / IN


def thermal(ci, co, NI_in, NI_out, area):
    I_in, I_out = NI_in / ci.turns, NI_out / co.turns
    th = solve_steady(I_in, area, coils=[ci.as_thermal_coil(), co.as_thermal_coil(I_out / I_in)])
    R_in = ci.resistance(th["T_coil_C"]["inner"])
    R_out = co.resistance(th["T_coil_C"]["outer"])          # one coil
    return dict(I_in=I_in, I_out=I_out, R_in_hot=R_in, R_out_hot=R_out,
                V_in=I_in * R_in, V_out=I_out * R_out * N_OUTER,
                P_in=I_in ** 2 * R_in, P_out=N_OUTER * I_out ** 2 * R_out,
                P_total=th["Q"], T_iron=th["T_iron_C"],
                T_in=th["T_coil_C"]["inner"], T_out=th["T_coil_C"]["outer"])


def main():
    m = final_metrics()
    NI_in, NI_out = m["ref_req_inner_At"], m["ref_req_outer_At"]
    area = m["envelope_area_m2"]
    print(f"FINAL geometry {FINAL_TAG}: req {NI_in:.1f} / {NI_out:.1f} At (committed ratio), "
          f"quality {m['ref_quality_score']:.2f}, shafts {m['ref_inner_shaft_B_max_T']:.2f} / "
          f"{m['ref_outer_shaft_B_max_T']:.2f} T, CAD {m['cad_microversion']}")

    ci, _, od_i = coil("inner", LAYERS["inner"])
    co, _, od_o = coil("outer", LAYERS["outer"])
    t = thermal(ci, co, NI_in, NI_out, area)
    a_cu = math.pi * (ci.d_bare / 2) ** 2
    print("\n=== SPEC (AWG 20 heavy-build polyimide, square packing, max OD) ===")
    for nm, c, od, I, R20, Rh, P, T, n in (
            ("inner", ci, od_i, t["I_in"], ci.resistance(20), t["R_in_hot"], t["P_in"], t["T_in"], 1),
            ("outer", co, od_o, t["I_out"], co.resistance(20), t["R_out_hot"], t["P_out"] / N_OUTER, t["T_out"], N_OUTER)):
        print(f"{nm}: {c.layers_used} layers x {c.turns // c.layers_used} = {c.turns} turns; "
              f"wound OD {od:.4f} in (CAD coil OD {0.75 + (FINAL['iem_coilt'] if nm == 'inner' else FINAL['oem_coilt']):.4f}); "
              f"I {I:.3f} A; R20 {R20:.3f} ohm, R_hot {Rh:.3f} ohm; P/coil {P:.2f} W; "
              f"T {T:.1f} C; wire {c.wire_length:.2f} m/coil ({c.wire_length * a_cu * RHO_CU_KG_M3 * 1e3:.0f} g)")
    print(f"channels: inner {t['I_in']:.3f} A @ {t['V_in']:.2f} V hot; outer (4 in series) "
          f"{t['I_out']:.3f} A @ {t['V_out']:.2f} V hot; total {t['P_total']:.2f} W; iron {t['T_iron']:.1f} C")

    # --- 1. insulation ------------------------------------------------------
    print("\n=== 1. insulation film: bare copper / heavy build / double heavy build ===")
    ins_rows = []
    for fs, lbl in ((0.0, "bare copper"), (1.0, "heavy build (spec)"), (2.0, "2x heavy film")):
        row = [lbl]
        for nm, c in (("inner", ci), ("outer", co)):
            L, pl, od = od_for_turns(nm, c.turns, fs)
            row.append((f"{L}x{pl}", od))
        ins_rows.append(row)
        print(f"{lbl:20s} inner {ci.turns} turns as {row[1][0]}: wound OD {row[1][1]:.4f} in | "
              f"outer {co.turns} turns as {row[2][0]}: wound OD {row[2][1]:.4f} in")
    print(f"limits: inner OD <= 1.385 in, outer OD <= 1.720 in")

    # --- 2. gauge -----------------------------------------------------------
    print("\n=== 2. gauge, in the specified CAD windows (fills whole layers) ===")
    print("AWG | inner turns  I_A    V    | outer turns  I_A    V(4 series) | P_total_W  T_iron_C")
    gauge_rows = []
    for awg in range(16, 29, 2):
        res = []
        for nm, ct in (("inner", FINAL["iem_coilt"]), ("outer", FINAL["oem_coilt"])):
            build = ct * IN / 2 - GROUND_INSULATION_MM * 1e-3
            n, pl, L, _ = turns_fit(build, height_m(nm), POLYIMIDE.od_mm(awg) * 1e-3, "square")
            res.append(coil(nm, L, awg=awg)[0])
        tg = thermal(res[0], res[1], NI_in, NI_out, area)
        gauge_rows.append((awg, tg))
        print(f"{awg:3d} | {res[0].turns:6d} {tg['I_in']:6.2f} {tg['V_in']:5.2f} | {res[1].turns:6d} "
              f"{tg['I_out']:6.2f} {tg['V_out']:6.2f}      | {tg['P_total']:8.2f} {tg['T_iron']:8.1f}")

    # --- 3. diameter vs power/temperature ----------------------------------
    print("\n=== 3. layers (= overall diameter) vs power and temperature ===")
    diam = {}
    for nm, other in (("inner", "outer"), ("outer", "inner")):
        diam[nm] = []
        print(f"{nm} (other at spec):  layers  turns  OD_in   I_A   P_total_W  T_coil_C")
        for L in range(1, 9):
            c, _, od = coil(nm, L)
            fixed = coil(other, LAYERS[other])[0]
            ci_, co_ = (c, fixed) if nm == "inner" else (fixed, c)
            tl = thermal(ci_, co_, NI_in, NI_out, area)
            fits = od <= (1.385 if nm == "inner" else 1.720)
            diam[nm].append((od, tl["P_total"], tl["T_in" if nm == "inner" else "T_out"], fits, L))
            print(f"{'':20s}{L:6d} {c.turns:6d} {od:6.3f} {tl['I_in' if nm == 'inner' else 'I_out']:6.2f} "
                  f"{tl['P_total']:9.2f} {tl['T_in' if nm == 'inner' else 'T_out']:8.1f}{'' if fits else '  (does not fit)'}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    for nm, mk in (("inner", "o"), ("outer", "s")):
        rows = diam[nm]
        axs[0].plot([r[0] for r in rows if r[3]], [r[1] for r in rows if r[3]], mk + "-", label=f"{nm} (fits)")
        axs[0].plot([r[0] for r in rows if not r[3]], [r[1] for r in rows if not r[3]], mk + ":", color="grey")
        sp = [r for r in rows if r[4] == LAYERS[nm]][0]
        axs[0].plot(sp[0], sp[1], "k*", ms=14)
    axs[0].set_xlabel("wound OD (in)"); axs[0].set_ylabel("total magnet power (W)")
    axs[0].set_title("diameter vs power (star = spec)")
    awgs = [r[0] for r in gauge_rows]
    axs[1].plot(awgs, [r[1]["P_total"] for r in gauge_rows], "o-")
    axs[1].axvline(AWG, color="k", ls=":")
    axs[1].set_xlabel("AWG (in the spec windows)"); axs[1].set_ylabel("total magnet power (W)")
    axs[1].set_title("gauge vs power")
    ax2 = axs[1].twinx()
    ax2.plot(awgs, [r[1]["I_in"] for r in gauge_rows], "s--", color="C1", label="inner A")
    ax2.plot(awgs, [r[1]["I_out"] for r in gauge_rows], "^--", color="C2", label="outer A")
    ax2.axhspan(1, 2, color="C1", alpha=0.08)
    ax2.set_ylabel("channel current (A)"); ax2.legend(fontsize=8)
    lbls = [r[0] for r in ins_rows]
    x = range(len(lbls))
    axs[2].bar([i - 0.2 for i in x], [r[1][1] for r in ins_rows], 0.4, label="inner")
    axs[2].bar([i + 0.2 for i in x], [r[2][1] for r in ins_rows], 0.4, label="outer")
    axs[2].axhline(1.385, color="C0", ls="--", lw=1); axs[2].axhline(1.720, color="C1", ls="--", lw=1)
    axs[2].set_xticks(list(x)); axs[2].set_xticklabels(lbls, fontsize=8)
    axs[2].set_ylim(0.9, 1.8); axs[2].set_ylabel("wound OD at spec turns (in)")
    axs[2].set_title("insulation film vs diameter (dashed = fit limit)")
    for a in (axs[0], axs[2]):
        a.legend(fontsize=8); a.grid(alpha=0.3)
    axs[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG, dpi=130)
    print(f"\nfigure -> {FIG}")


if __name__ == "__main__":
    main()
