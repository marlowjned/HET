"""
polecap.py

Inner pole cap shaping (user, 2026-10-15: "we can modify the inner pole cap
to be whatever"), on top of the tapered front plate (t 1.5 mm, w 15 mm, from
explore.py's taper study). The peak B_r sits between the two pole faces, so
moving it onto the exit plane needs BOTH poles concentrated there.

Stage A: cap variants, both plates and front plate only, no trim, each at its
minimum-power ratio. Stage B: the variant whose peak lands closest to the
exit, plus the T2 trim coil (trim.py), both configurations.

Cap parameters (explore.pole_cap): r radius (keep-out starts at 20.62 mm),
rim tip thickness t with a back-side taper of radial width w, and ext, a rim
lip protruding ext mm downstream of the exit plane (sits in the plume).

    python polecap.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from explore import (COLS, Q_SCAN, coil_power, evaluate, fieldline_fig, front_plate, pole_cap,  # noqa: E402
                     profile_fig, table)
from solve import solve_case  # noqa: E402
from trim import PLACEMENTS, best  # noqa: E402

FP = front_plate(t_tip=1.5, w=15.0)
CAPS = {
    "cap as built": dict(),
    "cap t3 w3": dict(t_tip=3.0, w=3.0),
    "cap t1.5 w3": dict(t_tip=1.5, w=3.0),
    "cap t1.5 w6": dict(t_tip=1.5, w=6.0),
    "cap r20 t3 w6": dict(r=20.0, t_tip=3.0, w=6.0),
    "cap r20 t1.5 w6": dict(r=20.0, t_tip=1.5, w=6.0),
    "cap r20 t1.5 w6 ext2": dict(r=20.0, t_tip=1.5, w=6.0, ext=2.0),
}
CONFIGS = {"both plates": {}, "front only": dict(back_plate=False)}


def min_power(res):
    ms = []
    for q in Q_SCAN:
        m, prof = evaluate(res, q)
        m["P_W"] = coil_power(m["NI_inner"], m["NI_outer"])
        ms.append((m, prof))
    return min(ms, key=lambda t: t[0]["P_W"])


def main():
    for c in ("P_W", "NI_trim", "rear_B_pct"):
        if c not in COLS:
            COLS.append(c)
    rows, fl = [], []
    for cname, cp in CONFIGS.items():
        curves, stage_a = [], []
        for capname, kw in CAPS.items():
            params = dict(cp, front_plate_poly=FP, pole_cap_poly=pole_cap(**kw) if kw else None)
            name = f"{cname} fp-taper {capname}"
            res = solve_case(name, params)
            m, prof = min_power(res)
            m["case"] = f"{cname}, {capname}"
            rows.append(m)
            stage_a.append((abs(m["peak_offset_mm"]), capname, kw))
            curves.append((f"{capname}, {m['P_W']:.1f} W", prof, m))
        # stage B: closest-to-exit cap without the plume-protruding lip, + T2 trim
        _, capname, kw = min((s for s in stage_a if "ext" not in s[2]), key=lambda s: (s[0], s[1]))
        win = PLACEMENTS["T2 behind anode"]
        params = dict(cp, front_plate_poly=FP, pole_cap_poly=pole_cap(**kw) if kw else None, trim=win)
        res = solve_case(f"{cname} fp-taper {capname} + T2", params, map_=True)
        m, prof = best(res, win, constrained=True)
        m["case"] = f"{cname}, {capname} + T2 trim"
        rows.append(m)
        curves.append((f"{capname} + T2 trim (x{m['NI_trim'] / m['NI_inner']:+.2f}), {m['P_W']:.1f} W", prof, m))
        fl.append((f"{cname}: tapered front plate, {capname}, T2 trim\n"
                   f"NI {m['NI_inner']:.0f}/{m['NI_outer']:.0f}/trim {m['NI_trim']:.0f}, {m['P_W']:.1f} W",
                   res, m["NI_inner"], m["NI_outer"], params, m["NI_trim"]))
        profile_fig(curves, f"polecap_{cname.replace(' ', '_')}",
                    f"Pole cap shaping + tapered front plate, {cname}")
    table(rows, "polecap")
    fieldline_fig(fl, "polecap")


if __name__ == "__main__":
    main()
