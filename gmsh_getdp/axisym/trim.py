"""
trim.py

Trim-coil study (user request 2026-10-15): a third, independently driven
winding near the back of the channel, normally run against the inner coil,
to cancel the residual field near the anode. Linear iron, so three basis
solves per geometry give every (inner, outer, trim) excitation:

    B = NI_in * (b_inner + q b_outer + t b_trim),  q = NI_out/NI_in, t = NI_trim/NI_in

For each plate configuration x trim placement, a grid over (q, t) is scored
on the channel profile scaled to 300 G peak B_r:
  rear_B_pct   mean |B| over the rear half of the channel (anode side), % of
               peak -- the "decay to as little as possible" objective
  reversal     whether B_r goes more than 1 G negative anywhere in the channel
and the best point is reported with and without the no-reversal constraint.

Placements (mm: r0, r1, y0, y1):
  T1 rear inner    ring outside the inner coil, anode plane to y 12.95
                   (r 14.8-19.8: clear of the coil OD 14.29 and the channel
                   wall keep-out 20.62)
  T2 behind anode  at channel radius, just behind the anode face

Power: inner 230 t, 0.643 ohm hot and outer 4 x 184 t, 0.494 ohm (spec
windings); the trim coil uses the inner coil's ohm per A-turn^2 scaled by
mean turn radius and inverse window area (same wire and packing).

    python trim.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from explore import COLS, coil_power, fieldline_fig, profile_fig, table  # noqa: E402
from geometry import CHANNEL  # noqa: E402
from metrics import TARGET_G, profile, shape_metrics  # noqa: E402
from solve import solve_case  # noqa: E402

PLACEMENTS = {"T1 rear inner": (14.8, 19.8, CHANNEL["y_anode"], 12.954),
              "T2 behind anode": (CHANNEL["r_in"], CHANNEL["r_out"], -2.0, 8.0)}
CONFIGS = {"both plates": {}, "front plate only": dict(back_plate=False), "no plates": dict(plates=False)}
QS = np.round(np.arange(0.0, 3.001, 0.05), 3)
TS = np.round(np.arange(-3.0, 3.001, 0.02), 3)

K_IN = 0.643 / 230 ** 2
R_MEAN_IN, A_IN = (9.525 + 14.2875) / 2, (14.2875 - 9.525) * (12.954 + 28.321)


def trim_power(ni, win):
    r0, r1, y0, y1 = win
    return ni ** 2 * K_IN * ((r0 + r1) / 2 / R_MEAN_IN) * (A_IN / ((r1 - r0) * (y1 - y0)))


def grid_search(res):
    y, bi, ai = profile(res["channel"]["inner"])
    _, bo, ao = profile(res["channel"]["outer"])
    has_trim = "trim" in res["channel"]
    _, bt, at = profile(res["channel"]["trim"]) if has_trim else (y, 0 * bi, 0 * ai)
    ts = TS if has_trim else np.array([0.0])
    Br = bi + QS[:, None, None] * bo + ts[None, :, None] * bt
    By = ai + QS[:, None, None] * ao + ts[None, :, None] * at
    ye, ya = CHANNEL["y_exit"], CHANNEL["y_anode"]
    sg = np.sign(Br[..., np.argmin(np.abs(y - ye))])
    sg[sg == 0] = 1
    Br, By = Br * sg[..., None], By * sg[..., None]
    win = (y >= ye - 15) & (y <= ya)
    pk = Br[..., win].max(axis=-1)
    ok = pk > 0
    s = np.where(ok, TARGET_G / np.where(ok, pk, 1), np.nan)
    inside = (y >= ye) & (y <= ya)
    rear = (y >= ye + (ya - ye) / 2) & (y <= ya)
    Bm = np.hypot(Br, By) * s[..., None]
    rear_B = Bm[..., rear].mean(axis=-1) / TARGET_G * 100
    rev = (Br[..., inside] * s[..., None]).min(axis=-1) < -1.0
    return y, Br, By, s, rear_B, rev, ts, sg


def best(res, win, constrained):
    y, Br, By, s, rear_B, rev, ts, sg = grid_search(res)
    obj = np.where(rev & constrained, np.inf, rear_B)
    obj = np.where(np.isnan(obj), np.inf, obj)
    i, j = np.unravel_index(np.argmin(obj), obj.shape)
    q, t, sc = QS[i], ts[j], s[i, j] * sg[i, j]
    m, _ = shape_metrics(y, Br[i, j], By[i, j])
    m.update(q=float(q), NI_inner=float(abs(sc)), NI_outer=float(abs(sc) * q), NI_trim=float(abs(sc) * t),
             rear_B_pct=float(rear_B[i, j]), reversal=bool(rev[i, j]))
    m["P_W"] = coil_power(m["NI_inner"], m["NI_outer"]) + (trim_power(m["NI_trim"], win) if win else 0.0)
    exc = {"inner": sc, "outer": sc * q, "trim": sc * t}
    Bfe = sum(exc[k] * res["iron"][k][:, 2:4] for k in res["iron"])
    m["iron_B_max_T"] = float(np.percentile(np.linalg.norm(Bfe, axis=1), 99.9))
    return m, (y, Br[i, j] * s[i, j], By[i, j] * s[i, j])


def main():
    for c in ("NI_trim", "rear_B_pct", "reversal", "P_W"):
        if c not in COLS:
            COLS.append(c)
    rows, fl = [], []
    for cname, cp in CONFIGS.items():
        curves = []
        for tname, win in [("no trim", None)] + list(PLACEMENTS.items()):
            name = f"{cname} + {tname}" if win else {"both plates": "baseline"}.get(cname, cname)
            params = dict(cp, trim=win) if win else cp
            res = solve_case(name, params, map_=True)
            for constrained in ((True, False) if win else (True,)):
                m, prof = best(res, win, constrained)
                tag = "" if not win else (" (no reversal)" if constrained else " (unconstrained)")
                m["case"] = f"{cname}, {tname}{tag}"
                rows.append(m)
                lab = (f"{tname}{tag}: q={m['q']:.2f}, trim/inner={m['NI_trim'] / m['NI_inner']:+.2f}, "
                       f"rear {m['rear_B_pct']:.0f}%, {m['P_W']:.0f} W")
                curves.append((lab, prof, m))
                if constrained:
                    fl.append((f"{cname}, {tname}\nq={m['q']:.2f}, NI {m['NI_inner']:.0f}/{m['NI_outer']:.0f}"
                               f"/trim {m['NI_trim']:.0f}, {m['P_W']:.0f} W", res, m["NI_inner"],
                               m["NI_outer"], params, m["NI_trim"]))
        profile_fig(curves, f"trim_{cname.replace(' ', '_')}", f"Trim coil, {cname}")
    table(rows, "trim")
    print(f"\n{'case':48s} {'NI_trim':>7s} {'rear%':>6s} {'rev':>5s} {'P_W':>6s}")
    for r in rows:
        print(f"{r['case']:48s} {r['NI_trim']:7.0f} {r['rear_B_pct']:6.1f} {str(r['reversal']):>5s} {r['P_W']:6.1f}")
    for k, cname in enumerate(CONFIGS):
        fieldline_fig(fl[3 * k:3 * k + 3], f"trim_fieldlines_{cname.replace(' ', '_')}")


if __name__ == "__main__":
    main()
