"""
permeability.py

Robustness of the recommended design (peak_forward best chamfers, built with
the stock DOM tubes) to the permeability of the steel, 2026-10-21.

Stock tubes (Speedy Metals, C1020 DOM, ASTM A513 Type 5, cold drawn):
  inner screen  1-1/2" OD x 0.120" wall  ->  r 16.00-19.05
  outer screen  3-1/2" OD x 0.188" wall  ->  r 39.67-44.45

Low-carbon steel spans roughly mu_r 100 (cold worked, at very low field) to
5000 (fully annealed, near the permeability peak). Cases:
  all       every iron part (cores, posts, plates, pole cap, screens) at one mu_r
  screens   the screens at their own mu_r, the rest held at 200 / 1000 / 5000

Scored against the reference profile (target_sigma.py) at 165 G peak.

    python permeability.py
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import target_sigma as ts  # noqa: E402
import peak_forward as pf  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import Y_EXIT  # noqa: E402

FIG = HERE / "figures" / "permeability"
RES = HERE / "results" / "permeability"
IN = 25.4
BEST = json.loads((pf.RES / "best_chamfers.json").read_text())
STOCK = dict(BEST, is_r0=round((0.75 - 0.120) * IN, 3), is_t=round(0.120 * IN, 3),
             os_r=round(3.124 / 2 * IN, 3), os_t=round(0.188 * IN, 3))
MU_ALL = [100.0, 200.0, 500.0, 1000.0, 1500.0, 3000.0, 5000.0]
MU_SCREEN = [100.0, 200.0, 500.0, 1000.0, 3000.0, 5000.0]
MU_REST = [200.0, 1000.0, 5000.0]
COLS = ["case", "mu_rest", "mu_screen", "ref_rms_pct", "ref_max_pct", "peak_offset_mm", "exit_Br_pct",
        "anode_Br_pct", "anode_B_pct", "sigma_mm", "half_depth_mm", "P165_W", "I165_inner_A", "I165_outer_A",
        "sat165", "sat_worst_part"]


def cases():
    out = [(f"all iron mu {m:g}", m, m, dict(STOCK, mu=m)) for m in MU_ALL]
    for mr in MU_REST:
        out += [(f"rest {mr:g}, screens {ms:g}", mr, ms, dict(STOCK, mu=mr, mu_screen=ms))
                for ms in MU_SCREEN if ms != mr]
    return out


def run():
    cs = cases()
    opt.solve_all([x for *_, x in cs])
    rows, profs = [], []
    for name, mr, ms, x in cs:
        m, prof = ts.score(x)
        m.update(case=name, mu_rest=mr, mu_screen=ms)
        rows.append(m)
        profs.append(prof)
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / "permeability.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    print(f"{'case':>28s} {'misfit':>6s} {'pk':>5s} {'exit%':>6s} {'anBr%':>6s} {'sig':>5s} {'P165':>5s} {'I_in':>5s}")
    for r in rows:
        print(f"{r['case']:>28s} {r['ref_rms_pct']:6.2f} {r['peak_offset_mm']:5.1f} {r['exit_Br_pct']:6.1f} "
              f"{r['anode_Br_pct']:6.2f} {r['sigma_mm']:5.1f} {r['P165_W']:5.2f} {r['I165_inner_A']:5.2f}")
    figure(rows, profs)
    return rows


def figure(rows, profs):
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.4), constrained_layout=True)
    ax = axs[0]
    d = np.linspace(0, ts.L_CH, 300)
    ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="k", lw=2.6, label="reference", zorder=5)
    cmap = plt.get_cmap("viridis")
    for r, (y, Br, By) in zip(rows, profs):
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        if r["case"].startswith("all"):
            c = cmap(np.log(r["mu_rest"] / 100) / np.log(50))
            ax.plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.4, label=f"all iron mu_r {r['mu_rest']:g}")
        else:
            ax.plot(1 - dd[sel] / ts.L_CH, Br[sel], color="0.6", lw=0.6, alpha=0.6)
    ax.set_xlim(0, 1.08)
    ax.set_ylim(0, 185)
    ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
    ax.set_ylabel("channel-averaged B_r (G)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title("All 25 cases (grey = screens at their own mu_r)", loc="left", fontsize=10)
    for ax, (k, lab) in zip(axs[1:], [("ref_rms_pct", "RMS misfit to reference (%)"),
                                      ("P165_W", "coil power at 165 G (W)")]):
        a = sorted([r for r in rows if r["case"].startswith("all")], key=lambda r: r["mu_rest"])
        ax.plot([r["mu_rest"] for r in a], [r[k] for r in a], "o-", color="k", lw=1.8, label="all iron")
        for i, mr in enumerate(MU_REST):
            s = sorted([r for r in rows if not r["case"].startswith("all") and r["mu_rest"] == mr]
                       + [r for r in a if r["mu_rest"] == mr], key=lambda r: r["mu_screen"])
            ax.plot([r["mu_screen"] for r in s], [r[k] for r in s], "o--", color=fd.PALETTE[i], ms=3.5,
                    label=f"rest at {mr:g}, x = screen mu_r")
        ax.set_xscale("log")
        ax.set_xlabel("mu_r")
        ax.set_ylabel(lab)
        ax.legend(frameon=False, fontsize=7.5)
    axs[1].set_ylim(0, None)
    axs[2].set_ylim(0, None)
    fig.suptitle("Recommended design with the stock DOM tubes: sensitivity to steel permeability (2-D, linear)",
                 x=0.01, ha="left", fontsize=10)
    fig.savefig(FIG / "permeability.png", dpi=140)
    plt.close(fig)




def lowend():
    """Reference vs the design with low-permeability screens (rest of the iron high)."""
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    picks = [("nominal: all iron mu_r 1500", dict(STOCK, mu=1500.0), "k", "-"),
             ("screens 500, rest 5000", dict(STOCK, mu=5000.0, mu_screen=500.0), fd.PALETTE[0], "-"),
             ("screens 200, rest 5000", dict(STOCK, mu=5000.0, mu_screen=200.0), fd.PALETTE[1], "-"),
             ("screens 100, rest 5000", dict(STOCK, mu=5000.0, mu_screen=100.0), fd.PALETTE[3], "-"),
             ("all iron mu_r 100", dict(STOCK, mu=100.0), "0.55", "--")]
    out = [(lab, *ts.score(x), c, ls) for lab, x, c, ls in picks]
    d = np.linspace(0, ts.L_CH, 300)
    xr, br = 1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    for ax in axs[:2]:
        ax.plot(xr, br, color="#888888", lw=6, alpha=0.45, label="reference (ideal)")
    for lab, m, (y, Br, By), c, ls in out:
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
        x = 1 - dd[sel] / ts.L_CH
        txt = f"{lab}: anode {m['anode_Br_pct']:.1f}%, misfit {m['ref_rms_pct']:.1f}%"
        axs[0].plot(x, Br[sel], color=c, ls=ls, lw=1.6, label=txt)
        axs[1].plot(x, np.clip(Br[sel], 0.1, None), color=c, ls=ls, lw=1.6)
        ins = (dd >= 0) & (dd <= ts.L_CH)
        axs[2].plot(1 - dd[ins] / ts.L_CH, Br[ins] - ts.B_REF * ts.ref_curve(dd[ins]), color=c, ls=ls, lw=1.6)
    axs[0].set_ylim(0, 180)
    axs[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    axs[0].set_title("B_r along the channel (165 G peak)", loc="left", fontsize=10)
    axs[1].set_yscale("log")
    axs[1].set_ylim(1, 250)
    axs[1].set_title("Same, log scale (shows the anode end)", loc="left", fontsize=10)
    axs[2].axhline(0, color="#888888", lw=4, alpha=0.45)
    axs[2].set_title("Difference from the reference (G)", loc="left", fontsize=10)
    axs[2].set_ylabel("B_r - reference (G)")
    for ax in axs:
        ax.set_xlim(0, 1.08 if ax is not axs[2] else 1.0)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.axvline(1, color="0.3", lw=0.7)
    axs[0].set_ylabel("channel-averaged B_r (G)")
    axs[1].set_ylabel("channel-averaged B_r (G)")
    fig.suptitle("Low screen permeability vs the reference (stock DOM tubes, 2-D linear)", x=0.01, ha="left",
                 fontsize=10)
    fig.savefig(FIG / "low_end.png", dpi=140)
    plt.close(fig)
    for lab, m, (y, Br, By), c, ls in out:
        dd = y - Y_EXIT
        at = lambda f: float(np.interp(f * ts.L_CH, dd, Br))
        print(f"{lab:30s} anode {m['anode_Br_pct']:.2f}% ({at(1.0):.1f} G)  x/L=0.25: {at(0.75):.1f} G  "
              f"x/L=0.5: {at(0.5):.1f} G  misfit {m['ref_rms_pct']:.2f}%  max {m['ref_max_pct']:.1f}%")
    print("reference", " ".join(f"{f}: {ts.B_REF*ts.ref_curve(f*ts.L_CH):.1f} G" for f in (1.0, 0.75, 0.5)))


if __name__ == "__main__":
    lowend() if sys.argv[1:] == ["lowend"] else run()
