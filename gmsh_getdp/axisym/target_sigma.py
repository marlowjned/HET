"""
target_sigma.py

Match the reference B_r profile supplied on 2026-10-21: a half-Gaussian
peaked at the exit plane, sigma = 0.353 L_ch (digitised; L_ch = anode-to-exit
= 33.02 mm, so sigma = 11.66 mm), Bmax ~165 G, ~2 % of peak at the anode.

The 10-21 final geometry (final_design.FINAL) only reaches sigma 7.2-9.9 mm
through the inner-screen gap, so this widens the sweep (gaps past 15 mm,
thinner screens, thicker pole rims, wider outer-screen gap) and runs a compass
search on the misfit to the reference curve itself.

Shape is scored at the spec ratio; drive, power and saturation are reported at
the 165 G target (linear iron: NI x 0.55, P x 0.3025, B_iron x 0.55).

    python target_sigma.py sweeps|search|figures|all
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import CHANNEL, DEFAULTS, Y_EXIT, ring_geometry  # noqa: E402

FIG = HERE / "figures" / "target"
RES = HERE / "results" / "target"
L_CH = CHANNEL["y_anode"] - Y_EXIT            # 33.02 mm
SIG_REF = 0.353 * L_CH                         # 11.66 mm
B_REF = 165.0                                  # G
K = B_REF / 300.0                              # field scale from the 300 G runs
X0 = dict(fd.FINAL)

SWEEPS = {
    "is_gap": [6.5, 10.0, 15.0, 18.0, 21.0, 24.0, 28.0, 32.0, 36.0, 40.0],
    "is_t": [1.5, 2.0, 2.5, 3.5],
    "cap_t": [1.5, 3.0, 5.0, 7.0, 9.525],
    "fp_t": [1.5, 3.0, 5.0, 7.0, 9.525],
    "os_gap": [5.5, 10.0, 15.0, 20.0, 30.0],
}
GRID_GAPS = [15.0, 21.0, 28.0]                 # is_t / rims / os_gap also swept at these gaps
FP_W_MAX = ring_geometry(DEFAULTS)[0] - DEFAULTS["fp_bore_r"] - 1.0   # taper must stop short of the posts (17.4)


def feasible(x):
    """optimise.feasible plus: the front-plate taper may not run under the outer
    posts (it would open an air gap at the post seat)."""
    if x["fp_w"] > FP_W_MAX + 1e-9:
        return "front-plate taper under the posts"
    return opt.feasible(x)


SIMPLE = dict(X0, cap_t=9.525, is_gap=21.0)  # the two changes that carry most of the match
SEARCH_VARS = ["is_gap", "is_t", "cap_t", "fp_t", "os_gap", "cap_w", "fp_w"]
COLS = ["case", "ref_rms_pct", "ref_max_pct", "objective", "sigma_mm", "gauss_rms_pct", "peak_offset_mm",
        "exit_Br_pct", "anode_B_pct", "half_depth_mm", "plume_B_10mm_pct", "P_W", "P165_W",
        "I165_inner_A", "I165_outer_A", "sat_worst", "sat165", "sat_worst_part", "inner_core_Bflux_max_T"] + list(opt.BASE)


def ref_curve(d):
    return np.exp(-d ** 2 / (2 * SIG_REF ** 2))


def score(x):
    m, (y, Br, By), _ = opt.evaluate_x(x)
    d = y - Y_EXIT
    inside = (d >= 0) & (d <= L_CH)
    pk = Br[inside].max()                      # 300 G by construction (peak window starts 15 mm out)
    err = Br[inside] / pk - ref_curve(d[inside])
    m["ref_rms_pct"] = float(np.sqrt(np.mean(err ** 2)) * 100)
    m["ref_max_pct"] = float(np.max(np.abs(err)) * 100)
    m["P165_W"] = m["P_W"] * K ** 2
    m["I165_inner_A"], m["I165_outer_A"] = m["I_inner_A"] * K, m["I_outer_A"] * K
    m["sat165"] = m["sat_worst"] * K
    m["objective"] = m["ref_rms_pct"] + 100 * max(0.0, m["sat165"] - opt.SAT_FLAG)
    return m, (y, Br * K, By * K)


def write(name, rows):
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / f"{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


def show(rows, title):
    print(f"\n{title}")
    print(f"{'case':>26s} {'refRMS':>6s} {'refMax':>6s} {'sig':>5s} {'an|B|%':>6s} {'exit%':>6s} {'pk':>5s} "
          f"{'P165':>5s} {'sat165':>6s}")
    for r in rows:
        print(f"{r['case']:>26s} {r['ref_rms_pct']:6.2f} {r['ref_max_pct']:6.2f} {r['sigma_mm']:5.1f} "
              f"{r['anode_B_pct']:6.2f} {r['exit_Br_pct']:6.1f} {r['peak_offset_mm']:5.1f} {r['P165_W']:5.2f} "
              f"{r['sat165']:6.2f}")


def sweeps():
    cases = [(f"is_gap={g:g}", dict(X0, is_gap=g)) for g in SWEEPS["is_gap"]]
    for g in GRID_GAPS:
        for v in ("is_t", "cap_t", "fp_t", "os_gap"):
            cases += [(f"gap{g:g} {v}={val:g}", dict(X0, is_gap=g, **{v: val})) for val in SWEEPS[v]]
    cases = [(n, x) for n, x in cases if feasible(x) is None]
    opt.solve_all([x for _, x in cases])
    rows = []
    for n, x in cases:
        m, _ = score(x)
        m["case"] = n
        rows.append(m)
    write("sweeps", rows)
    show(rows, "sweeps (shape vs the reference, drive at 165 G)")
    return rows


def search(x0=None, steps=(4.0, 2.0, 1.0, 0.5), max_iter=80):
    if x0 is None:
        rows = list(csv.DictReader(open(RES / "sweeps.csv")))
        best = min(rows, key=lambda r: float(r["objective"]))
        x0 = {k: float(best[k]) for k in opt.BASE}
        print(f"start from sweep best: {best['case']} ({float(best['objective']):.2f})")
    x = dict(x0)
    best, _ = score(x)
    best["case"] = "it0"
    trace = [best]
    it = 0
    for st in steps:
        while it < max_iter:
            it += 1
            cands = []
            for v in SEARCH_VARS:
                for sgn in (-1, 1):
                    val = x[v] + sgn * st
                    if v in ("cap_t", "fp_t"):
                        val = min(val, 9.525)
                    c = dict(x, **{v: round(val, 3)})
                    if c != x and feasible(c) is None:
                        cands.append(c)
            opt.solve_all(cands)
            ms = [score(c)[0] for c in cands]
            mb = min(ms, key=lambda m: m["objective"])
            if mb["objective"] < best["objective"] - 1e-3:
                best = dict(mb, case=f"it{it} step{st:g}")
                x = {k: mb[k] for k in opt.BASE}
                trace.append(best)
                print(f"it {it} step {st}: ref RMS {best['ref_rms_pct']:.3f}  " +
                      " ".join(f"{v}={x[v]:g}" for v in SEARCH_VARS), flush=True)
            else:
                break
    write("search_trace", trace)
    (RES / "optimum.json").write_text(json.dumps(x, indent=1))
    show(trace, "search trace")
    return x


def figures():
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    xo = json.loads((RES / "optimum.json").read_text())
    picks = [("10-21 final (gap 6.5)", X0), ("simple: flat cap + gap 21", SIMPLE), ("search optimum", xo)]
    rows = list(csv.DictReader(open(RES / "sweeps.csv")))
    out = [(lab, *score(x), x) for lab, x in picks]
    m0, (y0, B0, _) = fd.as_built()
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    ax = axs[0]
    d = np.linspace(0, L_CH, 200)
    ax.plot(1 - d / L_CH, B_REF * ref_curve(d), color="k", lw=2.4, label=f"reference (sigma {SIG_REF:.1f} mm)")
    for (lab, m, (y, Br, By), x), c in zip(out, [fd.PALETTE[0], fd.PALETTE[3], fd.PALETTE[1]]):
        dd = y - Y_EXIT
        sel = (dd >= -0.1 * L_CH) & (dd <= L_CH)
        ax.plot(1 - dd[sel] / L_CH, Br[sel], color=c, lw=1.6,
                label=f"{lab}: RMS {m['ref_rms_pct']:.1f}%, sigma {m['sigma_mm']:.1f}, anode B_r {m['anode_Br_pct']:.1f}%")
    dd = y0 - Y_EXIT
    sel = (dd >= -0.1 * L_CH) & (dd <= L_CH)
    ax.plot(1 - dd[sel] / L_CH, B0[sel] * K, color="0.55", ls="--", lw=1.2, label="as built")
    ax.set_xlabel("x / L_ch (0 = anode, 1 = exit)")
    ax.set_ylabel("channel-averaged B_r (G)")
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 200)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title("Profiles vs the reference (165 G peak)", loc="left", fontsize=10)
    ax = axs[1]
    for (lab, m, (y, Br, By), x), c in zip(out, [fd.PALETTE[0], fd.PALETTE[3], fd.PALETTE[1]]):
        dd = y - Y_EXIT
        sel = (dd >= 0) & (dd <= L_CH)
        ax.plot(1 - dd[sel] / L_CH, Br[sel] / B_REF * 100 - 100 * ref_curve(dd[sel]), color=c, lw=1.6, label=lab)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("x / L_ch")
    ax.set_ylabel("profile - reference (% of peak)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Residual", loc="left", fontsize=10)
    fig.savefig(FIG / "profiles_vs_reference.png", dpi=140)
    plt.close(fig)
    # sweep panels
    fig, axs = plt.subplots(1, 5, figsize=(17, 3.2), constrained_layout=True)
    g = [r for r in rows if r["case"].startswith("is_gap=")]
    axs[0].plot([float(r["is_gap"]) for r in g], [float(r["ref_rms_pct"]) for r in g], "o-", color=fd.PALETTE[0])
    axs[0].set_xlabel("inner-screen gap (mm)")
    axs[0].set_ylabel("RMS misfit to reference (%)")
    for ax, v in zip(axs[1:], ("is_t", "cap_t", "fp_t", "os_gap")):
        for i, gg in enumerate(GRID_GAPS):
            sub = sorted([r for r in rows if r["case"].startswith(f"gap{gg:g} {v}=")], key=lambda r: float(r[v]))
            ax.plot([float(r[v]) for r in sub], [float(r["ref_rms_pct"]) for r in sub], "o-", color=fd.PALETTE[i],
                    ms=3.5, label=f"gap {gg:g}")
        ax.set_xlabel(opt.LABEL[v])
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Misfit to the reference curve, one knob at a time", x=0.01, ha="left", fontsize=10)
    fig.savefig(FIG / "sweeps.png", dpi=130)
    plt.close(fig)
    fd.FIG = FIG
    fd.drawing(SIMPLE, fname="dimensioned_target",
               note=f"Simple match to the reference profile (sigma {SIG_REF:.1f} mm, 165 G): flat cap + inner-screen gap 21.")
    return out


if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in what:
        if w in ("sweeps", "all"):
            sweeps()
        if w in ("search", "all"):
            search()
        if w in ("figures", "all"):
            figures()
