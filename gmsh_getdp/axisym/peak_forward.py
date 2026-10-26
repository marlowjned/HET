"""
peak_forward.py

Pull the B_r peak toward the exit plane while keeping the wide reference
profile (sigma 11.66 mm, 165 G; target_sigma.py). Rule (user, 2026-10-21):
NOTHING past the exit plane. The flush pole cap is allowed.

Start: target_sigma.SIMPLE (flat pole cap, inner screen tip gap 21), whose
peak sits 2.0 mm inside the channel. Names follow figures/glossary.png.

Methods
  1 chamfers   grid over the pole cap edge face length (cap_t) x pole cap
               chamfer width (cap_w) x front plate bore face length (fp_t)
               x front plate chamfer width (fp_w), at inner screen tip gaps
               15 / 21 / 28, then a compass search on all knobs
  2 trim coil  a ring coil with its own current, either
                 TI  just upstream of the pole cap (r 14.69-20.2, between the
                     inner coil and the channel keep-out), length L
                 TO  just upstream of the front plate bore face (r 33.1-39.5,
                     between the channel outer keep-out and the outer screen)
               trim current per geometry chosen by grid (linear superposition)
  3 bore       front plate bore radius 33.1-34.925 (bore face moved radially
               toward the channel outer wall; keep-out 31.6 + 1.5)

Every case is scored on RMS misfit to the reference curve (which peaks on the
exit plane), with peak position, B_r at the exit, B_r at the anode and coil
power at 165 G reported.

    python peak_forward.py chamfers|search|trim|bore|figures|all
"""
import csv
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimise as opt  # noqa: E402
import target_sigma as ts  # noqa: E402
import final_design as fd  # noqa: E402
from geometry import CHANNEL, DEFAULTS, Y_EXIT, Y_FP_IN  # noqa: E402
from metrics import profile, shape_metrics  # noqa: E402
from shaping import coil_power  # noqa: E402
from solve import RUNS, solve_case  # noqa: E402

FIG = HERE / "figures" / "peak"
RES = HERE / "results" / "peak"
X0 = dict(ts.SIMPLE)
BORE_MIN = CHANNEL["r_out"] + 1.5
CAP_T = [1.5, 3.0, 5.0]
CAP_W = [0.5, 1.0, 2.0, 3.0, 6.0]
FP_T = [1.5, 3.0, 5.0]
FP_W = [1.0, 2.0, 4.0, 8.0, 14.0]
GAPS = [15.0, 21.0, 28.0]
SEARCH_VARS = ["cap_t", "cap_w", "fp_t", "fp_w", "is_gap", "is_t", "os_gap"]
TRIM_WIN = {"TI": (opt.IS_R0, 20.2), "TO": (BORE_MIN, 39.5)}
TRIM_L = [3.0, 6.0, 10.0, 15.0]
TRIM_T = np.round(np.arange(-2.0, 2.0001, 0.01), 3)        # NI_trim / NI_inner
COLS = ["family", "case", "ref_rms_pct", "ref_max_pct", "peak_offset_mm", "exit_Br_pct", "anode_Br_pct",
        "anode_B_pct", "sigma_mm", "half_depth_mm", "P165_W", "I165_inner_A", "I165_outer_A", "trim_NI165",
        "trim_P165_W", "sat165", "fp_bore"] + list(opt.BASE)


def feasible(x):
    if x.get("fp_bore", DEFAULTS["fp_bore_r"]) < BORE_MIN - 1e-9:
        return "bore inside the channel keep-out"
    return ts.feasible(x)


def write(name, rows):
    RES.mkdir(parents=True, exist_ok=True)
    with open(RES / f"{name}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


def read(name):
    out = []
    for r in csv.DictReader(open(RES / f"{name}.csv")):
        out.append({k: (v if k in ("family", "case") or v == "" else float(v)) for k, v in r.items()})
    return out


def show(rows, title, n=None):
    print(f"\n{title}")
    print(f"{'case':>40s} {'RMS':>5s} {'pk':>5s} {'exit%':>6s} {'anBr%':>6s} {'sig':>5s} {'P165':>5s}")
    for r in rows[:n]:
        print(f"{r['case']:>40s} {r['ref_rms_pct']:5.2f} {r['peak_offset_mm']:5.1f} {r['exit_Br_pct']:6.1f} "
              f"{r['anode_Br_pct']:6.2f} {r['sigma_mm']:5.1f} {r['P165_W']:5.2f}")


def run_geoms(cases, family):
    """cases: [(label, x)] of plain geometries; parallel solve, score."""
    cases = [(n, x) for n, x in cases if feasible(x) is None]
    opt.solve_all([x for _, x in cases])
    rows = []
    for n, x in cases:
        m, _ = ts.score(x)
        m.update(family=family, case=n, fp_bore=x.get("fp_bore", DEFAULTS["fp_bore_r"]))
        rows.append(m)
    return sorted(rows, key=lambda r: r["ref_rms_pct"])


# --------------------------------------------------------------------------
# 1 chamfers
def chamfers():
    cases = []
    for g in GAPS:
        caps = [(c, w) for c in CAP_T for w in CAP_W] + [(9.525, 6.0)]
        fps = [(c, w) for c in FP_T for w in FP_W] + [(9.525, 14.0)]
        for ct, cw in caps:
            for ft, fw in fps:
                x = dict(X0, is_gap=g, cap_t=ct, cap_w=cw, fp_t=ft, fp_w=fw)
                cases.append((f"gap{g:g} cap {ct:g}/{cw:g} fp {ft:g}/{fw:g}", x))
    rows = run_geoms(cases, "chamfers")
    write("chamfers", rows)
    show(rows, "chamfer grid, best 25 (cap = edge face length / chamfer width; fp likewise)", 25)
    return rows


def search():
    best0 = read("chamfers")[0]
    x = {k: best0[k] for k in opt.BASE}
    print(f"start: {best0['case']} ({best0['ref_rms_pct']:.3f})")
    best = run_geoms([("it0", x)], "search")[0]
    trace = [best]
    it = 0
    for st in (2.0, 1.0, 0.5, 0.25):
        while it < 80:
            it += 1
            cands = []
            for v in SEARCH_VARS:
                for sgn in (-1, 1):
                    val = round(x[v] + sgn * st, 3)
                    if v in ("cap_t", "fp_t"):
                        val = min(val, 9.525)
                    c = dict(x, **{v: val})
                    if c != x and feasible(c) is None:
                        cands.append((f"it{it}", c))
            rows = run_geoms(cands, "search")
            if rows and rows[0]["ref_rms_pct"] < best["ref_rms_pct"] - 1e-3:
                best = dict(rows[0], case=f"it{it} step{st:g}")
                x = {k: best[k] for k in opt.BASE}
                trace.append(best)
                print(f"it {it} step {st}: RMS {best['ref_rms_pct']:.3f} pk {best['peak_offset_mm']:.1f}  " +
                      " ".join(f"{v}={x[v]:g}" for v in SEARCH_VARS), flush=True)
            else:
                break
    write("search_trace", trace)
    (RES / "optimum.json").write_text(json.dumps(x, indent=1))
    show(trace, "search trace")
    return x


# --------------------------------------------------------------------------
# 2 trim coil
def trim_params(x, kind, L):
    r0, r1 = TRIM_WIN[kind]
    p = opt.to_params(x)
    p["trim"] = (r0, r1, Y_FP_IN, Y_FP_IN + L)
    return p


def trim_key(x, kind, L):
    s = json.dumps(dict({k: round(float(x[k]), 3) for k in sorted(x)}, kind=kind, L=L))
    return "peak/trim_" + hashlib.sha1(s.encode()).hexdigest()[:12]


def _trim_solve(a):
    x, kind, L = a
    solve_case(trim_key(x, kind, L), trim_params(x, kind, L))
    return True


def trim_score(x, kind, L):
    from trim import trim_power
    res = solve_case(trim_key(x, kind, L), trim_params(x, kind, L))
    y, bri, byi = profile(res["channel"]["inner"])
    _, bro, byo = profile(res["channel"]["outer"])
    _, brt, byt = profile(res["channel"]["trim"])
    d = y - Y_EXIT
    inside = (d >= 0) & (d <= ts.L_CH)
    ref = ts.ref_curve(d[inside])
    best = None
    for t in TRIM_T:
        Br = bri + opt.Q_SPEC * bro + t * brt
        sgn = np.sign(Br[np.argmin(np.abs(d))]) or 1.0
        b = sgn * Br[inside]
        if b.max() <= 0:
            continue
        e = np.sqrt(np.mean((b / b.max() - ref) ** 2))
        if best is None or e < best[0]:
            best = (e, t)
    e, t = best
    By = byi + opt.Q_SPEC * byo + t * byt
    Br = bri + opt.Q_SPEC * bro + t * brt
    m, s = shape_metrics(y, Br, By)
    ni_in = abs(s) * ts.K
    ni_out, ni_tr = ni_in * opt.Q_SPEC, ni_in * t
    p = trim_params(x, kind, L)
    pt = trim_power(abs(ni_tr), p["trim"])
    b = (Br * s)[inside]
    err = b / b.max() - ref
    m.update(family="trim", case=f"{kind} L{L:g} t{t:+.2f}", ref_rms_pct=float(np.sqrt(np.mean(err ** 2)) * 100),
             ref_max_pct=float(np.max(np.abs(err)) * 100), trim_NI165=ni_tr, trim_P165_W=pt,
             P165_W=coil_power(opt.to_params(x), ni_in, ni_out) + pt, I165_inner_A=ni_in / opt.N_IN,
             I165_outer_A=ni_out / opt.N_OUT, sat165=float("nan"), fp_bore=x.get("fp_bore", DEFAULTS["fp_bore_r"]))
    m.update(x)
    prof = (y, Br * s * ts.K, By * s * ts.K)
    return m, prof


def trim(bases=None):
    bases = bases or [("simple", X0)]
    jobs = [(x, k, L) for _, x in bases for k in TRIM_WIN for L in TRIM_L]
    new = [j for j in jobs if not (RUNS / trim_key(*j) / "basis.npz").exists()]
    if new:
        print(f"solving {len(new)} trim geometries", flush=True)
        with ProcessPoolExecutor(opt.WORKERS) as ex:
            list(ex.map(_trim_solve, new))
    rows = []
    for (bn, x) in bases:
        for k in TRIM_WIN:
            for L in TRIM_L:
                m, _ = trim_score(x, k, L)
                m["case"] = f"{bn}: {m['case']}"
                rows.append(m)
    rows.sort(key=lambda r: r["ref_rms_pct"])
    write("trim", rows)
    show(rows, "trim coil (TI behind the pole cap, TO behind the front plate bore face)")
    return rows


# --------------------------------------------------------------------------
# 3 bore radius
def bore(bases=None):
    bases = bases or [("simple", X0)]
    cases = [(f"{bn}: bore {b:g}", dict(x, fp_bore=b)) for bn, x in bases for b in (33.1, 33.5, 34.2, 34.925)]
    rows = run_geoms(cases, "bore")
    write("bore", rows)
    show(rows, "front plate bore radius")
    return rows


def figures():
    plt = fd._mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    best = read("chamfers")[0]
    xb = {k: best[k] for k in opt.BASE}
    (RES / "best_chamfers.json").write_text(json.dumps(xb, indent=1))
    picks = [("flat pole cap, gap 21 (previous match)", X0, fd.PALETTE[3]),
             ("steep chamfers, gap 21 (best here)", xb, fd.PALETTE[1])]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    d = np.linspace(0, ts.L_CH, 300)
    for ax, xl in zip(axs, ((0, 1.08), (0.8, 1.04))):
        ax.plot(1 - d / ts.L_CH, ts.B_REF * ts.ref_curve(d), color="k", lw=2.4, label="reference")
        for lab, x, c in picks:
            m, (y, Br, By) = ts.score(x)
            dd = y - Y_EXIT
            sel = (dd >= -0.1 * ts.L_CH) & (dd <= ts.L_CH)
            ax.plot(1 - dd[sel] / ts.L_CH, Br[sel], color=c, lw=1.7,
                    label=f"{lab}: misfit {m['ref_rms_pct']:.1f}%, peak {m['peak_offset_mm']:+.1f} mm")
            ax.plot(1 + m["peak_offset_mm"] / ts.L_CH, ts.B_REF, "v", color=c, ms=7)
        ax.axvline(1, color="0.3", lw=0.8)
        ax.set_xlim(*xl)
        ax.set_xlabel("x / L_ch (0 = anode, 1 = exit plane)")
        ax.set_ylabel("channel-averaged B_r (G)")
    axs[0].set_ylim(0, 190)
    axs[1].set_ylim(120, 172)
    axs[0].legend(frameon=False, fontsize=8, loc="upper left")
    axs[0].set_title("Profiles vs the reference (triangles = peak)", loc="left", fontsize=10)
    axs[1].set_title("Zoom on the exit plane", loc="left", fontsize=10)
    fig.savefig(FIG / "profiles.png", dpi=140)
    plt.close(fig)
    fd.FIG = FIG
    fd.drawing(xb, fname="dimensioned_best",
               note="Best within the rule (nothing past the exit plane): steep chamfers + inner screen tip gap 21.")


if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in what:
        if w in ("chamfers", "all"):
            chamfers()
        if w in ("search", "all"):
            search()
        if w in ("trim", "all"):
            trim()
        if w in ("bore", "all"):
            bore()
        if w in ("figures", "all"):
            figures()
