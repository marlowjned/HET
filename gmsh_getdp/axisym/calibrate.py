"""
calibrate.py

Checks the axisymmetric model against the 3-D solve of the same geometry
(final spec, frozen CAD): the saved 3-D basis fields
(../assembly/sweep_runs_bpl1000/ic0.750_.../basis.npz, channel elements only)
vs the 2-D basis fields, both binned into the same 20 axial stations that
field_quality.axial_profile uses. Both outer-ring smearing modes are tried.

Writes figures/calibration.png and prints per-basis amplitude ratios and
shape errors.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "assembly"))
from field_quality import axial_profile  # noqa: E402
from bpl1000_geometry import CHANNEL as CH3  # noqa: E402
from solve import solve_case  # noqa: E402

FIG = HERE / "figures"
BASIS_3D = (HERE.parent / "assembly" / "sweep_runs_bpl1000" /
            "ic0.750_it0.375_ip0.375_oc0.750_ot0.312_h1.625" / "basis.npz")
N = 20
SPEC_NI = (352.0, 223.4)


def prof3d():
    d = np.load(BASIS_3D)
    ax, az = d["axis"]
    out = {}
    for k in ("inner", "outer"):
        p = axial_profile(d["ch_cent"], d[f"ch_B_{k}"], ax, az, CH3["channel_inner_r"], CH3["channel_outer_r"],
                          CH3["channel_y_anode"], CH3["channel_y_exit"], n_stations=N)
        out[k] = (p["y"] * 1e3, p["B_r"] * 1e4, p["B_axial"] * 1e4)
    return out


def prof2d(res):
    """Same 20 stations from the 2-D channel grid (r-weighted, then averaged
    over each station's y-span)."""
    from metrics import profile
    edges = np.linspace(CH3["channel_y_anode"], CH3["channel_y_exit"], N + 1) * 1e3
    out = {}
    for k in ("inner", "outer"):
        y, br, by = profile(res["channel"][k])
        Br, By = np.empty(N), np.empty(N)
        for i in range(N):
            lo, hi = sorted((edges[i], edges[i + 1]))
            m = (y >= lo) & (y <= hi)
            Br[i], By[i] = br[m].mean(), by[m].mean()
        out[k] = ((edges[:-1] + edges[1:]) / 2, Br, By)
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(exist_ok=True)
    p3 = prof3d()
    modes = {}
    for mode in ("smeared", "thin"):
        res = solve_case(f"calib_{mode}", dict(ring_mode=mode))
        modes[mode] = prof2d(res)

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    for ax, k in zip(axs[:2], ("inner", "outer")):
        y3, br3, by3 = p3[k]
        s3 = np.sign(br3[-1])
        ax.plot(y3, s3 * br3, "ko-", ms=4, label="3-D B_r")
        ax.plot(y3, s3 * by3, "k^--", ms=4, label="3-D B_axial")
        for (mode, p2), c in zip(modes.items(), ("C0", "C3")):
            y2, br2, by2 = p2[k]
            s2 = np.sign(br2[-1])
            ax.plot(y2, s2 * br2, "-", color=c, label=f"2-D {mode} B_r")
            ax.plot(y2, s2 * by2, "--", color=c, label=f"2-D {mode} B_axial")
            amp = br2[-1] * s2 / (br3[-1] * s3)
            shape = np.sqrt(np.mean((s2 * br2 / (s2 * br2[-1]) - s3 * br3 / (s3 * br3[-1])) ** 2)) * 100
            print(f"{k:5s} basis, {mode:7s}: exit B_r 2-D/3-D = {amp:.3f}; "
                  f"shape RMS (normalised to exit) {shape:.1f}%; anode/exit 2-D {br2[0] / br2[-1] * 100:.1f}% "
                  f"vs 3-D {br3[0] / br3[-1] * 100:.1f}%")
        ax.set_title(f"{k} basis (1 A-turn{' per coil' if k == 'outer' else ''})")
        ax.set_xlabel("y (mm); exit -37.8, anode -4.8")
        ax.set_ylabel("G per A-turn")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    ax = axs[2]
    y3, bi3, ai3 = p3["inner"]
    _, bo3, ao3 = p3["outer"]
    br3 = SPEC_NI[0] * bi3 + SPEC_NI[1] * bo3
    s3 = np.sign(br3[-1])
    ax.plot(y3, s3 * br3, "ko-", ms=4, label="3-D")
    for (mode, p2), c in zip(modes.items(), ("C0", "C3")):
        y2, bi2, _ = p2["inner"]
        _, bo2, _ = p2["outer"]
        br2 = SPEC_NI[0] * bi2 + SPEC_NI[1] * bo2
        ax.plot(y2, np.sign(br2[-1]) * br2, "-", color=c, label=f"2-D {mode}")
        print(f"spec excitation {SPEC_NI}, {mode}: exit B_r 2-D {abs(br2[-1]):.0f} G vs 3-D {abs(br3[-1]):.0f} G")
    ax.set_title("Spec excitation 352 / 223.4 A-turns")
    ax.set_xlabel("y (mm)")
    ax.set_ylabel("B_r (G)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(FIG / "calibration.png", dpi=130)


if __name__ == "__main__":
    main()
