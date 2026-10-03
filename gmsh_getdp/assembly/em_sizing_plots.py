"""Figures for the BPL-1000 EM sizing sweep (em_sizing_sweep.py +
sizing/em_winding_trade.py). Writes PNGs into em_sizing_figures/.

Each figure varies one electromagnet with the other held at the Onshape
defaults, so inner and outer points are separated by which parameters
differ from DEFAULTS_IN."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from onshape.cache import DEFAULTS_IN  # noqa: E402

OUT = HERE / "em_sizing_figures"
KNEE_T = 1.5
D = DEFAULTS_IN


def family(df, em):
    """Rows varying `em` ('iem'/'oem') with everything else at defaults
    (the baseline belongs to both families)."""
    other = "oem" if em == "iem" else "iem"
    m = (df[f"{other}_coreod"] == D[f"{other}_coreod"]) & (df[f"{other}_coilt"] == D[f"{other}_coilt"]) \
        & (df["emag_height"] == D["emag_height"])
    if em == "iem":
        m &= df["iem_platet"] == D["iem_platet"]
    return df[m]


def per_coilt(ax, sub, em, ycol, label):
    for t, g in sorted(sub.groupby(f"{em}_coilt")):
        g = g.sort_values(f"{em}_coreod")
        ax.plot(g[f"{em}_coreod"], g[ycol], "o-", label=f"coilt {t:.3f} in")
    ax.set_xlabel(f"{em} core diameter (in)")
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)


def main():
    OUT.mkdir(exist_ok=True)
    geo = pd.read_csv(HERE / "em_sizing_results.csv")
    wt_path = HERE.parents[1] / "sizing" / "em_winding_trade.csv"
    wt = pd.read_csv(wt_path) if wt_path.exists() else None

    for em, name, at_col, sh_col in (("iem", "inner", "best_req_inner_At", "best_inner_shaft_B_max_T"),
                                      ("oem", "outer", "best_req_outer_At", "best_outer_shaft_B_max_T")):
        sub = family(geo, em)
        if sub.empty:
            continue
        fig, axs = plt.subplots(2, 2, figsize=(11, 8))
        per_coilt(axs[0, 0], sub, em, at_col, f"required {name} A-turns (best ratio)")
        per_coilt(axs[0, 1], sub, em, sh_col, f"max {name} shaft |B| (T)")
        axs[0, 1].axhline(KNEE_T, color="k", ls="--", lw=1, label="knee (linear invalid above)")
        per_coilt(axs[1, 0], sub, em, "best_quality_score", "field-quality score (best ratio)")
        per_coilt(axs[1, 1], sub, em, "best_ratio", "optimal outer/inner A-turn ratio")
        axs[1, 1].axhline(214.5 / 338, color="k", ls=":", lw=1, label="committed 214.5/338")
        for a in axs.flat:
            a.legend(fontsize=8)
        fig.suptitle(f"BPL-1000 {name} electromagnet: core diameter x build (other EM at defaults)")
        fig.tight_layout()
        fig.savefig(OUT / f"{name}_geometry.png", dpi=130)
        plt.close(fig)

        if wt is not None:
            w = wt[(wt.excitation == "best") & (wt.packing == "square") & wt.tag.isin(sub.tag)]
            fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
            for awg, g in w.groupby("awg"):
                g = g.sort_values(f"{name}_em_od_in")
                axs[0].scatter(g[f"{name}_em_od_in"], g["P20_W"], label=f"AWG {awg}", s=18)
                axs[1].scatter(g[f"{name}_em_od_in"], g[f"current_{name}_A"], label=f"AWG {awg}", s=18)
                axs[2].scatter(g[f"{name}_em_od_in"], g["T_coil_hot_C"], label=f"AWG {awg}", s=18)
            axs[0].set_ylabel("total magnet power at 20 C (W)")
            axs[1].set_ylabel(f"{name} channel current (A)")
            axs[2].set_ylabel("hottest coil (C)")
            for a in axs:
                a.set_xlabel(f"{name} EM overall OD (in)")
                a.grid(alpha=0.3)
                a.legend(fontsize=8)
            fig.suptitle(f"{name} EM winding trade, square packing, heavy-build polyimide")
            fig.tight_layout()
            fig.savefig(OUT / f"{name}_winding.png", dpi=130)
            plt.close(fig)

    plate = geo[(geo.iem_coreod == D["iem_coreod"]) & (geo.iem_coilt == D["iem_coilt"])
                & (geo.oem_coreod == D["oem_coreod"]) & (geo.oem_coilt == D["oem_coilt"])]
    if len(plate) > 1:
        plate = plate.sort_values("iem_platet")
        fig, axs = plt.subplots(1, 3, figsize=(14, 4))
        for a, c, l in zip(axs, ("best_req_inner_At", "best_quality_score", "best_inner_shaft_B_max_T"),
                           ("required inner A-turns", "field-quality score", "max inner shaft |B| (T)")):
            a.plot(plate.iem_platet, plate[c], "o-")
            a.set_xlabel("inner pole-plate thickness (in)")
            a.set_ylabel(l)
            a.grid(alpha=0.3)
        fig.suptitle("Inner pole-plate thickness (other parameters at defaults)")
        fig.tight_layout()
        fig.savefig(OUT / "plate_thickness.png", dpi=130)
        plt.close(fig)
    print(f"figures -> {OUT}")


if __name__ == "__main__":
    main()
