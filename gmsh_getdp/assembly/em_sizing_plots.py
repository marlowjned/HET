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
        & (df["emag_height"] == D["emag_height"]) & (df["iem_platet"] == D["iem_platet"])
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

    for em, name, at_col, sh_col in (("iem", "inner", "ref_req_inner_At", "ref_inner_shaft_B_max_T"),
                                      ("oem", "outer", "ref_req_outer_At", "ref_outer_shaft_B_max_T")):
        sub = family(geo, em)
        if sub.empty:
            continue
        fig, axs = plt.subplots(2, 2, figsize=(11, 8))
        per_coilt(axs[0, 0], sub, em, at_col, f"required {name} A-turns (committed 214.5/338 ratio)")
        per_coilt(axs[0, 1], sub, em, sh_col, f"max {name} shaft |B| (T)")
        axs[0, 1].axhline(KNEE_T, color="k", ls="--", lw=1, label="knee (linear invalid above)")
        per_coilt(axs[1, 0], sub, em, "ref_quality_score", "field-quality score (committed ratio)")
        per_coilt(axs[1, 1], sub, em, "best_ratio", "optimal outer/inner A-turn ratio")
        axs[1, 1].axhline(214.5 / 338, color="k", ls=":", lw=1, label="committed 214.5/338")
        for a in axs.flat:
            a.legend(fontsize=8)
        fig.suptitle(f"BPL-1000 {name} electromagnet: core diameter x build (other EM at defaults)")
        fig.tight_layout()
        fig.savefig(OUT / f"{name}_geometry.png", dpi=130)
        plt.close(fig)

        if wt is not None:
            w = wt[(wt.excitation == "ref") & (wt.packing == "square") & wt.tag.isin(sub.tag)]
            # This magnet's own coil power (I x V at 20 C), one panel per
            # gauge, one line per coil thickness: at a given overall OD the
            # thinner core with the thicker coil wins.
            w = w.assign(P_own=w[f"current_{name}_A"] * w[f"V20_{name}"])
            awgs = sorted(w.awg.unique())
            fig, axs = plt.subplots(1, len(awgs), figsize=(4 * len(awgs), 4.2), sharey=True)
            for ax, awg in zip(axs, awgs):
                for t, g in sorted(w[w.awg == awg].groupby(f"{em}_coilt")):
                    g = g.sort_values(f"{name}_em_od_in")
                    ax.plot(g[f"{name}_em_od_in"], g["P_own"], "o-", label=f"coil {t:.3f} in")
                    for _, r in g.iterrows():
                        ax.annotate(f"{r[f'{em}_coreod']:.3g}", (r[f"{name}_em_od_in"], r["P_own"]),
                                    textcoords="offset points", xytext=(3, 3), fontsize=6, color="0.4")
                ax.set_title(f"AWG {awg}")
                ax.set_xlabel(f"{name} EM overall OD (in)")
                ax.grid(alpha=0.3)
            axs[0].set_ylabel(f"{name} coil power at 20 C (W)" + (" (all 4)" if name == "outer" else ""))
            axs[0].legend(fontsize=8)
            fig.suptitle(f"{name} EM: coil power vs overall diameter (labels = core dia, in; committed ratio, square packing)")
            fig.tight_layout()
            fig.savefig(OUT / f"{name}_winding.png", dpi=130)
            plt.close(fig)

    plate = geo[(geo.iem_coreod == D["iem_coreod"]) & (geo.iem_coilt == D["iem_coilt"])
                & (geo.oem_coreod == D["oem_coreod"]) & (geo.oem_coilt == D["oem_coilt"])]
    if len(plate) > 1:
        plate = plate.sort_values("iem_platet")
        fig, axs = plt.subplots(1, 3, figsize=(14, 4))
        for a, c, l in zip(axs, ("ref_req_inner_At", "ref_quality_score", "ref_inner_shaft_B_max_T"),
                           ("required inner A-turns (committed ratio)", "field-quality score (committed ratio)", "max inner shaft |B| (T)")):
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
