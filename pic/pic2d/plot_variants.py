"""Compare the PIC runs: the 6 us example and the variants tried on the way (results/variants/*).

    python plot_variants.py        -> results/variants.png
"""
import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
ORDER = ["ignition_attempt", "lowcurrent_3us", "cutoff_walls_2us", "reflect_alpha64", "reflect_bohm_1us"]


def label(name, meta):
    a = meta["args"]
    walls = a.get("walls") or "cutoff"
    walls = f"cutoff {a['wall_cut_ev']:g} eV" if walls == "cutoff" else walls
    start = "seed 1e15" if a["n0"] < 1e16 else f"start {a['n0']:.0e}"
    return f"{name}: walls {walls}, alpha 1/{1 / a['alpha_in']:.0f} in, {start}"


def main():
    runs = [(n, RES / "variants" / n) for n in ORDER if (RES / "variants" / n / "monitor.csv").exists()]
    runs.append(("example (6 us)", RES))
    fig, axs = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    for i, (name, d) in enumerate(runs):
        m = pd.read_csv(d / "monitor.csv")
        meta = json.loads((d / "meta.json").read_text())
        lab = label(name, meta)
        kw = dict(lw=2.4, color="k") if name.startswith("example") else dict(lw=1.2, color=f"C{i}")
        roll = lambda s: s.rolling(8, min_periods=1).mean()          # noqa: E731
        axs[0, 0].plot(m["t_us"], roll(m["I_discharge"]), label=lab, **kw)
        axs[0, 1].plot(m["t_us"], roll(m["I_i_z_hi"]), **kw)
        axs[1, 0].plot(m["t_us"], m["n_e"], **kw)
        axs[1, 1].plot(m["t_us"], roll(m["phi_mid_exit"]), **kw)
    axs[0, 0].axhline(4.83, color="0.5", ls=":", lw=1)
    axs[0, 0].text(0.05, 4.83, "mass-flow equivalent 4.83 A", va="bottom", fontsize=8, color="0.4")
    axs[0, 0].set_yscale("log")
    axs[0, 0].set_ylim(1e-3, 8)
    axs[0, 0].set_title("discharge current (electrons reaching the anode) [A]", fontsize=10)
    axs[0, 1].set_title("ion beam current through the cathode plane [A]", fontsize=10)
    axs[1, 0].set_yscale("log")
    axs[1, 0].set_title("electron macroparticles (the plasma running down)", fontsize=10)
    axs[1, 1].set_title("potential at the exit plane, channel mid-radius [V]", fontsize=10)
    for ax in axs.ravel():
        ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel("t [us]")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle("2-D PIC runs of Design B (argon, 2 mg/s, 300 V): none reaches a self-sustained discharge yet",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.1, 1, 0.97))
    fig.savefig(RES / "variants.png", dpi=140)
    print(RES / "variants.png")


if __name__ == "__main__":
    main()
