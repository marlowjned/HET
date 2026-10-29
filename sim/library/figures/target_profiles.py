"""B_r profiles of a study's runs against the half-Gaussian target, on the normalised channel axis."""
import math

import matplotlib.pyplot as plt
import numpy as np

from bplsim.figures import figure


@figure(version="1.0.0")
def target_profiles(ctx, study: str, label_by: str | None = None, max_runs: int = 10,
                    sigma_per_length: float = 0.353, depth_path: str = "parts.channel.depth"):
    df = ctx.table(study)
    df = df[df["status"] == "done"].head(max_runs)
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    depth = None
    for _, r in df.iterrows():
        run = ctx.run(r["run_id"])
        depth = run.value(depth_path)
        p = run.probe("profile")
        z, br = np.ravel(p["z"]), np.ravel(p["B_r"])
        sel = (z >= -1.1 * depth) & (z <= 0.15 * depth)
        lab = f"{label_by} = {r[label_by]:.4g}" if label_by else r["run_id"][:8]
        ax.plot(1 + z[sel] / depth, br[sel] * 1e4, lw=1.3, label=lab)
    if depth:
        d = np.linspace(0, depth, 300)
        sigma = sigma_per_length * depth
        ax.plot(1 - d / depth, 165 * np.exp(-0.5 * (d / sigma) ** 2), color="0.5", lw=5, alpha=0.35,
                label=f"target (sigma {sigma * 1e3:.1f} mm)")
    ax.set_xlabel("x / L_ch (0 = anode face, 1 = exit plane)")
    ax.set_ylabel("B_r (G)")
    ax.legend(fontsize=7)
    return fig
