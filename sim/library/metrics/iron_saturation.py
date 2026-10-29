import numpy as np

from bplsim.metrics import Slot, metric

from het.part_types.post_array import PostArray


def _weighted_percentile(b, w, q):
    o = np.argsort(b)
    cw = np.cumsum(w[o]) / w.sum()
    return float(b[o][min(np.searchsorted(cw, q), len(b) - 1)])


def _flux_profile(z, bz, w, bin_m=1e-3):
    """Axial flux per unit length through 1 mm stations: sum(B_z dV) / dz (exact in 2-D axisymmetry)."""
    zb = np.round(z / bin_m)
    return np.array([np.sum(bz[zb == k] * w[zb == k]) / bin_m for k in np.unique(zb)])


@metric(version="1.0.0", inputs={"parts": "regions[B]"},
        outputs={"worst": "1", "worst_part": "text", "core_flux_B": "T", "posts_flux_B": "T"})
def iron_saturation(parts, design, core: Slot = "inner_core", percentile: float = 95.0):
    """The legacy saturation ratio (optimise.iron_metrics), worst over iron parts:
    - most parts: volume-weighted percentile of |B| over the part / the material's knee;
    - smeared post rings: the ring's axial flux / the real posts' area (mean B in the real posts) / knee;
    - the core: the larger of its percentile and the flux-mean B_z along the shaft / knee."""
    bmag = parts.magnitude("B")
    labels = parts["part"]
    w = parts["weight"]
    z, bz = parts["z"], parts["B_z"]
    worst, worst_part, core_b, posts_b = 0.0, "", 0.0, 0.0
    for slot in sorted(set(labels.tolist())):
        m = labels == slot
        part = design.part(slot)
        if not part.material.has("saturation_knee"):
            continue
        knee = float(part.material.saturation_knee)
        if isinstance(part._instance, PostArray):
            area = part.count * np.pi * part.post_radius ** 2
            b = float(np.max(np.abs(_flux_profile(z[m], bz[m], w[m])))) / area
            posts_b = b
        else:
            b = _weighted_percentile(bmag[m], w[m], percentile / 100)
            if slot == core:
                phi = _flux_profile(z[m], bz[m], w[m])
                core_b = float(np.max(np.abs(phi))) / (np.pi * part.radius ** 2)
                b = max(b, core_b)
        if b / knee > worst:
            worst, worst_part = b / knee, slot
    return {"worst": worst, "worst_part": worst_part, "core_flux_B": core_b, "posts_flux_B": posts_b}
