"""
field_quality.py

Quality metric for the discharge-channel B-field, for comparing geometries
in a sizing trade (emag_height / coil ID / coil OD sweep). Target: ~300 G
radially outward at the channel exit, weaker toward the anode -- the
standard Hall-thruster magnetic-lens shape for E x B electron confinement
(strong field near the exit, tapering off toward the anode to limit anode
electron loss/erosion).

Decoupled from any specific solve on purpose: takes raw (position, B)
samples plus explicit channel geometry, not a particular .pos/.vtu file --
so it can be unit-tested with synthetic fields (see `python
field_quality.py`) now, before the new (emag_height-parametrized) CAD has
been through the Gmsh/GetDP pipeline, and pointed at real solve output
later without changes. Axis/channel-annulus convention (axis_x, axis_z,
inner/outer radius, y-bounds) matches current_sweep.py's channel_mask, not
reinvented.

NOT YET VALIDATED AGAINST A REAL SOLVE. The target curve's B_floor/k and
the composite score's weights are first-pass defaults, explicitly meant to
be retuned once real swept-geometry field data exists to look at.
"""
import numpy as np


def cylindrical_components(positions: np.ndarray, B: np.ndarray, axis_x: float, axis_z: float):
    """
    positions, B: (N,3) arrays, global (x,y,z) -- thrust axis is global Y,
    passing through (axis_x, axis_z) in the X-Z plane (matches
    current_sweep.py's channel_mask convention).

    Returns (r, B_r, B_axial, B_theta, Bmag), each (N,).
    B_r: radially outward component. B_axial: along-thrust-axis component.
    B_theta: magnitude of whatever's left (out-of-plane "leakage") --
    sign is not physically meaningful here, only used as a leakage size.
    """
    dx = positions[:, 0] - axis_x
    dz = positions[:, 2] - axis_z
    r = np.hypot(dx, dz)
    r_safe = np.where(r > 1e-12, r, 1.0)
    rhat_x, rhat_z = dx / r_safe, dz / r_safe

    B_r = B[:, 0] * rhat_x + B[:, 2] * rhat_z
    B_axial = B[:, 1]
    Bmag = np.linalg.norm(B, axis=1)
    B_theta_sq = np.clip(Bmag**2 - B_r**2 - B_axial**2, 0.0, None)
    B_theta = np.sqrt(B_theta_sq)

    return r, B_r, B_axial, B_theta, Bmag


def axial_profile(
    positions: np.ndarray, B: np.ndarray,
    axis_x: float, axis_z: float,
    r_in: float, r_out: float,
    y_anode: float, y_exit: float,
    n_stations: int = 20,
):
    """
    Bins channel samples into n_stations evenly-spaced axial slices from
    y_anode to y_exit (need not have y_anode < y_exit numerically -- sign
    handled internally), averaging B_r/B_axial/B_theta/Bmag per slice.

    Returns dict with:
      s: (n_stations,) normalized station coordinate, 0=anode, 1=exit
      y: (n_stations,) actual y-coordinate per station
      B_r, B_axial, B_theta, Bmag: (n_stations,) slice-averaged values
      counts: (n_stations,) number of samples per slice (sanity check --
        a station with ~0 samples means the mesh didn't resolve that slice
        or the channel geometry bounds are wrong)
    """
    dx = positions[:, 0] - axis_x
    dz = positions[:, 2] - axis_z
    r = np.hypot(dx, dz)
    y = positions[:, 1]

    in_channel = (r >= r_in) & (r <= r_out)
    y_lo, y_hi = min(y_anode, y_exit), max(y_anode, y_exit)
    in_channel &= (y >= y_lo) & (y <= y_hi)

    _, B_r_all, B_axial_all, B_theta_all, Bmag_all = cylindrical_components(positions, B, axis_x, axis_z)

    edges = np.linspace(y_anode, y_exit, n_stations + 1)
    s = (np.arange(n_stations) + 0.5) / n_stations
    y_mid = (edges[:-1] + edges[1:]) / 2

    B_r = np.full(n_stations, np.nan)
    B_axial = np.full(n_stations, np.nan)
    B_theta = np.full(n_stations, np.nan)
    Bmag = np.full(n_stations, np.nan)
    counts = np.zeros(n_stations, dtype=int)

    lo_edges = np.minimum(edges[:-1], edges[1:])
    hi_edges = np.maximum(edges[:-1], edges[1:])

    for i in range(n_stations):
        mask = in_channel & (y >= lo_edges[i]) & (y < hi_edges[i] if i < n_stations - 1 else y <= hi_edges[i])
        counts[i] = mask.sum()
        if counts[i] > 0:
            B_r[i] = B_r_all[mask].mean()
            B_axial[i] = B_axial_all[mask].mean()
            B_theta[i] = B_theta_all[mask].mean()
            Bmag[i] = Bmag_all[mask].mean()

    return {"s": s, "y": y_mid, "B_r": B_r, "B_axial": B_axial, "B_theta": B_theta, "Bmag": Bmag, "counts": counts}


def target_profile(s: np.ndarray, B_peak: float, B_floor: float, k: float) -> np.ndarray:
    """B_target(s) = B_floor + (B_peak - B_floor) * exp(-k*(1-s)), s in [0,1] anode->exit."""
    return B_floor + (B_peak - B_floor) * np.exp(-k * (1.0 - s))


def evaluate_field_quality(
    profile: dict,
    B_peak: float = 300.0,   # Gauss
    B_floor: float = 30.0,   # Gauss, default 10% of peak
    k: float = 3.0,
    weights: dict | None = None,
) -> dict:
    """
    profile: output of axial_profile(), B_r/Bmag in the SAME units as
    B_peak/B_floor (Gauss by convention here -- convert before calling if
    your solve reports Tesla).

    Returns a dict of interpretable sub-metrics plus one composite score
    (0-100, higher = better). All defaults (B_floor, k, weights) are
    first-pass and meant to be retuned once real data exists.
    """
    if weights is None:
        weights = {"shape": 0.35, "purity": 0.35, "peak": 0.15, "monotonic": 0.15}

    s, B_r, Bmag = profile["s"], profile["B_r"], profile["Bmag"]
    valid = ~np.isnan(B_r)
    if not valid.all():
        n_bad = (~valid).sum()
        print(f"WARNING: {n_bad}/{len(valid)} station(s) had zero samples -- "
              f"check channel geometry bounds. Excluding from metric.")
    s, B_r, Bmag = s[valid], B_r[valid], Bmag[valid]

    target = target_profile(s, B_peak, B_floor, k)

    # 1. Peak error at exit (last valid station)
    peak_error_pct = (B_r[-1] - B_peak) / B_peak * 100.0

    # 2. RMS shape error vs. target curve, normalized by peak
    shape_rmse_pct = np.sqrt(np.mean((B_r - target) ** 2)) / B_peak * 100.0

    # 3. Radial purity: mean(B_r / |B|) -- unweighted across stations by default
    radial_purity_pct = float(np.mean(B_r / np.where(Bmag > 1e-12, Bmag, 1.0)) * 100.0)

    # 4. Monotonicity violation: B_r should rise moving anode -> exit (s increasing).
    #    Violation at station i->i+1 if B_r actually drops.
    drops = np.clip(B_r[:-1] - B_r[1:], 0.0, None)
    monotonicity_violation_pct = float(np.sum(drops) / B_peak * 100.0)

    goodness_peak = max(0.0, 100.0 - abs(peak_error_pct))
    goodness_shape = max(0.0, 100.0 - shape_rmse_pct)
    goodness_purity = radial_purity_pct
    goodness_monotonic = max(0.0, 100.0 - monotonicity_violation_pct)

    composite_score = (
        weights["peak"] * goodness_peak
        + weights["shape"] * goodness_shape
        + weights["purity"] * goodness_purity
        + weights["monotonic"] * goodness_monotonic
    )

    return {
        "peak_error_pct": float(peak_error_pct),
        "shape_rmse_pct": float(shape_rmse_pct),
        "radial_purity_pct": radial_purity_pct,
        "monotonicity_violation_pct": monotonicity_violation_pct,
        "composite_score": float(composite_score),
        "target_curve": target,
        "s": s,
        "B_r": B_r,
    }


if __name__ == "__main__":
    # Self-test with synthetic fields -- no solve required. Validates the
    # metric logic itself: a field built FROM the target curve should score
    # near-perfect; deliberately bad fields should score poorly on the
    # specific sub-metric they violate.
    rng = np.random.default_rng(0)
    axis_x, axis_z = 0.0, 0.0
    r_in, r_out = 0.030, 0.040   # 30-40mm channel annulus, arbitrary for the test
    y_anode, y_exit = 0.0, 0.080  # 80mm channel length, anode at y=0

    n_samples = 20000
    r = rng.uniform(r_in, r_out, n_samples)
    theta = rng.uniform(0, 2 * np.pi, n_samples)
    y = rng.uniform(y_anode, y_exit, n_samples)
    x = axis_x + r * np.cos(theta)
    z = axis_z + r * np.sin(theta)
    positions = np.stack([x, y, z], axis=1)
    s_true = (y - y_anode) / (y_exit - y_anode)

    def make_B(radial_mag, axial_mag=0.0, theta_leak=0.0):
        rhat_x, rhat_z = np.cos(theta), np.sin(theta)
        Bx = radial_mag * rhat_x
        Bz = radial_mag * rhat_z
        By = np.full_like(radial_mag, axial_mag)
        # add some out-of-plane leakage along the tangential direction
        that_x, that_z = -np.sin(theta), np.cos(theta)
        Bx = Bx + theta_leak * that_x
        Bz = Bz + theta_leak * that_z
        return np.stack([Bx, By, Bz], axis=1)

    print("=== Case 1: ideal field (built exactly from the target curve) ===")
    B_ideal_mag = target_profile(s_true, B_peak=300.0, B_floor=30.0, k=3.0)
    B = make_B(B_ideal_mag)
    profile = axial_profile(positions, B, axis_x, axis_z, r_in, r_out, y_anode, y_exit, n_stations=20)
    result = evaluate_field_quality(profile)
    for k_, v in result.items():
        if k_ not in ("target_curve", "s", "B_r"):
            print(f"  {k_}: {v:.2f}")

    print("\n=== Case 2: flat field (300G everywhere, no taper -- should hurt shape/monotonicity) ===")
    B = make_B(np.full(n_samples, 300.0))
    profile = axial_profile(positions, B, axis_x, axis_z, r_in, r_out, y_anode, y_exit, n_stations=20)
    result = evaluate_field_quality(profile)
    for k_, v in result.items():
        if k_ not in ("target_curve", "s", "B_r"):
            print(f"  {k_}: {v:.2f}")

    print("\n=== Case 3: correct shape but heavy axial leakage (should hurt purity) ===")
    B = make_B(B_ideal_mag, axial_mag=150.0)
    profile = axial_profile(positions, B, axis_x, axis_z, r_in, r_out, y_anode, y_exit, n_stations=20)
    result = evaluate_field_quality(profile)
    for k_, v in result.items():
        if k_ not in ("target_curve", "s", "B_r"):
            print(f"  {k_}: {v:.2f}")

    print("\n=== Case 4: correct shape but only 150G peak (should hurt peak error, not shape) ===")
    B = make_B(B_ideal_mag * 0.5)
    profile = axial_profile(positions, B, axis_x, axis_z, r_in, r_out, y_anode, y_exit, n_stations=20)
    result = evaluate_field_quality(profile)
    for k_, v in result.items():
        if k_ not in ("target_curve", "s", "B_r"):
            print(f"  {k_}: {v:.2f}")
