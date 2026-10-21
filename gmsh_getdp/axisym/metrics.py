"""
metrics.py

Channel B-field shape metrics for the target the user set on 2026-10-15:
B_r peaks EXACTLY at the exit plane and decays toward the anode to as little
as possible (a half-Gaussian, peak at the exit).

Inputs are the channel-grid samples from solve.py: (N, 4) arrays of r, y,
B_r, B_y (m, T) per basis. A profile is the volume (r-weighted) mean over the
channel's radial span at each y.

Reported per excitation, scaled so the PEAK B_r is TARGET_G:
  peak_offset_mm    peak position minus exit plane, + = downstream (outside
                    the channel), - = inside. Target 0.
  anode_Br_pct      B_r at the anode face, % of peak. Target 0.
  anode_B_pct       |B| at the anode face, % of peak (B_r alone can be
                    nulled while the axial component is not).
  exit_Br_pct       B_r at the exit plane, % of peak (100 if the peak is on it)
  half_depth_mm     distance from the exit into the channel where B_r first
                    falls to half the exit value -- how concentrated the
                    field is at the exit. Smaller = sharper.
  sigma_mm, gauss_rms_pct
                    least-squares half-Gaussian B_exit exp(-d^2 / 2 sigma^2)
                    over the channel (d = depth from the exit), and the RMS
                    misfit as % of peak
  NI_inner, NI_outer
                    ampere-turns for TARGET_G peak (outer per coil)
"""
import numpy as np

from geometry import CHANNEL

TARGET_G = 300.0
MM = 1e-3


def profile(ch):
    """(N, 4) r, y, Br, By samples -> y (mm, ascending), Br, By (G per unit
    excitation), r-weighted means across the channel."""
    r, y, br, by = ch.T
    ys = np.unique(np.round(y, 9))
    Br = np.empty(len(ys))
    By = np.empty(len(ys))
    for i, yy in enumerate(ys):
        m = np.abs(y - yy) < 1e-9
        w = r[m]
        Br[i] = np.sum(w * br[m]) / w.sum() * 1e4
        By[i] = np.sum(w * by[m]) / w.sum() * 1e4
    return ys / MM, Br, By


def combine(res, ni_in, ni_out):
    """Profile for a given excitation (A-turns; outer per coil)."""
    y, bri, byi = profile(res["channel"]["inner"])
    _, bro, byo = profile(res["channel"]["outer"])
    return y, ni_in * bri + ni_out * bro, ni_in * byi + ni_out * byo


def shape_metrics(y, Br, By, target_G=TARGET_G):
    """Metrics for one profile (any scale). Returns metrics and the scale
    factor applied to reach target_G at the peak."""
    ye, ya = CHANNEL["y_exit"], CHANNEL["y_anode"]
    sgn = np.sign(Br[np.argmin(np.abs(y - ye))]) or 1.0
    Br, By = sgn * Br, sgn * By
    # peak searched from 15 mm downstream of the exit to the anode
    win = (y >= ye - 15) & (y <= ya)
    i_pk = np.flatnonzero(win)[np.argmax(Br[win])]
    s = target_G / Br[i_pk]
    Br, By = Br * s, By * s
    pk = Br[i_pk]

    def at(yy, v):
        return float(np.interp(yy, y, v))

    B_exit = at(ye, Br)
    inside = (y >= ye) & (y <= ya)
    d = y[inside] - ye                      # depth into the channel, mm
    bi = Br[inside]
    below = np.flatnonzero(bi <= 0.5 * B_exit)
    half = float(d[below[0]]) if len(below) else float("nan")
    # half-Gaussian fit with the amplitude pinned at the exit value
    sig_grid = np.linspace(2, 80, 781)
    errs = [np.sqrt(np.mean((bi - B_exit * np.exp(-d ** 2 / (2 * sg ** 2))) ** 2)) for sg in sig_grid]
    j = int(np.argmin(errs))
    return dict(
        peak_offset_mm=float(ye - y[i_pk]),
        anode_Br_pct=at(ya, Br) / pk * 100,
        anode_B_pct=float(np.hypot(at(ya, Br), at(ya, By)) / pk * 100),
        exit_Br_pct=B_exit / pk * 100,
        half_depth_mm=half,
        sigma_mm=float(sig_grid[j]),
        gauss_rms_pct=float(errs[j] / pk * 100),
    ), s * sgn


def evaluate(res, q, target_G=TARGET_G):
    """Metrics at outer/inner A-turn ratio q (outer per coil), scaled to
    target_G peak, with the required A-turns and peak iron |B|."""
    y, Br, By = combine(res, 1.0, q)
    m, s = shape_metrics(y, Br, By, target_G)
    ni_in, ni_out = abs(s), abs(s) * q
    m.update(q=q, NI_inner=float(ni_in), NI_outer=float(ni_out))
    bi, bo = res["iron"]["inner"], res["iron"]["outer"]
    B = ni_in * bi[:, 2:4] + ni_out * bo[:, 2:4]
    m["iron_B_max_T"] = float(np.percentile(np.linalg.norm(B, axis=1), 99.9))
    return m, (y, Br * s, By * s)
