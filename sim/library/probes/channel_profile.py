import math

from bplsim.probes import Angle, Length, Slot, grid_rtz, probe, radial_mean
from bplsim.units import mm


@probe(version="1.1.0", kind="line")
def channel_profile(design, sol, channel: Slot = Slot("channel", type="het:part_types/channel_insulator"),
                    n_r: int = 11, dz: Length = 0.25 * mm, beyond_exit: Length = 30 * mm,
                    behind_anode: Length = 3 * mm, azimuth: Angle | str = "average", n_theta: int = 36):
    """The legacy channel profile: B averaged (r-weighted) over n_r radii spanning the channel, at axial
    positions every dz from beyond_exit downstream of the exit plane back to at most behind_anode behind the
    anode face. The grid is anchored at the downstream end (as the legacy grid), so it hits the exit plane.
    azimuth="average": also averaged over n_theta angles (identical in 2-D)."""
    ch = design.parts[channel]
    dr = (ch.r_out - ch.r_in) / (n_r - 1)
    n_z = math.floor((beyond_exit + ch.depth + behind_anode) / dz + 1e-9)
    z0 = beyond_exit - n_z * dz
    g = grid_rtz(r=(ch.r_in, ch.r_out, dr), z=(z0, beyond_exit, dz), azimuth=azimuth, n_theta=n_theta)
    return radial_mean(sol.sample(g, ["B"]), weight="r")
