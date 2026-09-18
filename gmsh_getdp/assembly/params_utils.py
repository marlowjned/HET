"""
params_utils.py

Parses assembly_params.txt (written once by build_assembly.py) so
downstream scripts -- generate_regions.py, current_sweep.py -- can get
winding geometry (pole type, cross-section area, placement) and channel
geometry without re-running the expensive geometry/meshing step.

Winding entries no longer carry zdir/xref: since build_assembly.py started
sourcing windings as real Onshape coil solids (not Gmsh-built primitives
placed via a hardcoded transform table), every winding's axis is
confirmed global Y directly from geometry, so only its (x,z) axis location
and cross-section area are needed -- see build_assembly.py's docstring.
"""
import re


def _parse_vec(s):
    return tuple(float(x) for x in s.strip("()").split(","))


def parse_params(path="assembly_params.txt"):
    tag = {}
    windings = []
    scalars = {}
    n_windings = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            m = re.match(r"^(\d+) (\w+)$", line)
            if m:
                tag[m.group(2)] = int(m.group(1))
                continue
            if line.startswith("n_windings "):
                n_windings = int(line.split()[1])
                continue
            m = re.match(r"^winding(\d+) pole=(\S+) A_cross=(\S+) loc=(\([^)]*\))$", line)
            if m:
                idx = int(m.group(1))
                windings.append(dict(
                    index=idx,
                    pole=m.group(2),
                    A_cross=float(m.group(3)),
                    loc=_parse_vec(m.group(4)),
                ))
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                scalars[parts[0]] = float(parts[1])

    windings.sort(key=lambda w: w["index"])
    assert len(windings) == n_windings, f"expected {n_windings} windings, parsed {len(windings)}"

    return dict(
        TAG=tag,
        n_windings=n_windings,
        windings=windings,
        axis_x=scalars["axis_x"],
        axis_z=scalars["axis_z"],
        channel_inner_r=scalars["channel_inner_r"],
        channel_outer_r=scalars["channel_outer_r"],
        channel_y_min=scalars["channel_y_min"],
        channel_y_max=scalars["channel_y_max"],
    )


if __name__ == "__main__":
    import pprint
    pprint.pprint(parse_params())
