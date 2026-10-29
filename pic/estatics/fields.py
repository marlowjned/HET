"""Read out/fields_rz.npz (written by run.py) and build the vacuum field of any case by superposition.

    import fields
    d = fields.load()
    phi, er, ez = fields.case(d, v_anode=300, iron="floating")    # or iron="grounded", or iron=<volts>
Arrays are indexed [r, z] on d["r"], d["z"] (metres); potential in V, field in V/m. d["region"]: 0 vacuum, 1 BN,
2 anode, 3 iron and coils.
"""
from pathlib import Path

import numpy as np

DEFAULT = Path(__file__).resolve().parent / "out" / "fields_rz.npz"


def load(path=DEFAULT):
    with np.load(path) as f:
        return {k: f[k] for k in f.files}


def iron_potential(d, v_anode, iron):
    if iron == "floating":
        return float(d["k_float"]) * v_anode                  # zero net charge on the iron
    if iron == "grounded":
        return 0.0
    return float(iron)


def case(d, v_anode=300.0, iron="floating"):
    vi = iron_potential(d, v_anode, iron)
    return tuple(v_anode * d[f"{q}_A"].astype(float) + vi * d[f"{q}_I"].astype(float) for q in ("phi", "Er", "Ez"))
