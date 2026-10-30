"""1-D reference model of Design B at the set point (HallThruster.jl): krypton 3.0 mg/s, 300 V.

Builds the HallThruster.jl JSON input from the shared geometry and the magnet map (../magnetics: channel-averaged
B_r at 165 G), runs it in Julia, and returns the output. Used to calibrate the anomalous transport to the
3.53 A set point and to predict the breathing mode for the PIC.

    python ht_case.py --c1 0.00625 --c2 0.0625 --duration 1e-3 --tag base
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from geometry import design_b as g  # noqa: E402

JULIA = Path.home() / ".juliaup/bin/julia"
BPROF = HERE.parent / "magnetics/out/bfield_profile_1d.csv"


def bfield_file(path):
    """HallThruster.jl reads a headerless CSV: z from the anode [m], B [T] (channel-averaged B_r here)."""
    d = np.loadtxt(BPROF, delimiter=",", skiprows=1)
    np.savetxt(path, d[:, [0, 2]], delimiter=",", fmt="%.6e")
    return path


def build(args, out):
    L = g.DEPTH * 1e-3
    return {
        "config": {
            "thruster": {
                "name": "BPL-1000 Design B (set point 2026-10-29)",
                "geometry": {"channel_length": L, "inner_radius": g.R_IN * 1e-3, "outer_radius": g.R_OUT * 1e-3},
                "magnetic_field": {"file": str(bfield_file(out / "bfield.csv"))},
            },
            "propellants": [{"gas": "Krypton", "flow_rate_kg_s": args.mdot, "max_charge": args.max_charge,
                             "velocity_m_s": args.v_neutral, "temperature_K": args.t_neutral}],
            "anode_mass_flow_rate": args.mdot,
            "discharge_voltage": args.voltage,
            "domain": [0.0, args.domain],
            "anom_model": {"type": "TwoZoneBohm", "c1": args.c1, "c2": args.c2},
            "wall_loss_model": {"type": "WallSheath", "material": "BoronNitride", "loss_scale": args.wall_scale},
        },
        "simulation": {"dt": 5e-9, "adaptive": True, "grid": {"type": "EvenGrid", "num_cells": args.cells},
                       "duration": args.duration, "num_save": args.num_save},
        "postprocess": {"output_file": str(out / "output.json"), "save_time_resolved": True,
                        "average_start_time": args.average_start or 0.5 * args.duration},
    }


def run(args):
    out = HERE / "runs" / args.tag
    out.mkdir(parents=True, exist_ok=True)
    inp = build(args, out)
    (out / "input.json").write_text(json.dumps(inp, indent=2))
    t = time.time()
    r = subprocess.run([str(JULIA), f"--project={HERE}", "-e",
                        "using HallThruster; HallThruster.run_simulation(ARGS[1]); nothing", str(out / "input.json")],
                       capture_output=True, text=True)
    (out / "julia.log").write_text(r.stdout + r.stderr)
    if r.returncode:
        sys.exit(f"HallThruster.jl failed:\n{(r.stdout + r.stderr)[-3000:]}")
    print(f"{args.tag}: {time.time() - t:.0f} s")
    return json.loads((out / "output.json").read_text())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    a = ap.add_argument
    a("--tag", default="base")
    a("--c1", type=float, default=1 / 160)
    a("--c2", type=float, default=1 / 16)
    a("--mdot", type=float, default=3.0e-6)
    a("--voltage", type=float, default=300.0)
    a("--max-charge", type=int, default=2)
    a("--v-neutral", type=float, default=195.0, help="neutral injection velocity [m/s] (Kr effusive at 600 K)")
    a("--t-neutral", type=float, default=600.0)
    a("--domain", type=float, default=0.08)
    a("--cells", type=int, default=200)
    a("--wall-scale", type=float, default=1.0)
    a("--duration", type=float, default=1e-3)
    a("--num-save", type=int, default=2000)
    a("--average-start", type=float, default=None)
    o = run(ap.parse_args())
    print(list(o.keys()), list(o.get("output", {}).keys()))
