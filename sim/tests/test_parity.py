"""Parity of bplsim's getdp_axisym + the het library with the legacy 2-D code (gmsh_getdp/axisym).

Legacy values were computed with the legacy code itself (AXISYM_CONFIG=bn, target_sigma.score on
results/bn/search_stock_stock.json and results/bn_noshield/best.json without the post shim), and the legacy
model's mesh sensitivity (+-20 % mesh size) was measured to set the tolerances (plan §17): see PARITY.md.

Run:  <bplsim venv>/python -m pytest sim/tests -m slow    (needs GetDP: BPLSIM_TOOL_GETDP or PATH)
"""
import os
import sqlite3
from pathlib import Path

import pytest

SIM = Path(__file__).resolve().parents[1]
GETDP = Path(__file__).resolve().parents[2] / "tools" / "getdp-3.5.0-Windows64" / "getdp.exe"

LEGACY = {
    "a": {"misfit.rms_pct": 1.1727446, "shape.anode_Br_pct": 1.7348204, "shape.exit_Br_pct": 94.954972,
          "shape.sigma": 15.2e-3, "shape.peak_z": -2.75e-3, "I_in": 0.91299775, "I_out": 0.72469196,
          "power.total": 1.5737348, "sat.worst": 0.28806643},
    "b": {"misfit.rms_pct": 1.9355217, "shape.anode_Br_pct": 2.7591922, "shape.exit_Br_pct": 95.615039,
          "shape.sigma": 15.4e-3, "shape.peak_z": -2.25e-3, "I_in": 0.92679805, "I_out": 0.73564595,
          "power.total": 1.8640766, "sat.worst": 0.22707838},   # B power: legacy fixed-turns formula
}
PCT_ABS = 0.1          # percentage points, for %-quantities
REL = 0.01             # relative, for physical quantities
GRID = 0.5e-3          # m: peak position is quantised to the 0.25 mm profile grid
SAT_REL = 0.03         # legacy p95 for the core includes the pole-cap zone (one region in legacy)


@pytest.mark.slow
@pytest.mark.parametrize("design", ["a", "b"])
def test_parity(design, tmp_path):
    if not os.environ.get("BPLSIM_TOOL_GETDP") and GETDP.exists():
        os.environ["BPLSIM_TOOL_GETDP"] = str(GETDP)
    from bplsim.execute.local import InlineExecutor
    from bplsim.study.driver import Driver
    d = Driver(SIM, SIM / "studies" / f"design_{design}_parity.yaml", tmp_path / "data", executor=InlineExecutor(),
               log=lambda *a: None)
    s = d.run()
    assert s.get("done") == 1, s
    con = sqlite3.connect(tmp_path / "data" / "trade.db")
    run = con.execute("SELECT run_id FROM row_runs").fetchone()[0]
    m = {r[0]: r[1] for r in con.execute("SELECT name, value_num FROM metrics WHERE run_id=?", [run])}
    p = {r[0]: r[1] for r in con.execute("SELECT path, value_num FROM params WHERE run_id=?", [run])}
    m["I_in"], m["I_out"] = p["normalized.op.circuits.inner.current"], p["normalized.op.circuits.outer.current"]
    for k, legacy in LEGACY[design].items():
        got = m[k]
        if k.endswith("_pct"):
            assert abs(got - legacy) <= PCT_ABS, (k, got, legacy)
        elif k == "shape.peak_z":
            assert abs(got - legacy) <= GRID + 1e-9, (k, got, legacy)
        elif k == "sat.worst":
            assert abs(got - legacy) <= SAT_REL * legacy, (k, got, legacy)
        else:
            assert abs(got - legacy) <= REL * abs(legacy), (k, got, legacy)
