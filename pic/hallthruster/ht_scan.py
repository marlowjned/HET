"""Scan the two-zone Bohm coefficients (c1 in the channel, c2 in the plume) at the set point and tabulate
discharge current, thrust, efficiencies and the breathing mode.

    python ht_scan.py        -> runs/scan_*, scan.csv
"""
import argparse
import csv
import itertools
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ht_analyse import summary  # noqa: E402

C1 = [1 / 400, 1 / 250, 1 / 160, 1 / 100, 1 / 64]
C2 = [1 / 32, 1 / 16, 1 / 8]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--duration", type=float, default=2e-3)
    a = ap.parse_args()
    rows = []
    for c1, c2 in itertools.product(C1, C2):
        tag = f"scan_c1_{1 / c1:.0f}_c2_{1 / c2:.0f}"
        r = subprocess.run([sys.executable, str(HERE / "ht_case.py"), "--tag", tag, "--c1", repr(c1), "--c2", repr(c2),
                            "--duration", repr(a.duration), "--num-save", "4000"], capture_output=True, text=True)
        if r.returncode:
            print(tag, "failed:", (r.stdout + r.stderr)[-300:])
            continue
        s = summary(HERE / "runs" / tag)
        s.update(c1_inv=round(1 / c1), c2_inv=round(1 / c2))
        rows.append(s)
        print(f"c1 1/{1 / c1:.0f} c2 1/{1 / c2:.0f}: I_d {s['I_d_mean']:.2f} A, T {s['thrust_mN']:.1f} mN, "
              f"eta_a {s['anode_eff']:.3f}, eta_m {s['mass_eff']:.3f}, f {s.get('f_peak_kHz', float('nan')):.1f} kHz, "
              f"osc {s.get('osc_rel', float('nan')):.0%}", flush=True)
    with open(HERE / "scan.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
