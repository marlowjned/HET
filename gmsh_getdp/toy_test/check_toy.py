"""
check_toy.py

Compares the toy solenoid solve against the exact analytic on-axis field
and FAILS if they disagree. Run after getdp.exe:

    python build_toy.py
    getdp.exe magnetostatics_toy.pro -msh toy.msh -solve Res_a -pos Map_b
    python check_toy.py

This exists because the validation it replaces was a number quoted in a
README, and the quoted "exact analytic" value had been computed with the
same wrong current-density area as the FEM it was checking -- so the two
agreed to 1.4% while both were 2.75x low. A comparison that can't fail
isn't a validation. Keep this runnable and keep it in the loop.
"""
import math
import sys

TOL = 0.05  # 5% -- the mesh is coarse by design; the real check is that
            # this is not off by a factor, which is how the bug presented


def read_params(path="toy_params.txt"):
    p = {}
    with open(path) as f:
        for line in f:
            k, _, v = line.strip().partition("=")
            if k:
                p[k] = float(v)
    return p


def read_b_center(path="b_center.txt"):
    """GetDP `Print[b, OnPoint, Format Table]` row: tags, coords, then the
    three field components."""
    with open(path) as f:
        row = f.read().split()
    bx, by, bz = (float(v) for v in row[-3:])
    return bx, by, bz, math.sqrt(bx * bx + by * by + bz * bz)


def main():
    p = read_params()
    bx, by, bz, bmag = read_b_center()
    b_ideal = p["B_ideal"]

    err = (bmag - b_ideal) / b_ideal
    print(f"FEM   |B| at coil center = {bmag*1e3:8.4f} mT   (Bx={bx*1e3:+.4f} By={by*1e3:+.4f} Bz={bz*1e3:+.4f})")
    print(f"exact thick-solenoid     = {b_ideal*1e3:8.4f} mT")
    print(f"error                    = {err*100:+.2f} %  (tolerance +-{TOL*100:.0f}%)")

    # Off-axis components should vanish by symmetry; a large one means the
    # gauge/formulation is misbehaving, which is what this caught originally.
    transverse = math.hypot(bx, by) / bmag
    print(f"transverse fraction      = {transverse*100:.3f} %  (should be ~0 by symmetry)")

    ok = abs(err) <= TOL and transverse < 0.02
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
