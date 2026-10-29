"""Check estatic.pro against concentric spheres (exact solution), in the same 2-D axisymmetric setting.

Inner sphere radius a (the 'anode', 1 V), BN-filled gap to a conducting shell b1..b2 (the 'iron'), vacuum to the
'chamber' sphere c (0 V). Shell floating: V_shell = C1 / (C1 + C2), Q_inner = C1 C2 / (C1 + C2);
shell grounded: Q_inner = C1, with C1 = 4 pi eps0 eps_r a b1 / (b1 - a), C2 = 4 pi eps0 b2 c / (c - b2).
Checks the BF_GroupOfNodes conductors, the floating constraint, the dielectric and the charge normalisation
(estatic.pro scales GetDP's per-radian axisymmetric charge by 2 pi)."""
import shutil
import subprocess
import sys
from math import pi
from pathlib import Path

import gmsh

HERE = Path(__file__).resolve().parent
EPS0 = 8.8541878128e-12
A, B1, B2, C, EPS_R = 0.1, 0.2, 0.25, 0.5, 4.0
GETDP = str(Path.home() / "tools/getdp/bin/getdp")


def mesh(path, h=0.004):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ

    def half(R):
        out, _ = occ.intersect([(2, occ.addDisk(0, 0, 0, R, R))], [(2, occ.addRectangle(0, -R, 0, R, 2 * R))])
        return out

    out, _ = occ.fragment(half(C), half(A) + half(B1) + half(B2))
    occ.synchronize()
    by_r = {}
    for _, s in gmsh.model.getEntities(2):
        by_r.setdefault(round(gmsh.model.getBoundingBox(2, s)[3], 6), []).append(s)
    arcs = {}
    for _, c in gmsh.model.getEntities(1):
        bb = gmsh.model.getBoundingBox(1, c)
        if bb[3] > 1e-9:
            arcs.setdefault(round(bb[3], 6), []).append(c)
    m = gmsh.model
    m.addPhysicalGroup(2, by_r[C], 1)       # vacuum
    m.addPhysicalGroup(2, by_r[B1], 2)      # dielectric ("BN")
    m.addPhysicalGroup(1, arcs[A], 11)
    m.addPhysicalGroup(1, arcs[B1] + arcs[B2], 12)
    m.addPhysicalGroup(1, arcs[C], 13)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    m.mesh.generate(2)
    gmsh.write(str(path))
    gmsh.finalize()


def solve(work, floating):
    shutil.copy(HERE / "estatic.pro", work / "estatic.pro")      # GetDP writes its outputs next to the .pro
    args = [GETDP, "estatic.pro", "-msh", "s.msh", "-solve", "ES", "-pos", "Post", "-v", "1",
            "-setnumber", "eps_r_bn", str(EPS_R), "-setnumber", "V_anode", "1",
            "-setnumber", "Iron_floating", str(int(floating)), "-setnumber", "V_iron", "0"]
    r = subprocess.run(args, cwd=work, capture_output=True, text=True)
    if r.returncode:
        sys.exit(r.stdout + r.stderr)
    last = lambda f: float((work / f).read_text().split()[-1])          # noqa: E731
    return last("Q_anode.txt"), last("V_iron.txt"), last("Q_iron.txt")


def main():
    work = HERE / "out" / "verify_spheres"
    work.mkdir(parents=True, exist_ok=True)
    mesh(work / "s.msh")
    c1 = 4 * pi * EPS0 * EPS_R * A * B1 / (B1 - A)
    c2 = 4 * pi * EPS0 * B2 * C / (C - B2)
    ok = True
    q, vs, qs = solve(work, floating=True)
    exp_v, exp_q = c1 / (c1 + c2), c1 * c2 / (c1 + c2)
    for name, got, exp in (("floating shell potential [V]", vs, exp_v), ("inner charge, shell floating [C]", q, exp_q),
                           ("shell net charge [C]", qs, 0.0)):
        err = abs(got - exp) / (abs(exp) if exp else 1e-11)
        ok &= err < 2e-3
        print(f"{name:34s} getdp {got: .6e}  exact {exp: .6e}  rel.err {err:.1e}")
    q, vs, qs = solve(work, floating=False)
    err = abs(q - c1) / c1
    ok &= err < 2e-3
    print(f"{'inner charge, shell grounded [C]':34s} getdp {q: .6e}  exact {c1: .6e}  rel.err {err:.1e}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
