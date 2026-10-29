"""Mesh Design B in 3-D for estatic.pro (Flag_Axi = 0), written in metres.

Differences from the 2-D model: the real 6 in square plates (with their bores), the four outer posts with one
coil each, and a cathode. Axisymmetric parts are the 2-D profiles revolved about the axis. The chamber wall is a
sphere centred on the exit plane.

Physical tags:
- dielectric regions: 1 vacuum, 2 BN
- conductor surfaces: 11 anode, 12 iron (iron and coils as one conductor), 13 chamber wall, 14 cathode
The conductor interiors are not in the mesh file.

PLACEHOLDER cathode: a cylinder of radius CATHODE_RADIUS whose axis is parallel to the thruster axis, at x = 110 mm
(beside the front plate); it runs from z = -40 mm to its orifice at +10 mm. The real mount position and angle are
still to come.
"""
import sys
from math import cos, pi, sin
from pathlib import Path

import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry import design_b as g  # noqa: E402

TAGS = {"vacuum": 1, "bn": 2, "anode": 11, "iron": 12, "chamber": 13, "cathode": 14}
CATHODE_X, CATHODE_RADIUS, CATHODE_Z0, CATHODE_Z1 = 110.0, 6.35, -40.0, 10.0


def build(msh_path, chamber_r=500.0, h_channel=1.5, h_body=5.0, h_gap=0.12, h_far=None, threads=4, verbose=False):
    h_far = h_far or chamber_r / 12
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.option.setNumber("General.NumThreads", threads)
        gmsh.model.add("design_b_3d")
        occ = gmsh.model.occ
        parts = []                                                  # (kind, [(3, tag), ...])

        def revolve(poly):
            pts = [occ.addPoint(r, 0, z) for r, z in poly]
            lines = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
            s = occ.addPlaneSurface([occ.addCurveLoop(lines)])
            out = occ.revolve([(2, s)], 0, 0, 0, 0, 0, 1, 2 * pi)
            return [e for e in out if e[0] == 3]

        for name, (kind, poly) in g.parts(posts="none").items():
            if name in ("front_plate", "back_plate"):
                continue
            parts.append((kind, revolve(poly)))
        side = 152.4                                                # the real square plates
        for z0, bore in ((-g.PLATE_T, g.FRONT_BORE_R), (g.BACK_Z - g.PLATE_T, g.BACK_BORE_R)):
            box = occ.addBox(-side / 2, -side / 2, z0, side, side, g.PLATE_T)
            hole = occ.addCylinder(0, 0, z0, 0, 0, g.PLATE_T, bore)
            cut, _ = occ.cut([(3, box)], [(3, hole)])
            parts.append(("iron", cut))
        length = g.BACK_Z - (-g.PLATE_T)
        for i in range(g.POST_COUNT):
            a = 2 * pi * i / g.POST_COUNT
            x, y = g.POST_RING_R * cos(a), g.POST_RING_R * sin(a)
            post = occ.addCylinder(x, y, g.BACK_Z, 0, 0, -length, g.POST_R)
            parts.append(("iron", [(3, post)]))
            outer = occ.addCylinder(x, y, -g.PLATE_T - g.OUTER_LEN, 0, 0, g.OUTER_LEN, g.POST_R + g.OUTER_BUILD)
            inner = occ.addCylinder(x, y, -g.PLATE_T - g.OUTER_LEN, 0, 0, g.OUTER_LEN, g.POST_R)
            coil, _ = occ.cut([(3, outer)], [(3, inner)])
            parts.append(("coil", coil))
        cath = occ.addCylinder(CATHODE_X, 0, CATHODE_Z0, 0, 0, CATHODE_Z1 - CATHODE_Z0, CATHODE_RADIUS)
        parts.append(("cathode", [(3, cath)]))
        sphere = occ.addSphere(0, 0, 0, chamber_r)

        tools = [e for _, ents in parts for e in ents]
        _, omap = occ.fragment([(3, sphere)], tools)
        occ.synchronize()
        frag = {k: [] for k in ("bn", "anode", "iron", "coil", "cathode")}
        used = set()
        idx = 1
        for kind, ents in parts:
            for _ in ents:
                frag[kind] += [t for _, t in omap[idx]]
                used |= {t for _, t in omap[idx]}
                idx += 1
        vacuum = sorted({t for _, t in omap[0]} - used)

        def boundary(vols):
            return sorted({abs(t) for _, t in gmsh.model.getBoundary([(3, v) for v in vols], combined=True,
                                                                      oriented=False)})

        everything = vacuum + frag["bn"] + frag["anode"] + frag["iron"] + frag["coil"] + frag["cathode"]
        m = gmsh.model
        m.addPhysicalGroup(3, vacuum, TAGS["vacuum"], "vacuum")
        m.addPhysicalGroup(3, sorted(set(frag["bn"])), TAGS["bn"], "bn")
        m.addPhysicalGroup(2, boundary(frag["anode"]), TAGS["anode"], "anode")
        m.addPhysicalGroup(2, boundary(frag["iron"] + frag["coil"]), TAGS["iron"], "iron")
        m.addPhysicalGroup(2, boundary(everything), TAGS["chamber"], "chamber")
        m.addPhysicalGroup(2, boundary(frag["cathode"]), TAGS["cathode"], "cathode")

        f = m.mesh.field
        fb = f.add("Box")
        for k, v in dict(VIn=h_channel, VOut=1e9, XMin=-50, XMax=50, YMin=-50, YMax=50, ZMin=-50, ZMax=30,
                         Thickness=10).items():
            f.setNumber(fb, k, v)
        ft = f.add("Box")
        for k, v in dict(VIn=h_body, VOut=1e9, XMin=-125, XMax=125, YMin=-100, YMax=100, ZMin=-75, ZMax=60,
                         Thickness=30).items():
            f.setNumber(ft, k, v)
        def is_edge_face(s):                                        # the pole cap edge face: r = CAP_R, -1.5 <= z <= 0
            bb = m.getBoundingBox(2, s)
            return abs(bb[3] - g.CAP_R) < 1e-6 and abs(bb[2] + g.CAP_EDGE) < 1e-6 and abs(bb[5]) < 1e-6
        edge = [s for _, s in m.getEntities(2) if is_edge_face(s)]
        fd = f.add("Distance")
        f.setNumbers(fd, "SurfacesList", edge)
        fg = f.add("Threshold")
        for k, v in dict(InField=fd, SizeMin=h_gap, SizeMax=h_channel, DistMin=0.15, DistMax=1.0,
                         StopAtDistMax=1).items():
            f.setNumber(fg, k, v)
        po = occ.addPoint(0, 0, 0)
        occ.synchronize()
        fo = f.add("Distance")
        f.setNumbers(fo, "PointsList", [po])
        ff = f.add("Threshold")
        for k, v in dict(InField=fo, SizeMin=h_body, SizeMax=h_far, DistMin=150, DistMax=chamber_r).items():
            f.setNumber(ff, k, v)
        fmin = f.add("Min")
        f.setNumbers(fmin, "FieldsList", [fb, ft, fg, ff])
        f.setAsBackgroundMesh(fmin)
        for k, v in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
                     "Mesh.MeshSizeFromCurvature": 0, "Mesh.Algorithm": 6, "Mesh.Algorithm3D": 10,
                     "Mesh.ScalingFactor": 1e-3}.items():
            gmsh.option.setNumber(k, v)
        m.mesh.generate(3)
        Path(msh_path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(msh_path))
        n_tet = sum(len(m.mesh.getElementsByType(4, v)[0]) for v in vacuum + sorted(set(frag["bn"])))
        return {"tetrahedra": n_tet, "nodes": len(m.mesh.getNodes()[0]), "chamber_r_mm": chamber_r,
                "h_channel_mm": h_channel, "h_gap_mm": h_gap, "edge_surfaces": len(edge)}
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    import time
    t = time.time()
    print(build(sys.argv[1] if len(sys.argv) > 1 else "out/work3d/design_b_3d.msh"), f"{time.time() - t:.0f} s")
