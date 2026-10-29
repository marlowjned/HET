"""Mesh the Design B r-z half-plane for estatic.pro (2-D axisymmetric, gmsh x = r, y = z, written in metres).

Domain: a half disc r >= 0 of radius `chamber_r` (the grounded chamber wall), minus the conductors. Physical tags:
  1 vacuum, 2 BN (dielectric regions)
  11 anode surface, 12 iron surface (iron and coils as one conductor), 13 chamber wall (conductor surfaces)
The conductor interiors are not in the mesh file; the axis needs no tag (natural boundary).
"""
import sys
from pathlib import Path

import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry import design_b as g  # noqa: E402

TAGS = {"vacuum": 1, "bn": 2, "anode": 11, "iron": 12, "chamber": 13}


def _bbox(dim, tag):
    return gmsh.model.getBoundingBox(dim, tag)          # xmin, ymin, zmin, xmax, ymax, zmax


def build(msh_path, chamber_r=1000.0, h_channel=0.4, h_thruster=1.5, h_gap=0.03, h_far=None,
          channel_box=(14.0, 50.0, -52.0, 40.0), verbose=False):
    """Write the mesh to `msh_path`; returns element and node counts. Lengths in mm (mesh scaled to m)."""
    h_far = h_far or max(h_thruster, chamber_r / 30.0)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.model.add("design_b_rz")
        occ = gmsh.model.occ
        R = chamber_r
        dom, _ = occ.intersect([(2, occ.addDisk(0, 0, 0, R, R))], [(2, occ.addRectangle(0, -R, 0, R, 2 * R))])
        names, kinds, tools = [], [], []
        for name, (kind, poly) in g.parts().items():
            pts = [occ.addPoint(r, z, 0) for r, z in poly]
            lines = [occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
            tools.append((2, occ.addPlaneSurface([occ.addCurveLoop(lines)])))
            names.append(name)
            kinds.append(kind)
        _, out_map = occ.fragment(dom, tools)
        occ.synchronize()

        frag = {k: [] for k in ("bn", "anode", "iron", "coil")}
        used = set()
        for kind, ents in zip(kinds, out_map[1:]):
            frag[kind] += [t for _, t in ents]
            used |= {t for _, t in ents}
        vacuum = sorted({t for _, t in out_map[0]} - used)
        conductor_iron = frag["iron"] + frag["coil"]

        def outer_boundary(surfs):
            curves = gmsh.model.getBoundary([(2, s) for s in surfs], combined=True, oriented=False)
            return sorted({abs(t) for _, t in curves if _bbox(1, abs(t))[3] > 1e-9})   # drop curves on the axis

        all_surfs = vacuum + frag["bn"] + frag["anode"] + conductor_iron
        chamber = [c for c in outer_boundary(all_surfs) if _bbox(1, c)[3] > 0.99 * R]
        anode_c = outer_boundary(frag["anode"])
        iron_c = outer_boundary(conductor_iron)
        assert chamber and anode_c and iron_c

        m = gmsh.model
        m.addPhysicalGroup(2, vacuum, TAGS["vacuum"], "vacuum")
        m.addPhysicalGroup(2, frag["bn"], TAGS["bn"], "bn")
        m.addPhysicalGroup(1, anode_c, TAGS["anode"], "anode")
        m.addPhysicalGroup(1, iron_c, TAGS["iron"], "iron")
        m.addPhysicalGroup(1, chamber, TAGS["chamber"], "chamber")

        # sizes: fine in the channel and near plume, finer still in the 0.12 mm gap between the pole cap edge
        # face and the BN bore, graded out to the chamber wall
        f = m.mesh.field
        r0, r1, z0, z1 = channel_box
        fb = f.add("Box")
        for k, v in dict(VIn=h_channel, VOut=1e9, XMin=r0, XMax=r1, YMin=z0, YMax=z1, Thickness=5.0).items():
            f.setNumber(fb, k, v)
        ft = f.add("Box")
        for k, v in dict(VIn=h_thruster, VOut=1e9, XMin=0, XMax=100, YMin=-70, YMax=80, Thickness=20.0).items():
            f.setNumber(ft, k, v)
        edge = [c for _, c in m.getEntities(1)
                if abs(_bbox(1, c)[0] - g.CAP_R) < 1e-6 and abs(_bbox(1, c)[3] - g.CAP_R) < 1e-6]
        fd = f.add("Distance")
        f.setNumbers(fd, "CurvesList", edge)
        f.setNumber(fd, "Sampling", 200)
        fg = f.add("Threshold")
        for k, v in dict(InField=fd, SizeMin=h_gap, SizeMax=h_channel, DistMin=0.15, DistMax=1.5,
                         StopAtDistMax=1).items():              # no constraint beyond DistMax (else 0.4 mm everywhere)
            f.setNumber(fg, k, v)
        po = occ.addPoint(0, 0, 0)
        occ.synchronize()
        fo = f.add("Distance")
        f.setNumbers(fo, "PointsList", [po])
        ff = f.add("Threshold")
        for k, v in dict(InField=fo, SizeMin=h_thruster, SizeMax=h_far, DistMin=150, DistMax=R).items():
            f.setNumber(ff, k, v)
        fmin = f.add("Min")
        f.setNumbers(fmin, "FieldsList", [fb, ft, fg, ff])
        f.setAsBackgroundMesh(fmin)
        for k, v in {"Mesh.MeshSizeExtendFromBoundary": 0, "Mesh.MeshSizeFromPoints": 0,
                     "Mesh.MeshSizeFromCurvature": 0, "Mesh.Algorithm": 6, "Mesh.ScalingFactor": 1e-3,
                     "Mesh.MshFileVersion": 4.1}.items():
            gmsh.option.setNumber(k, v)
        m.mesh.generate(2)
        Path(msh_path).parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(msh_path))
        n_tri = sum(len(t[0]) for t in [m.mesh.getElementsByType(2, s) for s in vacuum + frag["bn"]])
        n_nodes = len(m.mesh.getNodes()[0])
        return {"triangles": n_tri, "nodes": n_nodes, "chamber_r_mm": R, "h_channel_mm": h_channel,
                "h_gap_mm": h_gap, "h_thruster_mm": h_thruster, "h_far_mm": h_far}
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    print(build(sys.argv[1] if len(sys.argv) > 1 else "out/design_b_rz.msh", verbose=False))
