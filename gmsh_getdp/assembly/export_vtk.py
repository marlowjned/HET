"""
export_vtk.py

Exports the GetDP magnetostatics results (b_assembly.pos, b_iron.pos) to
VTK (.vtu) for use in ParaView, other simulation tools, or any other
VTK-standard consumer -- more interoperable than Gmsh's own .pos ASCII
format (see README.md's file-format notes) for layering multiple
simulation tools on the same result.

Gmsh's own gmsh.view.write(tag, "x.vtk") was tried first and rejected: it
silently writes a non-standard internal list-dump format when asked for
a "list-based" view (which is what a .pos file becomes once reloaded --
self-contained per-element geometry+value, not tied to the model mesh's
node IDs) -- no error, just not real VTK, confirmed by inspecting the
output (no "# vtk DataFile" header, no POINTS/CELLS sections).

Rather than fight that (or risk a silent element-ordering mismatch
correlating a separately-read assembly.msh against b_assembly.pos by
position in file), this builds the VTU directly and only from the .pos
data itself: each element already carries its own 4 corner-node
coordinates (see pos_utils.parse_pos_elements), so a valid unstructured
grid can be built without needing the mesh file's node numbering at all.
The one tradeoff: points aren't deduplicated/shared between neighboring
tets (each tet gets its own 4 points), so the file is larger than a
minimal mesh -- but this is actually a reasonable match for the physics
anyway, since B is genuinely piecewise-constant/discontinuous across
elements (see HOW_IT_WORKS.md sec. 9), not something that should be
smoothed across a shared node.
"""
import numpy as np
import meshio

from pos_utils import parse_pos_elements


def build_vtu(pos_path, out_path):
    nodes, B = parse_pos_elements(pos_path)
    n = len(nodes)
    print(f"{pos_path}: {n} elements")

    points = nodes.reshape(-1, 3)  # (4N, 3), NOT deduplicated -- see docstring
    cells = [("tetra", np.arange(4 * n).reshape(n, 4))]

    Bmag = np.linalg.norm(B, axis=1)
    cell_data = {"B": [B], "B_magnitude": [Bmag]}

    mesh = meshio.Mesh(points=points, cells=cells, cell_data=cell_data)
    mesh.write(out_path)
    print(f"Wrote {out_path} ({points.shape[0]} points, {n} cells)")


if __name__ == "__main__":
    build_vtu("b_assembly.pos", "assembly_field.vtu")
    build_vtu("b_iron.pos", "iron_field.vtu")
