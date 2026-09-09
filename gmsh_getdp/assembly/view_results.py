"""
view_results.py

Opens the real BPL-700 mesh and solved |B| field in Gmsh's own interactive
3-D GUI (bundled with the pip gmsh package -- no separate install needed).
Lets you rotate/zoom/clip the real geometry and the field together, unlike
post_process.py's flat cross-section plot.

Run: python view_results.py
(opens a native window -- run this from a real desktop session, not a
headless/SSH one)
"""
import gmsh

gmsh.initialize()
gmsh.open("assembly.msh")
gmsh.open("b_assembly.pos")   # |B| on iron+air+windings
gmsh.open("b_iron.pos")       # |B| on iron only

# A few sane defaults: hide the raw mesh edges, show the field as a color
# map with a clip plane through the assembly center so you can see inside.
gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
gmsh.option.setNumber("Mesh.VolumeEdges", 0)
gmsh.option.setNumber("View[0].RangeType", 2)   # custom range
gmsh.option.setNumber("View[0].CustomMax", 3.0)  # clip color scale at 3T (see README: peak is ~9.7T,
gmsh.option.setNumber("View[0].CustomMin", 0.0)  # an artifact of the linear no-saturation iron model)

print("Opening Gmsh GUI... View menu on the left lets you toggle b_assembly/b_iron,")
print("switch between vector/scalar/iso-surface display, and add clip planes.")
gmsh.fltk.run()
gmsh.finalize()
