"""
build_toy.py

Minimal validation case for the Gmsh + GetDP magnetostatics pipeline,
mirroring the center-coil parameters from ../../matlab/het_solenoid_bfield.m
(OD=2.5in, ID=1.0in, height=2.0in, 300 turns, 5A) so the FEM result can be
checked against the ideal-solenoid estimate B = mu0*N*I/L before trusting
this formulation on the real CAD assembly.

All geometry construction (Gmsh/OpenCASCADE) and meshing happens here in
Python only as orchestration -- the actual mesh generation is Gmsh's
compiled C++ kernel. The FEM solve itself happens later in getdp.exe
(separately, natively compiled, not Python).
"""
import gmsh
import math

in2m = 0.0254

OD = 2.5 * in2m
ID = 1.0 * in2m
height = 2.0 * in2m
turns = 300
current = 5.0  # A
mu0 = 4 * math.pi * 1e-7

A_cross = math.pi * ((OD / 2) ** 2 - (ID / 2) ** 2)
Jmag = turns * current / A_cross
B_ideal = mu0 * turns * current / height

print(f"A_cross = {A_cross:.6e} m^2")
print(f"Jmag    = {Jmag:.6e} A/m^2")
print(f"B_ideal (center, on-axis) = {B_ideal*1e3:.4f} mT")

air_radius = 8 * max(OD / 2, height / 2)

gmsh.initialize()
gmsh.model.add("toy_solenoid")
occ = gmsh.model.occ

# Air sphere, centered at origin
air = occ.addSphere(0, 0, 0, air_radius)

# Coil: hollow cylinder (outer cyl minus inner cyl), centered on z-axis,
# spanning z in [-height/2, height/2]
outer_cyl = occ.addCylinder(0, 0, -height / 2, 0, 0, height, OD / 2)
inner_cyl = occ.addCylinder(0, 0, -height / 2, 0, 0, height, ID / 2)
coil_cut = occ.cut([(3, outer_cyl)], [(3, inner_cyl)])
coil = coil_cut[0][0][1]

occ.synchronize()

# Fragment so the coil sits as a conformal sub-volume of the air sphere
out, out_map = occ.fragment([(3, air)], [(3, coil)])
occ.synchronize()

# out_map[0] = resulting tag(s) for the air input, out_map[1] = for coil input
air_tags = [t for (d, t) in out_map[0]]
coil_tags = [t for (d, t) in out_map[1]]
print("air fragments:", air_tags)
print("coil fragments:", coil_tags)

assert len(coil_tags) == 1, f"expected coil to stay one volume, got {coil_tags}"
# air fragment(s): should be the leftover air volume(s) after removing coil
air_only = [t for t in air_tags if t not in coil_tags]

# Physical groups (volumes). GetDP's Region[] only understands integer
# tags (no physical-name string lookup), so tags are pinned explicitly
# here and must match the numbers hardcoded in magnetostatics_toy.pro.
TAG_AIR, TAG_COIL, TAG_OUTERBND = 1, 2, 3
pg_air = gmsh.model.addPhysicalGroup(3, air_only, tag=TAG_AIR)
gmsh.model.setPhysicalName(3, pg_air, "Air")
pg_coil = gmsh.model.addPhysicalGroup(3, coil_tags, tag=TAG_COIL)
gmsh.model.setPhysicalName(3, pg_coil, "Coil")

# Outer boundary surface = boundary of the air sphere volume that is NOT
# shared with the coil (i.e. the true outer sphere surface)
air_boundary = gmsh.model.getBoundary([(3, t) for t in air_only], oriented=False)
coil_boundary = gmsh.model.getBoundary([(3, t) for t in coil_tags], oriented=False)
coil_surf_tags = set(t for (d, t) in coil_boundary)
outer_surf_tags = [t for (d, t) in air_boundary if t not in coil_surf_tags]
print("outer boundary surfaces:", outer_surf_tags)

pg_outer = gmsh.model.addPhysicalGroup(2, outer_surf_tags, tag=TAG_OUTERBND)
gmsh.model.setPhysicalName(2, pg_outer, "OuterBnd")

# Mesh sizing: coarse globally, finer on the coil wall thickness
wall = (OD - ID) / 2
gmsh.option.setNumber("Mesh.MeshSizeMin", wall / 3)
gmsh.option.setNumber("Mesh.MeshSizeMax", air_radius / 4)

# Refine near the coil using a distance+threshold field
coil_surfs = [t for (d, t) in coil_boundary]
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", coil_surfs)
gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", wall / 3)
gmsh.model.mesh.field.setNumber(2, "SizeMax", air_radius / 4)
gmsh.model.mesh.field.setNumber(2, "DistMin", wall)
gmsh.model.mesh.field.setNumber(2, "DistMax", air_radius / 3)
gmsh.model.mesh.field.setAsBackgroundMesh(2)

gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("toy.msh")

n_nodes = len(gmsh.model.mesh.getNodes()[0])
print(f"Mesh: {n_nodes} nodes")

gmsh.finalize()

with open("toy_params.txt", "w") as f:
    f.write(f"Jmag={Jmag!r}\nB_ideal={B_ideal!r}\n")
