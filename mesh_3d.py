import gmsh
import sys
import numpy as np
import math

gmsh.initialize(sys.argv)
gmsh.model.add("3d")

# geometry
cyl = gmsh.model.occ.addCylinder(0, 0, 5, 0, 0, -10, 3)
# sph = gmsh.model.occ.addSphere(0, 0, 0, 5)


###

# # Parameters
# height = 1.0
# half_angle_deg = 2.5  # half of 5 degree tip angle
# half_angle_rad = math.radians(half_angle_deg)
# radius = height * math.tan(half_angle_rad)

# # Add a cone: apex at (0, 0, 0), base centered at (0, 0, height)
# cone = gmsh.model.occ.addCone(
#     0, 0, 0,      # apex at origin
#     0, 0, height, # height vector
#     0,            # bottom radius (apex)
#     radius        # top radius (at height)
# )

# lc = 1
# angle = 5

# # Front bottom triangle (base at z=0)
# p0 = gmsh.model.occ.addPoint(0, 0, 0, lc)
# p1 = gmsh.model.occ.addPoint(1, 0, 0, lc)
# p2 = gmsh.model.occ.addPoint(math.cos(math.radians(angle)), math.sin(math.radians(angle)), 0, lc)

# # Top triangle (copy shifted in z)
# p3 = gmsh.model.occ.addPoint(0, 0, 1, lc)
# p4 = gmsh.model.occ.addPoint(1, 0, 1, lc)
# p5 = gmsh.model.occ.addPoint(math.cos(math.radians(angle)), math.sin(math.radians(angle)), 1, lc)

# # Lines bottom
# l1 = gmsh.model.occ.addLine(p0, p1)
# l2 = gmsh.model.occ.addLine(p1, p2)
# l3 = gmsh.model.occ.addLine(p2, p0)

# # Lines top
# l4 = gmsh.model.occ.addLine(p3, p4)
# l5 = gmsh.model.occ.addLine(p4, p5)
# l6 = gmsh.model.occ.addLine(p5, p3)

# # Vertical edges
# l7 = gmsh.model.occ.addLine(p0, p3)
# l8 = gmsh.model.occ.addLine(p1, p4)
# l9 = gmsh.model.occ.addLine(p2, p5)

# # Surface loops
# s1 = gmsh.model.occ.addCurveLoop([l1, l2, l3])
# s2 = gmsh.model.occ.addCurveLoop([l4, l5, l6])
# s3 = gmsh.model.occ.addCurveLoop([l1, l8, -l4, -l7])
# s4 = gmsh.model.occ.addCurveLoop([l2, l9, -l5, -l8])
# s5 = gmsh.model.occ.addCurveLoop([l3, l7, -l6, -l9])

# # Plane surfaces
# f1 = gmsh.model.occ.addPlaneSurface([s1])
# f2 = gmsh.model.occ.addPlaneSurface([s2])
# f3 = gmsh.model.occ.addPlaneSurface([s3])
# f4 = gmsh.model.occ.addPlaneSurface([s4])
# f5 = gmsh.model.occ.addPlaneSurface([s5])

# # Create surface loop and volume
# sl = gmsh.model.occ.addSurfaceLoop([f1, f2, f3, f4, f5])
# vol = gmsh.model.occ.addVolume([sl])

###


gmsh.model.occ.synchronize()

# Mesh options
gmsh.option.setNumber("Mesh.Algorithm", 6)     # frontal-delaunay
gmsh.option.setNumber("Mesh.Algorithm3D", 4)   # frontal
# num_runs = 5
# gmsh.option.setNumber("Mesh.Smoothing", num_runs)    # smoothing runs
msf = 3
gmsh.option.setNumber("Mesh.MeshSizeFactor", 1 / 2**(msf - 1))
# gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
gmsh.model.mesh.generate(3)
gmsh.model.mesh.setOrder(2)
# gmsh.model.mesh.optimize("Netgen")

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Save boundary node tags
tags_bdry, _, _ = gmsh.model.mesh.getNodes(2, -1, True, False)
np.savetxt(f"bdryn1{msf}.csv", tags_bdry, delimiter=",")

# Save all node coordinates
_, coords_all, _ = gmsh.model.mesh.getNodes()
coords_all = coords_all.reshape(-1, 3)
np.savetxt(f"coords1{msf}.csv", coords_all, delimiter=",")

# Save node tags for elements
dimTags = gmsh.model.getEntities(3)
eltTypes, _, nodeTags = gmsh.model.mesh.getElements(3, dimTags[0][1])
numNodesPerElt = gmsh.model.mesh.getElementProperties(eltTypes[0])[3]
print("num nodes / element", numNodesPerElt)
nodeTags = np.reshape(nodeTags[0], (-1, numNodesPerElt))
np.savetxt(f"nodes1{msf}.csv", nodeTags, delimiter=",")

gmsh.finalize()
