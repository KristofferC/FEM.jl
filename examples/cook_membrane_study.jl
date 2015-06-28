using FEM

# Load the VTK exporter
FEM.vtkexportmod()
using FEM.VTKExportMod

# Generates a mesh of the shape known as the "Cook membrane" which is a
# quadraterial with corners [0.0, 0.0], [0.0, 44.0], [48.0, 60.0], [48.0, 44.0]
# Possible mesh elements are
# - GeoTrig for 3 node triangles,
# - GeoQTrig for 6 node triangles
# - GeoQuad for 4 node quadraterials
n_ele_x = 100
n_ele_y = 100
geomesh = gencook(n_ele_x, n_ele_y, GeoQTrig)


# We create two node sets, one on the right edge and one on the left.
# This is done by giving an anonymous function that is satisfied by
# the edges.
push!(geomesh, gennodeset(n->n.coords.x>47.999999, "right", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords.x<0.00001, "left", geomesh.nodes))

# We create an element set containing all the elements
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))

# We create a material section of a linear isotropic material
# and assigns that section to all the elements.
mat_section = MaterialSection(FEM.linearisotropicmod().LinearIsotropic(1, 0.3))
push!(mat_section, geomesh.element_sets["all"])

# We create an element section and assign it to all the elements.
ele_section = ElementSection(FEM.quadtrigmod().QuadTrig)
push!(ele_section, geomesh.element_sets["all"])

# Apply Dirichlet BC to the left side in dofs for x-displacement (Du)
# and dofs in y-displacement (Dv)
bcs = Any[DirichletBC("0.0", [FEM.Du], geomesh.node_sets["left"]),
          DirichletBC("0.0", [FEM.Dv], geomesh.node_sets["left"])]

# Since we currently don't have edge load we give a nodal load
# on the node set on the right edge.
loads = Any[NodeLoad("1/($n_ele_y+1)", [FEM.Dv], geomesh.node_sets["right"])]

# Create the fe problem
fp = FEM.create_feproblem("cook_example_quad", geomesh, [ele_section], [mat_section], bcs, loads)


vtkexp = VTKExporter()
# Output fields are added by pushing them into the exporter.
# We want to export the stress and strain tensor as well as
# the von Mises stress so we push them.
push!(vtkexp, Stress)
push!(vtkexp, Strain)
push!(vtkexp, VonMises)
# We set the output of the exporter to be in ascii
set_binary!(vtkexp, false)

solver = NRSolver(rel_tol = 1e-7, max_iters = 4)

# Solve the fe problem using the solver and using the vtk exporter.
solve(solver, fp, vtkexp)
