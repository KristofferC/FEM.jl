using FEM
import FEM.write_data
# Nodes

# Generate geomesh and node / elementsets

n_ele = 2

geomesh = gencook(n_ele, n_ele, GeoTrig)

push!(geomesh, gennodeset(n->n.coords[1]>47.999999, "right", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords[1]<0.00001, "left", geomesh.nodes))
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))
# Material section
mat_section = MaterialSection(LinearIsotropic(1, 0.3))
push!(mat_section, geomesh.element_sets["all"])

# Element section
ele_section = ElementSection(LinTrig)
push!(ele_section, geomesh.element_sets["all"])

# Boundary conditions
bcs = [DirichletBC(0.0, [FEM.Du, FEM.Dv], geomesh.node_sets["left"])]

# Loads
loads = [NodeLoad(1/(n_ele+1), [FEM.Dv], geomesh.node_sets["right"])]

fp = create_feproblem("cook_example_quad", geomesh, [ele_section], [mat_section], bcs, loads)

vtkexp = VTKExporter()
# Output fields are added by pushing them into the exporter
push!(vtkexp, Stress)
push!(vtkexp, Strain)
set_compress!(vtkexp, false)
#push!(vtkexp, Displacement)

solver = NRSolver(1e-9, 2)

solve(solver, fp, vtkexp)

#write_data(fp, vtkexp)
#write_data(fp, vtkexp)

#write_VTKXML("test_bin.vtu", fp, false, false)


#exportVTK(fp, "cook.vtk", false)
