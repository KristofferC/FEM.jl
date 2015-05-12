import FEM.extload
import FEM.assembleK!
import FEM.intf
import FEM.assemble_intf
import FEM.updatedofs!

facts("FEM.FEProblem") do

# Setup a FEProblem
mesh = GeoMesh()

nodes = [GeoNode2(1, [0, 0]), GeoNode2(2, [1, 1]), GeoNode2(3, [1, 2]), GeoNode2(4, [0, 1])]


push!(mesh, nodes)

push!(mesh, GeoTrig(1, [1, 2, 3]))
push!(mesh, GeoTrig(2, [1, 2, 4]))


bottom_set = NodeSet("y0", [1])
top_set = NodeSet("x0", [2, 3])
push!(mesh, bottom_set)
push!(mesh, top_set)

element_set = ElementSet("all", [1, 2])
push!(mesh, element_set)

bcs =  Any[DirichletBC("0.1", [FEM.Du, FEM.Dv], mesh.node_sets["x0"])]
loads =  Any[NodeLoad("10e5", [FEM.Dv], mesh.node_sets["y0"])]




# Element set
mat_section = MaterialSection(LinearIsotropic(250e9, 0.3))
push!(mat_section, mesh.element_sets["all"])

ele_section = ElementSection(LinTrig)
push!(ele_section, mesh.element_sets["all"])


fp = create_feproblem("cook_example_quad", mesh, [ele_section], [mat_section], bcs, loads)

context("FEM.FEProblem") do
    K = FEM.create_sparse_structure(fp::FEProblem)
    colptrs = FEM.get_colptrs(K, fp::FEProblem)

    load = extload(fp)

    int_f = assemble_intf(fp)
    assembleK!(K, fp, colptrs)
    du = K \ (load - int_f)

    updatedofs!(fp, du)

    # Verify force balance
    load = extload(fp)
    int_f = assemble_intf(fp)
    @fact norm(load - int_f) / norm(load) => roughly(0.0, 10.0^(-10))


end # context

end # facts



