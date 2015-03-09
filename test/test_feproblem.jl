import FEM.createdofs
import FEM.extload
import FEM.assembleK
import FEM.intf
import FEM.assemble_intf
import FEM.updatedofs!

facts("FEM.FEProblem") do

# Setup a FEProblem
mesh = Mesh()

nodes = [Node([0, 0], 1), Node([1, 1], 2), Node([1, 2], 3), Node([0, 1], 4)]
addnodes!(mesh, nodes)

addelem!(mesh, LinTrig([1, 2, 3], 1))
addelem!(mesh, LinTrig([1, 2, 4], 2))

bottom_set = NodeSet("y0", [1])
top_set = NodeSet("x0", [2, 3])
addnodeset!(mesh, bottom_set)
addnodeset!(mesh, top_set)

element_set = ElementSet("all", [1, 2])

addelemset!(mesh, element_set)

bcs =  [DirichletBC(0.1, [Du(), Dv()], mesh.node_sets["x0"])]
loads =  [PointLoad(10e5, [Dv()], mesh.node_sets["y0"])]


# Element set
mat = LinearIsotropic(250e9, 0.3)
section = Section(mat)
addelemset!(section, mesh.element_sets["all"])


fp = FEProblem(mesh, bcs, loads, [section])

context("FEM.FEProblem.") do
    createdofs(fp) #TODO: Test this properly

    load = extload(fp)



    int_f = assemble_intf(fp)
    println(load)
    println(int_f)
    K = assembleK(fp)
    du = K \ (load - int_f)
    println(du1)

    updatedofs!(fp, du)

    load = extload(fp)
    int_f = assemble_intf(fp)

    du2 = K \ (load - int_f)

    u = du + du2
    norm()

    #@fact norm(u - [-2.21276596e-06, 1.30553191e-05, 1.54893617e-06, 1.15063830e-05]) => roughly(0.0, 10.0^(-13))




end # context

end # facts



