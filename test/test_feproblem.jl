import FEM.createdofs
import FEM.extload
import FEM.assembleK

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

bcs =  [DirichletBC(0.1, [Du, Dv], mesh.node_sets["x0"])]
loads =  [PointLoad(10e5, [Dv], mesh.node_sets["y0"])]

fp = FEProblem(mesh, bcs, loads)

context("FEM.FEProblem.") do
    createdofs(fp) #TODO: Test this properly

    f = extload(fp)

    K = assembleK(fp)

    println(K \ f)








end # context

end # facts



