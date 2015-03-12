# 4.4 sec, 1GB alloc

using FEM
import FEM.extload
import FEM.assembleK
import FEM.createdofs
# Nodes
n_ele = 128
mesh = gencook(n_ele, n_ele)

addnodeset!(mesh, gennodeset(x->x[1]>0.0479999, "right", mesh.nodes))
addnodeset!(mesh, gennodeset(x->x[1]<0.000001, "left", mesh.nodes))

mat = LinearIsotropic(250e9, 0.3)
section = Section(mat)


addelemset!(mesh, ElementSet("all", collect(1:2*n_ele*n_ele)))


bcs =  [DirichletBC(0.0, [Du(), Dv()], mesh.node_sets["left"])]
loads =  [NodeLoad(1/16/n_ele, [Dv()], mesh.node_sets["right"])]

# Element set
mat = LinearIsotropic(1e3, 0.3)
section = Section(mat)
addelemset!(section, mesh.element_sets["all"])

fp = FEProblem(mesh, bcs, loads, [section])

solver = NRSolver(1e-7, 5)

#createdofs(fp)
#K = assembleK(fp)

#println(K)

solve(solver, fp)


