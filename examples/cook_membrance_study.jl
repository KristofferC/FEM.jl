using FEM
import FEM.extload
import FEM.assembleK
import FEM.createdofs
# Nodes
n_ele = 128
mesh = gencook(n_ele, n_ele)

addnodeset!(mesh, gennodeset(x->x[1]>47.9999, "right", mesh.nodes))
addnodeset!(mesh, gennodeset(x->x[1]<0.000001, "left", mesh.nodes))



addelemset!(mesh, ElementSet("all", collect(1:2*n_ele*n_ele)))


bcs =  [DirichletBC(0.0, [Du, Dv], mesh.node_sets["left"])]
loads =  [NodeLoad(0.0078, [Dv], mesh.node_sets["right"])]

# Element set
mat = LinearIsotropic(1, 0.3)
section = Section(mat)
addelemset!(section, mesh.element_sets["all"])

fp = FEProblem(mesh, bcs, loads, [section])

solver = NRSolver(1e-7, 2)



solve(solver, fp)

#=
julia> @time include(".julia/v0.4/FEM/examples/cook_membrance_study.jl")
Starting Newton-Raphson solver..
        Iteration 1, relative residual 1.0
4.718447854656915e-15
[22.220644996854695]
        Iteration 2, relative residual 2.6125689174257096e-11
Converged!
elapsed time: 0.518623941 seconds (308 MB allocated, 8.98% gc time in 14 pauses with 1 full sweep)
=#


#=
Starting Newton-Raphson solver..
        Iteration 1, relative residual 1.0
        Iteration 2, relative residual 2.0046022544402292e-11
Converged!
elapsed time: 1.185315694 seconds (467 MB allocated, 5.20% gc time in 21 pauses with 1 full sweep)


=#

#=
julia> @time include(".julia/v0.4/FEM/examples/cook_membrance_study.jl")
Starting Newton-Raphson solver..
        Iteration 1, relative residual 1.0
        Iteration 2, relative residual 2.091766113241007e-11
Converged!
elapsed time: 3.132825849 seconds (1192 MB allocated, 4.90% gc time in 54 pauses with 1 full sweep)
=#


#createdofs(fp)
#K = assembleK(fp)

#println(K)

