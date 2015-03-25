using FEM


# Nodes

# Generate geomesh and node / elementsets
n_ele = 100
geomesh = gencook(n_ele, n_ele)
push!(geomesh, gennodeset(n->n.coords[1]>=48.0, "right", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords[1]<=0.0, "left", geomesh.nodes))
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))



# Material section
mat_section = MaterialSection(LinearIsotropic(1, 0.3))
push!(mat_section, geomesh.element_sets["all"])

# Element section
ele = LinTrig
ele_section = ElementSection(ele)
push!(ele_section, geomesh.element_sets["all"])

# Boundarý conditions
bcs = [DirichletBC(0.0, [Du, Dv], geomesh.node_sets["left"])]

# Loads
loads = [NodeLoad(0.0078, [Dv], geomesh.node_sets["right"])]

fp = create_feproblem(geomesh, [ele_section], [mat_section], bcs, loads)

solver = NRSolver(1e-7, 2)

solve(solver, fp)

#exportVTK(fp, "test.vtk")

#write_vtk_file(fp.FEMesh, "cook.jl", false)

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

