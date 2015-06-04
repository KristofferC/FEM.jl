using FEM

FEM.vtkexportmod()
using FEM.VTKExportMod

n_ele = 5

m = [[0.0; 0.0] [0.0; 10.0] [10.0; 10.0] [10.0; 0.0]]
geomesh = meshquad(n_ele, n_ele, m, GeoQTrig)


# Nodes


boundary_set = gennodeset(n->n.coords[2]>9.9999, "boundary", geomesh.nodes)
#boundary_set2 = gennodeset(n->n.coords[2]<0.0001, "boundary2", geomesh.nodes)
append!(boundary_set, gennodeset(n->n.coords[1]>9.9999, "boundary", geomesh.nodes))
append!(boundary_set, gennodeset(n->n.coords[2]<0.0001, "boundary", geomesh.nodes))
append!(boundary_set, gennodeset(n->n.coords[1]<0.0001, "boundary", geomesh.nodes))

push!(geomesh, boundary_set)
#push!(geomesh, boundary_set2)
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))
# Material section


E = 200000.e0
nu = 0.3e0
l = 1.e-2
Hg = 4.e7
Hl = 10_000.01
m = 2.0
sy = 1000.0
tstar = 1000.0
angles = [45.0, 105.0]
nslip = 2

 mat = FEM.gradmekhprimalmod().GradMekh(E, nu, l, Hg, Hl, m, sy, tstar, angles, nslip)

#=
faktor = 1
kinf = 1e+010
lambda_0 = 4.e-2
c_dam = 0.0
mat = FEM.gradmekhmod().GradMekh(E, nu, n, l, kinf, lambda_0,Hg, Hl, m, faktor, sy, tstar,c_dam, angles, nslip)
=#

mat_section = MaterialSection(mat)
push!(mat_section, geomesh.element_sets["all"])

# Element section
ele_section = ElementSection(FEM.gradtrigprimalmod().GradTrig)
push!(ele_section, geomesh.element_sets["all"])

# Boundary conditions
γ = 0.025
bcs = [DirichletBC("$(γ)*y*t", [FEM.Du], geomesh.node_sets["boundary"]),
       DirichletBC("0.0", [FEM.Dv], geomesh.node_sets["boundary"])]


fp = FEM.create_feproblem_grad("grad_jll", geomesh, [ele_section], [mat_section], bcs)



vtkexp = VTKExporter()
set_binary!(vtkexp, false)
push!(vtkexp, Stress)
push!(vtkexp, Strain)
push!(vtkexp, VonMises)
#push!(vtkexp, InvFp)
#push!(vtkexp, KAlpha)


solver = NRSolver(abs_tol = 1e-2, max_iters = 10)

solve(solver, fp, vtkexp)
