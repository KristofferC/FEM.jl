using FEM
n_ele = 2

m = [[0.0; 0.0] [0.0; 10.0] [10.0; 10.0] [10.0; 0.0]]
geomesh = meshquad(n_ele, n_ele, m, GeoQTrig)

using FEM
import FEM.write_data
# Nodes


push!(geomesh, gennodeset(n->n.coords[2]>9.9999, "top", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords[2]<0.0001, "bottom", geomesh.nodes))
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))
# Material section

E = 200000.e0
nu = 0.3e0
n = 1.e0
l = 1.e-2
kinf = 1e+010
lambda_0 = 4.e-2
Hg = 4.e7
Hl = 10_000.0
m = 2.0
faktor = 0.5
sy = 1000.0
tstar = 1000.0
c_dam = 1.0
angles = [45.0, 105.0]
nslip = 2

mat = GradMekh(E, nu, n, l, kinf, lambda_0,
               Hg, Hl, m, faktor, sy, tstar,
               c_dam, angles, nslip)

mat_section = MaterialSection(mat)
push!(mat_section, geomesh.element_sets["all"])

# Element section
ele_section = ElementSection(GradTrig)
push!(ele_section, geomesh.element_sets["all"])

# Boundary conditions
bcs = [DirichletBC(0.0, [FEM.Du, FEM.Dv], geomesh.node_sets["bottom"]),
       DirichletBC(0.0, [FEM.Dv], geomesh.node_sets["top"]),
        DirichletBC(0.125, [FEM.Du], geomesh.node_sets["top"])]


fp = FEM.create_feproblem_grad("grad", geomesh, [ele_section], [mat_section], bcs)

vtkexp = VTKExporter()
set_binary!(vtkexp, false)
push!(vtkexp, Stress)
push!(vtkexp, Strain)
push!(vtkexp, VonMises)
push!(vtkexp, InvFp)
push!(vtkexp, KappaVector)


solver = NRSolver(1e-4, 20)

solve(solver, fp, vtkexp)
#=
K = FEM.assembleK(fp)

#du = cholfact(Symmetric(K, :L)) \ force_imbalance

du = -K \ int_f

solve(solver, fp, vtkexp)

FEM.updatedofs!(fp, du)
FEM.update_feproblem(fp)


# Output fields are added by pushing them into the exporter

set_binary!(vtkexp, false)
write_data(fp, vtkexp)
=#