using FEM
n_ele = 60

m = [[0.0; 0.0] [0.0; 1.0] [1.0; 1.0] [1.0; 0.0]]
geomesh = meshquad(n_ele, n_ele, m, GeoQTrig)

using FEM
import FEM.write_data
# Nodes


push!(geomesh, gennodeset(n->n.coords[2]>0.9999, "top", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords[2]<0.0001, "bottom", geomesh.nodes))
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))
# Material section

mat_section = MaterialSection(LinearIsotropic(200000, 0.3))

push!(mat_section, geomesh.element_sets["all"])

# Element section
ele_section = ElementSection(QuadTrig)
push!(ele_section, geomesh.element_sets["all"])

# Boundary conditions
bcs = [DirichletBC(0.0, [FEM.Du, FEM.Dv], geomesh.node_sets["bottom"]),
       DirichletBC(0.0, [FEM.Dv], geomesh.node_sets["top"]),
        DirichletBC(0.025, [FEM.Du], geomesh.node_sets["top"])]


fp = create_feproblem("grad_iso", geomesh, [ele_section], [mat_section], bcs)


elem = fp.sections[1].elements[1]
matstat = elem.matstats[1]
temp_matstat = elem.temp_matstats[1]

int_f = FEM.assemble_intf(fp)

#=
for node in fp.nodes
    for dof in node.dofs
        if dof.active
            println("n: $(node.n) type: $(dof.dof_type) val: $(int_f[dof.eq_n])")
        end
    end
end
=#

vtkexp = VTKExporter()
set_binary!(vtkexp, false)
push!(vtkexp, Stress)
push!(vtkexp, Strain)
push!(vtkexp, VonMises)
solver = NRSolver(1e-1, 25)

solve(solver, fp, vtkexp)


#=
K = FEM.assembleK(fp)

#du = cholfact(Symmetric(K, :L)) \ force_imbalance

du = -K \ int_f

solve(solver, fp, vtkexp)




# Output fields are added by pushing them into the exporter

set_binary!(vtkexp, false)

=#