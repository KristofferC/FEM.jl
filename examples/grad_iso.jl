using FEM
n_ele = 2

m = [[0.0; 0.0] [0.0; 1.0] [1.0; 1.0] [1.0; 0.0]]
geomesh = meshquad(n_ele, n_ele, m, GeoTrig)



push!(geomesh, gennodeset(n->n.coords[2]>0.9999 , "top", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords[2]<0.0001, "bottom", geomesh.nodes))
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))
# Material section

mat_section = MaterialSection(LinearIsotropic(200000, 0.3))

push!(mat_section, geomesh.element_sets["all"])

# Element section
ele_section = ElementSection(LinTrig)
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

for elem in fp.sections[1].elements
    println(" ")
    for (dof, j) in activedofs(elem, fp.nodes)
        println("$j, $dof")
    end
end

println(" ")
for element in fp.sections[1].elements
    println(" ")
    i = 1
    for vertex in element.vertices
        for dof in fp.nodes[vertex].dofs
            if dof.active
                println("$i, $dof")
            end
            i+=1
        end
    end
end


#=
K = FEM.assembleK(fp)

#du = cholfact(Symmetric(K, :L)) \ force_imbalance

du = -K \ int_f

solve(solver, fp, vtkexp)




# Output fields are added by pushing them into the exporter

set_binary!(vtkexp, false)

=#
