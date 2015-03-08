import FEM.doftypes

type FEProblem
    mesh::Mesh
    bcs::Vector{DirichletBC}
    loads::Vector{PointLoad}
    dofs::Vector{Dof}
    node_doftypes::Dict{Int, Vector{DofType}}
    node_doftype_bc::Dict{(Int, DofType), DirichletBC}
    n_eqs::Int
    n_fixed::Int
end

function FEProblem(mesh::Mesh, bcs::Vector{DirichletBC},
                    loads::Vector{PointLoad})
    node_doftype_bc = Dict{Int, Vector{DofType}}()
    node_doftypes = Dict{Int, Vector{DofType}}()
    FEProblem(mesh, bcs, loads, Array(Dof, 0), node_doftypes, node_doftype_bc, 0, 0)
end


function createdofs(fp::FEProblem)

    # Create a dictionary between a node_id to
    # what dof types are in that node.
    eq_n = 1
    for element in fp.mesh.elements
        for vertex in element.vertices
            dof_types = doftypes(element, vertex)
            fp.node_doftypes[vertex] = dof_types
        end
    end

    # Create a dictionary between a tuple of node_id
    # and dof type to the BC for that tuple.
    dof_id = 1
    for bc in fp.bcs
        node_ids = bc.node_set.node_ids
        for node_id in node_ids
            node_number = fp.mesh.nodes[node_id].n
            for doftype in bc.dof_types
                fp.node_doftype_bc[(node_id, doftype)] = bc
            end
        end
    end

    # Create the dofs
    eq_n = 0
    pres_n = 0
    dofs = Array(Dof, 0)
    for node in fp.mesh.nodes
        for doftype in fp.node_doftypes[node.n]
            if haskey(fp.node_doftype_bc, (node.n, doftype))
                bc = fp.node_doftype_bc[(node.n, doftype)]
                pres_n += 1
                dof_id += 1
                push!(node.dofs, Dof(pres_n, dof_id, false, bc.value, doftype))
            else
                eq_n += 1
                dof_id += 1
                push!(node.dofs, Dof(eq_n, dof_id, true, 0.0, doftype))
            end
        end
    end
    #TODO: Make these have same name
    fp.n_eqs = eq_n
    fp.n_fixed = pres_n
end


function extload(fp::FEProblem)
    f = zeros(fp.n_eqs)
    for load in fp.loads
        for node_id in load.node_set.node_ids
            for dof in fp.mesh.nodes[node_id].dofs
                f[dof.dof_id] += load.value
            end
        end
    end
    return f
end

function assembleK(fp::FEProblem)
    dof_rows = Array(Int, 0)
    dof_cols = Array(Int, 0)
    k_values = Array(Float64, 0)
    mat = LinearIsotropic(200e9, 0.3)
    for element in fp.mesh.elements
        Ke = stiffness(element, fp.mesh.nodes, mat)
        for vertex1 in element.vertices
            for (i, dof1) in enumerate(fp.mesh.nodes[vertex1].dofs)
                for vertex2 in element.vertices
                    for (j, dof2) in enumerate(fp.mesh.nodes[vertex2].dofs)
                        if dof1.active && dof2.active
                            push!(dof_rows, dof1.eq_n)
                            push!(dof_cols, dof2.eq_n)
                            push!(k_values, Ke[i, j])
                        end
                    end
                end
            end
        end
    end
    return Base.sparse(I, J, V, fp.n_eqs, fp.n_eqs)
end
