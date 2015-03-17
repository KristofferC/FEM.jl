import FEM.doftypes

type FEProblem
    mesh::FEMesh
    bcs::Vector{DirichletBC}
    loads::Vector{NodeLoad}
    sections::Vector{Section}
    node_doftypes::Dict{Int, Vector{DofType}}
    node_doftype_bc::Dict{(Int, DofType), DirichletBC}
    n_eqs::Int
    n_fixed::Int
end

function FEProblem(mesh::FEMesh, bcs::Vector{DirichletBC}=Array(DirichletBC, 0),
                    loads::Vector{NodeLoad}=Array(NodeLoad, 0), sections=Array(Section, 0))
    node_doftype_bc = Dict{Int, Vector{DofType}}()
    node_doftypes = Dict{Int, Vector{DofType}}()
    FEProblem(mesh, bcs, loads, sections, node_doftypes, node_doftype_bc, 0, 0)
end

push!(fp::FEProblem, bc::DirichletBC) = push!(fp.bcs, bc)
push!(fp::FEProblem, load::NodeLoad) = push!(fp.loads, load)
push!(fp::FEProblem, sec::Section) = push!(fp.sections, section)


function createdofs(fp::FEProblem)

    # Create a dictionary between a node_id to
    # what dof types are in that node.
    eq_n = 1
    for (el_id, element) in fp.mesh.elements
        for vertex in element.vertices
            dof_types = doftypes(element, vertex)
            fp.node_doftypes[vertex] = dof_types
        end
    end

    # Create a dictionary between a tuple of node_id
    # and dof type to the BC for that tuple.

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
    id = 0
    dofs = Array(Dof, 0)
    for node in fp.mesh.nodes
        for doftype in fp.node_doftypes[node.n]
            if haskey(fp.node_doftype_bc, (node.n, doftype))
                bc = fp.node_doftype_bc[(node.n, doftype)]
                pres_n += 1
                id += 1
                push!(node.dofs, Dof(pres_n, id, false, bc.value, doftype))
            else
                eq_n += 1
                id += 1
                push!(node.dofs, Dof(eq_n, id, true, 0.0, doftype))
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
                if dof.dof_type in load.dof_types
                    f[dof.eq_n] += load.value
                end
            end
        end
    end
    return f
end


function assembleK(fp::FEProblem)
    dof_rows = Array(Int, 0)
    dof_cols = Array(Int, 0)
    k_values = Array(Float64, 0)
    for section in fp.sections
        mat = section.material
        for element_id in section.elements
            element = fp.mesh.elements[element_id]
            Ke = stiffness(element, fp.mesh.nodes, mat)
            dof1_n = 0
            for vertex1 in element.vertices
                for dof1 in fp.mesh.nodes[vertex1].dofs
                    dof1_n += 1
                    dof2_n = 0
                    for vertex2 in element.vertices
                        for dof2 in fp.mesh.nodes[vertex2].dofs
                            dof2_n += 1
                            if dof1.active && dof2.active
                                push!(dof_rows, dof1.eq_n)
                                push!(dof_cols, dof2.eq_n)
                                push!(k_values, Ke[dof1_n, dof2_n])
                            end
                        end
                    end
                end
            end
        end
    end
    K = Base.sparse(dof_rows, dof_cols, k_values, fp.n_eqs, fp.n_eqs)
    return K
end


function assemble_intf(fp::FEProblem)
    int_forces = intf(fp)
    int_forces_assem = zeros(fp.n_eqs)
    for node in fp.mesh.nodes
        for dof in node.dofs
            if dof.active
                int_forces_assem[dof.eq_n] += int_forces[dof.id]
            end
        end
    end
    return int_forces_assem
end

function intf(fp::FEProblem)
    fint = zeros(fp.n_eqs + fp.n_fixed)
    for section in fp.sections
        mat = section.material
        for element_id in section.elements
            element = fp.mesh.elements[element_id]
            finte = intf(element, fp.mesh.nodes, mat)
            i = 1
            for vertex in element.vertices
                for dof in fp.mesh.nodes[vertex].dofs
                    fint[dof.id] += finte[i]
                    i+=1
                end
            end
        end
    end
    return fint
end

function updatedofs!(fp::FEProblem, du::Vector{Float64})
    for node in fp.mesh.nodes
        for dof in node.dofs
            if dof.active
                dof.value += du[dof.eq_n]
            end
        end
    end
end



function assemble_dry(fp::FEProblem)
    dof_rows = Array(Int, 0)
    dof_cols = Array(Int, 0)
    k_values = Array(Float64, 0)
    for section in fp.sections
        mat = section.material
        for element_id in section.elements
            element = fp.mesh.elements[element_id]
            dof1_n = 0
            for vertex1 in element.vertices
                for dof1 in fp.mesh.nodes[vertex1].dofs
                    dof1_n += 1
                    dof2_n = 0
                    for vertex2 in element.vertices
                        for dof2 in fp.mesh.nodes[vertex2].dofs
                            dof2_n += 1
                            if dof1.active && dof2.active
                                push!(dof_rows, dof1.eq_n)
                                push!(dof_cols, dof2.eq_n)
                                push!(k_values, 0)
                            end
                        end
                    end
                end
            end
        end
    end
end
