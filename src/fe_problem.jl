type FEProblem
    nodes::Vector{FENode2}
    bcs::Vector{DirichletBC}
    loads::Vector{NodeLoad}
    sections::Vector{FESection}
    node_doftypes::Dict{Int, Vector{DofType}}
    node_doftype_bc::Dict{(Int, DofType), DirichletBC}
    n_eqs::Int
    n_fixed::Int
end

function FEProblem(nodes::Vector{FENode2}, bcs::Vector{DirichletBC}=Array(DirichletBC, 0),
                    loads::Vector{NodeLoad}=Array(NodeLoad, 0), sections=Array(FESection, 0))
    node_doftype_bc = Dict{Int, Vector{DofType}}()
    node_doftypes = Dict{Int, Vector{DofType}}()
    FEProblem(nodes, bcs, loads, sections, node_doftypes, node_doftype_bc, 0, 0)
end

push!(fp::FEProblem, bc::DirichletBC) = push!(fp.bcs, bc)
push!(fp::FEProblem, load::NodeLoad) = push!(fp.loads, load)
push!(fp::FEProblem, section::FESection) = push!(fp.sections, section)


function create_feproblem(geomesh, element_regions, material_regions, bcs, loads)

    gps = Dict{DataType, Vector{GaussPoint2}} ()
    interps = Dict{DataType, AbstractInterpolator} ()
    storage = Dict{DataType, ElemStorage} ()
    elem_types = Array(DataType, 0)

    for element_region in element_regions
        elem_type = element_region.element_type
        interps[elem_type] = get_interp(elem_type)
        gps[elem_type] = get_gps(elem_type)
        storage[elem_type] = get_storage(elem_type)
    end

    nodes = Array(FENode2, 0)
    for node in geomesh.nodes
        push!(nodes, FENode2(node.n, node.coords))
    end

    sections = Array(FESection, 0)
    for matregion in material_regions
        material = matregion.material
        for eleregion in element_regions
            ele_type = eleregion.element_type
            common = intersect(matregion.elements, eleregion.elements)
            gps_ele = gps[ele_type]
            elem_storage = storage[ele_type]
            interp = interps[ele_type]

            section = FESection(material, ele_type)

            for ele_id in common
                vertices = geomesh.elements[ele_id].vertices
               # mat_stat = create_matstat(typeof(material))

                element = ele_type(vertices, gps_ele, ele_id, interp,
                                   elem_storage)
                push!(section, element)
            end
            push!(sections, section)
        end
    end
    fe = FEProblem(nodes, bcs, loads, sections)
end

function createdofs(fp::FEProblem)

    # Create a dictionary between a node_id to
    # what dof types are in that node.
    eq_n = 1
    for section in fp.sections
        for element in values(section.elements)
            for vertex in element.vertices
                dof_types = doftypes(element, vertex)
                fp.node_doftypes[vertex] = dof_types
            end
        end
    end

    # Create a dictionary between a tuple of node_id
    # and dof type to the BC for that tuple.

    for bc in fp.bcs
        node_ids = bc.node_set.node_ids
        for node_id in node_ids
            node_number = fp.nodes[node_id].n
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
    for node in fp.nodes
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
            for dof in fp.nodes[node_id].dofs
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
        assemble_K_section(section, fp.nodes, dof_rows, dof_cols, k_values)
    end
    K = Base.sparse(dof_rows, dof_cols, k_values, fp.n_eqs, fp.n_eqs)
    return K
end

function assemble_K_section{T<:AbstractFElement, P <: AbstractMaterial}(section::FESection{T,P}, nodes::Vector{FENode2},
                                    dof_rows::Vector{Int}, dof_cols::Vector{Int}, k_values::Vector{Float64})
    mat = section.material
    for element in values(section.elements)
        Ke = stiffness(element, nodes, mat)
        dof1_n = 0
        for vertex1 in element.vertices
            for dof1 in nodes[vertex1].dofs
                dof1_n += 1
                dof2_n = 0
                for vertex2 in element.vertices
                    for dof2 in nodes[vertex2].dofs
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

function assemble_intf(fp::FEProblem)
    fint = zeros(fp.n_eqs)
    for section in fp.sections
        assemble_intf_section(section, fint, fp.nodes)
    end
    return fint
end

function assemble_intf_section{T<:AbstractFElement, P <: AbstractMaterial}(section::FESection{T, P},
                                                                           int_forces::Vector{Float64},
                                                                           nodes::Vector{FENode2})
    mat = section.material
    for element in values(section.elements)
        finte = intf(element, mat, nodes)
        i = 1
        for vertex in element.vertices
            for dof in nodes[vertex].dofs
                if dof.active
                    int_forces[dof.eq_n] += finte[i]
                end
                i+=1
            end
        end
    end
end

#=
function assemble_intf(fp::FEProblem)
    int_forces = intf(fp)
    int_forces_assem = zeros(fp.n_eqs)
    for node in fp.nodes
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
        assemble
        mat = section.material
        for element in values(section.elements)
            finte = intf(element, mat, fp.nodes)
            i = 1
            for vertex in element.vertices
                for dof in fp.nodes[vertex].dofs
                    fint[dof.id] += finte[i]
                    i+=1
                end
            end
        end
    end
    return fint
end
=#

function updatedofs!(fp::FEProblem, du::Vector{Float64})
    for node in fp.nodes
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
