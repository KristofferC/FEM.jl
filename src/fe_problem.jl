type FEProblem
    name::ASCIIString
    nodes::Vector{FENode2}
    bcs::Vector{DirichletBC}
    loads::Vector{NodeLoad}
    sections::Vector{FESection}
    node_doftypes::Dict{Int, Vector{DofType}}
    node_doftype_bc::Dict{(Int, DofType), DirichletBC}
    n_eqs::Int
    n_fixed::Int
end

function FEProblem(name::ASCIIString, nodes::Vector{FENode2}, bcs::Vector{DirichletBC},
                    loads, sections=Array(FESection, 0))
    node_doftype_bc = Dict{Int, Vector{DofType}}()
    node_doftypes = Dict{Int, Vector{DofType}}()
    FEProblem(name, nodes, bcs, loads, sections, node_doftypes, node_doftype_bc, 0, 0)
end

push!(fp::FEProblem, bc::DirichletBC) = push!(fp.bcs, bc)
push!(fp::FEProblem, load::NodeLoad) = push!(fp.loads, load)
push!(fp::FEProblem, section::FESection) = push!(fp.sections, section)


function create_feproblem(name, geomesh, element_regions, material_regions, bcs::Vector{DirichletBC}=Array(DirichletBC, 0), loads::Vector{NodeLoad}=Array(NodeLoad, 0))

    gps = Dict{DataType, Vector{GaussPoint2}} ()
    interps = Dict{DataType, AbstractInterpolator} ()
    storage = Dict{DataType, ElemStorage} ()
    elem_types = Array(DataType, 0)

    for element_region in element_regions
        elem_type = element_region.element_type
        interps[elem_type] = createinterp(elem_type)
        gps[elem_type] = creategps(elem_type)
        storage[elem_type] = createstorage(elem_type)
    end

    nodes = Array(FENode2, 0)
    for node in geomesh.nodes
        push!(nodes, FENode2(node.n, node.coords))
    end

    sections = Array(FESection, 0)
    for matregion in material_regions
        material = matregion.material
        matstat = create_matstat(typeof(material))
        for eleregion in element_regions
            ele_type = eleregion.element_type
            common = intersect(matregion.elements, eleregion.elements)
            common = collect(common)
            sort!(common) # TODO: Remove?
            gps_ele = gps[ele_type]
            elem_storage = storage[ele_type]
            interp = interps[ele_type]

            section = FESection(material, ele_type, typeof(matstat))

            for ele_id in common
                vertices = geomesh.elements[ele_id].vertices
                element = ele_type(vertices, ele_id, interp,
                                   elem_storage, gps_ele, matstat)
                push!(section, element)

            end
            push!(sections, section)
        end
    end
    fe = FEProblem(name, nodes, bcs, loads, sections)
    createdofs(fe)
    return fe
end


# Type stability sanitized
function set_dof_types_section!{T<:FESection}(section::T,
                                       node_doftypes::Dict{Int, Vector{DofType}})
    for element in section.elements
        for (i, v) in enumerate(element.vertices)
            dof_types = doftypes(element, i)
            node_doftypes[v] = dof_types
        end
    end
end


function createdofs(fp::FEProblem)

    # Create a dictionary between a node_id to
    # what dof types are in that node.
    eq_n = 1
    for section in fp.sections
        set_dof_types_section!(section, fp.node_doftypes)
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

    # TODO: Optimize
    for node in fp.nodes
        for doftype in fp.node_doftypes[node.n]
            if haskey(fp.node_doftype_bc, (node.n, doftype))
                bc = fp.node_doftype_bc[(node.n, doftype)]
                pres_n += 1
                id += 1
                push!(node.dofs, Dof(pres_n, id, false, 0.0, doftype))
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


function extload(fp::FEProblem, t::Number=0.0)
    f = zeros(fp.n_eqs)
    for load in fp.loads
        for node_id in load.node_set.node_ids
            for dof in fp.nodes[node_id].dofs
                if dof.dof_type in load.dof_types
                     extload(f, load, dof, fp.nodes[node_id], t)
                end
            end
        end
    end
    return f
end

@inline function extload{T}(f::Vector{Float64}, load::NodeLoad{T},
                            dof::Dof, node::FENode2, t::Number)
    f[dof.eq_n] += evaluate(load, node, t)
end



function assembleK!(K::SparseMatrixCSC, fp::FEProblem, colptrs::Vector{Int})
    fill!(K.nzval, 0.0)
    z = 1
    for section in fp.sections
        z = assembleK!(K, colptrs, z, section, fp.nodes)
    end
    return K
end


function assembleK!(K::SparseMatrixCSC, colptrs::Vector{Int}, z::Int,
                    section::FESection, nodes::Vector{FENode2})
    mat = section.material
    for element in section.elements
        Ke = stiffness(element, nodes, mat)
        dof1_n = 0
        for vertex1 in element.vertices
            for dof1 in nodes[vertex1].dofs
                dof1_n += 1
                if dof1.active
                    dof2_n = 0
                    for vertex2 in element.vertices
                        for dof2 in nodes[vertex2].dofs
                            dof2_n += 1
                            if dof2.active
                                K.nzval[colptrs[z]] += Ke[dof1_n, dof2_n]
                                z+=1
                            end
                        end
                    end
                end
            end
        end
    end
    return z
end

# This is nicer but the activedofs iterators are too slow right now
#=
function assembleK!(K::SparseMatrixCSC, colptrs::Vector{Int}, z::Int,
                    section::FESection, nodes::Vector{FENode2})
    mat = section.material
    for element in section.elements
        Ke = stiffness(element, nodes, mat)
        for (dof1, i) in activedofs(element, nodes)
            for (dof2, j) in activedofs(element, nodes)
                v = K[dof1.eq_n, dof2.eq_n]
                K.nzval[colptrs[z]] += Ke[i, j]
                z += 1
            end
        end
    end
    return z
end
=#




function create_sparse_structure(fp::FEProblem)
    dof_rows = Int[]
    dof_cols = Int[]
    for section in fp.sections
        create_sparse_structure(section, fp.nodes, dof_rows, dof_cols)
    end
    # Using ones until we get a sparse structure initializer in base
    return Base.sparse(dof_rows, dof_cols, ones(length(dof_rows)), fp.n_eqs, fp.n_eqs)
end

function create_sparse_structure(section::FESection, nodes::Vector{FENode2},
                                dof_rows::Vector{Int}, dof_cols::Vector{Int})
    for element in section.elements
        for (dof1, _) in activedofs(element, nodes)
            for (dof2, _) in activedofs(element, nodes)
                push!(dof_rows, dof1.eq_n)
                push!(dof_cols, dof2.eq_n)
            end
        end
    end
end

function assembleK(fp::FEProblem)
    dof_rows = Array(Int, 0)
    dof_cols = Array(Int, 0)
    k_values = Array(Float64, 0)
    for section in fp.sections
        assemble_K_section(section, fp.nodes, dof_rows, dof_cols, k_values)
    end

    return Base.sparse(dof_rows, dof_cols, k_values, fp.n_eqs, fp.n_eqs)
end

function assemble_K_section{T<:FESection}(section::T, nodes::Vector{FENode2},
                                          dof_rows::Vector{Int}, dof_cols::Vector{Int},
                                          k_values::Vector{Float64})
    mat = section.material
    for element in section.elements
        Ke = stiffness(element, nodes, mat)
        dof1_n = 0
        for vertex1 in element.vertices
            for dof1 in nodes[vertex1].dofs
                dof1_n += 1
                if dof1.active
                    dof2_n = 0
                    for vertex2 in element.vertices
                        for dof2 in nodes[vertex2].dofs
                            dof2_n += 1
                            if dof2.active
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
end

function assemble_intf(fp::FEProblem)
    fint = zeros(fp.n_eqs)
    for section in fp.sections
        assemble_intf_section(section, fint, fp.nodes)
    end
    return fint
end

function assemble_intf(fp::FEProblem)
    fint = zeros(fp.n_eqs)
    for section in fp.sections
        assemble_intf_section(section, fint, fp.nodes)
    end
    return fint
end

function assemble_intf_section{T<:FESection}(section::T,
                                            int_forces::Vector{Float64},
                                            nodes::Vector{FENode2})
    mat = section.material
    for element in section.elements
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

function updatedofs!(fp::FEProblem, du::Vector{Float64})
    for node in fp.nodes
        for dof in node.dofs
            if dof.active
                dof.value += du[dof.eq_n]
            end
        end
    end
end

function updatebcs!(fp::FEProblem, t::Number=0.0)
    for node in fp.nodes
        for dof in node.dofs
            if !dof.active
                bc = fp.node_doftype_bc[(node.n, dof.dof_type)]
                updatebc!(bc, dof, node, t)
            end
        end
    end
end

@inline function updatebc!{f}(bc::DirichletBC{f}, dof::Dof, node::FENode2, t::Number)
    dof.value = evaluate(bc, node, t)
end


function update_feproblem(fp::FEProblem)
    for section in fp.sections
        for element in section.elements
            for i in 1:length(element.matstats)
                element.matstats[i] = copy(element.temp_matstats[i])
            end
        end
    end
end


function FEProblem(name::ASCIIString, nodes::Vector{FENode2}, bcs,
                    loads, sections=Array(FESection, 0))
    node_doftype_bc = Dict{Int, Vector{DofType}}()
    node_doftypes = Dict{Int, Vector{DofType}}()
    FEProblem(name, nodes, bcs, loads, sections, node_doftypes, node_doftype_bc, 0, 0)
end

push!(fp::FEProblem, bc::DirichletBC) = push!(fp.bcs, bc)
push!(fp::FEProblem, load::NodeLoad) = push!(fp.loads, load)
push!(fp::FEProblem, section::FESection) = push!(fp.sections, section)


function create_feproblem_grad(name, geomesh, element_regions, material_regions, bcs::Vector{DirichletBC}=Array(DirichletBC, 0), loads::Vector{NodeLoad}=Array(NodeLoad, 0))

    gps = Dict{DataType, Vector{GaussPoint2}} ()
    interps = Dict{DataType, AbstractInterpolator} ()
    storage = Dict{DataType, AbstractElemStorage} ()
    elem_types = Array(DataType, 0)

    for element_region in element_regions
        elem_type = element_region.element_type
        interps[elem_type] = createinterp(elem_type)
        gps[elem_type] = creategps(elem_type)
        storage[elem_type] = createstorage(elem_type)
    end

    interp_grad = LinTrigInterp()

    nodes = Array(FENode2, 0)
    for node in geomesh.nodes
        push!(nodes, FENode2(node.n, node.coords))
    end

    sections = Array(FESection, 0)
    for matregion in material_regions
        material = matregion.material
        matstat = create_matstat(typeof(material))
        for eleregion in element_regions
            ele_type = eleregion.element_type
            common = intersect(matregion.elements, eleregion.elements)
            common = collect(common)
            sort!(common) # TODO: Remove?
            gps_ele = gps[ele_type]
            elem_storage = storage[ele_type]
            interp = interps[ele_type]

            section = FESection(material, ele_type, typeof(matstat))

            for ele_id in common
                vertices = geomesh.elements[ele_id].vertices
                element = ele_type(vertices, ele_id, interp, interp_grad,
                                   elem_storage, gps_ele, matstat)
                push!(section, element)

            end
            push!(sections, section)
        end
    end
    fe = FEProblem(name, nodes, bcs, loads, sections)
    createdofs(fe)
    return fe
end
