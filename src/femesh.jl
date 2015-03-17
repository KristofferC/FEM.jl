
##########
# FEMesh #
##########
abstract AbstractFEMesh

immutable FEMesh <: AbstractFEMesh
    nodes::Vector{FENode2}
    #elements::Dict{DataType, Vector}
    elements::Dict{Int, LinTrig}
    element_sets::Dict{ASCIIString, ElementSet}
    node_sets::Dict{ASCIIString, NodeSet}
end

function FEMesh(geomesh::GeoMesh)

    #elements = Dict{DataType, Vector} ()
    elements = Dict{Int, LinTrig} ()
    gps = Dict{DataType, Vector{GaussPoint2}} ()
    interps = Dict{DataType, Interpolator} ()
    storage = Dict{DataType, ElemStorage} ()
    elem_types = Array(DataType, 0)

    #geo_to_fe = Dict{Int, Int} ()

    for (i, elem2) in enumerate(geomesh.elements)
        elem = LinTrig
        elem_type = Type{elem}
        if !(elem_type in elem_types)
            push!(elem_types, elem_type)
            #elements[elem_type] = Array(elem, 0)
            interps[elem_type] = get_interp(elem)
            gps[elem_type] = get_gps(elem)
            storage[elem_type] = get_storage(elem)
        end
        elements[elem2.n] = elem(elem2.vertices, elem2.n, interps[elem_type],
                                        storage[elem_type],gps[elem_type])
    end

    nodes = Array(FENode2, 0)

    for node in geomesh.nodes
        push!(nodes, FENode2(node.n, node.coordinates))
    end
    element_sets = Dict{ASCIIString, ElementSet}()
    node_sets = Dict{ASCIIString, NodeSet}()
    FEMesh(nodes, elements, element_sets, node_sets)
end

function show(io::IO,mesh::FEMesh)
    print(io, string("FEMesh, ", length(mesh.elements)))
end

immutable Section
    material::AbstractMaterial
    elements::Set{Int} # The elements in the section
end


function Section(mat::AbstractMaterial)
    elements = Set{Int}()
    Section(mat, elements)
end

function push!(section::Section, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end
