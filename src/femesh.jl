
##########
# FEMesh #
##########
abstract AbstractFEMesh

immutable FEMesh <: AbstractFEMesh
    nodes::Vector{FENode2}
    elements::Dict{DataType, Vector}
    element_sets::Dict{ASCIIString, ElementSet}
    node_sets::Dict{ASCIIString, NodeSet}
end

function FEMesh(geomesh::GeoMesh)

    elements = Dict{DataType, Vector} ()
    gps = Dict{DataType, Vector{GaussPoint2}} ()
    interps = Dict{DataType, Interpolator} ()
    storage = Dict{DataType, ElemStorage} ()

    #geo_to_fe = Dict{Int, Int} ()



    for (i, elem2) in enumerate(geomesh.elements)
        elem = LinTrig
        elem_type = Type{elem}
        if !(elem_type in keys(elements))
            elements[elem_type] = Array(elem, 0)
            interps[elem_type] = get_interp(elem)
            gps[elem_type] = get_gps(elem)
            storage[elem_type] = get_storage(elem)
        end
        push!(elements[elem_type], elem(elem2.vertices, i, interps[elem_type],
                                        storage[elem_type],gps[elem_type]))
    end

    nodes = Array(FENode2, 0)
    element_sets = Dict{ASCIIString, ElementSet}()
    node_sets = Dict{ASCIIString, NodeSet}()
    FEMesh(nodes, elements, element_sets, node_sets)
end

function show(io::IO,mesh::FEMesh)
    print(io, string("FEMesh, ", length(mesh.elements)))
end

immutable Section
    material::Material
    element_set::ElementSet # The elements in the section
end


function Section(mat::Material)
    elements = Set(Int[])
    Section(mat, elements)
end

function addelemset!(section::Section, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end
