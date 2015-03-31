
immutable FESection{P <: AbstractFElement, T <: AbstractMaterial}
    elements::Vector{P}
    material::T
    ele_type::Type{P}
end

function FESection{P <: AbstractFElement, T <: AbstractMaterial}(mat::T, ele_type::Type{P})
    elements = Array(ele_type, 0)
    FESection(elements, mat, ele_type)
end


immutable MaterialSection{P <: AbstractMaterial}
    material::P
    elements::Set{Int}
end
MaterialSection{P <: AbstractMaterial}(mat::P) = MaterialSection(mat, Set{Int}())


immutable ElementSection{P <: AbstractFElement}
    element_type::Type{P}
    elements::Set{Int}
end

ElementSection{P <: AbstractFElement}(ele_type::Type{P}) = ElementSection{P}(ele_type, Set{Int}())


# TODO: Clean this up
function push!(section::FESection, elem::AbstractFElement)
    push!(section.elements, elem)
end


function push!{P <: AbstractMaterial}(section::MaterialSection{P}, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end

function push!{P <: AbstractFElement}(section::ElementSection{P}, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end

##########
# FEMesh #
##########
#=
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
            interps[elem_type] = get_interp(elem)
            gps[elem_type] = get_gps(elem)
            storage[elem_type] = get_storage(elem)
        end
        elements[elem2.n] = elem(elem2.vertices, elem2.n, interps[elem_type],
                                        storage[elem_type],gps[elem_type])
    end

    nodes = Array(FENode2, 0)

    for node in geomesh.nodes
        push!(nodes, FENode2(node.n, node.coords))
    end
    element_sets = Dict{ASCIIString, ElementSet}()
    node_sets = Dict{ASCIIString, NodeSet}()
    FEMesh(nodes, elements, element_sets, node_sets)
end

function show(io::IO,mesh::FEMesh)
    print(io, string("FEMesh, ", length(mesh.elements)))
end
=#