# Defines data types representing a purely geometrical mesh.

##########
# Points #
##########
immutable Point2<: FixedVector{Float64, 2}
    x::Float64
    y::Float64
end
typealias Vector2 Point2

immutable Point3 <: FixedVector{Float64, 3}
    x::Float64
    y::Float64
    z::Float64
end
typealias Vector3 Point3

#########
# Nodes #
#########
abstract AbstractGeoNode

immutable GeoNode2 <: AbstractGeoNode
    n::Int
    coords::Point2
end
# Check size here?
GeoNode2(n::Int, c::Vector{Float64}) = GeoNode2(n, Point2(c[1], c[2]))
GeoNode2(n::Int, c::Vector{Int}) = GeoNode2(n, convert(Vector{Float64}, c))
get_coord(node::GeoNode2) = Point3(node.coords.x, node.coords.y, 0.0)


immutable GeoNode3 <: AbstractGeoNode
    n::Int
    coords::Point3
end

# Check size here?
GeoNode3(n::Int, c::Vector{Float64}) = GeoNode3(n, Point3(c[1], c[2], c[3]))
GeoNode3(n::Int, c::Vector{Int}) = GeoNode3(n, convert(Vector{Float64}, c))


############
# Vertices #
############
abstract AbstractVertex{L} <: FixedVector{Int, L}

immutable Vertex2 <: AbstractVertex{2}
    v1::Int
    v2::Int
end

immutable Vertex3 <: AbstractVertex{3}
    v1::Int
    v2::Int
    v3::Int
end

immutable Vertex4 <: AbstractVertex{4}
    v1::Int
    v2::Int
    v3::Int
    v4::Int
end

immutable Vertex5 <: AbstractVertex{5}
    v1::Int
    v2::Int
    v3::Int
    v4::Int
end

immutable Vertex6 <: AbstractVertex{6}
    v1::Int
    v2::Int
    v3::Int
    v4::Int
    v5::Int
    v6::Int
end

###############
# GeoElements #
###############
abstract AbstractGeoElem
getindex(geoelem::AbstractGeoElem, i0::Real) = getindex(geoelem.vertices, i0)


abstract AbstractGeoElem2 <: AbstractGeoElem
abstract AbstractGeoElem3 <: AbstractGeoElem

immutable GeoTrig <: AbstractGeoElem2
    n::Int
    vertices::Vertex3
end
GeoTrig(n::Int, v::Vector{Int}) = GeoTrig(n, Vertex3(v[1], v[2], v[3]))


immutable GeoQTrig <: AbstractGeoElem2
    n::Int
    vertices::Vertex6
end
GeoQTrig(n::Int, v::Vector{Int}) = GeoQTrig(n, Vertex6(v[1], v[2], v[3], v[4], v[5], v[6]))


immutable GeoQuad <: AbstractGeoElem2
    n::Int
    vertices::Vertex4
end
GeoQuad(n::Int, v::Vector{Int}) = GeoQuad(n, Vertex4(v[1], v[2], v[3], v[4]))


immutable GeoTetra <: AbstractGeoElem3
    n::Int
    vertices::Vertex4
end
GeoTetra(n::Int, v::Vector{Int}) = GeoTetra(n, Vertex3(v[1], v[2], v[3], v[4]))

########
# Sets #
########
immutable ElementSet
    name::ASCIIString
    element_ids::Set{Int}
end
ElementSet(s::ASCIIString, elems) = ElementSet(s, Set{Int64}(elems))
Base.append!(ns1::ElementSet, elems) = union!(ns1.element_ids, Set{Int64}(element_ids))
Base.append!(es1::ElementSet, es2::ElementSet) = union!(es1.element_ids, ne2.element_ids)


immutable NodeSet
    name::ASCIIString
    node_ids::Set{Int}
end
NodeSet(s, nodes) = NodeSet(s, Set{Int64}(nodes))

Base.append!(ns1::NodeSet, nodes) = union!(ns1.node_ids, Set{Int64}(nodes))
Base.append!(ns1::NodeSet, ns2::NodeSet) = union!(ns1.node_ids, ns2.node_ids)



immutable EdgeSet
    name::String
    # (Element_id, edge_id) tuple
    edges::Vector{Tuple{Int, Int}}
end

immutable SurfaceSet
    name::String
    # (Element_id, surface_id) tuple
    surfaces::Vector{Tuple{Int, Int}}
end



###########
# GeoMesh #
###########
abstract AbstractGeoMesh

immutable GeoMesh <: AbstractGeoMesh
    nodes::Vector{GeoNode2}
    elements::Vector{AbstractGeoElem2}
    element_sets::Dict{ASCIIString, ElementSet}
    node_sets::Dict{ASCIIString, NodeSet}
end

function GeoMesh()
    nodes = Array(GeoNode2, 0)
    elements = Array(AbstractGeoElem2, 0)
    element_sets = Dict{ASCIIString, ElementSet}()
    node_sets = Dict{ASCIIString, NodeSet}()
    GeoMesh(nodes, elements, element_sets, node_sets)
end


# Push functions
push!(mesh::GeoMesh, node::AbstractGeoNode) = push!(mesh.nodes, node)

function push!(mesh::GeoMesh, nodes::Vector{GeoNode2})
    for n in nodes
        push!(mesh, n)
    end
end

push!(mesh::GeoMesh, elem::AbstractGeoElem) = push!(mesh.elements, elem)


function push!(mesh::GeoMesh, elem_set::ElementSet)
    mesh.element_sets[elem_set.name] = elem_set
end

function push!(mesh::GeoMesh, node_set::NodeSet)
     mesh.node_sets[node_set.name] = node_set
end


# Set generation functions
function gennodeset(f::Function, name::ASCIIString, nodes::Vector{GeoNode2})
    node_ids = Set{Int}()
    for node in nodes
        if f(node)
            push!(node_ids, node.n)
        end
    end
    return NodeSet(name, node_ids)
end
