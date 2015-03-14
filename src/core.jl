#typealias Vec3 Vector3{Float64}
#typealias Vec2 Vector2{Float64}
#typealias Mat2 Matrix2x2{Float64}
#typealias Mat3 Matrix3x3{Float64}

immutable Point2<: FixedVector{Float64, 2}
    x::Float64
    y::Float64
end

immutable Point3 <: FixedVector{Float64, 3}
    x::Float64
    y::Float64
    z::Float64
end

immutable GaussPoint
    local_coords::Point2
    weight::Float64
end


@enum DofType Du Dv Dw

#type Du <: DofType end
#type Dv <: DofType end
#type Dw <: DofType end

type Dof
    eq_n::Int
    id::Int
    active::Bool
    value::Float64
    dof_type::DofType # Not here
end


immutable Node
    coordinates::Point2
    n::Int
    dofs::Vector{Dof}
end

function Node(c::Vector{Float64}, n::Int)
    Node(Point2(c[1], c[2]), n, Array(Dof, 0))
end

Node(c::Vector{Int}, n::Int) = Node(convert(Vector{Float64}, c), n::Int)


immutable NodeSet
    name::String
    node_ids::Vector{Int}
end

function gennodeset(f::Function, name::ASCIIString, nodes::Vector{Node})
    node_ids = Int[]
    for node in nodes
        if f(node.coordinates)
            push!(node_ids, node.n)
        end
    end
    return NodeSet(name, node_ids)
end

immutable ElementSet
    name::String
    element_ids::Vector{Int}
end

immutable EdgeSet
    name::String
    # (Element_id, edge_id) tuple
    edges::Vector{(Int, Int)}
end

immutable SurfaceSet
    name::String
    # (Element_id, surface_id) tuple
    surfaces::Vector{(Int, Int)}
end


# TODO: Add more sets
# Only constant loads right now

immutable DirichletBC
    value::Float64
    dof_types::Vector{DofType}
    node_set::NodeSet
end

immutable NodeLoad
    value::Float64
    dof_types::Vector{DofType}
    node_set::NodeSet
end

immutable EdgeLoad
    value::Float64
    dof_types::Vector{DofType}
    edge_set::EdgeSet
end

# TODO Add more loads.
