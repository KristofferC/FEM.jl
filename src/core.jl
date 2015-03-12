#typealias Vec3 Vector3{Float64}
#typealias Vec2 Vector2{Float64}
#typealias Mat2 Matrix2x2{Float64}
#typealias Mat3 Matrix3x3{Float64}

typealias MatPool Dict{(Int, Int, ASCIIString), Array{Float64, 2}}
typealias VecPool Dict{(Int, ASCIIString), Array{Float64, 1}}

immutable GaussPoint
    local_coords::Vector{Float64}
    weight::Float64
end


abstract DofType

type Du <: DofType end
type Dv <: DofType end
type Dw <: DofType end

type Dof
    eq_n::Int
    id::Int
    active::Bool
    value::Float64
    dof_type::DofType # Not here
end


immutable Node
    coordinates::Vector{Float64}
    n::Int
    dofs::Vector{Dof}
end

function Node(c::Vector{Float64}, n::Int)
    Node(c, n, Array(Dof, 0))
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

function getmat(rows::Int, cols::Int, name::ASCIIString,
                matpool::MatPool)
    get!(() -> Array(Float64, rows, cols), matpool, (rows, cols, name))
end

function getvec(rows::Int, name::ASCIIString,
                vecpool::VecPool)
    get!(() -> Array(Float64, rows), vecpool, (rows, name))
end
