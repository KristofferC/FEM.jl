#typealias Vertex3 Vector3{Float64}
#typealias Vertex2 Vector2{Float64}

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
    dof_id::Int # TODO: rename dof_id to id
    active::Bool
    value::Float64
    dof_type::DofType
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

immutable ElementSet
    name::String
    element_ids::Vector{Int}
end


immutable DirichletBC
    value::Float64
    dof_types::Vector{DofType}
    node_set::NodeSet
end


immutable PointLoad
    value::Float64
    dof_types::Vector{DofType}
    node_set::NodeSet
end




