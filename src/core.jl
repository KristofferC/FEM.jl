#typealias Vec3 Vector3{Float64}
#typealias Vec2 Vector2{Float64}
#typealias Mat2 Matrix2x2{Float64}
#typealias Mat3 Matrix3x3{Float64}



abstract AbstractGaussPoint

immutable GaussPoint2 <: AbstractGaussPoint
    local_coords::Point2
    weight::Float64
end

immutable GaussPoint3 <: AbstractGaussPoint
    local_coords::Point3
    weight::Float64
end

getweight(gp::AbstractGaussPoint) = gp.weight

# Lint does not refognize enums
@lintpragma("Ignore use of undeclared variable DofType")
@lintpragma("Ignore use of undeclared variable Du")
@lintpragma("Ignore use of undeclared variable Dv")
@lintpragma("Ignore use of undeclared variable Dw")
@lintpragma("Ignore use of undeclared variable Pressure")

@enum DofType Du Dv Dw

type Dof
    eq_n::Int
    id::Int
    active::Bool
    value::Float64
    dof_type::DofType # Not here
end


#########
# Nodes #
#########
abstract AbstractFENode

immutable FENode2 <: AbstractFENode
    n::Int
    coords::Point2
    dofs::Vector{Dof}
end

#TODO: Dispatch on dimension
FENode2(n::Int, c::Point2) = FENode2(n, c, Array(Dof, 0))
# Check size here?
#FENode2(n::Int, c::Vector{Float64}) = FENode2(n, Point2(c[1], c[2]), Array(Dof, 0))
function FENode2{T <: Number}(n::Int, c::Vector{T})
    c_v = convert(Vector{Float64}, c)
    FENode2(n, Point2(c[1], c[2]), Array(Dof, 0))
end

# TODO: Move


get_coord(node::FENode2) = Point3(node.coords[1], node.coords[2], 0.0)
get_displacement(node::FENode2) = Point3(node.dofs[1].value, node.dofs[2].value, 0.0)


immutable FENode3 <: AbstractFENode
    n::Int
    coords::Point3
    dofs::Vector{Dof}
end

FENode3(n::Int, c::Point3) = FENode2(n, c, Array(Dof, 0))
# Check size here?
#FENode3(n::Int, c::Vector{Float64}) = FENode2(n, Point3(c[1], c[2], c[3]), Array(Dof, 0))
#FENode3(n::Int, c::Vector{Int}) = FENode2(n, convert(Vector{Float64}, c), Array(Dof, 0))
function FENode3{T <: Number}(n::Int, c::Vector{T})
    c_v = convert(Vector{Float64}, c)
    FENode2(n, Point3(c[1], c[2], c[3]), Array(Dof, 0))
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

function fill_from_start!{T}(v1::AbstractArray{T}, v2::AbstractArray{T})
    if length(v2) > length(v1)
        error("Array to be filled is shorter than array to fill with")
    end
    for i in 1:length(v2)
        v1[i] = v2[i]
    end
end
