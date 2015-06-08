#typealias Vec3 Vector3{Float64}
#typealias Vec2 Vector2{Float64}
#typealias Mat2 Matrix2x2{Float64}
#typealias Mat3 Matrix3x3{Float64}

import Base: Func, AddFun

immutable SubFun <: Func{2} end
call(::SubFun, x, y) = x - y

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
@lintpragma("Ignore use of undeclared variable Gu1")
@lintpragma("Ignore use of undeclared variable Gv1")
@lintpragma("Ignore use of undeclared variable Gu2")
@lintpragma("Ignore use of undeclared variable Gv2")

@enum DofType Du Dv Dw Gu1 Gv1 Gu2 Gv2 Gu3 Gv3 K1 K2


immutable Dof
    eq_n::Int
    dof_type::DofType # Not here
    active::Bool
end
@inline isactive(dof::Dof) = dof.active

immutable DofVals
    free_dof_values::Vector{Float64}
    presc_dof_values::Vector{Float64}
end
DofVals() = DofVals(zeros(0), zeros(0))

function get_value(dv::DofVals, dof::Dof)
    if isactive(dof)
        return dv.free_dof_values[dof.eq_n]
    else
        return dv.presc_dof_values[dof.eq_n]
    end
end

function mod_value(dv::DofVals, dof::Dof, f::Func, val)
    if isactive(dof)
        dv.free_dof_values[dof.eq_n] = f(dv.free_dof_values[dof.eq_n], val)
    else
        dv.presc_dof_values[dof.eq_n] =  f(dv.presc_dof_values[dof.eq_n] , val)
    end
end



type TimeStepper{T <: Range}
    timesteps::T
    t::Float64
    dt::Float64
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
@inline get_dof(node::AbstractFENode, i::Int) = node.dofs[i]
@inline get_dofs(node::AbstractFENode) = node.dofs



#TODO: Dispatch on dimension
FENode2(n::Int, c::Point2) = FENode2(n, c, Array(Dof, 0))
# Check size here?
#FENode2(n::Int, c::Vector{Float64}) = FENode2(n, Point2(c[1], c[2]), Array(Dof, 0))
function FENode2{T <: Number}(n::Int, c::Vector{T})
    c_v = convert(Vector{Float64}, c)
    FENode2(n, Point2(c[1], c[2]), Array(Dof, 0))
end

# TODO: Move
@inline get_coord(node::FENode2) = Point3(node.coords[1], node.coords[2], 0.0)
@inline get_displacement(node::FENode2, dv::DofVals) = Point3(get_value(dv, node.dofs[1]),
                                                              get_value(dv, node.dofs[2]),
                                                              0.0)


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
#=
immutable DirichletBC
    value::Float64
    dof_types::Vector{DofType}
    node_set::NodeSet
end
=#







immutable DirichletBC{f <: FastAnonymous.Fun}
    func::f
    dof_types::Vector{DofType}
    node_set::NodeSet
end
function DirichletBC(f::String, dof_types::Vector{DofType}, node_set::NodeSet)
    prog = string("ff = @anon (x,y,z,t) -> ", f)
    eval(parse(prog))
    DirichletBC(ff, dof_types,  node_set)
end

immutable NodeLoad{f <: FastAnonymous.Fun}
    func::f
    dof_types::Vector{DofType}
    node_set::NodeSet
end

function NodeLoad(f::String, dof_types::Vector{DofType}, node_set::NodeSet)
    prog = string("ff = @anon (x,y,z,t) -> ", f)
    eval(parse(prog))
    NodeLoad(ff, dof_types,  node_set)
end



function evaluate(bc::Union(DirichletBC, NodeLoad), node::FENode2, t::Number)
    return bc.func(node.coords.x, node.coords.y, 0.0, t)
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

function assemble!{T <: Real}(A::Array, X::AbstractArray, I::AbstractVector{T})
   count = 1
   for i in I
        A[count] = X[i]
        count += 1
   end
end

