abstract AbstractInterpolator

include("lin_trig_interp.jl")
include("bilin_quad_interp.jl")

function dNdxmatrix{T <: AbstractInterpolator}(interp::T, local_coords::Point2,
                    vertices::Vector{Int}, nodes::Vector{FENode2})

        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, local_coords, vertices, nodes, dN)
        dNdx = dN * (inv(J)')

        return dNdx
end


#=
function Jmatrix{T <: AbstractInterpolator}(::T, local_coords::Point2,
                 vertices::Vector{Int}, nodes::Vector{FENode2},
                 dN::Matrix{Float64})

    J = zeros(2, 2)


    for row in 1:size(dN, 1)
        x = nodes[vertices[row]].coords[1]
        y = nodes[vertices[row]].coords[2]


        J[1, 1] += dN[row, 1] * x
        J[1, 2] += dN[row, 1] * y
        J[2, 1] += dN[row, 2] * x
        J[2, 2] += dN[row, 2] * y
    end

    return J
end
=#