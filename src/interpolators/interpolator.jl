abstract Interpolator

include("lin_trig_interp.jl")
include("bilin_quad_interp.jl")

function dNdxmatrix(interp::Interpolator, local_coords::Vector{Float64},
                    vertices::Vector{Int}, nodes::Vector{Node}, mp::MatPool,
                    vp::VecPool)

        dN = dNmatrix(interp, local_coords, mp, vp)
        J = Jmatrix(interp, local_coords, vertices, nodes, dN, mp, vp)

        dNdx = dN * (inv(J)')

        return dNdx
end


function Jmatrix(interp::Interpolator, local_coords::Vector{Float64},
                 vertices::Vector{Int}, nodes::Vector{Node},
                 dN::Matrix{Float64}, mp::MatPool, vp::VecPool)

    J = getmat(2, 2, "J", mp)
    J[1, 1] = J[1, 2] = J[2, 1] = J[2,2] = 0.0

    for row in 1:size(dN, 1)
        x = nodes[vertices[row]].coordinates[1]
        y = nodes[vertices[row]].coordinates[2]


        J[1, 1] += dN[row, 1] * x
        J[1, 2] += dN[row, 1] * y
        J[2, 1] += dN[row, 2] * x
        J[2, 2] += dN[row, 2] * y
    end

    return J
end
