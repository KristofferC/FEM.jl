immutable LinTrigInterp <: Interpolator end

# Shape functions in local coordinates
function Nvec(interp::LinTrigInterp, loc_coords::Vector{Float64})

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = [ξ;
        η;
        1.0 - ξ - η]

    return N
end

function dNdxmatrix(interp::LinTrigInterp, local_coords::Vector{Float64},
                    vertices::Vector{Int}, nodes::Vector{Node})

        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, local_coords, vertices, nodes, dN)
        dNdx = dN * (inv(J)')

        return dNdx
end


function dNmatrix(interp::LinTrigInterp, loc_coords::Vector{Float64})
    # Derivative w.r.t ξ
    dN = [[ 1.0  0.0];
          [ 0.0  1.0];
          [-1.0 -1.0]]

    return dN
end


function Jmatrix(interp::LinTrigInterp, local_coords::Vector{Float64},
                 vertices::Vector{Int}, nodes::Vector{Node},
                 dN::Matrix{Float64})

    J = zeros(2, 2)

    for row in 1:3
        x = nodes[vertices[row]].coordinates[1]
        y = nodes[vertices[row]].coordinates[2]

        J[1, 1] += dN[row, 1] * x
        J[1, 2] += dN[row, 1] * y
        J[2, 1] += dN[row, 2] * x
        J[2, 2] += dN[row, 2] * y
    end

    return J
end
