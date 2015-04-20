immutable LinQuadInterp <: AbstractInterpolator
    N::Vector{Float64}
    dN::Matrix{Float64}
    dNdx::Matrix{Float64}
    J::Matrix{Float64}
end

function LinQuadInterp()
    N = zeros(4)
    dN = zeros(4,2)
    dNdx = Array(Float64, 4, 2)
    J = Array(Float64, 2, 2)
    LinQuadInterp(N, dN, dNdx, J)
end
get_allocated_N(i::LinQuadInterp) = i.N
get_allocated_dN(i::LinQuadInterp) = i.dN
get_allocated_dNdx(i::LinQuadInterp) = i.dNdx
get_allocated_J(i::LinQuadInterp) = i.J


# Shape functions in local coords
function Nvec(interp::LinQuadInterp, loc_coords::Point2)

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = get_allocated_N(interp)

    N[1] = (1 + ξ) * (1 + η) * 0.25
    N[2] = (1 - ξ) * (1 + η) * 0.25
    N[3] = (1 - ξ) * (1 - η) * 0.25
    N[4] = (1 + ξ) * (1 - η) * 0.25

    return N
end


function dNmatrix(interp::LinQuadInterp, loc_coords::Point2)
    ξ = loc_coords[1]
    η = loc_coords[2]

    dN = get_allocated_dN(interp)

   dN[1,1] =  (1 + η) * 0.25
   dN[1,2] =  (1 + ξ) * 0.25

   dN[2,1] = -(1 + η) * 0.25
   dN[2,2] =  (1 - ξ) * 0.25

   dN[3,1] = -(1 - η) * 0.25
   dN[3,2] = -(1 - ξ) * 0.25

   dN[4,1] =  (1 - η) * 0.25
   dN[4,2] = -(1 + ξ) * 0.25

    return dN
end

function Jmatrix(interp::LinQuadInterp,
                 vertices::Vertex4, nodes::Vector{FENode2},
                 dN::Matrix{Float64})

    J = get_allocated_J(interp)
    fill!(J, 0.0)
    for row in 1:size(dN, 1)
        x = nodes[vertices[row]].coords.x
        y = nodes[vertices[row]].coords.y

        J[1, 1] += dN[row, 1] * x
        J[1, 2] += dN[row, 1] * y
        J[2, 1] += dN[row, 2] * x
        J[2, 2] += dN[row, 2] * y
    end
    return J
end


function dNdxmatrix(interp::LinQuadInterp, local_coords::Point2,
                    vertices::Vertex4, nodes::Vector{FENode2})
        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, vertices, nodes, dN)
        A_mul_B!(interp.dNdx, dN, inv2x2t!(J))
        return interp.dNdx
end



