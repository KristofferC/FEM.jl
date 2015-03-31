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


# Shape functions in local coords
function Nvec(interp::LinQuadInterp, loc_coords::Point2)

    ξ = loc_coords[1]
    η = loc_coords[2]

    interp.N[1] = (1.0 + ξ) * (1.0 + η) * 0.25
    interp.N[2] = (1.0 - ξ) * (1.0 + η) * 0.25
    interp.N[3] = (1.0 - ξ) * (1.0 - η) * 0.25
    interp.N[4] = (1.0 + ξ) * (1.0 - η) * 0.25

    return interp.N
end


function dNmatrix(interp::LinQuadInterp, loc_coords::Point2)
    ξ = loc_coords[1]
    η = loc_coords[2]


    interp.dN[1,1] =  (1.0 + η) * 0.25
    interp.dN[1,2] =  (1.0 + ξ) * 0.25

    interp.dN[2,1] = -(1.0 + η) * 0.25
    interp.dN[2,2] =  (1.0 - ξ) * 0.25

    interp.dN[3,1] = -(1.0 - η) * 0.25
    interp.dN[3,2] = -(1.0 - ξ) * 0.25

    interp.dN[4,1] =  (1.0 - η) * 0.25
    interp.dN[4,2] = -(1.0 + ξ) * 0.25

    return interp.dN
end

function Jmatrix(interp::LinQuadInterp, local_coords::Point2,
                 vertices::Vertex4, nodes::Vector{FENode2},
                 dN::Matrix{Float64})

    fill!(interp.J, 0.0)
    for row in 1:size(dN, 1)
        x = nodes[vertices[row]].coords.x
        y = nodes[vertices[row]].coords.y


        interp.J[1, 1] += dN[row, 1] * x
        interp.J[1, 2] += dN[row, 1] * y
        interp.J[2, 1] += dN[row, 2] * x
        interp.J[2, 2] += dN[row, 2] * y
    end
    return interp.J
end


function dNdxmatrix(interp::LinQuadInterp, local_coords::Point2,
                    vertices::Vertex4, nodes::Vector{FENode2})
        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, local_coords, vertices, nodes, dN)
        A_mul_B!(interp.dNdx, dN, inv2x2t!(J))
        return interp.dNdx
end



