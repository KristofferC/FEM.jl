immutable LinTrigInterp <: Interpolator
    N::Vector{Float64}
    dN::Matrix{Float64}
    dNdx::Matrix{Float64}
    J::Matrix{Float64}
end

function LinTrigInterp()
    N = Array(Float64, 3)

    # Derivatives are constant
    #       ξ     η
    dN = [[ 1.0  0.0];
          [ 0.0  1.0];
          [-1.0 -1.0]]

    dNdx = Array(Float64, 3, 2)

    J = Array(Float64, 2, 2)
    LinTrigInterp(N, dN, dNdx, J)
end


# Shape functions in local coordinates
function Nvec(interp::LinTrigInterp, loc_coords::Point2)

    ξ = loc_coords[1]
    η = loc_coords[2]

    interp.N[1] = ξ
    interp.N[2] = η
    interp.N[3] = 1.0 - ξ - η

    return interp.N
end


function dNmatrix(interp::LinTrigInterp, loc_coords::Point2)
    return interp.dN
end

function Jmatrix(interp::LinTrigInterp, local_coords::Point2,
                 vertices::Vertex3, nodes::Vector{FENode2},
                 dN::Matrix{Float64})



    x1 =  nodes[vertices[1]].coordinates.x
    x2 =  nodes[vertices[2]].coordinates.x
    x3 =  nodes[vertices[3]].coordinates.x

    y1 =  nodes[vertices[1]].coordinates.y
    y2 =  nodes[vertices[2]].coordinates.y
    y3 =  nodes[vertices[3]].coordinates.y

    # Constant Jacobian
    interp.J[1, 1] = x1 - x3
    interp.J[2, 1] = x2 - x3
    interp.J[1, 2] = y1 - y3
    interp.J[2, 2] = y2 - y3

    return interp.J
end

function dNdxmatrix(interp::LinTrigInterp, local_coords::Point2,
                    vertices::Vertex3, nodes::Vector{FENode2})

        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, local_coords, vertices, nodes, dN)
        invJt = inv2x2t!(J)

        A_mul_B!(interp.dNdx, dN, invJt)
        #dNdx = dN * invJt

        return interp.dNdx
end

function inv2x2t!(J::Matrix{Float64})
    d = det2x2(J)
    tmp = J[1,1]
    J[1,1] = 1.0 / d * J[2,2]
    J[2,2] = tmp / d

    tmp = J[1,2]
    J[1,2] = -J[2,1]/ d
    J[2,1] = -tmp / d
    return J
end


det2x2(J::Matrix{Float64}) = J[1,1]*J[2,2] - J[1,2]*J[2,1]


