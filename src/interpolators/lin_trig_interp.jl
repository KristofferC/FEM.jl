module LinTrigInterpMod

using FEM

import FEM: AbstractInterpolator, Point2, FENode2, Vertex3, inv2x2t!, det2x2,
       dNmatrix, Jmatrix, dNdxmatrix, Nvec, get_area, mass_matrix, get_area

export LinTrigInterp

immutable LinTrigInterp <: AbstractInterpolator
    N::Vector{Float64}
    dN::Matrix{Float64}
    dNdx::Matrix{Float64}
    J::Matrix{Float64}
    ref_M::Matrix{Float64}
    M::Matrix{Float64}
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

    ref_M = 1/12 *
                     [[2.0  0.0  1.0  0.0  1.0  0.0];
                      [0.0  2.0  0.0  1.0  0.0  1.0];
                      [1.0  0.0  2.0  0.0  1.0  0.0];
                      [0.0  1.0  0.0  2.0  0.0  1.0];
                      [1.0  0.0  1.0  0.0  2.0  0.0];
                      [0.0  1.0  0.0  1.0  0.0  2.0]]

     ref_M2 = 1/12 *
                     [2.0 1.0 1.0;
                      1.0 2.0 1.0;
                      1.0 1.0 2.0]

    M = zeros(3,3)

    LinTrigInterp(N, dN, dNdx, J, ref_M2, M)
end
get_allocated_N(i::LinTrigInterp) = i.N
get_allocated_dN(i::LinTrigInterp) = i.dN
get_allocated_dNdx(i::LinTrigInterp) = i.dNdx
get_allocated_J(i::LinTrigInterp) = i.J


# Shape functions in local coords
function Nvec(interp::LinTrigInterp, loc_coords::Point2)

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = get_allocated_N(interp)

    N[1] = ξ
    N[2] = η
    N[3] = 1.0 - ξ - η

    return N
end


function dNmatrix(interp::LinTrigInterp, ::Point2)
    # dN is constant so just return it
    return get_allocated_dN(interp)
end

function Jmatrix(interp::LinTrigInterp,
                 vertices::Vertex3, nodes::Vector{FENode2},
                 ::Matrix{Float64})

    J = get_allocated_J(interp)

    x1 =  nodes[vertices[1]].coords.x
    x2 =  nodes[vertices[2]].coords.x
    x3 =  nodes[vertices[3]].coords.x

    y1 =  nodes[vertices[1]].coords.y
    y2 =  nodes[vertices[2]].coords.y
    y3 =  nodes[vertices[3]].coords.y

    # Constant Jacobian
    J[1, 1] = x1 - x3
    J[2, 1] = x2 - x3
    J[1, 2] = y1 - y3
    J[2, 2] = y2 - y3

    return J
end

function dNdxmatrix(interp::LinTrigInterp, local_coords::Point2,
                    vertices::Vertex3, nodes::Vector{FENode2})

        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, vertices, nodes, dN)
        J = inv2x2t!(J)
        A_mul_B!(interp.dNdx, dN, J)
        return interp.dNdx
end


function get_area(::LinTrigInterp, vertices::Vertex3, nodes::Vector{FENode2})

    x1 =  nodes[vertices[1]].coords.x
    x2 =  nodes[vertices[2]].coords.x
    x3 =  nodes[vertices[3]].coords.x

    y1 =  nodes[vertices[1]].coords.y
    y2 =  nodes[vertices[2]].coords.y
    y3 =  nodes[vertices[3]].coords.y

    area = 0.5 * (x1*(y2 - y3) + x2*(-y1 + y3) + x3*(y1 - y2))
    return area
end

function mass_matrix(interp::LinTrigInterp, vertices::Vertex3, nodes::Vector{FENode2})
    fill!(interp.M, 0.0)
    area = get_area(interp, vertices, nodes)
    scale!(interp.M, interp.ref_M, area)
    return interp.M
end

end # module