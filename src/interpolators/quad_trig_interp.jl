module QuadTrigInterpMod

import FEM: AbstractInterpolator, Point2, FENode2, Vertex6, inv2x2t!, det2x2,
       dNmatrix, Jmatrix, dNdxmatrix, Nvec, get_area

export QuadTrigInterp


immutable QuadTrigInterp <: AbstractInterpolator
    N::Vector{Float64}
    dN::Matrix{Float64}
    dNdx::Matrix{Float64}
    J::Matrix{Float64}
end

function QuadTrigInterp()
    N = zeros(6)
    dN = zeros(6,2)
    dNdx = Array(Float64, 6, 2)
    J = Array(Float64, 2, 2)
    QuadTrigInterp(N, dN, dNdx, J)
end
get_allocated_N(i::QuadTrigInterp) = i.N
get_allocated_dN(i::QuadTrigInterp) = i.dN
get_allocated_dNdx(i::QuadTrigInterp) = i.dNdx
get_allocated_J(i::QuadTrigInterp) = i.J


# Shape functions in local coords
# http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch16.d/IFEM.Ch16.pdf
function Nvec(interp::QuadTrigInterp, loc_coords::Point2)

    ξ = loc_coords[1]
    η = loc_coords[2]
    γ = 1 - ξ - η

    N = get_allocated_N(interp)

    N[1] = ξ * (2ξ - 1)
    N[2] = η * (2η - 1)
    N[3] = γ * (2γ - 1)
    N[4] = 4ξ * η
    N[5] = 4η * γ
    N[6] = 4ξ * γ

    return N
end


function dNmatrix(interp::QuadTrigInterp, loc_coords::Point2)
    ξ = loc_coords[1]
    η = loc_coords[2]
    γ = 1 - ξ - η


    dN = get_allocated_dN(interp)

    dN[1, 1] = 4ξ - 1
    dN[1, 2] = 0

    dN[2, 1] = 0
    dN[2, 2] = 4η - 1

    dN[3, 1] = -4γ + 1
    dN[3, 2] = -4γ + 1

    dN[4, 1] = 4η
    dN[4, 2] = 4ξ

    dN[5, 1] = -4η
    dN[5, 2] = 4(γ - η)

    dN[6, 1] = 4(γ - ξ)
    dN[6, 2] = -4ξ

    return dN
end

function Jmatrix(interp::QuadTrigInterp,
                 vertices::Vertex6, nodes::Vector{FENode2},
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

function dNdxmatrix(interp::QuadTrigInterp, local_coords::Point2,
                    vertices::Vertex6, nodes::Vector{FENode2})
        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, vertices, nodes, dN)
        A_mul_B!(interp.dNdx, dN, inv2x2t!(J))
        return interp.dNdx
end

end # module