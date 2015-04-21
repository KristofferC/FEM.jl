immutable QuadQuadInterp <: AbstractInterpolator
    N::Vector{Float64}
    dN::Matrix{Float64}
    dNdx::Matrix{Float64}
    J::Matrix{Float64}
end

function QuadQuadInterp()
    N = zeros(8)
    dN = zeros(8,2)
    dNdx = Array(Float64, 8, 2)
    J = Array(Float64, 2, 2)
    QuadQuadInterp(N, dN, dNdx, J)
end
get_allocated_N(i::QuadQuadInterp) = i.N
get_allocated_dN(i::QuadQuadInterp) = i.dN
get_allocated_dNdx(i::QuadQuadInterp) = i.dNdx
get_allocated_J(i::QuadQuadInterp) = i.J


# Shape functions in local coords
# http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch16.d/IFEM.Ch16.pdf
function Nvec(interp::QuadQuadInterp, loc_coords::Point2)

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = get_allocated_N(interp)

    N[1] = (1 + ξ) * (1 + η) * ( ξ + η - 1) * 0.25
    N[2] = (1 - ξ) * (1 + η) * (-ξ + η - 1) * 0.25
    N[3] = (1 - ξ) * (1 - η) * (-ξ - η - 1) * 0.25
    N[4] = (1 + ξ) * (1 - η) * ( ξ - η - 1) * 0.25
    N[5] = (1 - ξ*ξ) * (1 + η) * 0.5
    N[6] = (1 - η*η) * (1 - ξ) * 0.5
    N[7] = (1 - ξ*ξ) * (1 - η) * 0.5
    N[8] = (1 - η*η) * (1 + ξ) * 0.5

    return N
end


function dNmatrix(interp::QuadQuadInterp, loc_coords::Point2)
    ξ = loc_coords[1]
    η = loc_coords[2]

    dN = get_allocated_dN(interp)

   dN[1,1] =  (1 + η) * ( 2ξ + η) * 0.25
   dN[1,2] =  (1 + ξ) * ( 2η + η) * 0.25

   dN[2,1] = -(1 + η) * (-2ξ + η) * 0.25
   dN[2,2] =  (1 - ξ) * ( 2η + η) * 0.25

   dN[3,1] = -(1 - η) * (-2ξ + η) * 0.25
   dN[3,2] = -(1 - ξ) * (-2η + η) * 0.25

   dN[4,1] =  (1 - η) * ( 2ξ + η) * 0.25
   dN[4,2] = -(1 + ξ) * (-2η + η) * 0.25

   dN[5,1] = - ξ * (1 + η)
   dN[5,2] = (1 - ξ*ξ) * 0.5

   dN[6,1] = -(1 - η*η) * 0.5
   dN[6,2] =  -η * (1 - ξ)

   dN[7,1] = -ξ * (1 - η)
   dN[7,2] = -(1 - ξ*ξ) * 0.5

   dN[8,1] =  (1 - η*η) * 0.5
   dN[8,2] = -η * (1 + ξ)


    return dN
end

function Jmatrix(interp::QuadQuadInterp,
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


function dNdxmatrix(interp::QuadQuadInterp, local_coords::Point2,
                    vertices::Vertex6, nodes::Vector{FENode2})
        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, vertices, nodes, dN)
        A_mul_B!(interp.dNdx, dN, inv2x2t!(J))
        return interp.dNdx
end



