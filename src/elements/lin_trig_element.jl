immutable LinTrigElement <: Element
    vertices = Array{Int}
    gps = Array{GaussPoint}
    n_gp::Int
    n_dofs::Int
    N::Vector{Float64}
    dNdx::Vector{Float64}
    J::Matrix{Float64}
end

function LinTrigElement(vertices)
    n_gp = 1



Bmatrix(elem::LinTrig, gp::GaussPoint, nodes):
    dNdx = eval_dNdx(gp, elem.vertices, nodes)

    B[0, 0] = dNdx[0, 0]
    B[0, 2] = dNdx[1, 0]
    B[0, 4] = dNdx[2, 0]

    B[1, 1] = dNdx[0, 1]
    B[1, 3] = dNdx[1, 1]
    B[1, 5] = dNdx[2, 1]

    B[3, 0] = dNdx[0, 1]
    B[3, 1] = dNdx[0, 0]
    B[3, 2] = dNdx[1, 1]
    B[3, 3] = dNdx[1, 0]
    B[3, 4] = dNdx[2, 1]
    B[3, 5] = dNdx[2, 0]

    return B

