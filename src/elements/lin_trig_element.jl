immutable LinTrig <: Element
    vertices::Vector{Int}
    gps::Vector{GaussPoint}
    n_dofs::Int
    n::Int
    interp::LinTrigInterp
end

function LinTrig(vertices::Vector{Int}, n)

    n_dofs = 6
    gps = [GaussPoint([1/3; 1/3], 0.5)]
    interp = LinTrigInterp()

    LinTrig(vertices, gps, n_dofs, n, interp)
end

function doftypes(elem::LinTrig, vertex::Int)
    return [Du(), Dv()]
end


function Bmatrix(elem::LinTrig, gp::GaussPoint, nodes::Vector{Node})
    dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)
    B = zeros(4, 6)

    B[1, 1] = dNdx[1, 1]
    B[1, 3] = dNdx[2, 1]
    B[1, 5] = dNdx[3, 1]

    B[2, 2] = dNdx[1, 2]
    B[2, 4] = dNdx[2, 2]
    B[2, 6] = dNdx[3, 2]

    B[4, 1] = dNdx[1, 2]
    B[4, 2] = dNdx[1, 1]
    B[4, 3] = dNdx[2, 2]
    B[4, 4] = dNdx[2, 1]
    B[4, 5] = dNdx[3, 2]
    B[4, 6] = dNdx[3, 1]

    return B
end


