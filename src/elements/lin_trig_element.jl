

type LinTrigStorage <: ElemStorage
    B::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f_int::Vector{Float64}
end

LinTrigStorage() = LinTrigStorage(zeros(4, 6), zeros(4,6), zeros(6,6), zeros(4), zeros(6))

immutable LinTrig <: AbstractFElement
    vertices::Vertex3
    gps::Vector{GaussPoint2}
    n::Int
    interp::LinTrigInterp
    storage::LinTrigStorage
end

# Constructors
function LinTrig(vertices::Vertex3, n, interp::LinTrigInterp,
                 lts::LinTrigStorage, gps::Vector{GaussPoint2})

    LinTrig(vertices, gps, n, interp, lts)
end

function LinTrig(v::Vector{Int}, n, interp::LinTrigInterp,
                 lts::LinTrigStorage, gps::Vector{GaussPoint2})
    LinTrig(Vertex3(v[1], v[2], v[3]), gps, n, interp, lts)
end


get_geoelem(ele::LinTrig) = GeoTrig(ele.n, ele.vertices)
get_storage(::Type{LinTrig}) = LinTrigStorage()
get_interp(::Type{LinTrig}) = LinTrigInterp()
get_gps(::Type{LinTrig}) = [GaussPoint2(Point2(1/3, 1/3), 0.5)]


@inline function get_ndofs(::LinTrig)
    return 6
end


function doftypes(::LinTrig, ::Int)
    return [Du, Dv]
end


function Bmatrix(elem::LinTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)

    B = elem.storage.B

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

function weight(elem::LinTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end
