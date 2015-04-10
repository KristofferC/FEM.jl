type LinQuadStorage <: ElemStorage
    B::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f_int::Vector{Float64}
end

LinQuadStorage() = LinQuadStorage(zeros(4, 8), zeros(4,8),
                                  zeros(8,8), zeros(4), zeros(8))


type LinQuad{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex4
    gps::Vector{GaussPoint2}
    n::Int
    interp::LinQuadInterp
    storage::LinQuadStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end

gausspoints(elem::LinQuad) = elem.gps

# Constructors
function LinQuad{T <: AbstractMaterialStatus}(vertices::Vertex4, n, interp::LinQuadInterp,
                 storage::LinQuadStorage, gps::Vector{GaussPoint2}, matstat::T)
    matstats = T[]
    temp_matstats = T[]
    for i in 1:length(gps)
        push!(matstats, copy(matstat))
        push!(temp_matstats, copy(matstat))
    end
    LinQuad(vertices, gps, n, interp, storage, matstats, temp_matstats)
end

function LinQuad(v::Vector{Int}, n, interp::LinQuadInterp,
                 storage::LinQuadStorage, gps::Vector{GaussPoint2})
    LinQuad(Vertex4(v[1], v[2], v[3], v[4]), gps, n, interp, storage)
end


get_geoelem(ele::LinQuad) = GeoQuad(ele.n, ele.vertices)
get_geotype(::LinQuad) = GeoQuad

createstorage(::Type{LinQuad}) = LinQuadStorage()
createinterp(::Type{LinQuad}) = LinQuadInterp()

function creategps(::Type{LinQuad})
    p = 1 / 6
    [GaussPoint2(Point2(p, p), 1.0);
     GaussPoint2(Point2( p, -p), 1.0);
     GaussPoint2(Point2( p,  p), 1.0);
     GaussPoint2(Point2(-p,  p), 1.0)]
end

@inline function get_ndofs(::LinQuad)
    return 8
end

# Avoids allocation (remove?)
const LINQUAD_DOFTYPES = [Du, Dv]
function doftypes(::LinQuad, ::Int)
    return return LINQUAD_DOFTYPES
end

function Bmatrix(elem::LinQuad, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)
    B = elem.storage.B
    for i in 1:4
        B[1, 2*i - 1] = dNdx[i, 1]
        B[2, 2*i]     = dNdx[i, 2]

        B[4, 2*i - 1] = dNdx[i, 2]
        B[4, 2*i]     = dNdx[i, 1]
    end
    return B
end


function weight(elem::LinQuad, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end


# Get the stress in gausspoint i
get_field(elem::LinQuad, ::Type{Stress}, i::Int) = elem.matstats[i].stress
get_field(elem::LinQuad, ::Type{Strain}, i::Int) = elem.matstats[i].strain