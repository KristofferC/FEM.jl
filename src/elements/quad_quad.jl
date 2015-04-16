type QuadQuadStorage <: ElemStorage
    B::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f_int::Vector{Float64}
end

QuadQuadStorage() = QuadQuadStorage(zeros(4, 16), zeros(4,16),
                                  zeros(16,16), zeros(4), zeros(16))


type QuadQuad{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex6
    gps::Vector{GaussPoint2}
    n::Int
    interp::QuadQuadInterp
    storage::QuadQuadStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end
gausspoints(elem::QuadQuad) = elem.gps

# Constructor
function QuadQuad{T <: AbstractMaterialStatus}(vertices::Vertex6, n, interp::QuadQuadInterp,
                 storage::QuadQuadStorage, gps::Vector{GaussPoint2}, matstat::T)
    matstats = T[]
    temp_matstats = T[]
    for i in 1:length(gps)
        push!(matstats, copy(matstat))
        push!(temp_matstats, copy(matstat))
    end
    QuadQuad(vertices, gps, n, interp, storage, matstats, temp_matstats)
end


get_geoelem(ele::QuadQuad) = GeoQTrig(ele.n, ele.vertices)
get_geotype(::QuadQuad) = GeoQTrig

createstorage(::Type{QuadQuad}) = QuadQuadStorage()
createinterp(::Type{QuadQuad}) = QuadQuadInterp()

function creategps(::Type{QuadQuad})
end

@inline function get_ndofs(::QuadQuad)
    return 12
end

# Avoids allocation (remove this?)
const QUADTRIG_DOFTYPES = [Du, Dv]
function doftypes(::QuadQuad, ::Int)
    return return QUADTRIG_DOFTYPES
end

function Bmatrix(elem::QuadQuad, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)
    B = elem.storage.B
    for i in 1:6
        B[1, 2*i - 1] = dNdx[i, 1]
        B[2, 2*i]     = dNdx[i, 2]

        B[4, 2*i - 1] = dNdx[i, 2]
        B[4, 2*i]     = dNdx[i, 1]
    end
    return B
end


function weight(elem::QuadQuad, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end

# Get the stress/strain in gausspoint i
get_field(elem::QuadQuad, ::Type{Stress}, i::Int) = elem.matstats[i].stress
get_field(elem::QuadQuad, ::Type{Strain}, i::Int) = elem.matstats[i].strain


