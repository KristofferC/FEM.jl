FEM.linquadinterpmod()
using FEM.LinQuadInterpMod
println("Rand this shit")

module QuadTrigMod

using FEM
FEM.quadtriginterpmod()
using FEM.QuadTrigInterpMod

import FEM: AbstractMaterialStatus, Vertex6, Point2, GaussPoint2, Du, Dv, createstorage,
       createinterp, creategps, get_field, Bmatrix, doftypes, get_ndofs

export QuadTrig

type QuadTrigStorage <: AbstractElemStorage
    B::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f_int::Vector{Float64}
    u_field::Vector{Float64}
end

QuadTrigStorage() = QuadTrigStorage(zeros(4, 12), zeros(4,12),
                                  zeros(12,12), zeros(4), zeros(12), zeros(12))


type QuadTrig{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex6
    gps::Vector{GaussPoint2}
    n::Int
    interp::QuadTrigInterp
    storage::QuadTrigStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end
gausspoints(elem::QuadTrig) = elem.gps

# Constructor
function QuadTrig{T <: AbstractMaterialStatus}(vertices::Vertex6, n, interp::QuadTrigInterp,
                 storage::QuadTrigStorage, gps::Vector{GaussPoint2}, matstat::T)
    matstats = T[]
    temp_matstats = T[]
    for i in 1:length(gps)
        push!(matstats, copy(matstat))
        push!(temp_matstats, copy(matstat))
    end
    QuadTrig(vertices, gps, n, interp, storage, matstats, temp_matstats)
end

get_ref_area(::QuadTrig) = 0.5
get_geoelem(ele::QuadTrig) = GeoQTrig(ele.n, ele.vertices)
get_geotype(::QuadTrig) = GeoQTrig

createstorage(::Type{QuadTrig}) = QuadTrigStorage()
createinterp(::Type{QuadTrig}) = QuadTrigInterp()

function creategps(::Type{QuadTrig})
    p1 = 1/3
    p2 = 0.2
    p3 = 0.6
    [GaussPoint2(Point2(p1, p1), -0.281250000000000);
     GaussPoint2(Point2(p2, p3), 0.260416666666667);
     GaussPoint2(Point2(p2, p2), 0.260416666666667)
     GaussPoint2(Point2(p3, p2), 0.260416666666667)]
end

@inline function get_ndofs(::QuadTrig)
    return 12
end

# Avoids allocation (remove this?)
const QUADTRIG_DOFTYPES = [Du, Dv]
function doftypes(::QuadTrig, ::Int)
    return return QUADTRIG_DOFTYPES
end

function Bmatrix(elem::QuadTrig, gp::GaussPoint2, nodes::Vector{FENode2})
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



# Get the stress/strain in gausspoint i
get_field(elem::QuadTrig, ::Type{Stress}, i::Int) = elem.matstats[i].stress
get_field(elem::QuadTrig, ::Type{Strain}, i::Int) = elem.matstats[i].strain
function get_field(elem::AbstractFElement, ::Type{VonMises}, i::Int)
    σ = elem.matstats[i].stress
    m = (σ[1] + σ[2] + σ[3]) / 3
    return [sqrt(3/2) * sqrt((σ[1] - m)^2 + (σ[2] - m)^2 + (σ[3] - m)^2 +
                        2(σ[4]*σ[4]))]
end

end # module