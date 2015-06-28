# Load the Interpolator in the FEM module
FEM.lintriginterpmod()
using FEM.LinTrigInterpMod

module LinTrigMod

using FEM.LinTrigInterpMod

import FEM: createinterp, creategps, createstorage, get_field,
            Bmatrix, doftypes, get_ndofs, get_geotype, get_ref_area, dNdxmatrix

import FEM: AbstractMaterialStatus, AbstractElemStorage, AbstractFElement, FENode2
import FEM: Vertex3, Point2, GaussPoint2, Du, Dv, GeoTrig

export LinTrig

type LinTrigStorage <: AbstractElemStorage
    B::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f_int::Vector{Float64}
    u_field::Vector{Float64}
end

LinTrigStorage() = LinTrigStorage(zeros(4, 6), zeros(4,6), zeros(6,6),
                                  zeros(4), zeros(6), zeros(6))

type LinTrig{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex3
    gps::Vector{GaussPoint2}
    n::Int
    interp::LinTrigInterp
    storage::LinTrigStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end

# Constructors
function LinTrig{T <: AbstractMaterialStatus}(vertices::Vertex3, n, interp::LinTrigInterp,
                 lts::LinTrigStorage, gps::Vector{GaussPoint2}, matstat::T)
    matstats = Array(T, length(gps))
    temp_matstats = Array(T, length(gps))
    for i in 1:length(gps)
        matstats[i] = copy(matstat)
        temp_matstats[i] = copy(matstat)
    end

    LinTrig(vertices, gps, n, interp, lts, matstats, temp_matstats)
end

function LinTrig(v::Vector{Int}, n, interp::LinTrigInterp,
                 lts::LinTrigStorage, gps::Vector{GaussPoint2})
    LinTrig(Vertex3(v[1], v[2], v[3]), gps, n, interp, lts)
end

get_ref_area(::LinTrig) = 0.5
get_geoelem(ele::LinTrig) = GeoTrig(ele.n, ele.vertices)
get_geotype(::LinTrig) = GeoTrig
createstorage(::Type{LinTrig}) = LinTrigStorage()
createinterp(::Type{LinTrig}) = LinTrigInterp()
creategps(::Type{LinTrig}) = [GaussPoint2(Point2(1/3, 1/3), 0.5)]
get_ndofs(::LinTrig) = 6

const dofsa = [Du, Dv]
function doftypes(::LinTrig, ::Int)
    return dofsa
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

end # module
