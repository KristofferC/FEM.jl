type LinTrigUPStorage <: ElemStorage
    # Intermittient matrices
    B::Matrix{Float64}
    Bdiv::Vector{Float64}
    DB::Matrix{Float64}
    KE_uu::Matrix{Float64}
    KE_pp::Matrix{Float64}
    KE_up::Matrix{Float64}
    Ke::Matrix{Float64}
    #
    ɛ::Vector{Float64}
    f_int::Vector{Float64}
end

LinTrigUPStorage() = LinTrigStorage(zeros(4,6), zeros(6), zeros(4,6),
                                  zeros(6,6), zeros(3,3), zeros(6,3),
                                  zeros(9,9), zeros(4), zeros(9))

immutable LinTrigUP <: AbstractFElement
    vertices::Vertex3
    gps::Vector{GaussPoint2}
    n::Int
    interp::LinTrigInterp
    storage::LinTrigUPStorage
end

get_geoelem(ele::LinTrigUP) = GeoTrig(ele.n, ele.vertices)
createstorage(::Type{LinTrigUP}) = LinTrigUPStorage()
createinterp(::Type{LinTrigUP}) = LinTrigInterp()
creategps(::Type{LinTrigUP}) = [GaussPoint2(Point2(1/6, 1/6), 1/6),
                              GaussPoint2(Point2(2/3, 1/6), 1/6),
                              GaussPoint2(Point2(1/6, 2/3), 1/6)]


@inline function get_ndofs(::LinTrigUP)
    return 9
end

function LinTrigUP(vertices::Vertex3, n, interp::LinTrigInterp,
                 lts::LinTrigStorage, gps::Vector{GaussPoint2})
    LinTrigUP(vertices, gps, n, interp, lts)
end

function LinTrigUP(v::Vector{Int}, n, interp::LinTrigInterp,
                 lts::LinTrigStorage, gps::Vector{GaussPoint2})
    LinTrigUP(Vertex3(v[1], v[2], v[3]), gps, n, interp, lts)
end


function doftypes(::LinTrigUP, ::Int)
    return [Du, Dv, Pressure]
end

function stiffness{P <: AbstractMaterial}(elem::LinTrigUP,
                                          nodes::Vector{FENode2},
                                          material::P)
    Be = Bmatrix(elem, elem.gps[1], nodes) # Constant in mat, dummy gp
    Bdiv = Bdiv(elem, elem.gps[1], nodes) # Constant in mat, dummy gp
    De = stiffness(material, gps[1])
    for gp in elem.gps

        # Get deviatoric part

        A_mul_B!(elem.storage.DeBe, De, Be)
        At_mul_B!(elem.storage.Ke, Be, elem.storage.DeBe)

        dV = weight(elem, gp, nodes)
        scale!(elem.storage.Ke, dV)
    end
    return elem.storage.Ke
end

function KE_pp(elem::LinTrigUP, gp::GaussPoint2, nodes::Vector{FENode2})
     N = Nvec(

function Bdiv(elem::LinTrigUP, gp::GaussPoint2, nodes::Vector{FENode2})
     dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)
     Bdiv = elem.storage.Bdiv
     Bdiv[1] = dNdx[1,1]
     Bdiv[2] = dNdx[1,2]
     Bdiv[3] = dNdx[2,1]
     Bdiv[4] = dNdx[2,2]
     Bdiv[5] = dNdx[3,1]
     Bdiv[6] = dNdx[3,2]

function Bmatrix(elem::LinTrigUP, gp::GaussPoint2, nodes::Vector{FENode2})
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

