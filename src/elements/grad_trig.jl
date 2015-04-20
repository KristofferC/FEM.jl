const NSLIP = 2

function get_u_dof_idxs(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, nunodes * nuvars)
    count = 0;
    for i = 1:nunodes
        for j = 1:nuvars
            count++;
            if count < nuvars * ngradnodes + 1
                idxs[count] =  (i - 1) * (ngradvars + nuvars) + j
            else
                idxs[count] = count + ngradvars * ngradnodes
            end
        end
    end
    return idxs
end

function get_grad_dof_idxs(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, nunodes * nuvars)
    count = 0;
    for i = 1:ngradnodes;
        for j = 1:ngradvars
            count += 1
            idxs[count] = i * nuvars + (i - 1) * ngradvars + j;
        end
    end
    return idxs
end

function get_grad_idxs_plane(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int, alpha::Int)
    const NDIM = 2
    idxs = Array(Int, NDIM * ngradnodes)
    count = 0;

    for k=1:ngradnodes
        for j=1:ndim
            count += 1
            idxs[count] = (alpha - 1) * NDIM + (k - 1) * (NDIM * NSLIP) + j;
        end
    end
    return idxs
end


type GradTrigStorage <: ElemStorage
    B::Matrix{Float64}
    Bdiv::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f::Vector{Float64}
    f_u::Vector{Float64}
    f_grad::Vector{Float64}
end

function GradTrigStorage()
   B = zeros(4, 12)
   Bdiv = zeros(6)
   DeBe = zeros(4,12)
   Ke = zeros(12 + 2*3*NSLIP, 12 + 2*3*NSLIP)
   ɛ = zeros(4)
   f_u = zeros(12)
   f_grad = zeros(2 * 3 * NSLIP)
   f = zeros(12 + 2*3*NSLIP)
   GradTrigStorage(B, Bdiv, DeBe, Ke, ɛ, f_u, f_grad, f)
end


type GradTrig{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex6
    gps::Vector{GaussPoint2}
    n::Int
    interp_u::QuadTrigInterp
    interp_grad::LinTrigInterp
    storage::GradTrigStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end
gausspoints(elem::GradTrig) = elem.gps

# Constructor
function GradTrig{T <: AbstractMaterialStatus}(vertices::Vertex6, n, interp::GradTrigInterp,
                 storage::GradTrigStorage, gps::Vector{GaussPoint2}, matstat::T)
    matstats = T[]
    temp_matstats = T[]
    for i in 1:length(gps)
        push!(matstats, copy(matstat))
        push!(temp_matstats, copy(matstat))
    end
    GradTrig(vertices, gps, n, interp, storage, matstats, temp_matstats)
end

get_ndofs(::GradTrig) = 12 + 2*3*NSLIP
get_geoelem(ele::GradTrig) = GeoQTrig(ele.n, ele.vertices)
get_geotype(::GradTrig) = GeoQTrig

createstorage(::Type{GradTrig}) = GradTrigStorage()
createinterp(::Type{GradTrig}) = GradTrigInterp()

function creategps(::Type{GradTrig})
    p1 = 1/3
    p2 = 0.2
    p3 = 0.6
    [GaussPoint2(Point2(p1, p1), -0.281250000000000);
     GaussPoint2(Point2(p2, p3), 0.260416666666667);
     GaussPoint2(Point2(p2, p2), 0.260416666666667)
     GaussPoint2(Point2(p3, p2), 0.260416666666667)]
end

@inline function get_ndofs(::GradTrig)
    return 12 + 2*3*NSLIP
end

function doftypes(::GradTrig, v::Int)
    if v <= 3
        return [Du, Dv, Gu1, Gv1, Gu2, Gv2] # TODO: Fix from arbitrary NSLIP
    else
        return [Du, Dv]
end


function stiffness{P <: AbstractMaterial}(elem::GradTrig,
                                          nodes::Vector{FENode2},
                                          material::P)
    const H = 10e-7
    fill!(elem.storage.Ke, 0.0)
    f = intf(elem, mat, nodes)
    for node in nodes
        for dof in node.dofs
            if !dof.active
                continue
            end
            dof.value += H
            f_pert = intf(elem, mat, nodes)
            fdiff = (f_pert)

        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        A_mul_B!(elem.storage.DeBe, De, Be) # DeBe = De * Be
        # Ke += B' * DeBe * dV
        BLAS.gemm!('T', 'N' ,dV, Be, elem.storage.DeBe, 1.0, elem.storage.Ke)
    end
    return elem.storage.Ke
end

function intf_u{P <: AbstractMaterial}(elem::GradTrig, mat::P, nodes::Vector{FENode2})
    u = get_field(elem, nodes)
    u_f = u[get_u_dof_idxs....]
    ɛ = elem.storage.ɛ
    fill!(elem.storage.f_u, 0.0)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(ɛ, B, u)
        fill_from_start!(elem.temp_matstats[i].strain, ɛ)

        F_grad_mekh = [ɛ[1] + 1.0, ɛ[2] + 1.0, 1.0,
                       ɛ[4], 0.0, 0.0, 0.0, ɛ[4], 0.0]

        σ = stress(mat, F_grad_mekh, gp)
        fill_from_start!(elem.temp_matstats[i].stress, σ)

        dV = weight(elem, gp, nodes)

        # f_u += B' * σ * dV
        BLAS.gemv!('T', dV, B, σ, 1.0, elem.storage.f_u)
    end
    return elem.storage.f_u
end


function intf_grad{P <: AbstractMaterial}(elem::GradTrig, mat::P, nodes::Vector{FENode2})
    fill!(f_grad, 0.0)
    M = mass_matrix(elem.interp_grad, nodes)
    for slip in 1:NSLIP
        fe += M * field[u_grad_alpha]
    end


    for i in length(gausspoints)

    # Dummy gp, constant
    B = Bdiv(elem, elem.gps[1], nodes)
    dN = dNmatrix(elem.interp_u, gp.local_coords)
    J = Jmatrix(elem.interp_u, gp.local_coords, elem.vertices, nodes, dN)
    for slip in 1:NSLIP
        fe += k_alpha * B * A0 * abs(det2x2(J))
    end
    end


end


function intf{P <: AbstractMaterial}(elem::GradTrig, mat::P, nodes::Vector{FENode2})
    f_u = intf_u(elem, mat, nodes)
    f_grad = intf_grad(elem, matm nodes)
    f[] = f_u
    f[] = f_grad
    return f
end



function Bdiv(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp_grad, gp.local_coords, elem.vertices, nodes)
    Bdiv = elem.storage.Bdiv
    Bdiv[1] = dNdx[1, 1]
    Bdiv[2] = dNdx[1, 2]
    Bdiv[3] = dNdx[2, 1]
    Bdiv[4] = dNdx[2, 2]
    Bdiv[5] = dNdx[3, 1]
    Bdiv[6] = dNdx[3, 2]
    return Bdiv
end


function Bmatrix(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp_u, gp.local_coords, elem.vertices, nodes)
    B = elem.storage.B
    for i in 1:6
        B[1, 2*i - 1] = dNdx[i, 1]
        B[2, 2*i]     = dNdx[i, 2]

        B[4, 2*i - 1] = dNdx[i, 2]
        B[4, 2*i]     = dNdx[i, 1]
    end
    return B
end


function weight(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end

# Get the stress/strain in gausspoint i
get_field(elem::GradTrig, ::Type{Stress}, i::Int) = elem.matstats[i].stress
get_field(elem::GradTrig, ::Type{Strain}, i::Int) = elem.matstats[i].strain


