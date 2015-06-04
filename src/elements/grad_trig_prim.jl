# Load the Interpolator in the FEM module
FEM.lintriginterpmod()
using FEM.LinTrigInterpMod

FEM.quadtriginterpmod()
using FEM.QuadTrigInterpMod

module GradTrigPrimalMod

using Devectorize
using FEM
using FEM.LinTrigInterpMod
using FEM.QuadTrigInterpMod
import Base.AddFun

import FEM: SubFun
import FEM: createinterp, creategps, createstorage, get_field,
            Bmatrix, doftypes, get_ndofs, get_geotype, get_ref_area, get_kalpha, get_kalphas, stiffness, mod_value

import FEM: dNdxmatrix, intf, assemble!, fill_from_start!, stress, weight, mass_matrix, get_area, DofVals
import FEM: AbstractMaterialStatus, AbstractElemStorage, AbstractFElement, AbstractMaterial, FENode2
import FEM: Vertex3, Vertex6, Point2, GaussPoint2, Du, Dv, K1, K2, GeoQTrig, InvFp, KAlpha, Nvec

export GradTrig

const NSLIP = 2


function get_u_dof_idxs(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, nunodes * nuvars)
    count = 1
    for i = 1:nunodes
        for j = 1:nuvars
            if count < nuvars * ngradnodes + 1
                idxs[count] =  (i - 1) * (ngradvars + nuvars) + j
            else
                idxs[count] = count + ngradvars * ngradnodes
            end
            count += 1
        end
    end
    return idxs
end

function get_κ_α(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, ngradnodes * ngradvars)
    count = 1
    for i = 1:ngradnodes
        for j = 1:ngradvars
            idxs[count] = i * nuvars + (i - 1) * ngradvars + j
            count += 1
        end
    end
    return idxs
end

function get_κ(ngradnodes::Int, alpha::Int)
    ngradvars = 1
    idxs = Array(Int, ngradvars * ngradnodes)
    count = 1
    for i=1:ngradnodes
        for j=1:ngradvars
            idxs[count] = (alpha - 1) * ngradvars + (i - 1) * (ngradvars * NSLIP) + j;
            count += 1
        end
    end
    return idxs
end

get_κ_α
type GradTrigStorage <: AbstractElemStorage
    B::Matrix{Float64}
    Bdiv::Vector{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    u_field::Vector{Float64}
    u_u::Vector{Float64}
    u_grad::Vector{Float64}
    u_grad_plane::Vector{Float64}
    f_int::Vector{Float64}
    f_u::Vector{Float64}
    f_grad::Vector{Float64}
    f_grad_plane::Vector{Float64}
    dofs_idx_u::Vector{Int}
    dofs_idx_grad::Vector{Int}
    kappas::Vector{Float64}
    dofs_slip_plane::Vector{Vector{Int}}
end

function GradTrigStorage()

    B = zeros(4, 12)
    Bdiv = zeros(6)
    DeBe = zeros(4,12)
    Ke = zeros(12 + 3*NSLIP, 12 + 3*NSLIP)
    ɛ = zeros(4)
    u_field = zeros(12 + 3*NSLIP)
    u_u = zeros(12)
    u_grad = zeros(3 * NSLIP)
    u_grad_plane = zeros(3)
    f = zeros(12 + 3*NSLIP)
    f_u = zeros(12)
    f_grad = zeros(3 * NSLIP)
    f_grad_plane = zeros(3)
    dofs_idx_u = get_u_dof_idxs(6, 2, 3, NSLIP)
    dofs_idx_grad = get_κ_α(6, 2, 3, NSLIP)
    kappas = zeros(NSLIP)

    dofs_slip_plane = Array(Vector{Int}, NSLIP)
    for i = 1:NSLIP
        dofs_slip_plane[i] = get_κ(3, i)
    end

   GradTrigStorage(B, Bdiv, DeBe, Ke, ɛ, u_field, u_u, u_grad, u_grad_plane,
                   f, f_u, f_grad, f_grad_plane,
                   dofs_idx_u, dofs_idx_grad, kappas, dofs_slip_plane)
end


type GradTrig{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex6
    gps::Vector{GaussPoint2}
    n::Int
    interp::QuadTrigInterp
    interp_grad::LinTrigInterp
    storage::GradTrigStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end
gausspoints(elem::GradTrig) = elem.gps
get_dofs_slipplane(elem, i::Int) = elem.storage.dofs_slip_plane[i]

# Constructor
function GradTrig{T <: AbstractMaterialStatus}(vertices::Vertex6, n, interp_u::QuadTrigInterp,
                                                interp_grad::LinTrigInterp, storage::GradTrigStorage,
                                                gps::Vector{GaussPoint2}, matstat::T)
    matstats = T[]
    temp_matstats = T[]
    for i in 1:length(gps)
        push!(matstats, copy(matstat))
        push!(temp_matstats, copy(matstat))
    end

    GradTrig(vertices, gps, n, interp_u, interp_grad, storage,
             matstats, temp_matstats)
end

get_ndofs(::GradTrig) = 12 + 2*3*NSLIP
get_geoelem(ele::GradTrig) = GeoQTrig(ele.n, ele.vertices)
get_geotype(::GradTrig) = GeoQTrig
get_ref_area(::GradTrig) = 0.5
createstorage(::Type{GradTrig}) = GradTrigStorage()
createinterp(::Type{GradTrig}) = QuadTrigInterp()

function creategps(::Type{GradTrig})
    p1 = 2/3
    p2 = 1/6
    w = 1/6
    [GaussPoint2(Point2(p1, p2), w);
     GaussPoint2(Point2(p2, p1), w);
     GaussPoint2(Point2(p2, p2), w)]
end

@inline function get_ndofs(::GradTrig)
    return 12 + 2*3*NSLIP
end

function doftypes(::GradTrig, v::Int)
    if v <= 3
        return [Du, Dv, K1, K2] # TODO: Fix from arbitrary NSLIP
    else
        return [Du, Dv]
    end
end


function stiffness(elem::GradTrig,
                  nodes::Vector{FENode2},
                  material::AbstractMaterial,
                  dofvals::DofVals)
    H = 1e-7
    Ke = elem.storage.Ke
    fill!(Ke, 0.0)
    f = intf(elem, material, nodes, dofvals)

    # Need to copy because intf returns a reference
    # which will be overwritten on subseq call to intf
    ff = copy(f)
    p = 1
    col = 1
    for v in elem.vertices
        node = nodes[v]
        for dof in node.dofs
            mod_value(dofvals, dof, AddFun(), H)
            f_pert = intf(elem, material, nodes, dofvals)
            # Numeric derivative
            @inbounds @devec Ke[:, col] = (f_pert .- ff) ./ H
            mod_value(dofvals, dof, SubFun(), H)
            col += 1
        end
    end
    # Reset int variables
    intf(elem, material, nodes, dofvals)
    return Ke
end

M_2_V9(b) = M_2_V9!(zeros(9), b)
@inbounds function M_2_V9!(a, b::Matrix{Float64})
    a[1]=b[1,1]
    a[2]=b[2,2]
    a[3]=b[3,3]
    a[4]=b[1,2]
    a[9]=b[3,2]
    a[5]=b[2,3]
    a[8]=b[2,1]
    a[7]=b[1,3]
    a[6]=b[3,1]
    return a
end


sym_V9(b) = sym_V9!(zeros(9), b)
@inbounds function sym_V9!(a, b)

    a[1:3] =  b[1:3]
    a[4]   = (b[4]+b[8])/2
    a[8]   = (b[4]+b[8])/2
    a[5]   = (b[5]+b[9])/2
    a[9]   = (b[5]+b[9])/2
    a[6]   = (b[6]+b[7])/2
    a[7]   = (b[6]+b[7])/2
    return a
end

V9_d_V9(a, b) = zeros(9, a, b)
@inbounds function V9_d_V9!(c, a, b)
    c[1]=a[1]*b[1] + a[7]* b[6] + a[4]* b[8]
    c[2]=a[2]*b[2] + a[8]* b[4] + a[5]* b[9]
    c[3]=a[3]*b[3] + a[9]* b[5] + a[6]* b[7]
    c[4]=a[4]*b[2] + a[1]* b[4] + a[7]* b[9]
    c[5]=a[5]*b[3] + a[2]* b[5] + a[8]* b[7]
    c[6]=a[6]*b[1] + a[3]* b[6] + a[9]* b[8]
    c[7]=a[7]*b[3] + a[4]* b[5] + a[1]* b[7]
    c[8]=a[8]*b[1] + a[5]* b[6] + a[2]* b[8]
    c[9]=a[9]*b[2] + a[6]* b[4] + a[3]* b[9]
    return c
end

const I = Float64[1,1,1,0,0,0,0,0,0]
const STRESS_INDEX = [1,2,3,4,5,6]

function intf_u(elem::GradTrig, mat::AbstractMaterial, nodes::Vector{FENode2}, dof_vals::DofVals)
    E = mat.E
    ν = mat.ν
    sxm = mat.s_x_m

    u = get_field(elem, nodes, dof_vals)

    assemble!(elem.storage.u_u, u, elem.storage.dofs_idx_u) # 12 dofs
    assemble!(elem.storage.u_grad, u, elem.storage.dofs_idx_grad) # 6 dofs
    uu = elem.storage.u_u

    fill!(elem.storage.f_u, 0.0)
    ɛ = elem.storage.ɛ
    ε_p = zeros(9)
    γ = zeros(NSLIP)
    H = zeros(3,3)
    τ = zeros(NSLIP)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(ɛ, B, elem.storage.u_u) # B = 4 x 12
        fill_from_start!(elem.temp_matstats[i].strain, ɛ)

        dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)
        fill!(H, 0.0)
        H[1,1:2] = uu[1].*dNdx[1,1:2].+uu[3].*dNdx[2,1:2].+
                          uu[5].*dNdx[3,1:2].+uu[7].*dNdx[4,1:2].+
                          uu[9].*dNdx[5,1:2].+uu[11].*dNdx[6,1:2]

        H[2,1:2]=uu[2].*dNdx[1,1:2].+uu[4].*dNdx[2,1:2].+
                        uu[6].*dNdx[3,1:2].+uu[8].*dNdx[4,1:2].+
                        uu[10].*dNdx[5,1:2].+uu[12].*dNdx[6,1:2]


        ε_full = sym_V9(M_2_V9(H))
        N = Nvec(elem.interp_grad, gp.local_coords)
        for α = 1:NSLIP
            dofs_slip_plane = get_dofs_slipplane(elem, α)
            assemble!(elem.storage.u_grad_plane, elem.storage.u_grad, dofs_slip_plane)
            k = dot(N, elem.storage.u_grad_plane)
            γ[α] = -k
            ε_p += γ[α] * sxm[α] * sign(elem.matstats[i].n_τ[α])
        end
        ε_e = (ε_full - ε_p)
        G = E / (2*(1+ν))
        L = E*ν / ((1+ν)*(1-2*ν))
        # Calculate Mandel stress
        σ = L * dot(I, ε_e)*I + 2*G*ε_e
        fill_from_start!(elem.temp_matstats[i].stress, σ[STRESS_INDEX])

        for α = 1:NSLIP
            τ[α] = dot(σ, sxm[α])
        end


        κ = stress(mat, elem.matstats[i], elem.temp_matstats[i], τ, γ)
        fill_from_start!(elem.temp_matstats[i].κ, κ)


        dV = weight(elem, gp, nodes)
        fe = zeros(12)
        @inbounds for i = 1:6
            fe[2*i-1] += (σ[1]*dNdx[i,1]+σ[8]*dNdx[i,2] )
            fe[2*i]   += (σ[4]*dNdx[i,1]+σ[2]*dNdx[i,2] )
        end

        @devec elem.storage.f_u[:] += fe .* dV

        # f_u += B' * σ * dV
        #BLAS.gemv!('T', dV, B, σ[[1,2,3,4]], 1.0, elem.storage.f_u)
    end
    return elem.storage.f_u
end


function intf_grad(elem::GradTrig, mat::AbstractMaterial, nodes::Vector{FENode2}, dof_vals::DofVals)
    fill!(elem.storage.f_grad, 0.0)
    u = get_field(elem, nodes, dof_vals)
    assemble!(elem.storage.u_grad, u, elem.storage.dofs_idx_grad)

    vertslin = Vertex3(elem.vertices[1],elem.vertices[2], elem.vertices[3])
    M = mass_matrix(elem.interp_grad, vertslin, nodes)

    # Dummy gp, constant
    A = get_area(elem.interp_grad, vertslin, nodes)
    N = Nvec(elem.interp_grad, FEM.Point2(1/3,1/3))
    Bt = dNdxmatrix(elem.interp_grad, FEM.Point2(1/3,1/3),
                       vertslin, nodes)

    for α in 1:NSLIP
        dofs_slip_plane = get_dofs_slipplane(elem, α)
        assemble!(elem.storage.u_grad_plane, elem.storage.u_grad, dofs_slip_plane)
        k_α = elem.storage.u_grad_plane

        κ = 0.0
        for j in 1:length(elem.gps)
            κ += elem.temp_matstats[j].κ
        end
        κ /= length(elem.gps)

        Bts_α = Bt * mat.s[α]
        Q = mat.Hg * mat.l^2 * (Bts_α * Bts_α')
        term2 = (M + Q) * k_α

        @devec elem.storage.f_grad[dofs_slip_plane]  += N .* κ[α] .* A + term2
    end

    return elem.storage.f_grad
end


function intf(elem::GradTrig, mat::AbstractMaterial, nodes::Vector{FENode2}, dof_vals::DofVals)
  f_u = intf_u(elem, mat, nodes, dof_vals)
  f_grad = intf_grad(elem, mat, nodes, dof_vals)

  @devec elem.storage.f_int[elem.storage.dofs_idx_u] = f_u
  @devec elem.storage.f_int[elem.storage.dofs_idx_grad] = f_grad
  return elem.storage.f_int
end



function Bdiv(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    vertslin = Vertex3(elem.vertices[1],elem.vertices[2], elem.vertices[3])
    dNdx = dNdxmatrix(elem.interp_grad, gp.local_coords, vertslin, nodes)
    B = elem.storage.Bdiv
    B[1] = dNdx[1, 1]
    B[2] = dNdx[1, 2]
    B[3] = dNdx[2, 1]
    B[4] = dNdx[2, 2]
    B[5] = dNdx[3, 1]
    B[6] = dNdx[3, 2]
    return B
end


function Bmatrix(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
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
get_field(elem::GradTrig, ::Type{InvFp}, i::Int) = (invfp=elem.matstats[i].state[1:9]; invfp[1:3] -= 1.0; invfp)
get_field(elem::GradTrig, ::Type{KAlpha}, i::Int) = get_kalphas(elem.matstats[i])

end # module
