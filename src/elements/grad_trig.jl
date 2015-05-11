# Load the Interpolator in the FEM module
FEM.lintriginterpmod()
using FEM.LinTrigInterpMod

FEM.quadtriginterpmod()
using FEM.QuadTrigInterpMod

module GradTrigMod

using Devectorize

using FEM.LinTrigInterpMod
using FEM.QuadTrigInterpMod


import FEM: createinterp, creategps, createstorage, get_field,
            Bmatrix, doftypes, get_ndofs, get_geotype, get_ref_area, get_kalpha, stiffness

import FEM: dNdxmatrix, intf, assemble!, fill_from_start!, stress, weight, mass_matrix, get_area
import FEM: AbstractMaterialStatus, AbstractElemStorage, AbstractFElement, AbstractMaterial, FENode2
import FEM: Vertex3, Vertex6, Point2, GaussPoint2, Du, Dv, Gu1, Gv1, Gu2, Gv2, GeoQTrig, InvFp, KAlpha

export GradTrig

const NSLIP = 2
const STRESS_INDEX = [1,2,3,4,5,6]

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

function get_grad_dof_idxs(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, nunodes * nuvars)
    count = 1
    for i = 1:ngradnodes
        for j = 1:ngradvars
            idxs[count] = i * nuvars + (i - 1) * ngradvars + j
            count += 1
        end
    end
    return idxs
end

function get_grad_idxs_plane(ngradnodes::Int, alpha::Int)
    const NDIM = 2
    idxs = Array(Int, NDIM * ngradnodes)
    count = 1
    for i=1:ngradnodes
        for j=1:NDIM
            idxs[count] = (alpha - 1) * NDIM + (i - 1) * (NDIM * NSLIP) + j;
            count += 1
        end
    end
    return idxs
end


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
    Ke = zeros(12 + 2*3*NSLIP, 12 + 2*3*NSLIP)
    ɛ = zeros(4)
    u_field = zeros(12 + 2*3*NSLIP)
    u_u = zeros(12)
    u_grad = zeros(2 * 3 * NSLIP)
    u_grad_plane = zeros(6)
    f = zeros(12 + 2*3*NSLIP)
    f_u = zeros(12)
    f_grad = zeros(2 * 3 * NSLIP)
    f_grad_plane = zeros(6)
    dofs_idx_u = get_u_dof_idxs(6, 2, 3, 2*NSLIP)
    dofs_idx_grad = get_grad_dof_idxs(6, 2, 3, 2*NSLIP)
    kappas = zeros(NSLIP)

    dofs_slip_plane = Array(Vector{Int}, NSLIP)
    for i = 1:NSLIP
        dofs_slip_plane[i] = get_grad_idxs_plane(3, i)
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

#=
function creategps(::Type{GradTrig})
    p1 = 1/3
    p2 = 0.2
    p3 = 0.6
    [GaussPoint2(Point2(p1, p1), -0.281250000000000);
     GaussPoint2(Point2(p2, p3), 0.260416666666667);
     GaussPoint2(Point2(p2, p2), 0.260416666666667)
     GaussPoint2(Point2(p3, p2), 0.260416666666667)]
end
=#

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
        return [Du, Dv, Gu1, Gv1, Gu2, Gv2] # TODO: Fix from arbitrary NSLIP
    else
        return [Du, Dv]
    end
end


function compute_hardening(elem::GradTrig, nodes::Vector{FENode2})

    NDIM = 2

    l = 1e-2
    Hg = 4e7
    factor = 1

    u = get_field(elem, nodes)
    assemble!(elem.storage.u_grad, u, elem.storage.dofs_idx_grad)
    vertslin = Vertex3(elem.vertices[1],elem.vertices[2], elem.vertices[3])
    # Dummy gp

    dNdx = dNdxmatrix(elem.interp_grad, elem.gps[1].local_coords, vertslin, nodes)
    kappas = elem.storage.kappas
    B = Bdiv(elem, elem.gps[1], nodes)
    for i = 1:NSLIP
        slips = get_dofs_slipplane(elem, i)
        assemble!(elem.storage.u_grad_plane, elem.storage.u_grad, slips)
        kappas[i] = Hg * l * l * factor * dot(B, elem.storage.u_grad_plane)
    end

   # println("kappas $kappas")

    return kappas
end


function stiffness(elem::GradTrig,
                  nodes::Vector{FENode2},
                  material::AbstractMaterial)
   const H = 1e-7
    Ke = elem.storage.Ke
    fill!(Ke, 0.0)
    f = intf(elem, material, nodes)

    # Need to copy because intf returns a reference
    # which will be overwritten on subseq call to intf
    ff = copy(f)
    p = 1
    col = 1
    for v in elem.vertices
        node = nodes[v]
        for dof in node.dofs
            dof.value += H
            f_pert = intf(elem, material, nodes)
            # Numeric derivative
            @inbounds @devec Ke[:, col] = (f_pert .- ff) ./ H
            dof.value -= H
            col += 1
        end
    end
    intf(elem, material, nodes)
    return Ke
end


function intf_u(elem::GradTrig, mat::AbstractMaterial, nodes::Vector{FENode2})
    u = get_field(elem, nodes)

    assemble!(elem.storage.u_u, u, elem.storage.dofs_idx_u) # 12 dofs
    uu = elem.storage.u_u
    #println("u: $uu")
    fill!(elem.storage.f_u, 0.0)
    ɛ = elem.storage.ɛ
    F = zeros(3,3)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(ɛ, B, elem.storage.u_u) # B = 4 x 12
        fill_from_start!(elem.temp_matstats[i].strain, ɛ)

        # TODO:
        kappas = compute_hardening(elem, nodes)

        dNdx = dNdxmatrix(elem.interp, gp.local_coords, elem.vertices, nodes)
        fill!(F, 0.0)
        @devec F[1,1:2] = uu[1].*dNdx[1,1:2].+uu[3].*dNdx[2,1:2].+
                 uu[5].*dNdx[3,1:2].+uu[7].*dNdx[4,1:2].+
                 uu[9].*dNdx[5,1:2].+uu[11].*dNdx[6,1:2]

        @devec F[2,1:2]=uu[2].*dNdx[1,1:2].+uu[4].*dNdx[2,1:2].+
                 uu[6].*dNdx[3,1:2].+uu[8].*dNdx[4,1:2].+
                 uu[10].*dNdx[5,1:2].+uu[12].*dNdx[6,1:2]

        F[1,1] += 1
        F[2,2] += 1
        F[3,3] += 1


       # println(F)


        σ = stress(mat, elem.matstats[i], elem.temp_matstats[i], kappas, F)

        #@debug("σ = $σ")
        fill_from_start!(elem.temp_matstats[i].stress, σ[STRESS_INDEX])

        dV = weight(elem, gp, nodes)
        fe = zeros(12)
        @inbounds for i = 1:6
            fe[2*i-1] += (σ[1]*dNdx[i,1]+σ[8]*dNdx[i,2] )
            fe[2*i]   += (σ[4]*dNdx[i,1]+σ[2]*dNdx[i,2] )
        end
        #println("fe $fe")
        @devec elem.storage.f_u[:] += fe .* dV
       # println("fe_u contrib, $i: $fe")


        # f_u += B' * σ * dV
        #BLAS.gemv!('T', dV, B, σ[[1,2,3,4]], 1.0, elem.storage.f_u)
    end
    return elem.storage.f_u
end


function intf_grad(elem::GradTrig, mat::AbstractMaterial, nodes::Vector{FENode2})
    fill!(elem.storage.f_grad, 0.0)
    u = get_field(elem, nodes)

    assemble!(elem.storage.u_grad, u, elem.storage.dofs_idx_grad)

    vertslin = Vertex3(elem.vertices[1],elem.vertices[2], elem.vertices[3])
    M = mass_matrix(elem.interp_grad, vertslin, nodes)

    # Dummy gp, constant
    B = Bdiv(elem, elem.gps[1], nodes)
    A = get_area(elem.interp_grad, vertslin, nodes)

    buff = zeros(6)
    for i in 1:NSLIP
        dofs_slip_plane = get_dofs_slipplane(elem, i)
        assemble!(elem.storage.u_grad_plane, elem.storage.u_grad, dofs_slip_plane)

        k_alpha_tot = 0.0
        for j in 1:length(elem.gps)
            k_alpha = get_kalpha(elem.temp_matstats[j], i)
            k_alpha_tot += k_alpha
        end
        k_alpha_tot /= length(elem.gps)

        A_mul_B!(buff, M, elem.storage.u_grad_plane)
        @devec elem.storage.f_grad[dofs_slip_plane] += buff

      #  println("k_alph_tot $(k_alpha_tot)")

        # f_grad += k_alpha * B * A
        @devec elem.storage.f_grad[dofs_slip_plane] += k_alpha_tot .* B .* A

      #  println("f_grad: $(elem.storage.f_grad[dofs_slip_plane])")

       #  if abs(k_alpha_tot) > 0
       #     println(elem.storage.f_grad[dofs_slip_plane])
       #     println("UELEM")
       #     println(elem.storage.u_grad)
       # end
    end

    return elem.storage.f_grad
end


function intf(elem::GradTrig, mat::AbstractMaterial, nodes::Vector{FENode2})
  f_u = intf_u(elem, mat, nodes)
  f_grad = intf_grad(elem, mat, nodes)

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
get_field(elem::GradTrig, ::Type{KAlpha}, i::Int) = elem.matstats[i].state[10:11]

end # module