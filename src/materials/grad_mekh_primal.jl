
module GradMekhPrimalMod

using Devectorize
import Base.copy

# Function sto extend
import FEM: create_matstat, stiffness, stress, get_kalpha, get_kalphas

# Types needed
import FEM: AbstractMaterialStatus, AbstractMaterial, GaussPoint2

export GradMekhPrimalMS, GradMekh, compute_κ

const NSLIP = 2

type GradMekhPrimalMS <: AbstractMaterialStatus
    n_γ::Vector{Float64}
    κ::Vector{Float64}
    τ::Vector{Float64}

    strain::Vector{Float64}
    stress::Vector{Float64}
end

const NSLIP = 2

function GradMekhPrimalMS()
  GradMekhPrimalMS(zeros(NSLIP), zeros(NSLIP), 0.1*ones(NSLIP), zeros(6), zeros(6))
end

copy(matstat::GradMekhPrimalMS) = GradMekhPrimalMS(copy(matstat.n_γ),
                                                   copy(matstat.κ), copy(matstat.τ),
                                                   copy(matstat.strain), copy(matstat.stress))

immutable GradMekh <: AbstractMaterial
    E::Float64
    ν::Float64
    l::Float64
    Hg::Float64
    Hl::Float64
    m::Float64
    σy::Float64
    tstar::Float64
    angles::Vector{Float64}
    s_x_m::Vector{Vector{Float64}}
    s::Vector{Vector{Float64}}
    NSLIP::Int
end


function GradMekh(E, ν, l, Hg, Hl, m,
                  σy, tstar, angles, NSLIP)

    if length(angles) != NSLIP
        error("Need one angle per slip system")
    end

    s_x_m = Vector{Vector{Float64}}()
    sv = Vector{Vector{Float64}}()
    @devec angles[:] = angles[:] .* pi ./ 180.0
    for α = 1:NSLIP
        t = angles[α]
        s = [cos(t), sin(t), 0.0]
        mm = [cos(t + pi/2), sin(t + pi/2), 0.0]
        sxm_mat = s * mm'
        push!(s_x_m, M_2_V9(sxm_mat))
        push!(sv, [cos(t), sin(t)])
    end

    GradMekh(E, ν, l, Hg, Hl, m, σy, tstar, angles, s_x_m, sv, NSLIP)
end

create_matstat(::Type{GradMekh}) = GradMekhPrimalMS()


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

const STRESS_BUFFER = zeros(9)
const I = Float64[1,1,1,0,0,0,0,0,0]
const ε_R = 1e-3
const F_R = 1.0
function stress(mat::GradMekh, matstat::GradMekhPrimalMS, temp_matstat::GradMekhPrimalMS,
                τ, γ)

    ∆t = 2.000000000000000e-003 #TODO: Fix

    tstar = mat.tstar
    m = mat.m

    κ = zeros(NSLIP)
    μ = zeros(NSLIP)
    maxiter = 8



    for α = 1:NSLIP
        μ[α]  = (γ[α] - matstat.n_γ[α])
        niter = 0
        abstol = 1e-5
        κ[α] = resid(mat.σy, μ[α], ∆t, τ[α], mat.m, tstar)
        #=
        while true
            R = resid(mat.σy, ∆γ, ∆t, τ[α], κ[α], mat.m, tstar)
            if abs(R) < abstol*mat.σy
                break
            elseif niter > maxiter
                println(τ)
                println(γ)
                println(matstat)
                println(temp_matstat)
                error("Material failed to converge")
            else
                niter += 1
                H = 1.0
                κ[α] += H
                R_f = resid(mat.σy, ∆γ, ∆t, τ[α], κ[α], mat.m, tstar)
                κ[α] -= H
                J = (R_f - R) / H
                println("num")
                println(J)
                κ[α] = κ[α] - R / J
            end
        end
        =#
    end

    @devec temp_matstat.n_γ[:] = γ
    @devec temp_matstat.τ[:] = τ
    return κ
end

function resid(σy, μ, ∆t, τ, m, tstar)

    F = abs(τ) - σy


    #=
    if F <= 0
        ηR = ε_R + ε_R / pi * 2 * atan(F / F_R)
    else
        ηR = ε_R
    end
    ηf = max(0, F)^m + ηR
    =#

    if F <= 0
        κ = 0
    else
        κ = abs(τ) - σy - max(0,(tstar/∆t * μ))^(1/m)
        #println(κ)
        #println(τ)
        #println((tstar/∆t * μ))
        #println("---")
    end
    return κ
end

end # module
