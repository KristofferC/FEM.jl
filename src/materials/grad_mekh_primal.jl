
module GradMekhModPrimal

using Devectorize
import Base.copy

# Function sto extend
import FEM: create_matstat, stiffness, stress, get_kalpha, get_kalphas

# Types needed
import FEM: AbstractMaterialStatus, AbstractMaterial, GaussPoint2

export GradMekhPrimalMS, GradMekh

const NSLIP = 2

type GradMekhPrimalMS <: AbstractMaterialStatus
    n_γ::Vector{Float64}
    strain::Vector{Float64}
    stress::Vector{Float64}
end

const NSLIP = 2

function GradMekhPrimalMS()
  GradMekhPrimalMS(zeros(9), zeros(NSLIP), zeros(NSLIP), 0.1, zeros(6), zeros(6))
end

copy(matstat::GradMekhPrimalMS) = GradMekhPrimalMS(copy(matstat.n_ε_p),  copy(matstat.n_λ),
                                                   copy(matstat.n_k), matstat.n_τ,
                                                   copy(matstat.strain), copy(matstat.stress))

get_kalphas(ms::GradMekhPrimalMS) = ms.n_k
get_kalpha(ms::GradMekhPrimalMS, i::Int) = ms.n_k[i]

immutable GradMekh <: AbstractMaterial
    E::Float64
    ν::Float64
    n::Float64
    l::Float64
    Hg::Float64
    Hl::Float64
    m::Float64
    σy::Float64
    tstar::Float64
    angles::Vector{Float64}
    s_x_m::Vector{Vector{Float64}}
    NSLIP::Int
end


function GradMekh(E, ν, m, l, Hg, Hl, m,
                  σy, tstar, angles, NSLIP)

    if length(angles) != NSLIP
        error("Need one angle per slip system")
    end

    s_x_m = Vector{Vector{Float64}}()
    @devec angles[:] = angles[:] .* pi ./ 180.0
    for α = 1:NSLIP
        t = angles[α]
        s = [cos(t), sin(t), 0.0]
        mm = [cos(t + pi/2), sin(t + pi/2), 0.0]
        sxm_mat = s * mm'
        push!(s_x_m, M_2_V9(sxm_mat))
    end

    GradMekh(E, ν, m, l, Hg, Hl, m, σy, tstar, angles, s_x_m, NSLIP)
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

inv_V9(A) = inv_V9!(zeros(9), A)
@inbounds function inv_V9!(inv_A, A)
    det_A = (-A[7]*A[2]*A[6]+A[4]*A[5]*A[6]+
              A[7]*A[8]*A[9]-A[1]*A[5]*A[9]-
              A[4]*A[8]*A[3]+A[1]*A[2]*A[3])

    if abs(det_A) < eps(Float64)
        error("WARNING! Matrix singular to working precision")
    end

    inv_A[1]=(-A[5]*A[9]+A[2]*A[3])/det_A
    inv_A[4]=( A[7]*A[9]-A[4]*A[3])/det_A
    inv_A[7]=( A[5]*A[4]-A[2]*A[7])/det_A
    inv_A[8]=( A[5]*A[6]-A[8]*A[3])/det_A
    inv_A[2]=(-A[7]*A[6]+A[1]*A[3])/det_A

    inv_A[5]=(-A[5]*A[1]+A[8]*A[7])/det_A
    inv_A[6]=(-A[2]*A[6]+A[8]*A[9])/det_A
    inv_A[9]=( A[4]*A[6]-A[1]*A[9])/det_A
    inv_A[3]=( A[2]*A[1]-A[8]*A[4])/det_A
    return inv_A
end

@inbounds function trans_V9(b)
    a = zeros(9)
    a[1:3] = b[1:3]
    a[[4,5,6,7,8,9]] = b[[8,9,7,6,4,5]]
    return a
end


const STRESS_BUFFER = zeros(9)
const I = Float64[1,1,1,0,0,0,0,0,0]
function stress(mat::GradMekh, matstat::GradMekhPrimalMS, temp_matstat::GradMekhPrimalMS,
                τ, γ)

    ∆t = 2.000000000000000e-003 #TODO: Fix

    λ = -k
    ∆λ = λ - matstat.n_λ

    E = mat.E
    ν = mat.ν
    s_x_m = mat.s_x_m
    tstar = mat.tstar
    m = mat.m

    for α = 1:NSLIP
        ∆γ = γ[alpha] - matstat.n_γ[alpha]
        f = (γ[alpha] - matstat.n_γ[alpha])/dt * tstar
        if f <= 0.0
            κ[α] = 0.0
        else
            κ[α] = abs(τ[α]) - matstat.σy - f^(1/m)
        end
    end

    @devec temp_matstat.n_γ[:] = γ

end # module
