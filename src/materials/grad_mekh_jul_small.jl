
module GradMekhModJlSmall

using Devectorize
import Base.copy

# Function sto extend
import FEM: create_matstat, stiffness, stress, get_kalpha

# Types needed
import FEM: AbstractMaterialStatus, AbstractMaterial, GaussPoint2

export GradMekhMS, GradMekh

const NSLIP = 2

type GradMekhMS <: AbstractMaterialStatus
    n_ε_p::Vector{Float64}
    n_k::Vector{Float64}
    n_∆λ::Vector{Float64}
    strain::Vector{Float64}
    stress::Vector{Float64}
end

const NSLIP = 2

function GradMekhMS()
  GradMekhMS(zeros(9), zeros(NSLIP), zeros(NSLIP), zeros(6), zeros(6))
end

copy(matstat::GradMekhMS) = GradMekhMS(copy(matstat.n_ε_p), copy(matstat.n_k), copy(matstat.n_∆λ),
                                       copy(matstat.strain), copy(matstat.stress))

get_kalpha(ms::GradMekhMS, i::Int) = ms.n_k[i]

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


function GradMekh(E, ν, n, l, Hg, Hl, m,
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

    GradMekh(E, ν, n, l, Hg, Hl, m, σy, tstar, angles, s_x_m,
             NSLIP)
end

create_matstat(::Type{GradMekh}) = GradMekhMS()

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
function stress(mat::GradMekh, matstat::GradMekhMS, temp_matstat::GradMekhMS,
                κ_nl::Vector{Float64}, F1::Matrix{Float64})

    ε = sym_V9(M_2_V9(F1))
    ε[1] -= 1.0
    ε[2] -= 1.0
    ε[3] -= 1.0

    println("ε $ε")

    E = mat.E
    ν = mat.ν
    s_x_m = mat.s_x_m
    tstar = mat.tstar

    dt = 2.000000000000000e-003 #TODO: Fix

    n_ε_p = matstat.n_ε_p
      println("n_ε_p $n_ε_p ")

    ε_e = ε - n_ε_p

    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))
    σ = L * dot(I, ε_e)*I + 2*G*ε_e

    # determine the sign of s_alpha
    for α in 1:NSLIP
       sxm_α = s_x_m[α]
       t = dot(Mbar, sxm_α)
       if t < 0
          @devec sxm_α[:] = -sxm_α
       end
    end

    #Solve nonlinear system of equations
    #unknowns: invFp,dlambda(α),Φ(α)

    # Below could be used for guessing
    #n_k = matstat.n_k
    #n_∆λ = matstat.n_∆λ
    @devec ε_p = n_ε_p
    k = zeros(NSLIP)
    ∆λ = zeros(NSLIP)

    unbal = typemax(Float64)
    niter = 0
    σ = zeros(9)
    abstol = 1e-6
    H = 1e-7
    max_niter = 8
    J = zeros(NSLIP, NSLIP)
    while true
        # För ett givet ∆λ_α beräkning av obalanskrafter
        R, k = compute_imbalance(mat, matstat, ∆λ, ε, κ_nl, dt)
        println(R)
        if norm(R) < abstol * tstar
            println("Converged in $niter iterations")
            break
        else
            for α in 1:NSLIP
                ∆λ[α] += H
                Rdiff, _ = compute_imbalance(mat, matstat, ∆λ, ε, κ_nl, dt)
                J[:, α] = (Rdiff - R) / H
                ∆λ[α] -= H
            end
            if niter > max_niter
                 error("*** No convergence in const subroutine")
            end

            d∆λ = J \ R
             println("d∆λ $d∆λ")
            #Newton type of update, but dlambda never negative
            @devec ∆λ[:] = ∆λ - d∆λ
             println("∆λ: $∆λ")
            for α=1:NSLIP
                if ∆λ[α] < 0
                    ∆λ[α] = 1.e-10
                end
            end

            niter=niter+1
        end
    end  #end iteration on ∆λ

    @devec temp_matstat.n_ε_p[:] = ε_p
    @devec temp_matstat.n_k[:] = k
    @devec temp_matstat.n_∆λ[:] = ∆λ

    return σ
end

function compute_imbalance(mat, matstat, ∆λ, ε, κ_nl, dt)
    E = mat.E
    ν = mat.ν
    s_x_m = mat.s_x_m
    tstar = mat.tstar
    n_ε_p = matstat.n_ε_p
    n_k = matstat.n_k
    n = mat.m


    R = zeros(NSLIP)
    Φ = zeros(NSLIP)
    k = zeros(NSLIP)

    # εp = εp - Σ ∆λ_α s_α ⊗ m_α
    @devec ε_p = n_ε_p
    for α in 1:NSLIP
        ε_p -= ∆λ[α] * s_x_m[α]
    end

    ε_e = ε - ε_p

    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))
    # Calculate Mandel stress
    σ = L * dot(I, ε_e)*I + 2*G*ε_e

    for α=1:NSLIP
       k[α] = n_k[α] - ∆λ[α]
       Φ[α] = dot(σ, s_x_m[α]) - (-mat.Hl * k[α] + κ_nl[α]) - mat.σy
        # Computation of R_\Delta_\lambda_\alpha
       R[α] = tstar*∆λ[α] - dt * max(0, Φ[α])^n
    end
    return R, k
end


function  const_unbal_dam_se_lem(mat, matstat, ∆λ, ε, κ_nl, dt)

    E = mat.E
    ν = mat.nu
    s_x_m = mat.s_x_m
    tstar = mat.tstar
    m = mat.m
    Hl = mat.Hl

    n_ε_p = matstat.n_ε_p
    n_k = matstat.n_k

    J = zeros(NSLIP, NSLIP)
    R = zeros(NSLIP)
    Φ = zeros(NSLIP)
    k_α = zeros(NSLIP)

    # εp = εp - Σ ∆λ_α s_α ⊗ m_α
    @devec ε_p[:] = n_ε_p
    for α in 1:NSLIP
        ε_p -= ∆λ[α] * s_x_m[α]
    end

    ε_e = ε - ε_p

    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))
    # Calculate Mandel stress
    σ = L * dot(I, ε_e)*I + 2*G*ε_e


    for α=1:NSLIP
       k_α[α] = n_k[α] - λ[α]
       Φ[α] = dot(σ, s_x_m[α]) - (-Hl * k_α[α] + κ_nl[α]) - mat.σy
        # Computation of R_\Delta_\lambda_\alpha
       R[α] = tstar*λ[α] - dt*max(0, Φ[α])^n
    end

    for α = 1:NSLIP, β = 1:NSLIP
        if Φ[α] > 0
            J[α,β]= -2 * G * dot( s_x_m[α], dot(D, s_x_m[β]))
        else
            J[α,β] = 0
        end
    end

    for α =1:NSLIP, β = 1:NSLIP
        if Φ[α] > 0
             if α == β
                 J[α,β] -= Hl
                 J[α,β] = tstar - dt*n*Φ[α]^(n-1) * J[α,β]
             else
                 J[α,β] = -dt*n * Φ[α]^(n-1) * J[α,β]
             end
        else
            J[α,β] = 0
            if α == β
                J[α,β] = tstar
            end
        end
    end
    return J, R, σ, k_α, ε_p
end

end # module
