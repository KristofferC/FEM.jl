
module GradMekhModJl

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
    n_k_α::Vector{Float64}
    n_λ_α::Vector{Float64}
    strain::Vector{Float64}
    stress::Vector{Float64}
end

const NSLIP = 2

function GradMekhMS()
  GradMekhMS(zeros(9), zeros(NSLIP), zeros(NSLIP), zeros(6), zeros(6))
end

copy(matstat::GradMekhMS) = GradMekhMS(copy(matstat.n_ε_p), copy(matstat.n_k_α), copy(matstat.n_λ_α),
                                       copy(matstat.strain), copy(matstat.stress))

get_kalpha(ms::GradMekhMS, i::Int) = ms.n_k_α[i]

immutable GradMekh <: AbstractMaterial
    E::Float64
    ν ::Float64
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


function GradMekh(E, nu, n, l, Hg, Hl, m, factor,
                  sy, tstar, angles, NSLIP)

    if length(angles) != NSLIP
        error("Need one angle per slip system")
    end

    s_x_m = Vector{Vector{Float64}}()
    @devec angles[:] = angles[:] .* pi ./ 180.0
    for i = 1:NSLIP
        t = angles[i]
        s = [cos(t), sin(t), 0.0]
        mm = [cos(t + pi/2), sin(t + pi/2), 0.0]
        sxm_mat = s * mm'
        push!(s_x_m, M_2_V9(sxm_mat))
    end

    GradMekh(E, nu, n, l Hg, Hl, m,
             factor, sy, tstar, angles, s_x_m,
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

    ε = sym_V9!(M_2_V9(F1))
    ε[1] -= 1.0
    ε[2] -= 1.0
    ε[3] -= 1.0

    E = mat.E
    ν = mat.nu
    s_x_m = mat.s_x_m
    tstar = mat.tstar

    dt = 2.000000000000000e-003 #TODO: Fix

    n_ε_p = matstat.n_ε_p
    n_k_α = matstat.n_k_α
    n_λ_α = matstat.n_λ_α

    ε_e = ε - n_ε_p

    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))
    Mbar = L * dot(I, ε_e)*I + 2*G*ε_e

    # determine the sign of s_alpha
    for i in 1:NSLIP
       sxm = s_x_m[i]
       t = dot(Mbar, sxm)
       if t < 0
          @devec sxm[:] = -sxm
       end
    end

    #Solve nonlinear system of equations
    #unknowns: invFp,dlambda(ii),Φ(ii)
    # initial guess:
    ε_p = zeros(9)
    λ_α = zeros(NSLIP)
    k_α = zeros(NSLIP)

    unbal = typemax(Float64)
    niter = 0

    σ = zeros(9)
    abstol = 1e-6
    while true
        # För ett givet λ beräkning av obalanskrafter
        J, R, σ, k_α, ε_p = const_unbal_dam_se_lem(mat, matstat, λ, ε, κ_nl, dt)

        if niter > 8
             error("*** No convergence in const subroutine")
        end

        if norm(R) < abstol * tstar
            println("Converged in $niter iterations")
            break
        else
            ∆λ = J \ R
            #Newton type of update, but dlambda never negative
            @devec λ[:] = λ - ∆λ
            for ii=1:NSLIP
                if λ[ii] < 0
                    λ[ii]=1.e-10
                end
            end
            niter=niter+1
        end
    end  #end iteration on λ

    @devec temp_matstat.n_ε_p[:] = ε_p
    @devec temp_matstat.n_k_α[:] = k_α
    @devec temp_matstat.n_λ_α[:] = λ

    return σ
end

function  const_unbal_dam_se_lem(mat, matstat, λ, ε, κ_nl, dt)

    E = mat.E
    ν = mat.nu
    s_x_m = mat.s_x_m
    tstar = mat.tstar
    n = mat.m
    m = mat.m
    Hl = mat.Hl

    n_ε_p = matstat.n_ε_p
    n_k_α = matstat.n_k_α

    J = zeros(NSLIP, NSLIP)
    R = zeros(NSLIP)
    Φ = zeros(NSLIP)
    k_α = zeros(NSLIP)
    @devec ε_p[:] = n_ε_p

    for ii in 1:NSLIP
        ε_p -= λ[ii] * s_x_m[ii]
    end

    ε_e = ε - ε_p

    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))
    # Calculate Mandel stress
    σ = L * dot(I, ε_e)*I + 2*G*ε_e


    for ii=1:NSLIP
       k_α [ii] = n_k_α[ii] - λ[ii]
       Φ[ii] = dot(σ, s_x_m[ii]) - (-Hl * k_α[ii] + κ_nl[ii]) - mat.σy
        # Computation of R_\Delta_\lambda_\alpha
       R[ii] = tstar*λ[ii] - dt*max(0, Φ[ii])^n
    end

    for ii = 1:NSLIP, jj = 1:NSLIP
        if Φ[ii] > 0
            J[ii,jj]= -2 * G * dot( s_x_m[ii], dot(D, s_x_m[jj]))
        else
            J[ii,jj] = 0
        end
    end

    for ii =1:NSLIP, jj = 1:NSLIP
        if Φ[ii] > 0
             if ii == jj
                 J[ii,jj] -= Hl
                 J[ii,jj] = tstar - dt*n*Φ[ii]^(n-1) * J[ii,jj]
             else
                 J[ii,jj] = -dt*n * Φ[ii]^(n-1) * J[ii,jj]
             end
        else
            J[ii,jj] = 0
            if ii == jj
                J[ii,jj] = tstar
            end
        end
    end
    return J, R, σ, k_α, ε_p
end

end # module
