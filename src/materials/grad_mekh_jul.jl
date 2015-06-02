
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
    n_invFp::Vector{Float64}
    n_k_alpha::Vector{Float64}
    n_lambda_alpha::Vector{Float64}
    strain::Vector{Float64}
    stress::Vector{Float64}
end
get_kalpha(ms::GradMekhMS) = ms.n_k_alpha

const NSLIP = 2

function GradMekhMS()
  n_invFp = Float64[1,1,1,0,0,0,0,0,0]
  GradMekhMS(n_invFp , zeros(NSLIP), zeros(NSLIP), zeros(6), zeros(6))
end

copy(matstat::GradMekhMS) = GradMekhMS(copy(matstat.n_invFp), copy(matstat.n_k_alpha), copy(matstat.n_lambda_alpha),
                                       copy(matstat.strain), copy(matstat.stress))
get_kalpha(ms::GradMekhMS, i::Int) = ms.n_k_alpha[i]

immutable GradMekh <: AbstractMaterial
    E::Float64
    nu::Float64
    nn::Float64
    l::Float64
    Hg::Float64
    Hl::Float64
    m::Float64
    sy::Float64
    tstar::Float64
    angles::Vector{Float64}
    s_x_m::Vector{Vector{Float64}}
    NSLIP::Int
end


function GradMekh(E, nu, n, l, Hg, Hl, m,
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

    GradMekh(E, nu, n, l,Hg, Hl, m,
              sy, tstar,angles, s_x_m,
             NSLIP)
end

create_matstat(::Type{GradMekh}) = GradMekhMS()

@inbounds function M_2_V9(b::Matrix{Float64})
    a = zeros(9)
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


@inbounds function V9_d_V9(a, b)
    c = zeros(9)
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

@inbounds function sym_V9(b)
    a = zeros(9)
    a[1:3] =  b[1:3]
    a[4]   = (b[4]+b[8])/2
    a[8]   = (b[4]+b[8])/2
    a[5]   = (b[5]+b[9])/2
    a[9]   = (b[5]+b[9])/2
    a[6]   = (b[6]+b[7])/2
    a[7]   = (b[6]+b[7])/2
    return a
end

@inbounds function inv_V9(A)
    det_A = (-A[7]*A[2]*A[6]+A[4]*A[5]*A[6]+
              A[7]*A[8]*A[9]-A[1]*A[5]*A[9]-
              A[4]*A[8]*A[3]+A[1]*A[2]*A[3])

    inv_A = zeros(9)

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
                kappa_nl::Vector{Float64}, F1::Matrix{Float64})

    F = M_2_V9(F1)
   # println(F)
    E = mat.E
    ν = mat.nu
    s_x_m = mat.s_x_m
    tstar = mat.tstar

    dt = 2.000000000000000e-003 #TODO: Fix

    lkonv = 1

    n_invFp    = matstat.n_invFp
    n_k_alpha  = matstat.n_k_alpha
    dlambda    = matstat.n_lambda_alpha

    Ee = sym_V9(V9_d_V9(F, n_invFp) - I)


    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))

    # Calculate Mandel stress
    Mbar = L * dot(I,Ee)*I + 2*G*Ee
   # println("Mbar $Mbar")
    # determine the sign of s_alpha
    for i in 1:NSLIP
       sxm = s_x_m[i]
       t = dot(Mbar, sxm)
       if t < 0
          @devec sxm[:] = -sxm
       end
    end

    #Solve nonlinear system of equations
    #unknowns: invFp,dlambda(ii),phi(ii)
    # initial guess:
    @devec invFp = matstat.n_invFp
    dlambda = 0.0
    XX = zeros(NSLIP)

    unbal = typemax(Float64)
    niter = 0
    R = 0.0

    k_alpha = zeros(NSLIP)
    invFp = zeros(9)
    XX = zeros(2)
    Tau = zeros(9)
    while true
        # För ett givet XX beräkning av obalanskrafter
        Rold = R
        dR_dX = 0.0

        Jacob, R, Tau, k_alpha, invFp = const_unbal_dam_se_lem(mat, matstat, XX, F, kappa_nl, dt)


        unbal_old = unbal
        unbal = norm(R)
        #println(unbal)
        nint=0

#=
        while unbal > unbal_old && niter > 5 && nint < 20
            nint=nint+1
            @devec dXX[:] = dXX ./ 2
            @devec XX[:] = XXold .- dXX

             const_unbal_dam_se_lem(XX,n_invFp,n_k_alpha, kappa_nl,
                  invFp,R,Tau,k_alpha,ijacob)
             unbal=norm(R)
         end
=#

          if niter > 35
             error("*** No convergence in const subroutine")
          end

          #  if not convergence then update
          if unbal < 1e-6 * tstar
              break
          else
              XXold=XX
              dXX = Jacob \ R
             println(Jacob)

              #Newton type of update, but dlambda never negative
              @devec XX[:] = XXold - dXX
              for ii=1:NSLIP
                  if XX[ii] < 0
                      XX[ii]=1.d-10
                  end
             end
             niter=niter+1
       end
  end  #end iteration on XX


  @devec temp_matstat.n_invFp[:] = invFp
  @devec temp_matstat.n_k_alpha[:] = k_alpha
  @devec temp_matstat.n_lambda_alpha[:] = XX

  #Computation of 2nd Piola Kirchhoff
  S=V9_d_V9( V9_d_V9(inv_V9(F),Tau),inv_V9(trans_V9(F)) )
  #Computation of 1st Piola Kirchhoff
  Piola  = V9_d_V9(S,trans_V9(F))
  return Piola
end

function  const_unbal_dam_se_lem(mat, matstat, XX, F, kappa_nl, dt)

    E = mat.E
    ν = mat.nu
    s_x_m = mat.s_x_m
    tstar = mat.tstar
    n = mat.m

    n_invFp = matstat.n_invFp
    n_k_alpha = matstat.n_k_alpha


    n_Fp = inv_V9(n_invFp)
    Jakob = zeros(NSLIP, NSLIP)
    R = zeros(NSLIP)
    phi = zeros(NSLIP)
    k_alpha = zeros(NSLIP)

    m = mat.m
    Hl = mat.Hl

    invFp = I
    for ii in 1:NSLIP
       invFp -= XX[ii] * s_x_m[ii]
    end

    invFp = V9_d_V9(n_invFp,invFp)
    Fp = inv_V9(invFp)

    Fe = V9_d_V9(F, invFp)
    Ee = sym_V9(Fe - I)

    G = E / (2*(1+ν))
    L = E*ν / ((1+ν)*(1-2*ν))
    # Calculate Mandel stress
    Mbar = L * dot(I,Ee)*I + 2*G*Ee

    Tau = Mbar

    mclaur = 0
    for ii=1:NSLIP
       k_alpha[ii] = n_k_alpha[ii]-XX[ii]
       phi[ii] = dot(Mbar, s_x_m[ii]) - (-Hl * k_alpha[ii] + kappa_nl[ii]) - mat.sy
       mclaur = max(0, phi[ii])

        # Computation of R_\Delta_\lambda_\alpha
       R[ii] = tstar*XX[ii] - dt*mclaur^n
    end

    for ii = 1:NSLIP, jj = 1:NSLIP
        if phi[ii] > 0
            Jakob[ii,jj]= -2 * G * dot( s_x_m[ii], sym_V9(V9_d_V9(V9_d_V9(F,n_Fp) , s_x_m[jj])))
        else
           Jakob[ii,jj] = 0
        end
    end

    for ii =1:NSLIP, jj = 1:NSLIP
        if phi[ii] > 0
             if ii == jj
                 Jakob[ii,jj] -= Hl
                 Jakob[ii,jj] = tstar - dt*n*phi[ii]^(n-1) * Jakob[ii,jj]
             else
                 Jakob[ii,jj] = -dt*n * phi[ii]^(n-1) * Jakob[ii,jj]
             end
        else
            Jakob[ii,jj] = 0
            if ii == jj
                Jakob[ii,jj] = tstar
            end
        end
    end
    return Jakob, R, Tau, k_alpha, invFp
end

end # module
