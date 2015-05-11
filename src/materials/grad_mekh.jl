
module GradMekhMod

using Devectorize
import Base.copy

# Function sto extend
import FEM: create_matstat, stiffness, stress, get_kalpha

# Types needed
import FEM: AbstractMaterialStatus, AbstractMaterial, GaussPoint2

export GradMekhMS, GradMekh

type GradMekhMS <:AbstractMaterialStatus
    state::Vector{Float64}
    strain::Vector{Float64}
    stress::Vector{Float64}
end

const NSLIP = 2

function GradMekhMS()
  state = zeros(9 + 2*NSLIP + 1)
  state[1:3] = 1.0
  GradMekhMS(state, zeros(6), zeros(6))
end

copy(matstat::GradMekhMS) = GradMekhMS(copy(matstat.state), copy(matstat.strain), copy(matstat.stress))
get_kalpha(ms::GradMekhMS, i::Int) = ms.state[9 + i]

immutable GradMekh <: AbstractMaterial
    E::Float64
    nu::Float64
    nn::Float64
    l::Float64
    kinf::Float64
    lambda0::Float64
    Hg::Float64
    Hl::Float64
    m::Float64
    factor::Float64
    Sy::Float64
    tstar::Float64
    c_dam::Float64
    angles::Vector{Float64}
    s_x_m::Matrix{Float64}
    nslip::Int
    para::Vector{Float64}
    libmekh::Ptr{Void}
end

function M_2_V9(b::Matrix{Float64})
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

function GradMekh(E, nu, n, l , kinf, lambda_0, Hg, Hl, m, factor,
                  sy, tstar, c_dam, angles, nslip)

    if length(angles) != nslip
        error("Need one angle per slip system")
    end

    nmatpara = 13;
    npara = nmatpara + length(angles)
    nstate = 9 + 2 * nslip + 1


    para = Float64[E, nu, n, l, kinf, lambda_0, Hg, Hl, m, factor,
                        sy, tstar, c_dam]

    for angle in angles
        push!(para, angle)
    end


    s_x_m = zeros(9, nslip)
    @devec angles[:] = angles[:] .* pi ./ 180.0
    for i = 1:nslip
        t = angles[i]
        s = [cos(t), sin(t), 0.0]
        mm = [cos(t + pi/2), sin(t + pi/2), 0.0]
        sxm_mat = s * mm'
        s_x_m[:, i] = M_2_V9(sxm_mat)
    end

    mekh_lib = Libdl.dlopen("combined_material", Libdl.RTLD_GLOBAL)
    mekh_mat = Libdl.dlsym(mekh_lib, "__material_module_MOD_const_solve_dam_se_lem")

    GradMekh(E, nu, n, l, kinf, lambda_0, Hg, Hl, m,
             factor, sy, tstar, c_dam, angles, s_x_m,
             nslip, para, mekh_mat)
end

create_matstat(::Type{GradMekh}) = GradMekhMS()

const STRESS_BUFFER = zeros(9)
function stress(mat::GradMekh, matstat::GradMekhMS, temp_matstat::GradMekhMS, kappas::Vector{Float64}, F::Matrix{Float64})
    ɛ = temp_matstat.strain

  F1 = M_2_V9(F)

   dtime = 2.000000000000000e-003 #TODO: Fix
   LKONV = [1]
   npara = length(mat.para)

   state_old = matstat.state
   state_new = temp_matstat.state
   kappa_nl = kappas

    ccall(mat.libmekh,
               Void,
                (Ptr{Int}, Ptr{Float64}, Ptr{Int},
                 Ptr{Float64}, Ptr{Float64},
                 Ptr{Float64}, Ptr{Float64},
                 Ptr{Int}, Ptr{Float64},
                 Ptr{Float64},
                 Ptr{Float64},
                 Ptr{Float64},
                 Ptr{Int}),
                &length(mat.para), mat.para, &length(state_old),
                state_old, state_old, &dtime, F1,
                &NSLIP, kappa_nl, mat.s_x_m,
                state_new, STRESS_BUFFER, LKONV)

    if LKONV[1] != 1
      println("F: $F1")
      println("ɛ: $ɛ")
      println("u: $u")
      println("B: $B")
      println("invFP: $(state_old[1:9])")
      println("state_new $state_new")
      println("kappas: $kappas")
      println("stress: $STRESS_BUFFER")
      exit()
    end

   return STRESS_BUFFER
end

end
