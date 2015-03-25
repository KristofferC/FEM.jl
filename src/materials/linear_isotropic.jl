immutable LinearIsotropic <: AbstractMaterial
    E::Float64
    ν::Float64
    G::Float64
    k::Matrix{Float64}
    σ::Vector{Float64}
    matstats::Dict{Int, Vector{LinearIsotropicMS}}
    temp_matstats::Dict{Int, Vector{LinearIsotropicMS}}
end



function addmatstat!{T <: AbstractFElement}(mat::LinearIsotropic, ele::T)
    mat_stats = Array(LinearIsotropicMS, 0)
    # One material status per Gauss Point
    for i in 1:length(elem.gps)
        mat_stats[i] = create_matstat(mat)
    end
    mat.matstats[ele.n] = mat_stats
end


function LinearIsotropic(E, ν)
    G = E / (2 * (1 + ν))
    f = E / ((1.0 + ν) * (1.0 - 2.0 * ν))

    k = zeros(4, 4)

    k[1, 1] = k[2, 2] = k[3, 3] = f * (1.0 - ν)
    k[4, 4]  = G

    # Off diagonals
    k[1, 2] = k[1, 3] = k[2, 1] = k[2, 3] = k[3, 1] = k[3, 2] = f * ν
    LinearIsotropic(E, ν, G, k, zeros(4))
end
stiffness(mat::LinearIsotropic, ::GaussPoint2) = mat.k

function stress(mat::LinearIsotropic, ɛ::Vector{Float64}, gp::GaussPoint2)
    D = stiffness(mat, gp)
    A_mul_B!(mat.σ, D, ɛ)
    return mat.σ
end

create_matstat(::Type{LinearIsotropic}) = LinearIsotropicMS()

immutable LinearIsotropicMS <:AbstractMaterialStatus
    strain::Vector{Float64}
    stress::Vector{Float64}
end
LinearIsotropicMS() = LinearIsotropicMS(zeros(9), zeros(9))

update(mat::LinearIsotropic) = mat.matstat = mat.temp_matstat