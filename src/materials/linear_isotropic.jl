immutable LinearIsotropicMS <:AbstractMaterialStatus
    strain::Vector{Float64}
    stress::Vector{Float64}
end


LinearIsotropicMS() = LinearIsotropicMS(zeros(6), zeros(6))

copy(matstat::LinearIsotropicMS) = LinearIsotropicMS(copy(matstat.strain), copy(matstat.stress))


type LinearIsotropic <: AbstractMaterial
    E::Float64
    ν::Float64
    G::Float64
    k::Matrix{Float64}
    σ::Vector{Float64}
    matstats::Vector{Vector{LinearIsotropicMS}}
    temp_matstats::Vector{Vector{LinearIsotropicMS}}
end


function addmatstats!(mat::LinearIsotropic, n::Int)
    mat_stats = Array(LinearIsotropicMS, 0)
    temp_matstats = Array(LinearIsotropicMS, 0)
    for i in 1:n
        push!(mat_stats, create_matstat(typeof(mat)))
        push!(temp_matstats, create_matstat(typeof(mat)))
    end
    push!(mat.matstats, mat_stats)
    push!(mat.temp_matstats, temp_matstats)
end


#=
function addmatstat!(mat::LinearIsotropic, i::Int)
    push!(mat.matstats[i], create_matstat(typeof(mat)))
    push!(mat.temp_matstats[i], create_matstat(typeof(mat)))
end
=#


function LinearIsotropic(E, ν)
    G = E / (2 * (1 + ν))
    f = E / ((1.0 + ν) * (1.0 - 2.0 * ν))

    k = zeros(4, 4)

    k[1, 1] = k[2, 2] = k[3, 3] = f * (1.0 - ν)
    k[4, 4]  = G

    # Off diagonals
    k[1, 2] = k[1, 3] = k[2, 1] = k[2, 3] = k[3, 1] = k[3, 2] = f * ν

    matstats = Array(Vector{LinearIsotropicMS}, 0)
    temp_matstats = Array(Vector{LinearIsotropicMS}, 0)
    LinearIsotropic(E, ν, G, k, zeros(4), matstats, temp_matstats)
end

stiffness(mat::LinearIsotropic, ::GaussPoint2) = mat.k
create_matstat(::Type{LinearIsotropic}) = LinearIsotropicMS()


function stress(mat::LinearIsotropic, ɛ::Vector{Float64}, gp::GaussPoint2)
    D = stiffness(mat, gp)
    A_mul_B!(mat.σ, D, ɛ)
    return mat.σ
end


update!(mat::LinearIsotropic) = mat.matstats = mat.temp_matstats
