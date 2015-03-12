import FEM.GaussPoint

immutable LinearIsotropic <: Material
    E::Float64
    ν::Float64
    G::Float64
end

function LinearIsotropic(E, ν)
    G = E / (2 * (1 + ν))
    LinearIsotropic(E, ν, G)
end

function stiffness(mat::LinearIsotropic, gp::GaussPoint, mp::MatPool, vp::VecPool)
    ν = mat.ν
    f = mat.E / ((1.0 + ν) * (1.0 - 2.0 * ν))

    k = getmat(4, 4, "k", mp)

    k[1, 1] = k[2, 2] = k[3, 3] = f * (1.0 - ν)
    k[4, 4]  = mat.G

    # Off diagonals
    k[1, 2] = k[1, 3] = k[2, 1] = k[2, 3] = k[3, 1] = k[3, 2] = f * ν

    return k
end

function stress(mat::LinearIsotropic, ɛ::Vector{Float64}, gp::GaussPoint,
                mp::MatPool, vp::VecPool)
    D = stiffness(mat, gp, mp, vp)
    σ = getvec(4, "s", vp)
    #println(size(σ))
    #println(size(D*ɛ))
    #A_mul_B!(D, ɛ, σ)
    return D*ɛ
end
