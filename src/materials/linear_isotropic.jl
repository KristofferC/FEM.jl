immutable LinearIsotropicMaterial <: Material
    E::Float64
    ν::Float64
    G::Float64
end

function LinearIsotropicMaterial(E, ν)
    G = E / (2 * (1 + ν))
    LinearIsotropicMaterial(E, ν, G)
end

function stiffness(mat::LinearIsotropicMaterial, gp::GaussPoint)
    ν = mat.ν
    f = mat.E / ((1.0 + mat.ν) * (1.0 - 2.0 * ν))

    k[1, 1] = k[2, 2] = k[3, 3] = f * (1.0 - ν)
    k[4, 4] = self.G

    [k[i, j] = f * self.ν for i in

    # Off diagonals
    k[1, 2] = k[1, 3] = k[2, 1] = k[2, 3] = k[3, 1] = k[3, 2] = f * self.ν




