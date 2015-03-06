abstract Element

include("lin_trig_element.jl")

function stiffness(elem::Element, mesh::Mesh)
    Ke = zeros(elem.n_dofs, elem.n_dofs)

    for gp in elem.gps
        Be = Bmatrix(gp, mesh)
        De = constitutive(gp)
        dV = weight(gp, mesh)
        Ke += Be' * De * Be * dV
    end

    return Ke
end

function strain(elem::Element, gp::GaussPoint, mesh::Mesh)
    B = Bmatrix(elem, gp, mesh)
