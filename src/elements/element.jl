abstract Element

include("lin_trig_element.jl")

function stiffness(elem::Element, nodes::Vector{Node}, material::Material)
    Ke = zeros(elem.n_dofs, elem.n_dofs)

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        Ke += Be' * De * Be * dV
    end
    return Ke
end

#function strain(elem::Element, gp::GaussPoint, nodes::Vector{Node})
#    B = Bmatrix(elem, gp, nodes)'
#end
