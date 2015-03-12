abstract Element

include("lin_trig_element.jl")
include("bilin_quad_element.jl")

function stiffness(elem::Element, nodes::Vector{Node}, material::Material,
                   mp::MatPool, vp::VecPool)

    Ke = getmat(elem.n_dofs, elem.n_dofs, "Ke", mp)
    for i in size(Ke, 1)
        for j in size(Ke, 2)
            Ke[i,j] = 0.0
        end
    end

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes, mp, vp)
        De = stiffness(material, gp, mp, vp)
        dV = weight(elem, gp, nodes, mp, vp)
        Ke += Be' * De * Be * dV
    end

    return Ke
end

function get_field(elem::Element, nodes::Vector{Node})
    u = zeros(elem.n_dofs)
    i = 1
    for vert in elem.vertices
        for dof in nodes[vert].dofs
            u[i] = dof.value
            i += 1
        end
    end
    return u
end

function intf(elem::Element, nodes::Vector{Node}, mat::Material, mp::MatPool, vp::VecPool)
    f_int = zeros(elem.n_dofs)
    u = get_field(elem, nodes)
    for gp in elem.gps
        B = Bmatrix(elem, gp, nodes, mp, vp)
        ɛ = B * u
        σ = stress(mat, ɛ, gp, mp, vp)
        dV = weight(elem, gp, nodes, mp, vp)
        f_int += B' * σ * dV
    end
    return f_int
end

function strain(elem::Element, gp::GaussPoint, nodes::Vector{Node},
                u::Vector{Float64}, mp::MatPool, vp::VecPool)
    B = Bmatrix(elem, gp, nodes, mp, vp)
    strain = getvec(4, "e", vp)
    A_mul_B!(B, u, strain)
    return strain
end

function weight(elem::Element, gp::GaussPoint, nodes::Vector{Node}, mp::MatPool, vp::VecPool)
    dN = dNmatrix(elem.interp, gp.local_coords, mp, vp)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN, mp, vp)
    return det(J) * gp.weight
end
