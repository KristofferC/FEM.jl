abstract AbstractFElement

include("lin_trig_element.jl")
include("bilin_quad_element.jl")


getindex(elem::AbstractFElement, i0::Real) = getindex(elem.vertices, i0)

function stiffness(elem::AbstractFElement, nodes::Vector{FENode2}, material::Material)
    Ke = zeros(elem.n_dofs, elem.n_dofs)

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        Ke += Be' * De * Be * dV
    end

    return Ke
end

function get_field(elem::AbstractFElement, nodes::Vector{FENode2})
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


function intf(elem::AbstractFElement, nodes::Vector{FENode2}, mat::Material)
    f_int = zeros(elem.n_dofs)
    u = get_field(elem, nodes)
    for gp in elem.gps
        B = Bmatrix(elem, gp, nodes)
        ɛ = B * u
        σ = stress(mat, ɛ, gp)
        dV = weight(elem, gp, nodes)
        f_int += B' * σ * dV
    end
    return f_int
end


function strain(elem::AbstractFElement, gp::AbstractGaussPoint, nodes::Vector{FENode2}, u::Vector{Float64})
    B = Bmatrix(elem, gp, nodes)
    return B * u
end

function weight(elem::AbstractFElement, gp::AbstractGaussPoint, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return det(J) * gp.weight
end
