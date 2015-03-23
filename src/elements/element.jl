abstract AbstractFElement #{T <: AbstractMaterialStatus}

include("lin_trig_element.jl")
#include("bilin_quad_element.jl")


getindex{T <: AbstractFElement}(elem::T, i0::Real) = getindex(elem.vertices, i0)

#=
function stiffness{T <: AbstractFElement, P <: AbstractMaterial}
                  (elem::T, mat::P, nodes::Vector{FENode2})
    Ke = zeros(elem.n_dofs, elem.n_dofs)

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(mat, gp)
        dV = weight(elem, gp, nodes)
        Ke += Be.' * De * Be * dV
    end

    return Ke
end
=#

function stiffness{T <: AbstractFElement,  P <: AbstractMaterial}(elem::T,
                                                                  nodes::Vector{FENode2},
                                                                  material::P)
    Ke = zeros(elem.n_dofs, elem.n_dofs)

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        A_mul_B!(elem.lts.DeBe, De, Be)
        transpose!(elem.lts.Bet, Be)
        A_mul_B!(elem.lts.Ke, elem.lts.Bet, elem.lts.DeBe)
        for i in 1:elem.n_dofs
            for j in 1:elem.n_dofs
                Ke[i,j] += elem.lts.Ke[i,j] * dV
            end
        end
    end
    return Ke
end

function get_field{T <: AbstractFElement}(elem::T, nodes::Vector{FENode2})
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


function intf{T <: AbstractFElement, P <: AbstractMaterial}(elem::T, mat::P, nodes::Vector{FENode2})
    f_int = zeros(elem.n_dofs)
    u = get_field(elem, nodes)
    for gp in elem.gps
        B = Bmatrix(elem, gp, nodes)
        ɛ = B * u
        σ = stress(mat, ɛ, gp)
        dV = weight(elem, gp, nodes)
        f_int += B.' * σ * dV
    end
    return f_int
end


function strain{T <: AbstractFElement}(elem::T, gp::GaussPoint2,
                nodes::Vector{FENode2}, u::Vector{Float64})
    B = Bmatrix(elem, gp, nodes)
    return B * u
end

function weight{T <: AbstractFElement}(elem::T, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return det(J) * gp.weight
end
