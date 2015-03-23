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
    n_dofs = get_ndofs(elem)

    for gp in elem.gps
        #Be = Bmatrix(elem, gp, nodes)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        A_mul_B!(elem.lts.DeBe, De, elem.lts.B)
        #transpose!(elem.lts.Bet, Be)
        A_mul_B!(elem.lts.Ke, elem.lts.Bet, elem.lts.DeBe)
        for i in 1:n_dofs
            for j in 1:n_dofs
                elem.lts.Ke[i,j] *= dV
            end
        end
    end
    return elem.lts.Ke
end

function get_field{T <: AbstractFElement}(elem::T, nodes::Vector{FENode2})
    u = zeros(get_ndofs(elem))
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
    u = get_field(elem, nodes)
    for gp in elem.gps
        B = Bmatrix(elem, gp, nodes)
        transpose!(elem.lts.Bet, B)
        A_mul_B!(elem.lts.ɛ, B, u)
        #ɛ = B * u
        σ = stress(mat, elem.lts.ɛ, gp)
        dV = weight(elem, gp, nodes)
        A_mul_B!(elem.lts.f_int, elem.lts.Bet, σ)
        for i in 1:get_ndofs(elem)
            elem.lts.f_int[i] *= dV
        end
    end
    return elem.lts.f_int
end


#function strain{T <: AbstractFElement}(elem::T, gp::GaussPoint2,
#                nodes::Vector{FENode2}, u::Vector{Float64})
#    B = Bmatrix(elem, gp, nodes)
#    return B * u
#end

function weight{T <: AbstractFElement}(elem::T, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return det(J) * gp.weight
end
