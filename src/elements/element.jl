abstract AbstractFElement #{T <: AbstractMaterialStatus}

include("lin_trig_element.jl")
#include("bilin_quad_element.jl")


getindex{T <: AbstractFElement}(elem::T, i0::Real) = getindex(elem.vertices, i0)


function stiffness{T <: AbstractFElement, P <: AbstractMaterial}(elem::T,
                                                                 nodes::Vector{FENode2},
                                                                 mat::P)
    Ke = zeros(get_ndofs(elem), get_ndofs(elem))

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(mat, gp)
        dV = weight(elem, gp, nodes)
        Ke += Be.' * De * Be * dV
    end

    return Ke
end

#=
function stiffness{T <: AbstractFElement,  P <: AbstractMaterial}(elem::T,
                                                                  nodes::Vector{FENode2},
                                                                  material::P)
    n_dofs = get_ndofs(elem)

    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        println(Be)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        A_mul_B!(elem.lts.DeBe, De, Be)
        At_mul_B!(elem.lts.Ke, Be, elem.lts.DeBe)
        scale!(elem.lts.Ke, dV)
    end
    return elem.lts.Ke
end
=#

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

function intf{T <: AbstractFElement, P <: AbstractMaterial}(elem::T, mat::P,
                          nodes::Vector{FENode2})
    f_int = zeros(get_ndofs(elem))
    u = get_field(elem, nodes)
    for gp in elem.gps
        B = Bmatrix(elem, gp, nodes)
        println(B)
        ɛ = B * u
        σ = stress(mat, ɛ, gp)
        dV = weight(elem, gp, nodes)
        f_int += B.' * σ * dV
    end
    return f_int
end


#=

function intf{T <: AbstractFElement, P <: AbstractMaterial}(elem::T, mat::P, nodes::Vector{FENode2})
    u = get_field(elem, nodes)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(elem.lts.ɛ, B, u)

        σ = stress(mat, elem.lts.ɛ, gp)
        dV = weight(elem, gp, nodes)
        At_mul_B!(elem.lts.f_int, B, σ)
        scale!(elem.lts.f_int, dV)
            # Store in matstat

        copy!(mat.temp_matstats[elem.n][i].strain, elem.lts.ɛ)
        copy!(mat.temp_matstats[elem.n][i].stress, σ)
    end
    return elem.lts.f_int
end
=#

function weight{T <: AbstractFElement}(elem::T, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return det(J) * gp.weight
end
