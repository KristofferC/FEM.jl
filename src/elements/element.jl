abstract AbstractFElement{T <: AbstractMaterialStatus}

abstract ElemStorage

include("lin_trig_element.jl")
include("lin_quad_element.jl")


getindex{T <: AbstractFElement}(elem::T, i0::Real) = getindex(elem.vertices, i0)

function show{T <: AbstractFElement}(io::IO,elem::T)
    print(io, string(typeof(elem), ":", elem.vertices))
end



function stiffness{T <: AbstractFElement,  P <: AbstractMaterial}(elem::T,
                                                                  nodes::Vector{FENode2},
                                                                  material::P)
    fill!(elem.storage.Ke, 0.0)
    for gp in elem.gps
        Be = Bmatrix(elem, gp, nodes)
        De = stiffness(material, gp)
        dV = weight(elem, gp, nodes)
        A_mul_B!(elem.storage.DeBe, De, Be) # DeBe = De * Be
        # Ke += B' * DeBe * dV
        BLAS.gemm!('T', 'N' ,dV, Be, elem.storage.DeBe, 1.0, elem.storage.Ke)
    end
    return elem.storage.Ke
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
    fill!(elem.storage.f_int, 0.0)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(elem.storage.ɛ, B, u)

        σ = stress(mat, elem.storage.ɛ, gp)
        dV = weight(elem, gp, nodes)
        # f_int += B' * σ * dV
        BLAS.gemv!('T', dV, B, σ, 1.0, elem.storage.f_int)

        copy!(elem.temp_matstats[i].strain, elem.storage.ɛ)
        copy!(elem.temp_matstats[i].stress, σ)
    end
    return elem.storage.f_int
end


#=
function weight{T <: AbstractFElement}(elem::T, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)

    return det(J) * gp.weight
end
=#