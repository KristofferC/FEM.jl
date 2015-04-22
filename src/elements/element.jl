abstract AbstractFElement{T <: AbstractMaterialStatus}

abstract ElemStorage

include("lin_trig.jl")
include("lin_quad.jl")
include("quad_trig.jl")
include("grad_trig.jl")

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
    @inbounds for vert in elem.vertices
        for dof in nodes[vert].dofs
            u[i] = dof.value
            i += 1
        end
    end
    return u
end


function intf{T <: AbstractFElement, P <: AbstractMaterial}(elem::T, mat::P, nodes::Vector{FENode2})
    u = get_field(elem, nodes)
    ɛ = elem.storage.ɛ
    fill!(elem.storage.f_int, 0.0)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(ɛ, B, u)
        fill_from_start!(elem.temp_matstats[i].strain, ɛ)

        σ = stress(mat, ɛ, gp)
        fill_from_start!(elem.temp_matstats[i].stress, σ)

        dV = weight(elem, gp, nodes)

        # f_int += B' * σ * dV
        BLAS.gemv!('T', dV, B, σ, 1.0, elem.storage.f_int)
    end
    return elem.storage.f_int
end

#get_cell_data{T <: AbstractScalar}(::AbstractFElement, ::Type{T}) = 0.0
#get_point_data{T <: AbstractScalar}(::AbstractFElement, ::Type{T}) = 0.0

#get_cell_data{T <: AbstractTensor}(::AbstractFElement, ::Type{T}) = zeros(6)
#get_point_data{T <: AbstractTensor}(::AbstractFElement, ::Type{T}) = zeros(6)

#get_cell_data{T <: AbstractVector}(::AbstractFElement, ::Type{T}) = zeros(3)
#get_point_data{T <: AbstractVector}(::AbstractFElement, ::Type{T}) = zeros(3)


function get_cell_data{T <: AbstractTensor}(elem::AbstractFElement, field::Type{T})
    cellfield = zeros(get_ncomponents(field))
    for (i, gp) in enumerate(elem.gps)
        gpfield = get_field(elem, field, i)
        axpy!(getweight(gp), gpfield, cellfield)
    end
    return cellfield
end

function weight(elem::AbstractFElement, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end