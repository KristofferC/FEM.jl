abstract AbstractFElement{T <: AbstractMaterialStatus}

abstract AbstractElemStorage

@lazymod LinTrigMod "lin_trig.jl"
@lazymod LinQuadMod "lin_quad.jl"
@lazymod QuadTrigMod "quad_trig.jl"
@lazymod GradTrigMod "grad_trig.jl"
@lazymod GradTrigPrimalMod "grad_trig_prim.jl"

# Interface for a FE-element
createstorage() = error("Not Implemented")
createinterp() = error("Not Implemented")
Bmatrix() = error("Not Implemented")
creategps() = error("Not Implemented")
doftypes() = error("Not Implemented")
get_ndofs() = error("Not Implemented")
get_geotype() = error("Not Implemented")
get_ref_area() = error("Not Implemented")


getindex{T <: AbstractFElement}(elem::T, i0::Int) = getindex(elem.vertices, i0)
vertices(elem::AbstractFElement) = elem.vertices
gausspoints(elem::AbstractFElement) = elem.gps


function show{T <: AbstractFElement}(io::IO,elem::T)
    print(io, string(typeof(elem), ":", elem.vertices))
end


function stiffness(elem::AbstractFElement,
                  nodes::Vector{FENode2},
                  material::AbstractMaterial,
                  dof_vals)
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

# Extracts the values of the dofs for an element.
function get_field(elem::AbstractFElement, nodes::Vector{FENode2}, dof_vals::DofVals)
    u = elem.storage.u_field
    i = 1
    @inbounds for vert in elem.vertices
        for dof in nodes[vert].dofs
            u[i] = get_value(dof_vals, dof)
            i += 1
        end
    end
    return u
end

# Computes the internal forces for an element.
function intf(elem::AbstractFElement, mat::AbstractMaterial, nodes::Vector{FENode2}, dof_vals::DofVals)
    u = get_field(elem, nodes, dof_vals)
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


function weight(elem::AbstractFElement, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end

const ZERO_SCALAR = zeros(1)
const ZERO_VECTOR = zeros(3)
const ZERO_TENSOR = zeros(9)
get_field(elem::AbstractFElement, field::Type{AbstractScalar}, ::Int) = ZERO_SCALAR
get_field(elem::AbstractFElement, field::Type{AbstractVector}, ::Int) = ZERO_VECTOR
get_field(elem::AbstractFElement, field::Type{AbstractTensor}, ::Int) = ZERO_TENSOR

function get_cell_data{T <: AbstractField}(elem::AbstractFElement, field::Type{T})
    cellfield = zeros(get_ncomponents(field))
    for (i, gp) in enumerate(elem.gps)
        gpfield = get_field(elem, field, i)
        axpy!(getweight(gp)/get_ref_area(elem), gpfield, cellfield)
    end
    return cellfield
end


# Iterator for active dofs, too slow right now.
immutable ActiveDofsIterator{T <: AbstractFElement}
    nodes::Vector{FENode2}
    ele::T
end

get_field(elem::AbstractFElement, ::Type{Stress}, i::Int) = elem.matstats[i].stress
get_field(elem::AbstractFElement, ::Type{Strain}, i::Int) = elem.matstats[i].strain

function get_field(elem::AbstractFElement, ::Type{VonMises}, i::Int)
    σ = elem.matstats[i].stress
    m = (σ[1] + σ[2] + σ[3]) / 3
    return [sqrt(3/2) * sqrt((σ[1] - m)^2 + (σ[2] - m)^2 + (σ[3] - m)^2 +
                        2(σ[4]*σ[4]))]
end

function activedofs(element::AbstractFElement, nodes::Vector{FENode2})
    ActiveDofsIterator(nodes, element)
end


function Base.start(adi::ActiveDofsIterator)
    # We need to find first active dof because
    # we cant start the iterator loop unless
    # we actually have active dofs
    skip_inactive(adi, (1, 0, 0))
end


# Iterator for
function Base.done(adi::ActiveDofsIterator, state)
    state[1] > length(vertices(adi.ele))
end

function Base.next(adi::ActiveDofsIterator, state)
    v = state[1]
    i = state[2]
    n = state[3]

    node = adi.nodes[vertices(adi.ele)[v]]
    dof = get_dof(node, i)
    return((dof, n), skip_inactive(adi, state))
end


function skip_inactive(adi::ActiveDofsIterator, state)
    v = state[1]
    i = state[2] + 1
    n = state[3]

    node = adi.nodes[vertices(adi.ele)[v]]

    # Finish dofs for this node
    for i in i:length(get_dofs(node))
        n+= 1
        if isactive(get_dof(node, i))
            return (v, i, n)
        end
    end

    # Start looking for new active dof
    v+=1
    for v in v:length(vertices(adi.ele))
        node = adi.nodes[vertices(adi.ele)[v]]
        for i in 1:length(get_dofs(node))
            n += 1
            dof = get_dof(node, i)
            if isactive(dof)
                return (v, i, n)
            end
        end
    end
    return (v+1, i, n)
end
