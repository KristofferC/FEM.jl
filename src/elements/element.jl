abstract Element

include("lin_trig_element.jl")
include("bilin_quad_element.jl")


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

# Calculates the internal forces
function intf(elem::Element, nodes::Vector{Node}, mat::Material)
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

# TODO: Test this
# TODO: Not working, fix
#=
function bodyf(elem::Element, nodes::Vector{Node}, mat::Material, f::Float64)
    fb = zeros(elem.n_dofs)
    for gp in elem.gps
        N = Nvec(elem.interp, gp.local_coords)
        dV = weight(elem, gp, nodes)
        fb += N' * f * dV
    end
    return fb
end
=#

function strain(elem::Element, gp::GaussPoint, nodes::Vector{Node}, u::Vector{Float64})
    B = Bmatrix(elem, gp, nodes)
    return B * u
end


function weight(elem::Element, gp::GaussPoint, nodes::Vector{Node})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return det(J) * gp.weight
end
