abstract AbstractInterpolator

export AbstractInterpolator

@lazymod LinTrigInterpMod "lin_trig_interp.jl"
@lazymod LinQuadInterpMod "lin_quad_interp.jl"
@lazymod QuadTrigInterpMod "quad_trig_interp.jl"

dNmatrix() = error("Not implemented")
dNdxmatrix() = error("Not implemented")
Jmatrix() = error("Not implemented")
Nvec() = error("Not implemented")
get_area() = error("Not implemented")
mass_matrix() = error("Not implemented")
mass_matrix_big() = error("Not implemented")
@inline function inv2x2t!(J::Matrix{Float64})
    d = det2x2(J)
    J[1,1], J[2,2] = J[2,2]/d, J[1,1]/d
    J[1,2], J[2,1] = -J[2,1]/d, -J[1,2]/d
    return J
end

@inline function det2x2(J::Matrix{Float64})
    d = J[1,1]*J[2,2] - J[1,2]*J[2,1]
    return d
end


#=
function dNdxmatrix{T <: AbstractInterpolator}(interp::T, local_coords::Point2,
                    vertices::Vector{Int}, nodes::Vector{FENode2})

        dN = dNmatrix(interp, local_coords)
        J = Jmatrix(interp, vertices, nodes, dN)
        dNdx = dN * (inv(J)')

        return dNdx
end



function Jmatrix{T <: AbstractInterpolator}(::T,
                 vertices::Vector{Int}, nodes::Vector{FENode2},
                 dN::Matrix{Float64})

    J = zeros(2, 2)

    for row in 1:size(dN, 1)
        x = nodes[vertices[row]].coords[1]
        y = nodes[vertices[row]].coords[2]

        J[1, 1] += dN[row, 1] * x
        J[1, 2] += dN[row, 1] * y
        J[2, 1] += dN[row, 2] * x
        J[2, 2] += dN[row, 2] * y
    end

    return J
end
=#
