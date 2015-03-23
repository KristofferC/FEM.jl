immutable BilinQuadInterp <: AbstractInterpolator end

# Shape functions in local coords
function Nvec(::BilinQuadInterp, loc_coords::Vector{Float64})

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = [(1.0 + ξ) * (1.0 + η) * 0.25;
         (1.0 - ξ) * (1.0 + η) * 0.25;
         (1.0 - ξ) * (1.0 - η) * 0.25;
         (1.0 + ξ) * (1.0 - η) * 0.25]

    return N
end


function dNmatrix(::BilinQuadInterp, loc_coords::Vector{Float64})
    ξ = loc_coords[1]
    η = loc_coords[2]

    # Derivative w.r.t
    #              ξ                     η
    dN = [[ 0.25 * (1.0 + η)     0.25 * (1.0 + ξ)];
          [-0.25 * (1.0 + η)     0.25 * (1.0 - ξ)];
          [-0.25 * (1.0 - η)    -0.25 * (1.0 - ξ)];
          [ 0.25 * (1.0 - η)    -0.25 * (1.0 + ξ)]]

    return dN
end



