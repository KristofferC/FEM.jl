immutable LinTrigInterp <: Interpolator end

# Shape functions in local coordinates
function Nvec(interp::LinTrigInterp, loc_coords::Vector{Float64})

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = [ξ;
         η;
         1.0 - ξ - η]

    return N
end


function dNmatrix(interp::LinTrigInterp, loc_coords::Vector{Float64})
    # Derivative w.r.t
    #       ξ     η
    dN = [[ 1.0  0.0];
          [ 0.0  1.0];
          [-1.0 -1.0]]

    return dN
end


