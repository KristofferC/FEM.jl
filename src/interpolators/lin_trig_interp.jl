immutable LinTrigInterp <: Interpolator end

# Shape functions in local coordinates
function Nvec(interp::LinTrigInterp, loc_coords::Vector{Float64}, mp::MatPool,
              vp::VecPool)

    ξ = loc_coords[1]
    η = loc_coords[2]

    N = getvec(3, "N", vp)

    N[1] = ξ
    N[2] = η
    N[3] = 1.0 - ξ - η

    return N
end


function dNmatrix(interp::LinTrigInterp, loc_coords::Vector{Float64}, mp::MatPool,
              vp::VecPool)
    # Derivative w.r.t
    #       ξ     η

    dN = getmat(3, 2, "dN", mp)

    dN[1, 1] = 1.0
    dN[1, 2] = 0.0
    dN[2, 1] = 0.0
    dN[2, 2] = 1.0
    dN[3, 1] = -1.0
    dN[3, 2] = -1.0

    return dN
end


