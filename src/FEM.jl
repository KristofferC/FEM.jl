module FEM

    include("core\\dof.jl")
    include("core\\gauss_quadrature.jl")
    include("core\\mesh.jl")

    include("elements\\element.jl")
    include("interpolators\\interpolator.jl")
    include("materials\\material.jl")


end