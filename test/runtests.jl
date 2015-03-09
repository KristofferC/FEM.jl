module TestFEM
    using FactCheck
    using FEM

    include("test_mesh.jl")
    include("test_materials.jl")
    include("test_elements.jl")
    include("test_interpolators.jl")
    include("test_feproblem.jl")
end
