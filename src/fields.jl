abstract AbstractField

type ScalarField <: AbstractField
    name::UTF8String
    data::Vector{Float64}
    n::Int
end

type VectorField <: AbstractField
    name::UTF8String
    data::Vector{Vector2}
    n::Int
end

type TensorField <: AbstractField
    name::UTF8String
    data::Vector{Vector6}
    n::Int
end

type ExportData
    scalarfields::Dict{Symbol, ScalarField}
    vectorfields::Dict{Symbol, VectorField}
    tensorfields::Dict{Symbol, TensorField}
end


push!(strain, export_data.vectorfields[:strain][node.n] =