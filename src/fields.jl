



abstract AbstractField

# Tensor fields (right now 6 components)
abstract AbstractTensor <: AbstractField
get_ncomponents{T <: AbstractTensor}(::Type{T}) = 6

# Vector fields
abstract AbstractVector <: AbstractField
get_ncomponents{T <: AbstractVector}(::Type{T}) = 3

# Scalar fields
abstract AbstractScalar <: AbstractField
#get_ncomponents(::Type{AbstractTensor}) = 1 # Not needed for scalars?

# TODO: Maybe we need a distinction between primary and secondary fields?
# Exporting a primary field as cell_data makes no sense? Or does it?

# Tensors
type Strain <: AbstractTensor end
type Stress <: AbstractTensor end

# Vectors
type Displacement <: AbstractVector end

# Scalars
type Pressure <: AbstractScalar end
