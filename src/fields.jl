abstract AbstractField

# Tensor fields (right now 6 components)
abstract AbstractTensor <: AbstractField
get_ncomponents{T <: AbstractTensor}(::Type{T}) = 6

# Vector fields
abstract AbstractVector <: AbstractField
get_ncomponents{T <: AbstractVector}(::Type{T}) = 3

# Scalar fields
abstract AbstractScalar <: AbstractField
get_ncomponents{T <: AbstractScalar}(::Type{T}) = 1
#get_ncomponents(::Type{AbstractTensor}) = 1 # Not needed for scalars?

abstract AbstractFullTensor <: AbstractField
get_ncomponents{T <: AbstractFullTensor}(::Type{T}) = 9

type KAlpha <: AbstractField end
get_ncomponents(::Type{KAlpha}) = 2

type Kappa <: AbstractField end
get_ncomponents(::Type{Kappa}) = 2

type τ <: AbstractField end
get_ncomponents(::Type{τ}) = 2
# TODO: Maybe we need a distinction between primary and secondary fields?
# Exporting a primary field as cell_data makes no sense? Or does it?

# Tensors
type Strain <: AbstractTensor end
type Stress <: AbstractTensor end

type VonMises <: AbstractScalar end

# Vectors
type Displacement <: AbstractVector end

# Scalars
type Pressure <: AbstractScalar end

type FullStress <: AbstractFullTensor end
type FullStrain <: AbstractFullTensor end
type InvFp <: AbstractFullTensor end
type FullStrain <: AbstractFullTensor end

