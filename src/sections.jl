
immutable FESection{P <: AbstractFElement, T <: AbstractMaterial}
    elements::Vector{P}
    material::T
end

function FESection{P <: AbstractFElement, T <: AbstractMaterial,
                   Q <: AbstractMaterialStatus}(mat::T, ele_type::Type{P}, matstat_type::Type{Q})
    elements = ele_type{matstat_type}[]
    FESection(elements, mat)
end


immutable MaterialSection{P <: AbstractMaterial}
    material::P
    elements::Set{Int}
end
MaterialSection{P <: AbstractMaterial}(mat::P) = MaterialSection(mat, Set{Int}())


immutable ElementSection{P <: AbstractFElement}
    element_type::Type{P}
    elements::Set{Int}
end

ElementSection{P <: AbstractFElement}(ele_type::Type{P}) = ElementSection{P}(ele_type, Set{Int}())


# TODO: Clean this up
function push!(section::FESection, elem::AbstractFElement)
    push!(section.elements, elem)
end


function push!{P <: AbstractMaterial}(section::MaterialSection{P}, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end

function push!{P <: AbstractFElement}(section::ElementSection{P}, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end
