using Base.Order
function get_colptr_LUT{T,Ti}(A::SparseMatrixCSC{T,Ti}, i0::Integer, i1::Integer)
    i0 = convert(Ti, i0)
    i1 = convert(Ti, i1)
    if !(1 <= i0 <= A.m && 1 <= i1 <= A.n); throw(BoundsError()); end
    r1 = Int(A.colptr[i1])
    r2 = Int(A.colptr[i1+1]-1)

    i = (r1 > r2) ? r1 : searchsortedfirst(A.rowval, i0, r1, r2, Forward)

    return i
end

function get_colptrs(K::SparseMatrixCSC, fp::FEProblem)
    #colptrs = zeros(Int, length(K.nzval))
    colptrs = Int[]
    z = 1
    for section in fp.sections
        z = get_colptrs(K, z, colptrs, section, fp.nodes)
    end
    return colptrs
end

function get_colptrs(K::SparseMatrixCSC, z::Int, colptrs::Vector{Int},
                     section::FESection, nodes::Vector{FENode2})
    for element in section.elements
        for (dof1, i) in activedofs(element, nodes)
            for (dof2, j) in activedofs(element, nodes)
                push!(colptrs, get_colptr_LUT(K, dof1.eq_n, dof2.eq_n))
                #z += 1
            end
        end
    end
    return z
end

