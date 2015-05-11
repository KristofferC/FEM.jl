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
               dof1_n = 0
        for vertex1 in element.vertices
            for dof1 in nodes[vertex1].dofs
                dof1_n += 1
                if dof1.active
                    dof2_n = 0
                    for vertex2 in element.vertices
                        for dof2 in nodes[vertex2].dofs
                            dof2_n += 1
                            if dof2.active
                                push!(colptrs, get_colptr_LUT(K, dof1.eq_n, dof2.eq_n))
                            end
                        end
                    end
                end
            end
        end
    end
    return z
end

#Ke = stiffness(element, nodes, mat)


#=using Base.Order

immutable FEMSParseMatrix{T, Ti}
    K::SparseMatrixCSC{T,Ti}
    eqn2nzval::Vector{Int}
end

@inline get_mat(K::FEMSParseMatrix) = K.K
@inline addval(K::FEMSParseMatrix, eqn::Int, val) = K.K.nzval[K.eqn2nzval[eqn]] += val

function gen_sparse(fp::FEProblem)
    dof_rows = Int[]
    dof_cols = Int[]
    for section in fp.sections
        gen_sparse(section, fp.nodes, dof_rows, dof_cols)
    end

    K = Base.sparse(dof_rows, dof_cols, ones(length(dof_rows)), fp.n_eqs, fp.n_eqs)

    for section in fp.sections
        gen_eq2nzval(section, fp.nodes, dof_rows, dof_cols)
    end

    eqn2nzval =
    for element in section.elements
        for (dof1, i) in activedofs(element, nodes)
            for (dof2, j) in activedofs(element, nodes)
                idx = get_colptr_LUT(K, dof1.eq_n, dof2.eq_n)
                push!(colptrs, g)
                #z += 1
            end
        end
    end
    return z
end

function
    #colptrs = zeros(Int, length(K.nzval))
    colptrs = Int[]
    z = 1
    for section in fp.sections
        z = get_colptrs(K, z, colptrs, section, fp.nodes)
    end
    return colptrs
end


function gen_sparse(section::FESection, nodes::Vector{FENode2},
                   dof_rows::Vector{Int}, dof_cols::Vector{Int})
    for element in section.elements
        for (dof1, _) in activedofs(element, nodes)
            for (dof2, _) in activedofs(element, nodes)
                push!(dof_rows, dof1.eq_n)
                push!(dof_cols, dof2.eq_n)
            end
        end
    end
end

function get_nzval{T,Ti}(A::SparseMatrixCSC{T,Ti}, i0::Integer, i1::Integer)
    i0 = convert(Ti, i0)
    i1 = convert(Ti, i1)
    if !(1 <= i0 <= A.m && 1 <= i1 <= A.n); throw(BoundsError()); end
    r1 = Int(A.colptr[i1])
    r2 = Int(A.colptr[i1+1]-1)

    i = (r1 > r2) ? r1 : searchsortedfirst(A.rowval, i0, r1, r2, Forward)

    return i
end

=#