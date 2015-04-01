#https://github.com/JuliaLang/julia/issues/10694#issuecomment-88168907
function sparse_plus{Tv,Ti<:Integer}(I::AbstractVector{Ti}, J::AbstractVector{Ti},
                                V::AbstractVector{Tv},
                                nrow::Integer, ncol::Integer)

    if length(I) == 0; return spzeros(eltype(V),nrow,ncol); end
    N = length(I)
    ((N == length(J)) && (N == length(V))) || throw(BoundsError())

    # Work array
    Wj = Array(Ti, max(nrow,ncol)+1)

    # Allocate sparse matrix data structure
    # Count entries in each row
    Rnz = zeros(Ti, nrow+1)
    Rnz[1] = 1
    nz = 0
    for k=1:N
        if V[k] != 0
            Rnz[I[k]+1] += 1
            nz += 1
        end
    end
    Rp = cumsum(Rnz)
    Ri = Array(Ti, nz)
    Rx = Array(Tv, nz)

    # Construct row form
    # place triplet (i,j,x) in column i of R
    # Use work array for temporary row pointers
    @simd for i=1:nrow; @inbounds Wj[i] = Rp[i]; end

    @inbounds for k=1:N
        iind = I[k]
        jind = J[k]
        ((iind > 0) && (jind > 0)) || throw(BoundsError())
        ((iind <= nrow) && (jind <= ncol)) || throw(BoundsError())
        p = Wj[iind]
        Vk = V[k]
        if Vk != 0
            Wj[iind] += 1
            Rx[p] = Vk
            Ri[p] = jind
        end
    end

    # Reset work array for use in counting duplicates
    @simd for j=1:ncol; @inbounds Wj[j] = 0; end

    # Sum up duplicates and squeeze
    anz = 0
    @inbounds for i=1:nrow
        p1 = Rp[i]
        p2 = Rp[i+1] - 1
        pdest = p1

        for p = p1:p2
            j = Ri[p]
            pj = Wj[j]
            if pj >= p1
                Rx[pj] = Rx[pj] + Rx[p]
            else
                Wj[j] = pdest
                if pdest != p
                    Ri[pdest] = j
                    Rx[pdest] = Rx[p]
                end
                pdest += 1
            end
        end

        Rnz[i] = pdest - p1
        anz += (pdest - p1)
    end

    # Transpose from row format to get the CSC format
    RiT = Array(Ti, anz)
    RxT = Array(Tv, anz)

    # Reset work array to build the final colptr
    Wj[1] = 1
    @simd for i=2:(ncol+1); @inbounds Wj[i] = 0; end
    @inbounds for j = 1:nrow
        p1 = Rp[j]
        p2 = p1 + Rnz[j] - 1
        for p = p1:p2
            Wj[Ri[p]+1] += 1
        end
    end
    RpT = cumsum(Wj[1:(ncol+1)])

    # Transpose
    @simd for i=1:length(RpT); @inbounds Wj[i] = RpT[i]; end
    @inbounds for j = 1:nrow
        p1 = Rp[j]
        p2 = p1 + Rnz[j] - 1
        for p = p1:p2
            ind = Ri[p]
            q = Wj[ind]
            Wj[ind] += 1
            RiT[q] = j
            RxT[q] = Rx[p]
        end
    end

    return SparseMatrixCSC(nrow, ncol, RpT, RiT, RxT)
end