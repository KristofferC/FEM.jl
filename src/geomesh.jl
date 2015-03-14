abstract AbstractVertex{L} <: FixedVector{Int, L}

immutable Vertex2 <: AbstractVertex{2}
    v1::Int
    v2::Int
end


immutable Vertex3 <: AbstractVertex{3}
    v1::Int
    v2::Int
    v3::Int
end


immutable Vertex4 <: AbstractVertex{4}
    v1::Int
    v2::Int
    v3::Int
    v4::Int
end
