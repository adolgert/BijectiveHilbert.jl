
abstract type HilbertAlgorithm{T} end


index_type(::HilbertAlgorithm{T}) where {T} = T


"""
    encode_hilbert(ha::HilbertAlgorithm{T}, X::Vector{A})

A 1-based Hilbert encoding. Both the Hilbert index and the axes start counting
at 1 instead of 0.
"""
function encode_hilbert(gg::HilbertAlgorithm{T}, X::Vector{A}) where {A, T}
    encode_hilbert_zero(gg, X .- one(A)) + one(T)
end


"""
    decode_hilbert!(ha::HilbertAlgorithm{T}, X::Vector{A})

A 1-based Hilbert decode, from [`decode_hilbert_zero!`](@ref). Both the Hilbert
index and the axes start counting at 1 instead of 0.
"""
function decode_hilbert!(gg::HilbertAlgorithm{T}, X::Vector{A}, h::T) where {A,T}
    decode_hilbert_zero!(gg, X, h - one(T))
    X .+= one(A)
end
