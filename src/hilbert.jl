module Bijective
import ..BijectiveHilbert: encode_hilbert_zero, decode_hilbert_zero!, HilbertAlgorithm, axis_type

struct Simple2D{T} <: HilbertAlgorithm{T}
end


Simple2D(::Type{T}) where {T} = Simple2D{T}()


axis_type(::Simple2D{T}) where {T} = Int


# Columns are odd resolution, even resolution (rmin)
# Columns are quadrants in x and y like you'd graph it.
# 1 2
# 0 3
hilbert_variant_encode_table = Function[
    ((x, y, w) -> (x, y))          ((x, y, w) -> (x, y))
    ((x, y, w) -> (y-w, x))        ((x, y, w) -> (y, x-w))
    ((x, y, w) -> (y-w, x-w))      ((x, y, w) -> (y-w, x-w))
    ((x, y, w) -> ((w<<1)-x-one(w), w-y-one(w))) ((x, y, w) -> (w-x-one(w), (w<<1)-y-one(w)))
]


"""
    encode_hilbert_zero(::Simple2D{T}, X::Vector{A})

Computes an integer Hilbert index for x and y using a variant algorithm.

Given two integer indices for a 2-dimensional plane, return a single index.
This index is designed to increase locality for 1-dimensional access.
It does this by keeping nearby points in 2 dimensions also nearby in
1 dimension.

`x` and `y` need to be integers that have bit-shifting operations.

The variant algorithm used differs from the usual Hilbert code because it
doesn't need to know the size of the whole grid before computing the code [^1].
It looks like a slightly-rotated version of the Hilbert curve, but it
has the benefit that it is 1-1 between `(x, y)` and `z`, so you can translate
back and forth.

This function is zero-based. `0 <= x < 2^n`, `0 <= y < 2^n`, and the result
is `0 <= z < 4^n`.

See also: [`decode_hilbert_zero`](@ref), [`encode_hilbert`](@ref).

[^1]: N. Chen, N. Wang, B. Shi, A new algorithm for encoding and decoding the Hilbert order. Software—Practice and Experience 2007; 37(8): 897–908.
"""
function encode_hilbert_zero(::Simple2D{T}, X::Vector{A}) where {A, T}
    x = X[1]
    y = X[2]
    z = zero(T)
    if x == zero(A) && y == zero(A)
        return z
    end
    rmin = convert(Int, floor(log2(max(x, y))) + 1)
    w = one(A)<<(rmin - 1)
    while rmin > 0
        if rmin&1 == 1  # odd
            quadrant = x < w ? (y < w ? 0 : 1) : (y < w ? 3 : 2)
            parity = 1
        else  # even
            quadrant = x < w ? (y < w ? 0 : 3) : (y < w ? 1 : 2)
            parity = 2
        end
        z = (z<<2) + T(quadrant)
        x, y = hilbert_variant_encode_table[quadrant+1, parity](x, y, w)
        rmin -= 1
        w >>= 1
    end
    z
end


hilbert_variant_decode_table = Function[
    ((x, y, w) -> (x, y))          ((x, y, w) -> (x, y))
    ((x, y, w) -> (y, x+w))        ((x, y, w) -> (y+w, x))
    ((x, y, w) -> (y+w, x+w))      ((x, y, w) -> (y+w, x+w))
    ((x, y, w) -> ((w<<1)-x-one(w), w-y-one(w))) ((x, y, w) -> (w-x-one(2), (w<<1)-y-one(w)))
]


"""
    decode_hilbert_zero(z::Integer) -> (x, y)

Computes the (x, y) from a Hilbert code.

This function is zero-based. `0 <= x < 2^n`, `0 <= y < 2^n`, and
`0 <= z < 4^n`.

See also: [`encode_hilbert_zero`](@ref), [`decode_hilbert`](@ref).
"""
function decode_hilbert_zero!(::Simple2D{T}, X::Vector{A}, z::T) where {A,T}
    r = z & T(3)
    x, y = A.([(0, 0), (0, 1), (1, 1), (1, 0)][r + 1])
    z >>= 2
    rmin = 2
    w = one(z) << 1
    while z > zero(T)
        r = z & T(3)
        parity = 2 - rmin&1
        x, y = hilbert_variant_decode_table[r+1, parity](x, y, w)
        z >>= 2
        rmin += 1
        w <<= 1
    end
    X[1] = x
    X[2] = y
end

end
