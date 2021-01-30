# The GlobalGray algorithm. This follows an article,
# "Programming the Hilbert Curve," by John Skilling.
# 707 (2004), http://dx.doi.org/10.1063/1.1751381.
# I call it "global Gray" because the insight of the article
# is that a single, global Gray code can be applied to all
# np bits of a Hilbert length.

"""
    transpose_to_axes(coord, b, n)

Convert from Hilbert index to a coordinatefor `b` bits
in `n` dimensions.
"""
function transpose_to_axes!(X::Vector{T}, b, n) where {T <: Unsigned}
    N = T(2) << (b - 1)
    # Gray decode by H^(H/2)
    t = X[n] >> 1
    for i = n:-1:2
        X[i] ⊻= X[i - 1]
    end
    X[1] ⊻= t
    # Undo excess work
    Q = T(2)
    while Q != N
        P = Q - one(T)
        for i = n:-1:1
            if (X[i] & Q) != zero(T)  # invert
                X[1] ⊻= P
            else  # exchange
                t = (X[1] ⊻ X[i]) & P
                X[1] ⊻= t
                X[i] ⊻= t
            end
        end
        Q <<= one(T)
    end
end


function axes_to_transpose!(X::Vector{T}, b, n) where {T <: Unsigned}
    M = one(T) << (b - 1)
    # Inverse undo
    Q = M
    while Q > one(T)
        P = Q - one(T)
        for i = 1:n
            if (X[i] & Q) != zero(T)
                X[1] ⊻= P
            else
                t = (X[1] ⊻ X[i]) & P
                X[1] ⊻= t
                X[i] ⊻= t
            end
        end
        Q >>= one(T)
    end
    # Gray encode
    for i = 2:n
        X[i] ⊻= X[i - 1]
    end
    t2 = zero(T)
    Q = M
    while Q > one(T)
        if (X[n] & Q) != 0
            t2 ⊻= (Q - one(T))
        end
        Q >>= one(T)
    end
    for i = 1:n
        X[i] ⊻= t2
    end
end


function interleave_transpose(X::Vector{T}, b, n) where {T <: Unsigned}
    h = zero(UInt64)
    for i in 0:(b - 1)
        for d in 1:n
            h |= ((X[d] & (1<<i))) << (i*(n - 1) + d - 1)
        end
    end
    h
end


function outerleave_transpose!(X::Vector{T}, h, b, n) where {T <: Unsigned}
    X .= zero(T)
    for i in 0:(b-1)
        for d in 1:n
            X[d] |= (h & (1 << (i * n + d - 1))) >> (i * (n - 1) + d - 1)
        end
    end
end


struct GlobalGray{A,T} <: HilbertAlgorithm{A,T}
    b::Int
    n::Int
end


function GlobalGray(b, n)
    atype = large_enough_unsigned(b)
    ttype = large_enough_unsigned(b * n)
    GlobalGray{atype, ttype}(b, n)
end


function encode_hilbert_zero!(g::GlobalGray{A,T}, X::Vector{A})::T where {A,T}
    axes_to_transpose!(X, g.b, g.n)
    interleave_transpose(X, g.b, g.n)
end


function encode_hilbert_zero(g::GlobalGray{A,T}, X::Vector{A})::T where {A,T}
    Y = copy(X)
    encode_hilbert_zero(g, Y)
end


function decode_hilbert_zero!(g::GlobalGray{A,T}, X::Vector{A}, h::T) where {A,T}
    outerleave_transpose!(X, h, g.b, g.n)
    transpose_to_axes!(X, g.b, g.n)
end
