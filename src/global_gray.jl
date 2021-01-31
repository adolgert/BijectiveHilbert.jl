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
    for io = n:-1:2
        X[io] ⊻= X[io - 1]
    end
    X[1] ⊻= t
    # Undo excess work
    Q = T(2)
    while Q != N
        P = Q - one(T)
        for jo = n:-1:1
            if (X[jo] & Q) != zero(T)
                X[1] ⊻= P  # invert
            else  # exchange
                t = (X[1] ⊻ X[jo]) & P
                X[1] ⊻= t
                X[jo] ⊻= t
            end
        end
        Q <<= 1
    end
end


function axes_to_transpose!(X::Vector{T}, b, n) where {T <: Unsigned}
    M = one(T) << (b - 1)
    # Inverse undo
    Q = M
    while Q > one(T)
        P = Q - one(T)
        for io = 1:n
            if (X[io] & Q) != zero(T)
                X[1] ⊻= P
            else
                t = (X[1] ⊻ X[io]) & P
                X[1] ⊻= t
                X[io] ⊻= t
            end
        end
        Q >>= 1
    end
    # Gray encode
    for jo = 2:n
        X[jo] ⊻= X[jo - 1]
    end
    t2 = zero(T)
    Q = M
    while Q > one(T)
        if (X[n] & Q) != 0
            t2 ⊻= (Q - one(T))
        end
        Q >>= one(T)
    end
    for ko = 1:n
        X[ko] ⊻= t2
    end
end

"""
Takes a vector of length `n` and places the bits of all `n` integers
into a single integer. The vector's 1st component is the most significant bit.
"""
function interleave_transpose(X::Vector{T}, b, n) where {T <: Integer}
    h = zero(UInt64)
    for i in 0:(b - 1)
        for d in 1:n
            ith_bit = (X[d] & (1<<i)) >> i
            h |= ith_bit << (i * n + n - d)
        end
    end
    h
end


"""
Takes a single integer and places its values into components of a vector,
bit-by-bit.
"""
function outerleave_transpose!(X::Vector{T}, h, b, n) where {T <: Integer}
    X .= zero(T)
    for i in 0:(b-1)
        for d in 1:n
            ith_bit = (h & (1 << (i * n + n - d))) >> (i * n + n - d)
            X[d] |= ith_bit << i
        end
    end
end


"""
Takes a vector of length `n` and places the bits of all `n` integers
into a single integer. The vector's 1st component is the least significant bit.
"""
function interleave_transpose_low(X::Vector{T}, b, n) where {T <: Integer}
    h = zero(UInt64)
    for i in 0:(b - 1)
        for d in 1:n
            h |= ((X[d] & (1<<i))) << (i*(n - 1) + d - 1)
        end
    end
    h
end


function outerleave_transpose_low!(X::Vector{T}, h, b, n) where {T <: Integer}
    X .= zero(T)
    for i in 0:(b-1)
        for d in 1:n
            X[d] |= (h & (1 << (i * n + d - 1))) >> (i * (n - 1) + d - 1)
        end
    end
end


struct GlobalGray{T} <: HilbertAlgorithm{T}
    b::Int
    n::Int
end


axis_type(gg::GlobalGray) = large_enough_unsigned(gg.b)


function GlobalGray(b, n)
    ttype = large_enough_unsigned(b * n)
    GlobalGray{ttype}(b, n)
end


function GlobalGray(::Type{T}, b, n) where {T}
    GlobalGray{T}(b, n)
end


function encode_hilbert_zero!(g::GlobalGray{T}, X::Vector)::T where {T}
    axes_to_transpose!(X, g.b, g.n)
    interleave_transpose(X, g.b, g.n)
end


function encode_hilbert_zero(g::GlobalGray{T}, X::Vector)::T where {T}
    Y = copy(X)
    encode_hilbert_zero!(g, Y)
end


function decode_hilbert_zero!(g::GlobalGray{T}, X::Vector, h::T) where {T}
    outerleave_transpose!(X, h, g.b, g.n)
    transpose_to_axes!(X, g.b, g.n)
end
