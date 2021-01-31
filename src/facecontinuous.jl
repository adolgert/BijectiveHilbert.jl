"""
This is the Butz algorithm, as presented by Lawder. Haverkort2017
says it is face continuous. The code is in lawder.c.
The original paper had an error, and Lawder put a correction on his website.
http://www.dcs.bbk.ac.uk/~jkl/publications.html
"""
struct FaceContinuous{T} <: HilbertAlgorithm{T}
    b::Int
    n::Int
end


function FaceContinuous(b::Int, n::Int)
    T = large_enough_unsigned(b * n)
    FaceContinuous{T}(b, n)
end


FaceContinuous(::Type{T}, b::Int, n::Int) where {T} = FaceContinuous{T}(b, n)

"""
Used during testing to pick a type for the xyz coordinate.
"""
axis_type(gg::FaceContinuous) = large_enough_unsigned(gg.b)

# Mark one-based variables with a trailing o.
# Mark zero-based variables with a trailing z.
# The C code takes place in vectors, without the packed Hilbert index.
# ORDER -> b  The number of bits.
# DIM -> n  The number of dimensions.
# U_int -> A, the type for the Axis coordinate.
# C translation rules:
# * / -> ÷
# * ^ -> ⊻
# * % -> %  This is the same choice. mod() is wrong for negative numbers.
# Julia has different operator precedence among addition, subtraction, bit shifts, and and xor.
# Lawder knew C operator precedence well, so there are few parentheses. Add them liberally.
# https://en.cppreference.com/w/c/language/operator_precedence. Relevant ones:
# (++, --), (* / %), (+ -), (<< >>), (< >), (== !=), (&), (^), (= += <<= &= ^= |=)

# zero-based i.
g_mask(::Type{A}, n, iz) where {A} = one(A) << (n - iz - 1)


function calc_P(gg::FaceContinuous, i, H::Vector{A}) where {A}
    elementz = i ÷ gg.b
    P = H[elementz + 1]
    if (i % gg.b) > (gg.b - gg.n)
        temp1 = H[elementz + 2]
        P >>= (i % gg.b)
        temp1 <<= (gg.b - (i % gg.b))
        P |= temp1
    else
        P >>= (i % gg.b)
    end
    if gg.n < gg.b
        P &= (one(A) << gg.n) - one(A)
    end
    P
end


function calc_P2(gg::FaceContinuous, S::A) where {A}
    P = S & g_mask(A, gg.n, 0)
    for i = 1:(gg.n - 1)
        if ((S & g_mask(A, gg.n, i)) ⊻ ((P >> 1) & g_mask(A, gg.n, i))) != zero(A)
            P |= g_mask(A, gg.n, i)
        end
    end
    P
end


function calc_J(gg::FaceContinuous, P::A) where {A}
    J = A(gg.n)
    i = one(A)
    while i < gg.n
        if ((P >> i) & one(A)) != (P & one(A))
            break
        end
        i += one(A)
    end
    if i != gg.n
        J -= i
    end
    J
end


function calc_T(P::A) where {A}
    if P < A(3)
        0
    elseif (P % 2) == 0
        (P - one(A)) ⊻ ((P - one(A)) << 1)
    else
        (P - A(2)) ⊻ ((P - A(2)) << 1)
    end
end


function calc_tS_tT(gg::FaceContinuous, xJ, val::A) where {A}
    if (xJ % gg.n) != zero(A)
        temp1 = val >> (xJ % gg.n)
        temp2 = val << (gg.n - (xJ % gg.n))
        retval = temp1 + temp2
        retval & ((one(A) << gg.n) - one(A))
    else
        val
    end
end


function H_decode!(gg::FaceContinuous, H::Vector{A}, pt::Vector{A}) where {A}
    wordbits = 8 * sizeof(A)
    mask = one(A) << (gg.b - 1)
    W = zero(A)
    pt .= zero(A)
    i = gg.b * gg.n - gg.n
    P = calc_P(gg, i, H)
    J = calc_J(gg, P)
    xJ = J - one(A)
    tS = P ⊻ (P << 1)  # brgc(P)
    S = tS
    A1 = tS
    T = calc_T(P)
    tT = T
    jo = gg.n
    while P > 0
        if (P & one(A)) != 0
            pt[jo] |= mask
        end
        P >>= 1
        jo -= 1
    end
    i -= gg.n
    mask >>= 1
    while i >= 0
        P = calc_P(gg, i, H)
        S = P ⊻ (P << 1)  # brgc(P)
        tS = calc_tS_tT(gg, xJ, S)
        W ⊻= tT
        A1 = W ⊻ tS
        jo = gg.n
        while A1 > 0
            if (A1 & one(A)) != 0
                pt[jo] |= mask
            end
            A1 >>= 1
            jo -= 1
        end
        if i > 0
            T = calc_T(P)
            tT = calc_tS_tT(gg, xJ, T)
            J = calc_J(gg, P)
            xJ += J + one(A)
        end
        i -= gg.n
        mask >>= 1
    end
end


function H_encode!(gg::FaceContinuous, pt::Vector{A}, h::Vector{A}) where {A}
    wordbits = 8 * sizeof(A)
    mask = one(A) << (wordbits - 1)
    W = zero(A)
    P = zero(A)
    h .= zero(A)
    i = gg.n * gg.b - gg.b
    A1 = zero(A)
    for jo = 1:gg.n
        if (pt[jo] & mask) != zero(A)
            A1 |= g_mask(A, gg.n, jo - 1)
        end
    end
    tS = A1
    S = A1
    P |= (S & g_mask(A, gg.n, 0))
    for jz = 1:(gg.n - 1)
        gm = g_mask(A, gg.n, j)
        if ((S & gm) ⊻ ((P >> 1) & gm)) != 0
            P |= gm
        end
    end

    elementz = i ÷ wordbits
    if (i % wordbits) > (wordbits - gg.n)
        h[elementz + 1] |= (P << (i % wordbits))
        h[elementz + 2] |= P >> (wordbits - (i % wordbits))
    else
        h[elementz + 1] |= (P << (i - elementz * wordbits))
    end
    J = gg.n
    j = 1
    while j < gg.n
        if ((P >> j) & one(A)) != (P & one(A))
            break
        end
        j += 1
    end
    if j != gg.n
        J -= A(j)
    end
    if P < A(3)
        T = 0
    else
        if (P % 2) != 0
            T = (P - one(A)) ⊻ ((P - one(A)) << 1)
        else
            T = (P - A(2)) ⊻ ((P - A(2)) << 1)
        end
    end
    tT = T
    i -= gg.n
    mask >>= 1
    while i >= 0
        A1 = zero(A)
        for jz = 0:(gg.n - 1)
            if pt[jz + 1] & mask != zero(A)
                A1 |= g_mask(A, gg.n, jz)
            end
        end
        W ⊻= tT
        tS = A1 ⊻ W
        if xJ % gg.n != 0
            temp1 = tS << (xJ % gg.n)
            temp2 = tS >> (gg.n - (xJ % gg.n))
            S = temp1 | temp2
            S &= (one(A) << gg.n) - one(A)
        else
            S = tS
        end

        P = S & g_mask(A, gg.n, 0)
        for jz = 1:(gg.n - 1)
            gn = g_mask(A, gg.n, j)
            if ((S & gn) ⊻ ((P >> 1) & gn)) != 0
                P |= gn
            end
        end
        elementz = i ÷ wordbits
        if (i % wordbits) > (wordbits - gg.n)
            h[elementz + 1] |= (P << (i % wordbits))
            h[elementz + 2] |= (P >> (wordbits - (i % wordbits)))
        else
            h[elementz + 1] |= P << (i - elementz * wordbits)
        end
        if i > 0
            if P < 3
                T = 0
            else
                if P % 2 != 0
                    T = (P - one(A)) ⊻ ((P - 1) << 1)
                else
                    T = (P - A(2)) ⊻ ((P - A(2)) << 1)
                end
            end
            if xJ % gg.n != 0
                temp1 = T << (xJ % gg.n)
                temp2 = T >> (gg.n - (xJ % gg.n))
                S = temp1 | temp2
                S &= (one(A) << gg.n) - one(A)
            else
                tT = T
            end

            J = gg.n
            j = 1
            while j < gg.n
                if ((P >> j) & one(A)) != (P & one(A))
                    break
                end
                j += 1
            end
            if j != gg.n
                J -= A(j)
            end
            xJ += J - one(A)
        end

        i -= gg.n
        mask >>= 1
    end
end


function encode_hilbert_zero(fc::FaceContinuous{T}, X::Vector{A})::T where {A, T}
    h = zeros(A, fc.n)
    H_encode!(fc::FaceContinuous, X, h)
    interleave_transpose(T, X, fc.b, fc.n)
end


function decode_hilbert_zero!(fc::FaceContinuous{T}, X::Vector{A}, h::T) where {A, T}
    H = zeros(A, fc.n)
    outerleave_transpose!(T, H, h, fc.b, fc.n)
    H_decode!(fc::FaceContinuous, H, X)
end
