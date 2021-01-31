"""
This is the Butz algorithm, as presented by Lawder. Haverkort2017
says it is face continuous. The code is in lawder.c.
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


axis_type(gg::FaceContinuous) = large_enough_unsigned(gg.b)
g_mask(::Type{A}, n) where {A} = [one(A) << i for i = (n-1):-1:0]
g_mask(::Type{A}, n, i) where {A} = one(A) << (n-i)


function calc_P(gg::FaceContinuous, i, H::Vector{A}) where {A}
    element = i ÷ gg.b
    P = H[element + 1]
    if (i % gg.b) > (gg.b - gg.n)
        temp1 = H[element + 2]
        P >>= i % gg.b
        temp1 <<= gg.b - i % gg.b
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
    P = S & g_mask(A, gg.n, 1)
    for i = 2:gg.n
        if (S & g_mask(A, gg.n, i) ⊻ (P >> 1) & g_mask(A, gg.n, i)) != zero(A)
            P |= g_mask(A, gg.n, i)
        end
    end
    P
end


function calc_J(gg::FaceContinuous, P::A) where {A}
    J = A(gg.n)
    i = one(A)
    while i < gg.n
        if ((P >> i) & one(P)) != (P & one(P))
            break
        end
        i += 1
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
        (P - one(A)) ⊻ (P - one(A)) ÷ A(2)
    else
        (P - A(2)) ⊻ (P - A(2)) ÷ A(2)
    end
end


function calc_tS_tT(gg::FaceContinuous, xJ::A, val::A) where {A}
    if (xJ % gg.n) != 0
        temp1 = val >> (xJ % gg.n)
        temp2 = val << (gg.n - (xJ % gg.n))
        retval = temp1 + temp2
        retval & (one(A) << gg.n) - one(A)
    else
        val
    end
end


function H_decode!(gg::FaceContinuous, H::Vector{A}, pt::Vector{A}) where {A}
    mask = one(A) << (gg.b - 1)
    W = zero(A)
    pt .= zero(A)
    i = gg.b * gg.n - gg.n
    P = calc_P(gg, i, H)
    J = calc_J(gg, P)
    xJ = J - one(A)
    tS = brgc(P)
    S = tS
    A1 = S
    T = calc_T(P)
    tT = T
    j = gg.n
    while P > 0
        if P & one(A)
            pt[j] |= mask
        end
        P >>= 1
        j -= 1
    end
    i -= gg.n
    mask >>= 1
    while i >= 0
        P = calc_P(gg, i, H)
        S = brgc(P)
        tS = calc_tS_tT(gg, xJ, S)
        W ⊻= tT
        A1 = W ⊻ tS
        j = gg.n  # one-based
        while A1 > 0
            if (A1 & one(A)) != 0
                pt[j] |= mask
            end
            A1 >>= 1
            j -= 1
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
    mask = one(A) << (gg.b - 1)
    W = zero(A)
    P = zero(A)
    h .= zero(A)
    i = gg.n * gg.b - gg.b
    A1 = zero(A)
    for j = 1:gg.n
        if (pt[j] & mask) != 0
            A1 |= g_mask(A, gg.n, i)
        end
    end
    tS = A1
    S = A1
    P = calc_P2(gg, S)
    element = i ÷ gg.b
    if (i % gg.b) > (gg.b - gg.n)
        h[element + 1] |= (P << (i % gg.b))
        h[element + 2] |= P >> (gg.b - (i % gg.b))
    else
        h[element + 1] |= (P << (i - element * gg.b))
    end
    J = calc_J(gg, P)
    xJ = J - one(A)
    T = calc_T(P)
    tT = T
    i -= gg.n
    mask >>= 1
    while i >= 0
        A1 = 0
        for j = 1:gg.n
            if pt[j] & mask != 0
                A1 |= g_mask(A, gg.n, i)
            end
            W ⊻= tT
            tS = A1 ⊻ W
            S = calc_tS_tT(gg, xJ, tS)
            P = calc_P2(gg, S)
            element = i ÷ gg.b
            if (i % gg.b) > (gg.b - gg.n)
                h[element + 1] |= (P << (i % gg.b))
                h[element + 2] |= (P >> (gg.b - (i % gg.b)))
            else
                h[element + 1] |= P << (i - element * order)
            end
            if i > 0
                T = calc_T(P)
                tT = calc_tS_tT(gg, xJ, T)
                J = calc_J(gg, P)
                xJ += J - one(A)
            end
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
