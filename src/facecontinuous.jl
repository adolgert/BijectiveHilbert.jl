"""
This is the Butz algorithm, as presented by Lawder. Haverkort2017
says it is face continuous.
"""
struct FaceContinuous{T} <: HilbertAlgorithm{T}
    b::Int
    n::Int
end


axis_type(gg::FaceContinuous{T}) = large_enough_unsigned(gg.b)
g_mask(::Type{T}, n) = [one(T) << i for i = (n-1):-1:0]
g_mask(::Type{T}, n, i) = one(T) << (n-i)


function calc_P(gg::FaceContinuous, i, H::Vector{A}) where {A}
    element = i / gg.b
    P = h[element + 1]
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
    for i = 1:(gg.n - 1)
        if ((P >> i) & one(P)) != (P & one(P))
            break
        end
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
        (P - one(A)) ⊻ (P - one(A)) / A(2)
    else
        (P - A(2)) ⊻ (P - A(2)) / A(2)
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


function H_decode!(gg::FaceContinuous, H::Vector{A}, P::Vector{A}) where {A}
    mask = one(A) << (gg.n - 1)
    W = zero(A)
    P .= zero(A)
    i = gg.b * gg.n - gg.n
    P = calc_P(gg, i, H)
    J = calc_J(gg, P)
    xJ = J - one(A)
    tS = brgc(P)
    S = tS
    A = S
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
        P = calc_P(i, H)
        S = brgc(P)
        tS = calc_tS_tT(gg, xJ, S)
        W ⊻= tT
        A1 = W ⊻ tS
        j = gg.n  # one-based
        while A1 > 0
            if A1 & one(A)
                pt[j] |= mask
            end
            A1 >>= 1
            j -= 1
        end
        if i > 0
            T = calc_T(P)
            tT = calc_tS_tT(gg, xJ, T)
            J = calc_J(P)
            xJ += J + one(A)
        end
        i -= gg.n
        mask >>= 1
    end
end


function encode_hilbert_zero(g::FaceContinuous{T}, X::Vector)::T where {T}
    
end


function decode_hilbert_zero!(g::FaceContinuous{T}, X::Vector, h::T) where {T}
    
end
