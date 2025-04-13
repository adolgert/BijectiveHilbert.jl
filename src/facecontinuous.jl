"""
    FaceContinuous(b::Int, n::Int)
    FaceContinuous(::Type{T}, b::Int, n::Int)

For `n` dimensions, use `b` bits of precision in this Hilbert curve.
If you specify a type `T`, this will be used as the type of the Hilbert encoding.
If not, the smallest unsigned integer that can hold `n*b` bits will be used
for the Hilbert index data type.

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

# g_mask[x] can be replaced by (1 << DIM - 1 - x)
g_mask(::Type{A}, n, iz) where {A} = one(A) << (n - iz - 1)


function H_encode!(gg::FaceContinuous{F}, pt::Vector{K}, h::Vector{K}) where {K,F}
    wordbits = 8 * sizeof(K)
    W = zero(K)
    P = zero(K)
    h .= zero(K)
    # Start from the high bit in each element and work to the low bit.
    mask = one(K) << (gg.b - 1)  # lawder puts wordbits here.
    # The H-index needs b * n bits of storage. i points to the base of the current level.
    i = gg.b * gg.n - gg.n
    # A will hold bits from all pt elements, at the current level.
    A = zero(K)
    for j = 0:(gg.n - 1)
        if (pt[j + 1] & mask) != zero(K)
            A |= g_mask(K, gg.n, j)
        end
    end
    S = tS = A
    P |= (S & g_mask(K, gg.n, 0))
    for j = 1:(gg.n - 1)
        gm = g_mask(K, gg.n, j)
        if ((S & gm) ⊻ ((P >> 1) & gm)) != 0
            P |= gm
        end
    end

    # P is the answer, but it has to be packed into a vector.
    element = i ÷ wordbits
    if (i % wordbits) > (wordbits - gg.n)
        h[element + 1] |= (P << (i % wordbits))
        h[element + 2] |= P >> (wordbits - (i % wordbits))
	else
        h[element + 1] |= (P << (i - element * wordbits))
    end

    J = gg.n
    j = 1
    while j < gg.n
        if ((P >> j) & one(K)) == (P & one(K))
            j += 1
			continue
		else
            break
        end
    end
    if j != gg.n
        J -= j
    end
    xJ = J - 1

    if P < K(3)
        T = zero(K)
    else
        if (P % K(2)) != zero(K)
            T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
        else
            T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
        end
    end
    tT = T
    i -= gg.n
    mask >>= 1
    while i >= 0
        # println("i=$i mask=$mask T=$T P=$P")
		# println("h[0]=$(h[1]) h[1]=$(h[2]) h[2]=$(h[3])")
        A = zero(K)
        for j = 0:(gg.n - 1)
            if pt[j + 1] & mask != zero(K)
                A |= g_mask(K, gg.n, j)
            end
        end
        W ⊻= tT
        tS = A ⊻ W
        if xJ % gg.n != 0
            temp1 = tS << (xJ % gg.n)
            temp2 = tS >> (gg.n - (xJ % gg.n))
            S = temp1 | temp2
            S &= (one(K) << gg.n) - one(K)
        else
            S = tS
        end

        P = S & g_mask(K, gg.n, 0)
        for j = 1:(gg.n - 1)
            gn = g_mask(K, gg.n, j)
            if ((S & gn) ⊻ ((P >> 1) & gn)) != 0
                P |= gn
            end
        end
        element = i ÷ wordbits
        if (i % wordbits) > (wordbits - gg.n)
            h[element + 1] |= (P << (i % wordbits))
            h[element + 2] |= (P >> (wordbits - (i % wordbits)))
        else
            h[element + 1] |= P << (i - element * wordbits)
        end
        if i > 0
            if P < K(3)
                T = 0
            else
                if P % K(2) != zero(K)
                    T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
                else
                    T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
                end
            end
            if xJ % gg.n != 0
                temp1 = T >> (xJ % gg.n)
                temp2 = T << (gg.n - (xJ % gg.n))
                tT = temp1 | temp2
                tT &= (one(K) << gg.n) - one(K)
            else
                tT = T
            end

            J = gg.n
            j = 1
            while j < gg.n
                if ((P >> j) & one(K)) == (P & one(K))
                    j += 1
                    continue
                else
                    break
                end
            end
            if j != gg.n
                J -= j
            end
            xJ += J - 1
        end

        i -= gg.n
        mask >>= 1
    end
    # println("h[0]=$(h[1]) h[1]=$(h[2]) h[2]=$(h[3])")
end


mutable struct FCLevel{K}
    mask::K
    i::Int
    n::Int
end


function FCLevel(fc, ::Type{K}) where {K}
    # XXX lawder uses wordbits instead of fc.b
    FCLevel{K}(one(K) << (fc.b - 1), fc.b * fc.n - fc.n, fc.n)
end


function downlevel!(l::FCLevel)
    l.mask >>= 1
    l.i -= l.n
end


function index_at_level(H::Vector{K}, l::FCLevel{K}) where {K}
    wordbits = 8 * sizeof(K)
    element = l.i ÷ wordbits
    P = H[element + 1]
    if (l.i % wordbits) > (wordbits - l.n)
        temp1 = H[element + 2]  # one-based
        P >>= l.i % wordbits
        temp1 <<= wordbits - (l.i % wordbits)
        P |= temp1
    else
        P >>= (l.i % wordbits)
    end
    
    # 	/* the & masks out spurious highbit values */
    if l.n < wordbits
        P &= (one(K) << l.n) - one(K)
    end
    P
end


function distribute_to_coords!(bits::K, axes::Vector{K}, l::FCLevel{K}) where K
    j = l.n - 1
    while bits > 0
        if bits & one(K) != 0
            axes[j + 1] |= l.mask
        end
        bits >>= 1
        j -= 1
    end
end


function fc_parity_match(P, n)
    J = n
    j = 1
    parity = P & one(P)
    while j < n
        if ((P >> j) & one(P)) != parity
            break
        end
        j += 1
    end
    j
end


function fc_rotation(P, n)
    j = fc_parity_match(P, n)
    if j == n
        xJ = n - 1
    else
        xJ = n - j - 1
    end
    xJ
end


function fc_flip(P::K) where {K}
    if P < K(3)
        T = zero(K)
    else
        if P % 2 != 0
            T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
        else
            T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
        end
    end
    T
end


function H_decode!(gg::FaceContinuous{F}, H::Vector{K}, pt::Vector{K}) where {K,F}
    pt .= zero(K)
    l = FCLevel(gg, K)
    P = index_at_level(H, l)
    xJ = fc_rotation(P, gg.n)
    A = S = tS = brgc(P)
    tT = T = fc_flip(P)
    # XXX Lawder's code puts P here instead of A.
    distribute_to_coords!(A, pt, l)

    nmask = (one(K) << gg.n) - one(K)  # mask of n bits.

    W = zero(K)
    downlevel!(l)
    while l.i >= 0
        P = index_at_level(H, l)
        S = brgc(P)
        @assert S & ~nmask == 0
        tS = rotateright(S, xJ % gg.n, gg.n)
    
        W ⊻= tT
        A = W ⊻ tS
        distribute_to_coords!(A, pt, l)
    
        if l.i > 0
            T = fc_flip(P)
            @assert T & ~nmask == 0
            tT = rotateright(T, xJ % gg.n, gg.n)
        
            xJ += fc_rotation(P, gg.n)
        end
        downlevel!(l)
    end
end


function encode_hilbert_zero(fc::FaceContinuous{T}, X::Vector{A})::T where {A,T}
    hvec = zeros(A, fc.n)
    H_encode!(fc, X, hvec)
    # Encoding packs H into a vector of A, using all bits in the A type.
    h = zero(T)
    for i in fc.n:-1:1
        h <<= 8 * sizeof(A)
        h |= hvec[i]
    end
    h
end


function decode_hilbert_zero!(fc::FaceContinuous{T}, X::Vector{A}, h::T) where {A,T}
    # H is in a larger type T but algorithm expects it packed into a vector of A.
    hvec = zeros(A, fc.n)
    for i in 1:fc.n
        hvec[i] |= h & ~zero(A)
        h >>= 8 * sizeof(A)
    end
    H_decode!(fc, hvec, X)
end
