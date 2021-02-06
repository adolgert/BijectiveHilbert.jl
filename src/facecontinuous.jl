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

# g_mask[x] can be replaced by (1 << DIM - 1 - x)
g_mask(::Type{A}, n, iz) where {A} = one(A) << (n - iz - 1)


function H_encode!(gg::FaceContinuous, pt::Vector{K}, h::Vector{K}) where {K}
    wordbits = 8 * sizeof(K)
    W = zero(K)
    P = zero(K)
    h .= zero(K)
    # Start from the high bit in each element and work to the low bit.
    mask = one(K) << (wordbits - 1)
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
        println("i=$i mask=$mask T=$T P=$P")
		println("h[0]=$(h[1]) h[1]=$(h[2]) h[2]=$(h[3])")
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
    println("h[0]=$(h[1]) h[1]=$(h[2]) h[2]=$(h[3])")
end


function H_decode!(gg::FaceContinuous, H::Vector{K}, pt::Vector{K}) where {K}
    wordbits = 8 * sizeof(K)
    mask = one(K) << (wordbits - 1)
    i = gg.b * gg.n - gg.n
    
    # 	/*--- P ---*/
    element = i ÷ wordbits
    P = H[element + 1]
    if (i % wordbits) > (wordbits - gg.n)
        temp1 = H[element + 2]  # one-based
        P >>= i % wordbits
        temp1 <<= wordbits - (i % wordbits)
        P |= temp1
    else
        P >>= (i % wordbits)
    end
    
    # 	/* the & masks out spurious highbit values */
    if gg.n < wordbits
        P &= (one(K) << gg.n) - one(K)
    end
    
    # 	/*--- xJ ---*/
    J = gg.n
    j = 1
    while j < gg.n
        if ((P >> j) & one(K)) != (P & one(K))
            break
        end
        j += 1
    end
    if j != gg.n
        J -= j
    end
    xJ = J - 1
    
    # 	/*--- S, tS, A ---*/
    A = S = tS = P ⊻ (P >> 1)    
    
    # 	/*--- T ---*/
    if P < K(3)
        T = 0
    else
        if P % 2 != 0
            T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
        else
            T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
        end
    end
    
    # 	/*--- tT ---*/
    tT = T
    
    # 	/*--- distrib bits to coords ---*/
    j = gg.n - 1
    while P > zero(K)
        if P & one(K) != 0
            pt[j + 1] |= mask
        end
        P >>= 1
        j -= 1
    end
    
    W = zero(K)
    i -= gg.n
    mask >>= 1
    while i >= 0
        println("i=$i mask=$mask T=$T P=$P, xJ=$xJ, tT=$tT, W=$W")
		println("pt[0]=$(pt[1]) pt[1]=$(pt[2]) pt[2]=$(pt[3])")
        # 		/*--- P ---*/
        element = i ÷ wordbits
        P = H[element + 1]
        if (i % wordbits) > (wordbits - gg.n)
            temp1 = H[element + 2]
            P >>= i % wordbits
            temp1 <<= wordbits - (i % wordbits)
            P |= temp1
        else
            P >>= i % wordbits
        end
    
    # 		/* the & masks out spurious highbit values */
        if gg.n < wordbits
            P &= (one(K) << gg.n) - one(K)
        end 
    
        # 		/*--- S ---*/
        S = P ⊻ (P >> 1)
    
        # 		/*--- tS ---*/
        if xJ % gg.n != 0
            temp1 = S >> (xJ % gg.n)
            temp2 = S << (gg.n - (xJ % gg.n))
            tS = temp1 | temp2
            tS &= (one(K) << gg.n) - 1
        else
            tS = S
        end
    
        # 		/*--- W ---*/
        W ⊻= tT
    
        # 		/*--- A ---*/
        A = W ⊻ tS
    
        # 		/*--- distrib bits to coords ---*/
        j = gg.n - 1
        while A > 0
            if A & one(K) != 0
                pt[j + 1] |= mask
            end
            A >>= 1
            j -= 1
        end
    
        if i > 0
            # 			/*--- T ---*/
            if P < K(3)
                T = 0
            else
                if P % 2 != 0
                    T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
                else
                    T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
                end
            end
    
            # 			/*--- tT ---*/
            if xJ % gg.n != 0
                temp1 = T >> (xJ % gg.n)  # C precedence
                temp2 = T << (gg.n - (xJ % gg.n))
                tT = temp1 | temp2
                tT &= (one(K) << gg.n) - one(K)
            else
                tT = T
            end
    
            # 			/*--- xJ ---*/
            J = gg.n
            j = 1
            while j < gg.n
                if ((P >> j) & one(K)) != P & one(K)
                    break
                end
                j += 1
            end
            println("j=$j P=$P J=$J")
            if j != gg.n
                J -= j
            end
            xJ += J - 1
        end
        i -= gg.n
        mask >>= 1
    end
end


function encode_hilbert_zero(fc::FaceContinuous{T}, X::Vector{A})::T where {A,T}
    hvec = zeros(A, fc.n)
    H_encode!(fc::FaceContinuous, X, hvec)
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
    H_decode!(fc::FaceContinuous, hvec, X)
end
