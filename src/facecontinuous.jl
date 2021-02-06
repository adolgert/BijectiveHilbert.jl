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
    # U_int	mask = (U_int)1 << WORDBITS - 1, element, temp1, temp2,
    # 	A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    mask = one(K) << (wordbits - 1)
    W = zero(K)
    P = zero(K)
	# Hcode	h = {0};
    h .= zero(K)
	# int	i = NUMBITS * DIM - DIM, j;
    i = gg.b * gg.n - gg.n
    A = zero(K)
	# for (j = A = 0; j < DIM; j++)
    for j = 0:(gg.n - 1)
        # if (p.hcode[j] & mask)
        if (pt[j + 1] & mask) != zero(K)
			# A |= g_mask[j];
            A |= g_mask(K, gg.n, j)
        end
    end
	# S = tS = A;
    S = tS = A
	# P |= S & g_mask[0];
    P |= (S & g_mask(K, gg.n, 0))
    # for (j = 1; j < DIM; j++)
    for j = 1:(gg.n - 1)
        gm = g_mask(K, gg.n, j)
		# if( S & g_mask[j] ^ (P >> 1) & g_mask[j])
        if ((S & gm) ⊻ ((P >> 1) & gm)) != 0
			# P |= g_mask[j];
            P |= gm
        end
    end

    # element = i / WORDBITS;
    elementz = i ÷ wordbits
	# if (i % WORDBITS > WORDBITS - DIM)
    if (i % wordbits) > (wordbits - gg.n)
		# h.hcode[element] |= P << i % WORDBITS;
        h[elementz + 1] |= (P << (i % wordbits))
		# h.hcode[element + 1] |= P >> WORDBITS - i % WORDBITS;
        h[elementz + 2] |= P >> (wordbits - (i % wordbits))
    # else
	else
		# h.hcode[element] |= P << i - element * WORDBITS;
        h[elementz + 1] |= (P << (i - elementz * wordbits))
    end
	# J = DIM;
    J = gg.n
	# for (j = 1; j < DIM; j++)
    j = 1
    while j < gg.n
		# if ((P >> j & 1) == (P & 1))
        if ((P >> j) & one(K)) == (P & one(K))
            # continue;
            j += 1
			continue
		else
		# 	break;
            break
        end
    end
	# if (j != DIM)
    if j != gg.n
		# J -= j;
        J -= j
    end
	# xJ = J - 1;
    xJ = J - 1

	# if (P < 3)
    if P < K(3)
		# T = 0;
        T = zero(K)
	# else
    else
		# if (P % 2)
        if (P % K(2)) != zero(K)
			# T = (P - 1) ^ (P - 1) / 2;
            T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
		# else
        else
			# T = (P - 2) ^ (P - 2) / 2;
            T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
        end
    end
	# tT = T;
    tT = T
	# for (i -= DIM, mask >>= 1; i >=0; i -= DIM, mask >>= 1)
    #     {
    i -= gg.n
    mask >>= 1
    while i >= 0
        println("i=$i mask=$mask T=$T P=$P")
		println("h[0]=$(h[1]) h[1]=$(h[2]) h[2]=$(h[3])")
		# for (j = A = 0; j < DIM; j++)
        A = zero(K)
        for j = 0:(gg.n - 1)
			# if (p.hcode[j] & mask)
            if pt[j + 1] & mask != zero(K)
				# A |= g_mask[j];
                A |= g_mask(K, gg.n, j)
            end
        end
		# W ^= tT;
        W ⊻= tT
		# tS = A ^ W;
        tS = A ⊻ W
		# if (xJ % DIM != 0)
        #     {
        if xJ % gg.n != 0
			# temp1 = tS << xJ % DIM;
            temp1 = tS << (xJ % gg.n)
			# temp2 = tS >> DIM - xJ % DIM;
            temp2 = tS >> (gg.n - (xJ % gg.n))
			# S = temp1 | temp2;
            S = temp1 | temp2
			# S &= ((U_int)1 << DIM) - 1;
            S &= (one(K) << gg.n) - one(K)
        #     }
		# else
        else
			# S = tS;
            S = tS
        end

		# P = S & g_mask[0];
        P = S & g_mask(K, gg.n, 0)
		# for (j = 1; j < DIM; j++)
        for j = 1:(gg.n - 1)
			# if( S & g_mask[j] ^ (P >> 1) & g_mask[j])
            gn = g_mask(K, gg.n, j)
            if ((S & gn) ⊻ ((P >> 1) & gn)) != 0
				# P |= g_mask[j];
                P |= gn
            end
        end
		# element = i / WORDBITS;
        element = i ÷ wordbits
		# if (i % WORDBITS > WORDBITS - DIM)
        #     {
        if (i % wordbits) > (wordbits - gg.n)
			# h.hcode[element] |= P << i % WORDBITS;
            h[element + 1] |= (P << (i % wordbits))
			# h.hcode[element + 1] |= P >> WORDBITS - i % WORDBITS;
            h[element + 2] |= (P >> (wordbits - (i % wordbits)))
        #     }
		# else
        else
			# h.hcode[element] |= P << i - element * WORDBITS;
            h[element + 1] |= P << (i - element * wordbits)
        end
		# if (i > 0)
        #     {
        if i > 0
			# if (P < 3)
            if P < K(3)
				# T = 0;
                T = 0
			# else
            else
				# if (P % 2)
                if P % K(2) != zero(K)
					# T = (P - 1) ^ (P - 1) / 2;
                    T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
				# else
                else
					# T = (P - 2) ^ (P - 2) / 2;
                    T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
                end
            end
			# if (xJ % DIM != 0)
            #     {
            if xJ % gg.n != 0
				# temp1 = T >> xJ % DIM;
                temp1 = T >> (xJ % gg.n)
				# temp2 = T << DIM - xJ % DIM;
                temp2 = T << (gg.n - (xJ % gg.n))
				# tT = temp1 | temp2;
                tT = temp1 | temp2
				# tT &= ((U_int)1 << DIM) - 1;
                tT &= (one(K) << gg.n) - one(K)
            #     }
			# else
            else
				# tT = T;
                tT = T
            end

			# J = DIM;
            J = gg.n
			# for (j = 1; j < DIM; j++)
            j = 1
            while j < gg.n
				# if ((P >> j & 1) == (P & 1))
                if ((P >> j) & one(K)) == (P & one(K))
                #     continue;
                    j += 1
                    continue
                # else
                else
				# 	break;
                    break
                end
            end
			# if (j != DIM)
            if j != gg.n
				# J -= j;
                J -= j
            end
			# xJ += J - 1;
            xJ += J - 1
        end

        i -= gg.n
        mask >>= 1
    end
    println("h[0]=$(h[1]) h[1]=$(h[2]) h[2]=$(h[3])")
end


# Point H_decode (Hcode H)
# {
function H_decode!(gg::FaceContinuous, H::Vector{K}, pt::Vector{K}) where {K}
    wordbits = 8 * sizeof(K)
    # 	U_int	mask = (U_int)1 << WORDBITS - 1, element, temp1, temp2,
    # 		A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    mask = one(K) << (wordbits - 1)
    # 	Point	pt = {0};
    # 	int	i = NUMBITS * DIM - DIM, j;
    i = gg.b * gg.n - gg.n
    
    # 	/*--- P ---*/
    # 	element = i / WORDBITS;
    element = i ÷ wordbits
    # 	P = H.hcode[element];
    P = H[element + 1]
    # 	if (i % WORDBITS > WORDBITS - DIM)
    # 	{
    if (i % wordbits) > (wordbits - gg.n)
    # 		temp1 = H.hcode[element + 1];
        temp1 = H[element + 2]  # one-based
    # 		P >>= i % WORDBITS;
        P >>= i % wordbits
    # 		temp1 <<= WORDBITS - i % WORDBITS;
        temp1 <<= wordbits - (i % wordbits)
    # 		P |= temp1;
        P |= temp1
    # 	}
    else
    # 	else
    # 		P >>= i % WORDBITS;	/* P is a DIM bit hcode */
        P >>= (i % wordbits)
    end
    
    # 	/* the & masks out spurious highbit values */
    # 	#if DIM < WORDBITS
    if gg.n < wordbits
    # 		P &= (1 << DIM) -1;
        P &= (one(K) << gg.n) - one(K)
    # 	#endif
    end
    
    # 	/*--- xJ ---*/
    # 	J = DIM;
    J = gg.n
    # 	for (j = 1; j < DIM; j++)
    j = 1
    while j < gg.n
    # 		if ((P >> j & 1) == (P & 1))
        if ((P >> j) & one(K)) != (P & one(K))
            break
        end
        j += 1
    end
    # 	if (j != DIM)
    if j != gg.n
    # 		J -= j;
        J -= j
    end
    # 	xJ = J - 1;
    xJ = J - 1
    
    # 	/*--- S, tS, A ---*/
    # 	A = S = tS = P ^ P / 2;
    A = S = tS = P ⊻ (P >> 1)
    
    
    # 	/*--- T ---*/
    # 	if (P < 3)
    if P < K(3)
    # 		T = 0;
        T = 0
    # 	else
    else
    # 		if (P % 2)
        if P % 2 != 0
    # 			T = (P - 1) ^ (P - 1) / 2;
                # precedence changes for Julia
            T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
    # 		else
        else
    # 			T = (P - 2) ^ (P - 2) / 2;
            T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
        end
    end
    
    # 	/*--- tT ---*/
    # 	tT = T;
    tT = T
    
    # 	/*--- distrib bits to coords ---*/
    # 	for (j = DIM - 1; P > 0; P >>=1, j--)
    j = gg.n - 1
    while P > zero(K)
    # 		if (P & 1)
        if P & one(K) != 0
    # 			pt.hcode[j] |= mask;
            pt[j + 1] |= mask  # one-based
        end
        P >>= 1
        j -= 1
    end
    
    W = zero(K)
    # 	for (i -= DIM, mask >>= 1; i >=0; i -= DIM, mask >>= 1)
    i -= gg.n
    mask >>= 1
    while i >= 0
        println("i=$i mask=$mask T=$T P=$P, xJ=$xJ, tT=$tT, W=$W")
		println("pt[0]=$(pt[1]) pt[1]=$(pt[2]) pt[2]=$(pt[3])")
    # 	{
    # 		/*--- P ---*/
    # 		element = i / WORDBITS;
        element = i ÷ wordbits
    # 		P = H.hcode[element];
        P = H[element + 1]  # one-based
    # 		if (i % WORDBITS > WORDBITS - DIM)
    # 		{
        if (i % wordbits) > (wordbits - gg.n)
    # 			temp1 = H.hcode[element + 1];
            temp1 = H[element + 2]  # one-based
    # 			P >>= i % WORDBITS;
            P >>= i % wordbits
    # 			temp1 <<= WORDBITS - i % WORDBITS;
            temp1 <<= wordbits - (i % wordbits)  # C precedence
    # 			P |= temp1;
            P |= temp1
    # 		}
    # 		else
        else
    # 			P >>= i % WORDBITS;	/* P is a DIM bit hcode */
            P >>= i % wordbits
        end
    
    # 		/* the & masks out spurious highbit values */
    # 		#if DIM < WORDBITS
        if gg.n < wordbits
    # 			P &= (1 << DIM) -1;
            P &= (one(K) << gg.n) - one(K)
    # 		#endif
        end 
    
    # 		/*--- S ---*/
    # 		S = P ^ P / 2;
        S = P ⊻ (P >> 1)
    
    # 		/*--- tS ---*/
    # 		if (xJ % DIM != 0)
    # 		{
        if xJ % gg.n != 0
    # 			temp1 = S >> xJ % DIM;
            temp1 = S >> (xJ % gg.n)  # C precedence
    # 			temp2 = S << DIM - xJ % DIM;
            temp2 = S << (gg.n - (xJ % gg.n))
    # 			tS = temp1 | temp2;
            tS = temp1 | temp2
    # 			tS &= ((U_int)1 << DIM) - 1;
            tS &= (one(K) << gg.n) - 1
    # 		}
    # 		else
        else
    # 			tS = S;
            tS = S
        end
    
    # 		/*--- W ---*/
    # 		W ^= tT;
        W ⊻= tT
    
    # 		/*--- A ---*/
    # 		A = W ^ tS;
        A = W ⊻ tS
    
    # 		/*--- distrib bits to coords ---*/
    # 		for (j = DIM - 1; A > 0; A >>=1, j--)
        j = gg.n - 1
        while A > 0
    # 			if (A & 1)
            if A & one(K) != 0
    # 				pt.hcode[j] |= mask;
                pt[j + 1] |= mask  # one-based
            end
            A >>= 1
            j -= 1
        end
    
    # 		if (i > 0)
    # 		{
        if i > 0
    # 			/*--- T ---*/
    # 			if (P < 3)
            if P < K(3)
    # 				T = 0;
                T = 0
    # 			else
            else
    # 				if (P % 2)
                if P % 2 != 0
    # 					T = (P - 1) ^ (P - 1) / 2;
                    T = (P - one(K)) ⊻ ((P - one(K)) >> 1)
    # 				else
                else
    # 					T = (P - 2) ^ (P - 2) / 2;
                    T = (P - K(2)) ⊻ ((P - K(2)) >> 1)
                end
            end
    
    # 			/*--- tT ---*/
    # 			if (xJ % DIM != 0)
    # 			{
            if xJ % gg.n != 0
    # 				temp1 = T >> xJ % DIM;
                temp1 = T >> (xJ % gg.n)  # C precedence
    # 				temp2 = T << DIM - xJ % DIM;
                temp2 = T << (gg.n - (xJ % gg.n))
    # 				tT = temp1 | temp2;
                tT = temp1 | temp2
    # 				tT &= ((U_int)1 << DIM) - 1;
                tT &= (one(K) << gg.n) - one(K)
    # 			}
    # 			else
            else
    # 				tT = T;
                tT = T
            end
    
    # 			/*--- xJ ---*/
    # 			J = DIM;
            J = gg.n
    # 			for (j = 1; j < DIM; j++)
            j = 1
            while j < gg.n
    # 				if ((P >> j & 1) == (P & 1))
                if ((P >> j) & one(K)) != P & one(K)
    # 					continue;
                    break
                end
                j += 1
    # 				else
    # 					break;
            end
            println("j=$j P=$P J=$J")
    # 			if (j != DIM)
            if j != gg.n
    # 				J -= j;
                J -= j
            end
    # 			xJ += J - 1;
            xJ += J - 1
    # 		}
        end
        i -= gg.n
        mask >>= 1
    end
    # 	}
    # 	return pt;
    # }
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
