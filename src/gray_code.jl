
"""
Binary-reflected Gray code.
"""
brgc(i) = i ⊻ (i >> 1)


"""
This is how the paper describes the inverse of the
binary-reflected gray code.
"""
function brgc_inv_naive(g::Integer)
    i = g
    m = log_base2(g)
    for j = 1:m
        i = i ⊻ (g >> j)
    end
    i
end


"""
Inverse of the binary-reflected Gray code.
This takes you from the Gray code to its index.
"""
function brgc_inv(i::UInt128)
    i ⊻= (i >> 1)
    i ⊻= (i >> 2)
    i ⊻= (i >> 4)
    i ⊻= (i >> 8)
    i ⊻= (i >> 16)
    i ⊻= (i >> 32)
    i ⊻= (i >> 64)
end


function brgc_inv(i::UInt64)
    i ⊻= (i >> 1)
    i ⊻= (i >> 2)
    i ⊻= (i >> 4)
    i ⊻= (i >> 8)
    i ⊻= (i >> 16)
    i ⊻= (i >> 32)
end


function brgc_inv(i::UInt32)
    i ⊻= (i >> 1)
    i ⊻= (i >> 2)
    i ⊻= (i >> 4)
    i ⊻= (i >> 8)
    i ⊻= (i >> 16)
end


function brgc_inv(i::UInt16)
    i ⊻= (i >> 1)
    i ⊻= (i >> 2)
    i ⊻= (i >> 4)
    i ⊻= (i >> 8)
end


function brgc_inv(i::UInt8)
    i ⊻= (i >> 1)
    i ⊻= (i >> 2)
    i ⊻= (i >> 4)
end


"""
GrayCodeRank, Algorithm 3 from the paper
mask μ
pattern π
"""
function brgc_rank2(mask::T, w::T, n) where {T}
    r = zero(w)
    for j = (n-1):-1:0
        if (mask >>j ) & one(w) != 0
            r = (r<<1) | ((w>>j) & one(w))
        end
    end
    r
end


"""
Gray code rank from Hamilton's code.
"""
function brgc_rank(mask::T, w::T, n) where {T}
    r = zero(T)
    im = one(T)
    jm = one(T)
    for i = 0:(n-1)
        if mask & im != 0
            if w & im != 0
                r |= jm
            end
            jm <<= 1
        end
        im <<= 1
    end
    r
end


function brgc_rank_inv(mask::T, ptrn, r, n, b) where {T}
    g = zero(T)
    gi = zero(T)
    i = n - 1
    ir = 0  # rack
    im = one(T) << i
    j = b - 1
    jr = 0  # rack
    jm = one(T) << j
    gi0 = zero(T)
    gi1 = zero(T)
    g0 = zero(T)
    while i >= 0
        if mask & im != 0
            gi1 = gi0
            gi0 = ((r & jm) > 0) ? one(T) : zero(T)
            g0 = gi0 ⊻ gi1
            if gi0 != 0
                gi |= im
            end
            if g0 != 0
                g |= im
            end
            jm >>= 1
            if jm == 0
                jm = one(T) << (8 * sizeof(T) - 1)
                jr -= 1
            end
        else
            g0 = ((ptrn & im) > 0) ? one(T) : zero(T)
            gi1 = gi0
            gi0 = g0 ⊻ gi1
            if gi0 != 0
                gi |= im
            end
            if g0 != 0
                g |= im
            end
        end
        im >>= 1
        if im == 0
            im = one(T) << (8 * sizeof(T) - 1)
            ir -= 1
        end
        i -= 1
    end
    g, gi
end
