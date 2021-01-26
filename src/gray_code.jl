
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
