
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
