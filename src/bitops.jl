function large_enough_unsigned(bit_cnt)
    unsigned_types = [UInt8, UInt16, UInt32, UInt64, UInt128]
    atype = nothing
    for xtype in unsigned_types
        if sizeof(xtype) * 8 >= bit_cnt
            atype = xtype
            break
        end
    end
    return atype
end


"""
    is_power_of_two(v::Base.BitInteger)

This integer is a power of two.
"""
function is_power_of_two(v::Base.BitInteger)::Bool
    return v != 0 && ((v & (v - 1)) == 0)
end



"""
    log_base2(v::Integer)

Count the number of bits in an integer, rounding down.
"""
log_base2(v::Integer)::Int = 8 * sizeof(v) - leading_zeros(v) - 1



"""
    trailing_zero_bits(v::Integer)

The number of zero bits after the last one-bit.
"""
trailing_zero_bits(v::Integer)::Int = trailing_zeros(v)


"""
    trailing_set_bits(v::Integer)

The number of one-bits after the last zero-bit.
"""
trailing_set_bits(v::Integer)::Int = trailing_ones(v)


"""
Treat `x` as an `n`-bit unsigned integer. Rotate the bits
`k` places to the left, wrapping them around the right side.
"""
function rotateleft_naive(x::Base.BitInteger, k::Integer, n::Integer)
    y = zero(x)
    for i = 0:(n - 1)
        ind = (i - k) % n
        ind = (ind >= 0) ? ind : ind + n
        if (x & (one(x) << ind)) > 0
            y |= (one(x) << i)
        end  # else leave unset.
    end
    y
end


function rotateleft(x::T, k::Integer, n::Integer) where {T <: Base.BitInteger}
    @assert k >= 0
    @assert n > 0
    @assert k < n
    x &= fbvn1s(T, n)
    x = (x << k) | (x >> (n - k))
    x &= fbvn1s(T, n)
    x
end


function rotateright_naive(x::Base.BitInteger, k::Integer, n::Integer)
    y = zero(x)
    for i = 0:(n - 1)
        ind = (i + k) % n
        ind = (ind >= 0) ? ind : ind + n
        if (x & (one(x) << ind)) > 0
            y |= (one(x) << i)
        end  # else leave unset.
    end
    y
end


function rotateright(x::T, k::Integer, n::Integer) where {T <: Base.BitInteger}
    x &= fbvn1s(T, n)
    x = (x >> k) | (x << (n - k))
    x &= fbvn1s(T, n)
    x
end


"""
2^n - 1, computed for an unsigned integer.
"""
function fbvn1s(T::DataType, n)
    if n == T(sizeof(T) * 8)
        ~zero(T)
    else
        (one(T) << n) - one(T)
    end
end


function fbvmod(i, m)
    if i >= m
        i -= m * i / m
    else
        i
    end
end


function reverse_bits(w, n)
    r = zero(w)
    for i in 0:(n - 1)
        r |= ((w & (one(w)<<i)) >> i) << (n - 1 - i)
    end
    r
end


function set_bits(h, n, i, w)
    h |= w << (i * n)
    h
end


function set_bits_naive(h, n, i, w)
    for j in 0:(n - 1)
        b = (w & (1 << j)) >> j
        h |= (b << (n * i + j))
    end
    h
end


"""
Get b bits from v, starting at place i.
"""
function get_bits(v::T, b, i) where {T}
    n = 8 * sizeof(T)
    (v << (n - b - i)) >> (n - b)
end


"""
    first_set_bit(i::Integer)

Returns the 1-indexed position of the least significant set bit, or 0 if i == 0.
"""
first_set_bit(i::Integer)::Int = i == 0 ? 0 : trailing_zeros(i) + 1
