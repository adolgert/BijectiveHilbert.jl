bshow(i) = bitstring(i)[end-7:end]


function large_enough_unsigned(bit_cnt)
    unsigned_types = [UInt8, UInt16, UInt32, UInt64, UInt128]
    atype = nothing
    for xtype in unsigned_types
        if sizeof(xtype) * 8 >= bit_cnt
            atype = xtype
            break
        end
    end
    atype
end


"""
    is_power_of_two(v::Base.BitInteger)

This integer is a power of two.
"""
function is_power_of_two(v::Base.BitInteger)
    v != 0 && ((v & (v - 1)) == 0)
end



"""
    log_base2(v::Integer)

Count the number of bits in an integer, rounding down.
"""
function log_base2(v::Integer)::Integer
    r = zero(v)
    while (v = v >> one(v)) > 0
      r += one(v)
    end
    r
end


"""
    count_set_bits(v::Integer)

A naive algorithm to count bits that are set.
Look here to improve: https://graphics.stanford.edu/~seander/bithacks.html.
"""
function count_set_bits(v::Integer)
    c = zero(v)
    while v != zero(v)
        c += v & one(v)
        v = v>>one(v)
    end
    c
end


"""
    msb(v::Integer)

Most-significant bit, zero-based count. So `0b1` is 0,
`0b1010` is 3.
"""
function msb(v::Integer)
    r = 0
    while (v >>= 1) != 0
        r += 1
    end
    r
end


"""
    trailing_zero_bits(v::Integer)

The number of zero bits after the last one-bit.
"""
function trailing_zero_bits(v::Integer)
    c = zero(v)
    if v != zero(v)
        # Set v's trailing 0s to 1s and zero rest
        v = (v ⊻ (v - one(v))) >> one(v)
        while v != zero(v)
            v = v >> 1
            c += one(v)
        end
    else
        c = 8 * sizeof(v)
    end
    c
end


"""
    trailing_set_bits(v::Integer)

The number of one-bits after the last zero-bit.
"""
trailing_set_bits(v) = trailing_zero_bits(~v)


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


"""
2^n - 1, computed for an unsigned integer.
"""
function fbvn1s(_::T, n) where {T}
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


function trailing_set_bits_hamilton(i::UInt64)
    T = UInt64
    c = zero(Int)
    if i == ~zero(T)
        return T(8) * sizeof(T)
    end
    if (i & fbvn1s(T, 32)) == fbvn1s(T, 32)
        i >>= T(32)
        c ⊻= T(32)
    end
    if (i & fbvn1s(T, 16)) == fbvn1s(T, 16)
        i >>= T(16)
        c ⊻= T(16)
    end
    if (i & fbvn1s(T, 8)) == fbvn1s(T, 8)
        i >>= T(8)
        c ⊻= T(8)
    end
    if (i & fbvn1s(T, 4)) == fbvn1s(T, 4)
        i >>= T(4)
        c ⊻= T(4)
    end
    if (i & fbvn1s(T, 2)) == fbvn1s(T, 2)
        i >>= T(2)
        c ⊻= T(2)
    end
    if (i & fbvn1s(T, 1)) == fbvn1s(T, 1)
        i >>= T(1)
        c ⊻= T(1)
    end
    c
end


function trailing_set_bits_hamilton(i::UInt32)
    T = UInt32
    c = zero(Int)
    if i == ~zero(T)
        return T(8) * sizeof(T)
    end
    if (i & fbvn1s(T, 16)) == fbvn1s(T, 16)
        i >>= T(16)
        c ⊻= T(16)
    end
    if (i & fbvn1s(T, 8)) == fbvn1s(T, 8)
        i >>= T(8)
        c ⊻= T(8)
    end
    if (i & fbvn1s(T, 4)) == fbvn1s(T, 4)
        i >>= T(4)
        c ⊻= T(4)
    end
    if (i & fbvn1s(T, 2)) == fbvn1s(T, 2)
        i >>= T(2)
        c ⊻= T(2)
    end
    if (i & fbvn1s(T, 1)) == fbvn1s(T, 1)
        i >>= T(1)
        c ⊻= T(1)
    end
    c
end


function trailing_set_bits_hamilton(i::UInt16)
    T = UInt16
    c = zero(Int)
    if i == ~zero(T)
        return T(8) * sizeof(T)
    end
    if (i & fbvn1s(T, 8)) == fbvn1s(T, 8)
        i >>= T(8)
        c ⊻= T(8)
    end
    if (i & fbvn1s(T, 4)) == fbvn1s(T, 4)
        i >>= T(4)
        c ⊻= T(4)
    end
    if (i & fbvn1s(T, 2)) == fbvn1s(T, 2)
        i >>= T(2)
        c ⊻= T(2)
    end
    if (i & fbvn1s(T, 1)) == fbvn1s(T, 1)
        i >>= T(1)
        c ⊻= T(1)
    end
    c
end


function trailing_set_bits_hamilton(i::UInt8)
    T = UInt8
    c = zero(Int)
    if i == ~zero(T)
        return T(8) * sizeof(T)
    end
    if (i & fbvn1s(T, 4)) == fbvn1s(T, 4)
        i >>= T(4)
        c ⊻= T(4)
    end
    if (i & fbvn1s(T, 2)) == fbvn1s(T, 2)
        i >>= T(2)
        c ⊻= T(2)
    end
    if (i & fbvn1s(T, 1)) == fbvn1s(T, 1)
        i >>= T(1)
        c ⊻= T(1)
    end
    c
end


function first_set_bit(i::T)::Int where {T <: Integer}
    if i == 0
        return 0
    end
    for j = 0:(8 * sizeof(T)  - 1)
        if (i & (T(1) << j)) != 0
            return j + 1
        end
    end
    return 0 # for obvious type stability
end


function first_set_bit(i::UInt128)::Int
    T = UInt128
    c = zero(Int)
    if i == zero(T)
        return c
    end
    if (i & fbvn1s(T, 64)) == zero(T)
        i >>= T(64)
        c ⊻= 64
    end
    if (i & fbvn1s(T, 32)) == zero(T)
        i >>= T(32)
        c ⊻= 32
    end
    if (i & fbvn1s(T, 16)) == zero(T)
        i >>= T(16)
        c ⊻= 16
    end
    if (i & fbvn1s(T, 8)) == zero(T)
        i >>= T(8)
        c ⊻= 8
    end
    if (i & fbvn1s(T, 4)) == zero(T)
        i >>= T(4)
        c ⊻= 4
    end
    if (i & fbvn1s(T, 2)) == zero(T)
        i >>= T(2)
        c ⊻= 2
    end
    if (i & fbvn1s(T, 1)) == zero(T)
        i >>= T(1)
        c ⊻= 1
    end
    c + one(Int)
end


function first_set_bit(i::UInt64)::Int
    T = UInt64
    c = zero(Int)
    if i == zero(T)
        return c
    end
    if (i & fbvn1s(T, 32)) == zero(T)
        i >>= T(32)
        c ⊻= 32
    end
    if (i & fbvn1s(T, 16)) == zero(T)
        i >>= T(16)
        c ⊻= 16
    end
    if (i & fbvn1s(T, 8)) == zero(T)
        i >>= T(8)
        c ⊻= 8
    end
    if (i & fbvn1s(T, 4)) == zero(T)
        i >>= T(4)
        c ⊻= 4
    end
    if (i & fbvn1s(T, 2)) == zero(T)
        i >>= T(2)
        c ⊻= 2
    end
    if (i & fbvn1s(T, 1)) == zero(T)
        i >>= T(1)
        c ⊻= 1
    end
    c + one(Int)
end


function first_set_bit(i::UInt32)::Int
    T = UInt32
    c = zero(Int)
    if i == zero(T)
        return c
    end
    if (i & fbvn1s(T, 16)) == zero(T)
        i >>= T(16)
        c ⊻= 16
    end
    if (i & fbvn1s(T, 8)) == zero(T)
        i >>= T(8)
        c ⊻= 8
    end
    if (i & fbvn1s(T, 4)) == zero(T)
        i >>= T(4)
        c ⊻= 4
    end
    if (i & fbvn1s(T, 2)) == zero(T)
        i >>= T(2)
        c ⊻= 2
    end
    if (i & fbvn1s(T, 1)) == zero(T)
        i >>= T(1)
        c ⊻= 1
    end
    c + one(Int)
end


function first_set_bit(i::UInt16)::Int
    T = UInt16
    c = zero(Int)
    if i == zero(T)
        return c
    end
    if (i & fbvn1s(T, 8)) == zero(T)
        i >>= T(8)
        c ⊻= 8
    end
    if (i & fbvn1s(T, 4)) == zero(T)
        i >>= T(4)
        c ⊻= 4
    end
    if (i & fbvn1s(T, 2)) == zero(T)
        i >>= T(2)
        c ⊻= 2
    end
    if (i & fbvn1s(T, 1)) == zero(T)
        i >>= T(1)
        c ⊻= 1
    end
    c + one(Int)
end


function first_set_bit(i::UInt8)::Int
    T = UInt8
    c = zero(Int)
    if i == zero(T)
        return c
    end
    if (i & fbvn1s(T, 4)) == zero(T)
        i >>= T(4)
        c ⊻= 4
    end
    if (i & fbvn1s(T, 2)) == zero(T)
        i >>= T(2)
        c ⊻= 2
    end
    if (i & fbvn1s(T, 1)) == zero(T)
        i >>= T(1)
        c ⊻= 1
    end
    c + one(Int)
end
