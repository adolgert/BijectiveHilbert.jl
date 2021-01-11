bshow(i) = bitstring(i)[end-7:end]

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
Rotate bits around an n-bit window. This extends bitrotate,
which operates on the number of bits in a type.
"""
function bitrotate_a(x::Base.BitInteger, k::Integer, n::Integer)
    k == 0 && return x
    k = k % n
    if k > 0
        keep = one(x)<<(n - k) - one(x)
        shiftmask = ~zero(x) ⊻ keep
        x << k | (x & shiftmask) >> (n - k)
    else  # shift right for k < )
        x << k | ((x & (one(x)<< -k - one(x))) << (n + k))
    end
end


"""
For a right rotation (-k)
Given [b_{n-1}...b_0] return [b_{n-1+i%n}...b_{i%n}]
"""
function bitrotaten(x::Base.BitInteger, k::Integer, n::Integer)
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


function rotateleft(x::Base.BitInteger, k::Integer, n::Integer)
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

function rotateright(x::Base.BitInteger, k::Integer, n::Integer)
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
