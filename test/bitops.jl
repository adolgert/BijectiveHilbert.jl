import Base: bitrotate

bshow(i) = bitstring(i)[end-7:end]

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


for i in 1:100
    m = log_base2(i) + 1
    println(i, " ", string(i, base = 2), " ", m, " ", ceil(Int, log(2, i)))
    @assert((1<<(m - 1)) & i != 0)
    @assert(1<<m & i == 0)
end

"""
Rotate bits around an n-bit window. This extends bitrotate,
which operates on the number of bits in a type.
"""
function bitrotate(x::Base.BitInteger, k::Integer, n::Integer)
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

tries = [
    Any[0b1, 0b10, 1, 8],
    Any[0b1, 0b10000000, -1, 8],
    Any[0b1, 0b01000000, -1, 7],
    Any[0b1, 0b01000000, -2, 8],
    Any[0b1, 0b00001000, -1, 4],
    Any[0b1, 0b00000100, -1, 3],
    Any[0b1, 0b00000010, -2, 3],
    Any[0b1, 0b00000001, -3, 3],
    Any[0b1, 0b00000100, -4, 3],
    Any[0b1, 0b00000010, -5, 3],
    Any[0b1, 0b100, 2, 8],
    Any[0b1, 0b1000, 3, 8],
    Any[0b1, 0b1000, 3, 4],
    Any[0b1, 0b0001, 4, 4],
    Any[0b1, 0b0010, 5, 4],
    Any[0b1, 0b1000, 7, 4],
    Any[0b1, 0b0001, 8, 4]
]
for (try_idx, trial) in enumerate(tries)
    x, y, k, n = trial
    y1 = bitrotate(x, k, n)
    if y != y1
        println("k ", k, " n ", n, " ", bitstring(y), " ", bitstring(y1))
    end
end
