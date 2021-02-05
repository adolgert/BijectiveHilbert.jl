module BijectivePaper

# Columns are odd resolution, even resolution (rmin)
# Columns are quadrants in x and y like you'd graph it.
# 1 2
# 0 3
hilbert_variant_encode_table = Function[
    ((x, y, w) -> (x, y))          ((x, y, w) -> (x, y))
    ((x, y, w) -> (y-w, x))        ((x, y, w) -> (y, x-w))
    ((x, y, w) -> (y-w, x-w))      ((x, y, w) -> (y-w, x-w))
    ((x, y, w) -> (2w-x-1, w-y-1)) ((x, y, w) -> (w-x-1, 2w-y-1))
]


function encode_hilbert_zero(x::Integer, y::Integer)
    z = zero(x)
    if x == 0 && y == 0
        return z
    end
    rmin = convert(typeof(x), floor(log2(max(x, y))) + 1)
    w = 1<<(rmin - 1)
    while rmin > 0
        if rmin&1 == 1  # odd
            quadrant = x < w ? (y < w ? 0 : 1) : (y < w ? 3 : 2)
            parity = 1
        else  # even
            quadrant = x < w ? (y < w ? 0 : 3) : (y < w ? 1 : 2)
            parity = 2
        end
        z = (z<<2) + quadrant
        x, y = hilbert_variant_encode_table[quadrant+1, parity](x, y, w)
        rmin -= 1
        w >>= 1
    end
    z
end


hilbert_variant_decode_table = Function[
    ((x, y, w) -> (x, y))          ((x, y, w) -> (x, y))
    ((x, y, w) -> (y, x+w))        ((x, y, w) -> (y+w, x))
    ((x, y, w) -> (y+w, x+w))      ((x, y, w) -> (y+w, x+w))
    ((x, y, w) -> (2w-x-1, w-y-1)) ((x, y, w) -> (w-x-1, 2w-y-1))
]


function decode_hilbert_zero(z)
    r = z & 3
    x, y = typeof(z).([(0, 0), (0, 1), (1, 1), (1, 0)][r + 1])
    z >>= 2
    rmin = 2
    w = convert(typeof(z), 2)
    while z > 0
        r = z & 3
        parity = 2 - rmin&1
        x, y = hilbert_variant_decode_table[r+1, parity](x, y, w)
        z >>= 2
        rmin += 1
        w <<= 1
    end
    x, y
end


encode_hilbert(x::Integer, y::Integer) = encode_hilbert_zero(x - 1, y - 1) + 1


function decode_hilbert(z::Integer)
    x, y = decode_hilbert_zero(z - 1)
    x + 1, y + 1
end


"""
    hilbert_order(v::AbstractArray, subdivisions)

Calculates the permutation of `v` that would order a 2D array by Hilbert curve,
where `v` is an array of real numbers, dimension (2, N) and subdivisions is
a real number that specifies how many boxes to place the `v` into, per side.

This tells you how to order an array of 2-dimensional numbers so that
they have more memory locality in 1 dimension.
Given an array of real numbers of dimension `(2, n)`, subdivide them in each
dimension by `subdivisions`, and assign each point a Hilbert code. Return
the permutation that would sort the given array by that Hilbert code.

See also: [`encode_hilbert_zero`](@ref).

# Example
```julia
rng = MersenneTwister(984720987)
points_in_space = zeros(2, 100)
rand!(rng, points_in_space)
points_reordered = points_in_space[:, hilbert_order(points_in_space, 50)]
```
"""
function hilbert_order(v::AbstractArray, subdivisions)
    lowx = minimum(v[1, :])
    lowy = minimum(v[2, :])
    highx = maximum(v[1, :])
    highy = maximum(v[2, :])
    iv = zeros(Int64, size(v, 2))
    for i in 1:size(v, 2)
        iv[i] = encode_hilbert_zero(
                    round(Int64, (v[1, i] - lowx) * subdivisions / (highx - lowx)),
                    round(Int64, (v[2, i] - lowy) * subdivisions / (highy - lowy))
                    )
    end
    sortperm(iv)
end

export encode_hilbert_zero, decode_hilbert_zero
export encode_hilbert, decode_hilbert, hilbert_order

end
