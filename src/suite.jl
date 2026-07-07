using TestItems

# This is a test suite that should work for every Hilbert curve.
@testmodule HilbertTestSuite begin

using BijectiveHilbert: HilbertAlgorithm, decode_hilbert_zero!, encode_hilbert_zero
using BijectiveHilbert: index_type, axis_type

export check_own_inverse, check_complete_set

"""
    check_own_inverse(algorithm, extent, dimensions)

Assert that encoding and decoding are inverses of each other. The extent is
log2(side length). Dimensions are the number of dimensions for the curve.
This shows that encode(decode(z)) == z for the Hilbert coordinate.
It also checks that the Cartesian coordinate is within the given dimensions.
"""
function check_own_inverse(gg::HilbertAlgorithm, b::Int, n)
    success = true
    p = zeros(axis_type(gg), n)
    H = index_type(gg)
    extent = 1<<b
    for h in 0:(1<<(n*b) - 1)
        decode_hilbert_zero!(gg, p, H(h))
        for idx in 1:n
            if p[idx] < 0 || p[idx] >= extent
                println("for b $(b) value is $(p)")
                return false
            end
        end
        h2 = encode_hilbert_zero(gg, p)
        if h2 != H(h)
            return false
        end
    end
    success
end


"""
    check_own_inverse(algorithm, ms::Vector{Int})

Assert that encoding and decoding are inverses of each other.
Here, `ms` contains the number of bits per axis (not extents).
This shows that encode(decode(z)) == z for every value of the Hilbert coordinate
and that the decoded value is within the given dimensions.
"""
function check_own_inverse(gg::HilbertAlgorithm, ms::Vector{Int})
    n = length(ms)
    A = axis_type(gg)
    p = zeros(A, n)
    H = index_type(gg)
    total_bits = sum(ms)
    for _ in 0:(1<<total_bits - 1)
        h = encode_hilbert_zero(gg, p)
        for idx in 1:n
            if p[idx] < 0 || p[idx] >= one(A)<<ms[idx]
                println("for bits $(ms) value is $(p)")
                return false
            end
        end
        q = similar(p)
        decode_hilbert_zero!(gg, q, h)
        if p != q
            return false
        end
        for inc in 1:n
            p[inc] += one(A)
            if p[inc] == one(A)<<ms[inc]
                p[inc] = zero(A)
            else
                break
            end
        end
    end
    true
end


function check_complete_set(gg::HilbertAlgorithm{H}, b, n) where {H}
    dim_cnt = n
    m = b
    A = axis_type(gg)
    success = true
    indices = Base.IteratorsMD.CartesianIndices(tuple(collect(1<<m for i in 1:dim_cnt)...))
    seen2 = Dict{H, Vector{A}}()
    for idx in indices
        # -1 for zero-based.
        vidx = convert(Vector{A}, [Tuple(idx)...] .- 1)
        h = encode_hilbert_zero(gg, vidx)
        seen2[h] = vidx
        if h < zero(H)
            success = false
        end
        if h >= one(H)<<(dim_cnt * m)
            success = false
        end
    end
    if length(seen2) != 1<<(dim_cnt * m)
        success = false
    end
    for ihidx in 0:(1<<(dim_cnt*m) - 2)  # compare with next, so stop one early.
        hidx = H(ihidx)
        differ = seen2[hidx] .!= seen2[hidx + one(H)]
        if sum(differ) != 1
            @show ihidx, seen2[hidx], seen2[hidx + one(H)]
            success = false
        end
        if sum(differ) == 1
            a = seen2[hidx][differ][1]
            b = seen2[hidx + 1][differ][1]
            dx = (a > b) ? a - b : b - a
            if H(dx) != one(H)
                success = false
                break
            end
        end
    end
    success
end


"""
    check_complete_set(gg::HilbertAlgorithm, ms::Vector{Int})

Anisotropic complete-set check. `ms` holds the number of bits per axis. Walks
the whole anisotropic box `[0, 2^ms[1]) x ... x [0, 2^ms[n])` with an odometer,
asserts that `encode_hilbert_zero` is a bijection onto `0:2^sum(ms)-1`, and that
consecutive Hilbert indices decode to points that differ in exactly one axis by
exactly 1.
"""
function check_complete_set(gg::HilbertAlgorithm{H}, ms::Vector{Int}) where {H}
    n = length(ms)
    A = axis_type(gg)
    total_bits = sum(ms)
    success = true

    # Odometer over the anisotropic box; record each point by its Hilbert index.
    seen = Dict{H, Vector{A}}()
    p = zeros(A, n)
    # Guard the upper-bound check against index-type overflow when total_bits
    # equals the index-type width (then any representable h is in range).
    for _ in 0:(1 << total_bits - 1)
        h = encode_hilbert_zero(gg, p)
        if h < zero(H) || (total_bits < 8 * sizeof(H) && h >= one(H) << total_bits)
            success = false
        end
        seen[h] = copy(p)
        for inc in 1:n
            p[inc] += one(A)
            if p[inc] == one(A) << ms[inc]
                p[inc] = zero(A)
            else
                break
            end
        end
    end

    if length(seen) != 1 << total_bits
        success = false
    end

    for ihidx in 0:(1 << total_bits - 2)  # compare with next, so stop one early.
        hidx = H(ihidx)
        (haskey(seen, hidx) && haskey(seen, hidx + one(H))) || (success = false; break)
        differ = seen[hidx] .!= seen[hidx + one(H)]
        if sum(differ) != 1
            @show ihidx, seen[hidx], seen[hidx + one(H)]
            success = false
            break
        end
        a = seen[hidx][differ][1]
        b = seen[hidx + one(H)][differ][1]
        dx = (a > b) ? a - b : b - a
        if A(dx) != one(A)
            success = false
            break
        end
    end
    success
end

end
