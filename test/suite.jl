module HilbertTestSuite
using BijectiveHilbert: HilbertAlgorithm, decode_hilbert_zero!, encode_hilbert_zero
using BijectiveHilbert: index_type, axis_type


function check_own_inverse(gg::HilbertAlgorithm, b::Int, n)
    success = true
    p = zeros(axis_type(gg), n)
    H = index_type(gg)
    for h in 0:(1<<(n*b) - 1)
        decode_hilbert_zero!(gg, p, H(h))
        h2 = encode_hilbert_zero(gg, p)
        if h2 != H(h)
            success = false
        end
    end
    success
end


function check_own_inverse(gg::HilbertAlgorithm, ms::Vector{Int}, n)
    success = true
    A = axis_type(gg)
    p = zeros(A, n)
    H = index_type(gg)
    total = prod(ms)
    for i in 0:(1<<total - 1)
        h = encode_hilbert_zero(gg, p)
        q = similar(p)
        decode_hilbert_zero!(gg, q, h)
        if p != q
            success = false
        end
        for inc in 1:n
            p[inc] += one(A)
            if p[inc] == one(A)<<ms[inc]
                p[inc] = zero(A)
            else
                break
            end
        end
        @assert p != q
    end
    success
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

end
