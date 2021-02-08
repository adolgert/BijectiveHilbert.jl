using BijectiveHilbert
using BenchmarkTools
using Base.Iterators: take

function encode(gg::HilbertAlgorithm{H}, ms::Vector) where {H}
    A = axis_type(gg)
    success = true
    indices = Base.IteratorsMD.CartesianIndices(tuple(ms...))
    seen2 = Dict{H, Vector{A}}()

    for idx in take(indices, 200)
        # -1 for zero-based.
        vidx = convert(Vector{A}, [Tuple(idx)...] .- 1)
        h = encode_hilbert_zero(gg, vidx)
        seen2[h] = vidx
    end
    seen2
end


function decode(gg::HilbertAlgorithm{H}, ms::Vector) where {H}
    A = axis_type(gg)
    success = true
    indices = Base.IteratorsMD.CartesianIndices(tuple(ms...))
    seen2 = Dict{H, Vector{A}}()
    mvec = zeros(A, length(ms))
    total = H(prod(ms) - 1)
    mid = total >> 1
    cnt = H(200)
    hs = vcat(H(0):minimum([mid, cnt]), maximum([mid, total - cnt]):total)
    for h in hs
        # -1 for zero-based.
        decode_hilbert_zero!(gg, mvec, h)
        seen2[h] = mvec
    end
    seen2
end

HA = GlobalGray
trials = [(5, 2), (31, 2), (5, 3), (20, 3), (5, 10), (11, 10)]
results = zeros(2, length(trials))
for (trial_idx, (bits, dims)) in enumerate(trials)
    ha = HA(bits, dims)
    ms = [2^bits for i in 1:dims]
    @show trial_idx, index_type(ha)
    benchen = @benchmark encode(ha, ms)
    benchde = @benchmark decode(ha, ms)
    results[:, trial_idx] = median.([benchen.times, benchde.times])
end
b = 5
n = 2
ms = [2^b for i in 1:n]
H = UInt64
gg_u64 = GlobalGray(H, b, n)
sg_u64 = SpaceGray(H, b, n)
cc_u64 = Compact(H, ms)
fc_u64 = FaceContinuous(H, b, n)
if n == 2
    ss_u64 = Simple2D(H)
    @benchmark encode(ss_u64, ms)
end
@benchmark encode(gg_u64, ms)
@benchmark encode(sg_u64, ms)
@benchmark encode(cc_u64, ms)
@benchmark encode(fc_u64, ms)
if n == 2
    ss_u64 = Simple2D(H)
    @benchmark decode(ss_u64, ms)
end
@benchmark decode(gg_u64, ms)
@benchmark decode(sg_u64, ms)
@benchmark decode(cc_u64, ms)
b = @benchmark decode(fc_u64, ms)
