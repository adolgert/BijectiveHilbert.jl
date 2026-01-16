
using StaticArraysCore: SVector, MVector

function normalizer_to_12d(xv::Vector{SVector{D,Float64}}) where{D}
  extreme_axes = [extrema(x->x[j], xv) for j in 1:D]
  divisor = maximum(x->x[2]-x[1], extreme_axes)*(1 + 10*eps()) 
  minx    = SVector{D,Float64}([x[1] for x in extreme_axes])
  (minx, divisor)
end

# send the float in [1,2) to a 52 bit UInt, and then add 12 zeros in front of the binary to make a UInt64.
function float_to_int(x::Float64)
  as_uint = reinterpret(UInt64, x)
  UInt64(as_uint & ((UInt64(1) << 52) - 1))
end
float_to_int(x::SVector{D,Float64}) where{D} = float_to_int.(x)

struct NormalizingHilbertEncoder{D,T,B}
  minx::SVector{D,Float64}
  divisor::Float64
  encoder::SpaceGray{T,B}
end

function (he::NormalizingHilbertEncoder{D,T,B})(x) where{D,T,B} 
  x_12 = SVector{D,Float64}(ntuple(j->(x[j] - he.minx[j])/he.divisor + 1, D))
  encode_hilbert(he.encoder, float_to_int(x_12))
end

# Maximum bits per axis that fits in UInt128 for a given dimension
max_bits_per_axis(D::Int) = 128 ÷ D

# In benchmarking, it is probably faster to actually allocate and fill a Vector{UIntXXX} and then
# sort on that instead of using the in-place with by=... But I can grind optimizations later once
# I am convinced of correctness.
function hilbertsort!(pts::Vector{SVector{D,Float64}}; 
                      bits_per_axis::Int=min(52, max_bits_per_axis(D))) where{D}
  (minx, divisor) = normalizer_to_12d(pts)
  encoder = NormalizingHilbertEncoder(minx, divisor, SpaceGray(bits_per_axis, D))
  sort!(pts, by = x -> encoder(x), alg=QuickSort)
end


