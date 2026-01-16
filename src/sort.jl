
using StaticArraysCore: SVector

function normalizer_to_12d(xv::AbstractVector{<:AbstractVector})
  D = length(first(xv))
  extreme_axes = [extrema(x->x[j], xv) for j in 1:D]
  divisor = maximum(x->x[2]-x[1], extreme_axes)*(1 + 10*eps()) 
  minx    = SVector{D,Float64}([x[1] for x in extreme_axes])
  (minx, divisor)
end

# send the float in [1,2) to a 52 bit UInt, and then add 12 zeros in front of
# the binary to make a UInt64.
function float_to_int(x::Float64)
  as_uint = reinterpret(UInt64, x)
  UInt64(as_uint & ((UInt64(1) << 52) - 1))
end
float_to_int(x::SVector{D,Float64}) where{D} = float_to_int.(x)

struct NormalizingHilbertEncoder{D,E}
  minx::SVector{D,Float64}
  divisor::Float64
  encoder::E
end

function NormalizingHilbertEncoder(pts, encoder::Type{SpaceGray},
                                   bits_per_axis::Int)
  (minx, divisor) = normalizer_to_12d(pts)
  encoder = SpaceGray(bits_per_axis, length(minx))
  NormalizingHilbertEncoder(minx, divisor, encoder)
end

function NormalizingHilbertEncoder(pts, encoder::Type{Simple2D},
                                   bits_per_axis::Int)
  (minx, divisor) = normalizer_to_12d(pts)
  encoder = Simple2D(UInt128) 
  NormalizingHilbertEncoder(minx, divisor, encoder)
end

function (he::NormalizingHilbertEncoder{D,E})(x) where{D,E} 
  x_12 = SVector{D,Float64}(ntuple(j->(x[j] - he.minx[j])/he.divisor + 1, D))
  encode_hilbert(he.encoder, float_to_int(x_12))
end

function default_encoder(pts::AbstractVector{<:AbstractVector})
  p1 = first(pts)
  length(p1) == 2 ? Simple2D : SpaceGray
end
default_encoder(pts::Vector{SVector{2,Float64}}) = Simple2D

# Maximum bits per axis that fits in UInt128 for a given dimension
max_bits_per_axis(D::Int) = min(52, 128 ÷ D)

"""
    hilbertsort!(pts::AbstractVector{<:AbstractVector};
                 encoder=default_encoder(pts),
                 with_buffer=false,
                 bits_per_axis=max_bits_per_axis(length(first(pts))))

Sorts `pts` using a Hilbert space-filling curve in-place.
"""
function hilbertsort!(pts::AbstractVector{<:AbstractVector}; 
                      encoder=default_encoder(pts),
                      with_buffer=false,
                      bits_per_axis::Int=max_bits_per_axis(length(first(pts))))
  norm_encoder = NormalizingHilbertEncoder(pts, encoder, bits_per_axis)
  if with_buffer
    buffer = norm_encoder.(pts)
    sp     = sortperm(buffer)
    permute!(pts, sp)
  else
    sort!(pts, by = norm_encoder, alg=QuickSort)
  end
  pts
end

"""
    hilbertsort(pts; kwargs...)

See docstrings for `hilbertsort!`.
"""
function hilbertsort(pts::AbstractVector{<:AbstractVector}; kwargs...)
  _pts = copy(pts)
  hilbertsort!(_pts; kwargs...)
end

