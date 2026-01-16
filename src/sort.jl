
using StaticArraysCore: SVector

function normalizer_to_12d(xv::AbstractVector{<:AbstractVector})
  D = length(first(xv))
  extreme_axes = [extrema(x->x[j], xv) for j in 1:D]
  divisor = maximum(x->x[2]-x[1], extreme_axes)*(1 + 10*eps()) 
  minx    = SVector{D,Float64}([x[1] for x in extreme_axes])
  (minx, divisor)
end

# send the float in [1,2) to its mantissa bits as a UInt
# For Float64: extracts the 52 explicit mantissa bits
# Works generically for Float16, Float32, Float64, etc.
function float_to_int(x::T) where {T <: AbstractFloat}
  UI = Base.uinttype(T)
  nbits = Base.significand_bits(T)
  # For floats in [1,2), subtract 1 to get [0,1), then scale to [0, 2^nbits)
  UI((x - one(T)) * (one(UI) << nbits))
end
float_to_int(x::SVector{D,T}) where {D,T<:AbstractFloat} = float_to_int.(x)

struct NormalizingHilbertEncoder{D,E}
  minx::SVector{D,Float64}
  divisor::Float64
  encoder::E
  shift::Int  # Right shift amount to extract top bits_per_axis from mantissa
  function NormalizingHilbertEncoder(minx::SVector{D}, divisor::Float64,
        encoder::E, bits_per_axis) where {D,E <: HilbertAlgorithm}
      new{D,E}(minx, divisor, encoder, Base.significand_bits(Float64) - bits_per_axis)
  end
end

function NormalizingHilbertEncoder(minx, divisor, encoder::Type{SpaceGray},
                                   bits_per_axis::Int)
  encoder = SpaceGray(bits_per_axis, length(minx))
  NormalizingHilbertEncoder(minx, divisor, encoder, bits_per_axis)
end

function NormalizingHilbertEncoder(minx, divisor, encoder::Type{Simple2D},
                                   bits_per_axis::Int)
  encoder = Simple2D(UInt128)
  NormalizingHilbertEncoder(minx, divisor, encoder, bits_per_axis)
end

function (he::NormalizingHilbertEncoder{D,E})(x) where{D,E}
  x_12 = SVector{D,Float64}(ntuple(j->(x[j] - he.minx[j])/he.divisor + 1, D))
  # Extract mantissa bits and shift to use the top bits_per_axis bits
  scaled_coords = float_to_int(x_12) .>> he.shift
  encode_hilbert(he.encoder, scaled_coords)
end

function default_encoder(pts::AbstractVector{<:AbstractVector})
  p1 = first(pts)
  length(p1) == 2 ? Simple2D : SpaceGray
end
default_encoder(pts::Vector{SVector{2,Float64}}) = Simple2D

# Maximum bits per axis that fits in UInt128 for a given dimension
max_bits_per_axis(D::Int) = min(Base.significand_bits(Float64), (8*sizeof(UInt128)) ÷ D)

"""
    hilbertsort!(pts::AbstractVector{<:AbstractVector};
                 encoder=default_encoder(pts),
                 with_buffer=false,
                 bits_per_axis=max_bits_per_axis(length(first(pts))))

Sorts `pts` using a Hilbert space-filling curve in-place.
"""
function hilbertsort!(pts::AbstractVector{<:AbstractVector}; 
                      encoder=default_encoder(pts),
                      bits_per_axis::Int=max_bits_per_axis(length(first(pts))))
  norm_encoder = NormalizingHilbertEncoder(normalizer_to_12d(pts)..., encoder, bits_per_axis)
  permute!(pts, sortperm(norm_encoder.(pts)))
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
