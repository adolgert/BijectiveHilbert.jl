module BijectiveHilbert

include("bitops.jl")
include("gray_code.jl")
include("hilbert_algorithm.jl")
export index_type
export encode_hilbert
export decode_hilbert!

include("global_gray.jl")
export GlobalGray
export axis_type
export encode_hilbert_zero
export decode_hilbert_zero!

include("hamilton.jl")
export SpaceGray
export Compact

# include("bijective.jl")
include("simple2d.jl")
import .Bijective: Simple2D

@doc """
Simple2D(::Type{T})

The type is a data type to hold the Hilbert index. It should have
enough bits to hold all bits of integer axes to encode.

The variant algorithm used differs from the usual Hilbert code because it
doesn't need to know the size of the whole grid before computing the code.
It looks like a slightly-rotated version of the Hilbert curve, but it
has the benefit that it is 1-1 between `(x, y)` and `z`, so you can translate
back and forth.

It comes from a paper: 
N. Chen, N. Wang, B. Shi, A new algorithm for encoding and decoding
the Hilbert order. Software—Practice and Experience 2007; 37(8): 897–908.
"""
Simple2D

export Simple2D

include("suite.jl")
end
