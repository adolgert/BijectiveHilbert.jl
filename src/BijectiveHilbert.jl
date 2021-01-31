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
export Simple2D

include("facecontinuous.jl")
export FaceContinuous

include("suite.jl")
end
