module BijectiveHilbert

include("hilbert.jl")
export encode_hilbert_zero
export decode_hilbert_zero
export encode_hilbert
export decode_hilbert
export hilbert_order

include("global_gray.jl")
export GlobalGray
export axis_type
export transpose_type
export encode_hilbert_zero!
export decode_hilbert_zero!
end
