module BijectiveHilbert

include("bitops.jl")
include("gray_code.jl")
include("hilbert.jl")
include("hamilton.jl")

export encode_hilbert_zero
export decode_hilbert_zero
export encode_hilbert
export decode_hilbert
export hilbert_order

export brgc, brgc_inv
export hi_d, hi_e, hi_f, hi_g, hi_T, hi_t_inv
include("global_gray.jl")
export GlobalGray
export axis_type
export transpose_type
export encode_hilbert_zero!
export decode_hilbert_zero!
end
