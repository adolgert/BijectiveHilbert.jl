
abstract type HilbertAlgorithm{T} end

index_type(::HilbertAlgorithm{T}) where {T} = T
