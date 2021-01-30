
abstract type HilbertAlgorithm{A,T} end


axis_type(::HilbertAlgorithm{A,T}) where {A,T} = A
index_type(::HilbertAlgorithm{A,T}) where {A,T} = T
