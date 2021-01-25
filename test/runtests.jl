using BijectiveHilbert
using SafeTestsets
using Test

@testset "BijectiveHilbert.jl" begin
include("test_bitops.jl")
incluce("test_gray_code.jl")
include("test_hilbert.jl")
include("test_hamilton.jl")
include("test_global_gray.jl")
end
