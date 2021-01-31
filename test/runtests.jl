using BijectiveHilbert
using SafeTestsets
using Test

@testset "BijectiveHilbert.jl" begin
include("test_bitops.jl")
include("test_gray_code.jl")
include("test_simple2d.jl")
include("test_hamilton.jl")
include("test_global_gray.jl")
include("test_facecontinuous.jl")
end
