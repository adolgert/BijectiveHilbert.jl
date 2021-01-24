using BijectiveHilbert
using SafeTestsets
using Test

@testset "BijectiveHilbert.jl" begin
    include("test_hilbert.jl")
    include("test_global_gray.jl")
end
