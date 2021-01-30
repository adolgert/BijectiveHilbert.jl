@safetestset gray_code_symmetry = "Gray code symmetry" begin
using BijectiveHilbert
# for i in 0:10
#     println(lpad(i, 3, " "), " ", lpad(string(BijectiveHilbert.brgc(i), base = 2), 8, "0"))
# end
# Here's a test of the gray code.
n = 5
for i in 0x0:(0x1 << n - 0x1)
    # println(BijectiveHilbert.brgc(1<<n - 1 - i), " ", BijectiveHilbert.brgc(i) ⊻ (1 << (n-1)))
    @test(BijectiveHilbert.brgc(1<<n - 1 - i) == BijectiveHilbert.brgc(i) ⊻ (1 << (n - 1)))
end
end

@safetestset gray_code_single = "Gray code changes one bit" begin
using BijectiveHilbert: brgc, is_power_of_two
for i in 0x1:0x1000
    @test is_power_of_two(brgc(i) ⊻ brgc(i + 1))
end
end


@safetestset brgc_own_inverse = "Gray code is own inverse" begin
using BijectiveHilbert
for i in 0x0:0x1000
    @test(BijectiveHilbert.brgc_inv(BijectiveHilbert.brgc(i)) == i)
end
end


@safetestset brgc_equals_naive = "Gray code matches paper description" begin
using BijectiveHilbert: brgc_inv
using Random
rng = MersenneTwister(9719742)
for T in [UInt8, UInt16, UInt32, UInt64]
    for trial in 1:10000
        v = rand(rng, T)
        n = brgc_inv(Int128(v))  # The paper version is used for integers.
        i = brgc_inv(v)
        @test n == i
    end
end
end


@safetestset brgc_ranks_equal = "the rank calculation is the same as the paper" begin
    using BijectiveHilbert: brgc_rank, brgc_rank2
    using Random
    rng = MersenneTwister(974073242)
    for trial in 1:1000
        n = rand(rng, 2:7)
        w = UInt8(rand(rng, 0:(1<<n - 1)))
        mask = UInt8(rand(rng, 0:(1<<n - 1)))
        a = brgc_rank(mask, w, n)
        b = brgc_rank2(mask, w, n)
        @test a == b
    end
end
