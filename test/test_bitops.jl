@safetestset log_base_two_float = "log base two agrees with float version" begin
using BijectiveHilbert
for i in 1:100
    m = BijectiveHilbert.log_base2(i) + 1
    # println(i, " ", string(i, base = 2), " ", m, " ", ceil(Int, log(2, i)))
    @test((1<<(m - 1)) & i != 0)
    @test(1<<m & i == 0)
end
end


@safetestset bitrotate_cases = "bitrotate rotates particular cases" begin
using BijectiveHilbert    
tries = [
    Any[0b1, 0b10, 1, 8],
    Any[0b1, 0b10000000, -1, 8],
    Any[0b1, 0b01000000, -1, 7],
    Any[0b1, 0b01000000, -2, 8],
    Any[0b1, 0b00001000, -1, 4],
    Any[0b1, 0b00000100, -1, 3],
    Any[0b1, 0b00000010, -2, 3],
    Any[0b1, 0b00000001, -3, 3],
    Any[0b1, 0b00000100, -4, 3],
    Any[0b1, 0b00000010, -5, 3],
    Any[0b1, 0b100, 2, 8],
    Any[0b1, 0b1000, 3, 8],
    Any[0b1, 0b1000, 3, 4],
    Any[0b1, 0b0001, 4, 4],
    Any[0b1, 0b0010, 5, 4],
    Any[0b1, 0b1000, 7, 4],
    Any[0b1, 0b0001, 8, 4]
]
for (try_idx, trial) in enumerate(tries)
    x, y, k, n = trial
    y1 = BijectiveHilbert.bitrotaten(x, k, n)
    if y != y1
        println("k ", k, " n ", n, " ", bitstring(y), " ", bitstring(y1))
    end
    @test y == y1
end

end


@safetestset bitrotate_matches_naive = "bitrotate matches naive" begin
using BijectiveHilbert
for n = 3:3
    for k = -2n:2n
        for x = 0b0:(0b1<<n - 0b1)
            @test(typeof(x) == UInt8)
            bd = bitstring(BijectiveHilbert.bitrotaten(x, k, n))
            b = bitstring(BijectiveHilbert.bitrotate_a(x, k, n))
            bx = bitstring(x)
            if b != bd
                @show k, n, bx, b, bd
            end
        end
    end
end

end

@safetestset bitrotate_matches_eight = "bitrotate matches eight bit version" begin
using BijectiveHilbert    
# It agrees with the 8-bit version when n=8.
n = 8
for k = -3n:3n
    for x = 0b0:(0b1<<n - 0b1)
        @test(BijectiveHilbert.bitrotaten(x, k, n) == bitrotate(x, k))
    end
end

end

@safetestset rotating_is_inverse = "rotating is its own inverse" begin
using BijectiveHilbert
n = 7
for k in 0:n
    for i in 0x1:(0x1<<n - 0x1)
        j = BijectiveHilbert.bitrotaten(BijectiveHilbert.bitrotaten(i, k, n), -k, n)
        @test(j == i)
    end
end
end


@safetestset bitrotate_left = "rotateleft cases" begin
using BijectiveHilbert    
tries = [
    Any[0b1, 0b100, 2, 8],
    Any[0b1, 0b1000, 3, 8],
    Any[0b1, 0b1000, 3, 4],
    Any[0b10, 0b0001, 3, 4],
    Any[0b1000, 0b0001, 1, 4],
    Any[0b1000, 0b0010, 2, 4],
    Any[0b1110, 0b1011, 2, 4]
]
for (try_idx, trial) in enumerate(tries)
    x, y, k, n = trial
    y1 = BijectiveHilbert.rotateleft(x, k, n)
    if y != y1
        println("k ", k, " n ", n, " ", bitstring(y), " ", bitstring(y1))
    end
    @test y == y1
end

end


@safetestset bitrotate_right = "rotateright cases" begin
using BijectiveHilbert
tries = [
    Any[0b1, 0b10000000, 1, 8],
    Any[0b1, 0b01000000, 1, 7],
    Any[0b1, 0b01000000, 2, 8],
    Any[0b1, 0b00001000, 1, 4],
    Any[0b1, 0b00000100, 1, 3],
    Any[0b1, 0b00000010, 2, 3],
]
for (try_idx, trial) in enumerate(tries)
    x, y, k, n = trial
    y1 = BijectiveHilbert.rotateright(x, k, n)
    if y != y1
        println("k ", k, " n ", n, " ", bitstring(y), " ", bitstring(y1))
    end
    @test y == y1
end

end


@safetestset rotateright_random = "rotateright random" begin
    using BijectiveHilbert: fbvn1s, rotateright, rotateright_naive
    using Random
    rng = MersenneTwister(9871879147)
    for T in [UInt8, UInt16, UInt32, UInt64]
        for trial in 1:1000
            bit_max = 8 * sizeof(T)
            width = rand(rng, 2:bit_max)
            shift = rand(rng, 0:(width - 1))
            val = rand(rng, zero(T):fbvn1s(T, width))
            naive = rotateright_naive(val, shift, width)
            fancy = rotateright(val, shift, width)
            @test fancy == naive
        end
    end
end


@safetestset rotateleft_random = "rotateleft random" begin
    using BijectiveHilbert: fbvn1s, rotateleft, rotateleft_naive
    using Random
    rng = MersenneTwister(9871879147)
    for T in [UInt8, UInt16, UInt32, UInt64]
        for trial in 1:1000
            bit_max = 8 * sizeof(T)
            width = rand(rng, 2:bit_max)
            shift = rand(rng, 0:(width - 1))
            val = rand(rng, zero(T):fbvn1s(T, width))
            naive = rotateleft_naive(val, shift, width)
            fancy = rotateleft(val, shift, width)
            @test fancy == naive
        end
    end
end


@safetestset reverse_trials = "reversebits trials" begin
using BijectiveHilbert: reverse_bits
trials = [
    (0b11001, 5, 0b10011),
    (0b01001, 5, 0b10010),
    (0b00011, 5, 0b11000)
]
for (a, n, b) in trials
    r = reverse_bits(a, n)
    @test bitstring(r) == bitstring(b)
end
end


@safetestset setbits_trials = "setbits does some sets" begin
    using BijectiveHilbert: set_bits, set_bits_naive
    trials = [
        [0b0, 2, 0, 0b11, 0b11],
        [0b0, 2, 1, 0b11, 0b1100],
        [0b0, 2, 1, 0b01, 0b0100],
        [0b0, 2, 1, 0b10, 0b1000],
        [0b0, 3, 1, 0b111, 0b111000],
        [0b0, 3, 1, 0b101, 0b101000],
        [0b0, 3, 1, 0b100, 0b100000],
        [0b1, 3, 1, 0b100, 0b100001]
    ]
    for (h, n, i, w, r) in trials
        a = set_bits_naive(h, n, i, w)
        @test bitstring(a) == bitstring(r)
    end
end


@safetestset setbits_agrees = "setbits is same as hamilton" begin
    using BijectiveHilbert: set_bits, set_bits_naive
    using Random
    rng = MersenneTwister(82147985)
    for i in 1:10000
        a = rand(rng, UInt64)
        n = rand(rng, 1:7)
        i = rand(rng, 1:8)
        w = rand(rng, 0:(1<<n - 1))
        h1 = set_bits(a, n, i, w)
        h2 = set_bits_naive(a, n, i, w)
        @test h1 == h2
    end
end
