using TestItemRunner

@testitem "log base two agrees with float version" begin
using BijectiveHilbert
for i in 1:100
    m = BijectiveHilbert.log_base2(i) + 1
    # println(i, " ", string(i, base = 2), " ", m, " ", ceil(Int, log(2, i)))
    @test((1<<(m - 1)) & i != 0)
    @test(1<<m & i == 0)
end
end


@testitem "rotateleft cases" begin
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


@testitem "rotateright cases" begin
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


@testitem "rotateright random" begin
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


@testitem "rotateleft random" begin
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


@testitem "reversebits trials" begin
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


@testitem "setbits does some sets" begin
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


@testitem "setbits is same as hamilton" begin
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


@testitem "getbits trials" begin
using BijectiveHilbert: get_bits
    trials = [
        (0b1, 1, 0, 0b1),
        (0b1, 1, 1, 0b0),
        (0b1, 1, 2, 0b0),
        (0b1, 1, 7, 0b0),
        (0b1, 2, 0, 0b1),
        (0b1, 2, 1, 0b0),
        (0b100, 1, 0, 0b0),
        (0b100, 1, 1, 0b0),
        (0b100, 1, 2, 0b1),
        (0b100, 1, 3, 0b0),
        (0b100, 1, 4, 0b0),
        (0b10000000, 1, 7, 0b1),
        (0b10000000, 1, 6, 0b0),
        (0b10000000, 1, 5, 0b0),
        (0b111000, 1, 2, 0b0),
        (0b111000, 1, 3, 0b1),
        (0b111000, 1, 4, 0b1),
        (0b111000, 1, 5, 0b1),
        (0b111000, 1, 6, 0b0),
        (0b111000, 1, 7, 0b0),
        (0b111000, 2, 0, 0b0),
        (0b111000, 2, 1, 0b0),
        (0b111000, 2, 2, 0b10),
        (0b111000, 2, 3, 0b11),
        (0b111000, 2, 4, 0b11),
        (0b111000, 2, 5, 0b01),
        (0b111000, 2, 6, 0b0),
        (0b111000, 3, 0, 0b0),
        (0b111000, 3, 1, 0b100),
        (0b111000, 3, 2, 0b110),
        (0b111000, 3, 3, 0b111),
        (0b111000, 3, 4, 0b011),
        (0b111000, 3, 5, 0b001),
        (0b111000, 3, 6, 0b000)
    ]
    for (v, b, i, r) in trials
        a = get_bits(v, b, i)
        @test a == r
    end
end


@testitem "first_set_bit unsigned trials" begin
using BijectiveHilbert: first_set_bit
    trials = [
        (0b0, 0),
        (UInt64(0), 0),
        (UInt32(0), 0),
        (UInt16(0), 0),
        (UInt8(0), 0),
        (UInt64(1), 1),
        (UInt32(1), 1),
        (UInt16(1), 1),
        (UInt8(1), 1),
        (UInt128(0b1100), 3),
        (UInt64(0b1100), 3),
        (UInt32(0b1100), 3),
        (UInt16(0b1100), 3),
        (UInt8(0b1100), 3),
        (0x80, 8),
        (0x81, 1),
        (0x82, 2),
        (0x84, 3)
    ]
    for (v, i) in trials
        r = first_set_bit(v)
        @test r == i
    end
end


@testitem "first_set_bit signed matches signed" begin
using BijectiveHilbert: first_set_bit
using Random
trials = [Int8, Int16, Int32, Int64]
rng = MersenneTwister(349720)
for T in trials
    for i in 1:100
        v = rand(rng, T)
        a = first_set_bit(v)
        ui = parse(UInt64, "0b" * bitstring(v))
        b = first_set_bit(ui)
        @test a == b
    end
end
end
