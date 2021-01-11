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
