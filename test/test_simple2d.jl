using TestItemRunner


@testitem "hilbert decode zero" begin
  using BijectiveHilbert
  table2 = [
    15 12 11 10
    14 13  8  9
    1  2  7  6
    0  3  4  5
  ]
  table2 = table2[end:-1:1, :]

  table2ind = CartesianIndices(table2)
  gg = Simple2D(Int)
  X = zeros(Int, 2)
  for z in 0:15
    BijectiveHilbert.decode_hilbert_zero!(gg, X, z)
    found = table2ind[table2 .== z][1]
    @test found[2] - 1 == X[1]
    @test found[1] - 1 == X[2]
  end
end


@testitem "simple2d is its own inverse" setup=[HilbertTestSuite] begin
  using BijectiveHilbert: Simple2D
  gg = Simple2D(Int)
  for b in 2:7
    @test HilbertTestSuite.check_own_inverse(gg, b, 2)
  end
end


@testitem "simple2d is a complete set" setup=[HilbertTestSuite] begin
  using BijectiveHilbert: Simple2D
  gg = Simple2D(Int)
  for b in 2:7
    @test HilbertTestSuite.check_complete_set(gg, b, 2)
  end
end


@testitem "Simple2D type interactions" begin
    using BijectiveHilbert
    using UnitTestDesign
    using Random
    rng = Random.MersenneTwister(9790323)
    for retrial in 1:5
        AxisTypes = shuffle(rng, [Int8, Int, UInt, Int128, UInt8, UInt128])
        IndexTypes = shuffle(rng, [Int8, UInt8, Int, UInt, Int128, UInt128])
        Count= shuffle(rng, [0, 1])
        Dims = shuffle(rng, [2, 3, 4])
        Bits = shuffle(rng, [2, 3, 4, 5])
        test_set = all_pairs(
            AxisTypes, IndexTypes, Count, Dims, Bits;
        )
        for (A, I, C, D, B) in test_set
            gg = Simple2D(I)
            if B * D > log2(typemax(I))
                continue
            end
            # Add this because these values aren't in Hilbert order because
            # It's an asymmetrical space.
            if B * D > log2(typemax(A))
              continue
            end
          last = (one(I) << (B * D)) - one(I) + I(C)
            mid = one(I) << (B * D - 1)
            few = 5
            X = zeros(A, D)
            hlarr = vcat(C:min(mid, few), max(mid + 1, last - few):last)
            for hl in hlarr
                local hli = I(hl)
                if C == 0
                    decode_hilbert_zero!(gg, X, hli)
                    hl2 = encode_hilbert_zero(gg, X)
                    if hl2 != hli
                        @show A, I, C, D, B, X
                        @test hl2 == hli
                    end
                    @test typeof(hl2) == typeof(hli)
                else
                    decode_hilbert!(gg, X, hli)
                    hl3 = encode_hilbert(gg, X)
                    @test hl3 == hli
                    @test typeof(hl3) == typeof(hli)
                end
            end
        end
    end
end
