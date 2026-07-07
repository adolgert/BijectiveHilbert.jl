using TestItemRunner


@testitem "Compact constructor validation" begin
    using BijectiveHilbert: Compact

    # Valid construction (default BRGC / standard).
    g = Compact{UInt64, UInt32}([3, 2, 4])
    @test g.max_m == 4
    @test g.m_sum == 9
    @test g.n == 3
    @test g.m_bits == [4, 3, 2]        # sorted non-increasing

    # Edge case: uniform dimensions.
    g2 = Compact{UInt64, UInt32}([3, 3, 3])
    @test g2.max_m == 3
    @test g2.m_sum == 9

    # Edge case: single dimension.
    g3 = Compact{UInt64, UInt32}([5])
    @test g3.max_m == 5

    # Edge case: all zeros (degenerate).
    g4 = Compact{UInt64, UInt32}([0, 0])
    @test g4.max_m == 0
    @test g4.m_sum == 0

    # Error: empty vector.
    @test_throws ArgumentError Compact{UInt64, UInt32}(Int[])

    # Error: negative values.
    @test_throws ArgumentError Compact{UInt64, UInt32}([-1, 2])

    # Error: too many bits for index type.
    @test_throws ArgumentError Compact{UInt8, UInt32}([5, 5])   # 10 bits > 8

    # Error: too many bits for coordinate type.
    @test_throws ArgumentError Compact{UInt64, UInt8}([10])     # 10 > 8

    # Error: illegal (gray, path) combinations.
    @test_throws ArgumentError Compact([1, 1, 1, 1]; gray = :brgc, path = :catalog)
    @test_throws ArgumentError Compact([1, 1, 1, 1]; gray = :random, path = :standard)

    # Error: gray_index out of range at k_max (k_max = 4, gray_count(4) = 10).
    @test_throws ArgumentError Compact([1, 1, 1, 1]; gray = :random, gray_index = 99)

    # Error: path_index out of range at k_max.
    @test_throws ArgumentError Compact([1, 1, 1, 1]; gray = :random, path = :catalog, path_index = 99)

    # Error: :random with k_max > 10 has no catalog entry.
    @test_throws ArgumentError Compact(fill(1, 11); gray = :random)
end


@testitem "Compact convenience constructor type selection" begin
    using BijectiveHilbert: Compact, index_type, axis_type

    @test index_type(Compact([3, 2, 4])) == UInt16   # sum 9 -> UInt16
    @test axis_type(Compact([3, 2, 4])) == UInt8      # max 4 -> UInt8
    @test index_type(Compact([2, 2, 2])) == UInt8     # sum 6 -> UInt8
    @test index_type(Compact([9])) == UInt16          # sum 9 -> UInt16
    @test index_type(Compact(fill(9, 8))) == UInt128  # sum 72 -> UInt128
    @test axis_type(Compact(fill(9, 8))) == UInt16    # max 9 -> UInt16

    # Signed index type argument is converted to unsigned (like CompactHamilton).
    g = Compact(Int64, [3, 2, 4])
    @test index_type(g) == UInt64
end


@testitem "Compact basic encode/decode" begin
    using BijectiveHilbert: Compact, encode_hilbert_zero, decode_hilbert_zero!

    g = Compact{UInt64, UInt32}([3, 3, 3])
    X = UInt32[5, 2, 7]
    Y = zeros(UInt32, 3)
    decode_hilbert_zero!(g, Y, encode_hilbert_zero(g, X))
    @test X == Y

    g2 = Compact{UInt64, UInt32}([3, 2, 4])
    X2 = UInt32[5, 2, 11]
    Y2 = zeros(UInt32, 3)
    decode_hilbert_zero!(g2, Y2, encode_hilbert_zero(g2, X2))
    @test X2 == Y2

    # Degenerate all-zero domain.
    g3 = Compact{UInt64, UInt32}([0, 0])
    @test encode_hilbert_zero(g3, UInt32[0, 0]) == 0
end


@testitem "Compact is its own inverse and complete set (uniform)" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact
    for b in 1:3, n in 2:4
        g = Compact{UInt64, UInt32}(fill(b, n))
        @test HilbertTestSuite.check_own_inverse(g, b, n)
        @test HilbertTestSuite.check_complete_set(g, b, n)
    end
end


@testitem "Compact is its own inverse and complete set (anisotropic)" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact

    test_cases = [
        [2, 5, 3],
        [3, 4, 2, 3],
        [4, 3, 2],
        [1, 2, 3, 4],
        [2, 0, 3],
    ]
    for ms in test_cases
        g = Compact(ms)
        @test HilbertTestSuite.check_own_inverse(g, ms)
        @test HilbertTestSuite.check_complete_set(g, ms)
    end
end


@testitem "Compact curve-family selection" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact

    combos = [(:brgc, :standard), (:brgc, :hub), (:random, :hub), (:random, :catalog)]
    for (gray, path) in combos
        # Uniform b = 1..2, n = 3..4 (random needs k_max = n >= 3). Explicit
        # wide types (like test_compact) so the complete-set bound check does
        # not overflow when the total bit count equals the index-type width.
        for b in 1:2, n in 3:4
            g = Compact{UInt64, UInt32}(fill(b, n); gray = gray, path = path)
            @test HilbertTestSuite.check_own_inverse(g, b, n)
            @test HilbertTestSuite.check_complete_set(g, b, n)
        end
        # One anisotropic case (k_max = 3, with a smaller k = 2 level).
        ms = [3, 2, 2]
        g = Compact{UInt64, UInt32}(ms; gray = gray, path = path)
        @test HilbertTestSuite.check_own_inverse(g, ms)
        @test HilbertTestSuite.check_complete_set(g, ms)
    end
end


@testitem "Compact exhaustive catalog selection at k_max" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact, gray_count, path_count

    # k_max = 3: every catalog gray_index (just one) and every path_index.
    ms3 = [1, 1, 1]
    for gi in 0:(gray_count(3) - 1)
        gh = Compact(ms3; gray = :random, path = :hub, gray_index = gi)
        @test HilbertTestSuite.check_complete_set(gh, ms3)
        for pi in 0:(path_count(3, gi) - 1)
            gc = Compact(ms3; gray = :random, path = :catalog, gray_index = gi, path_index = pi)
            @test HilbertTestSuite.check_complete_set(gc, ms3)
        end
    end

    # k_max = 4: every catalog gray_index (ten) and a few path_index values.
    ms4 = [1, 1, 1, 1]
    for gi in 0:(gray_count(4) - 1)
        gh = Compact(ms4; gray = :random, path = :hub, gray_index = gi)
        @test HilbertTestSuite.check_complete_set(gh, ms4)
        pcount = path_count(4, gi)
        for pi in unique([0, 3, pcount - 1])
            gc = Compact(ms4; gray = :random, path = :catalog, gray_index = gi, path_index = pi)
            @test HilbertTestSuite.check_complete_set(gc, ms4)
        end
    end
end


@testitem "Compact uniform-m matches CompactHamilton" begin
    using BijectiveHilbert: Compact, CompactHamilton, encode_hilbert_zero
    # Empirically, the uniform (isotropic) case produces identical indices to
    # CompactHamilton even though the two use different axis-order / rotation
    # conventions. (Anisotropic indices generally differ - see the docstring.)
    for b in 1:2, n in 2:3
        g = Compact(fill(b, n))
        c = CompactHamilton(fill(b, n))
        for idx in CartesianIndices(ntuple(_ -> 1 << b, n))
            X = UInt16[Tuple(idx)...] .- 1
            @test encode_hilbert_zero(g, X) == encode_hilbert_zero(c, X)
        end
    end
end


@testitem "Compact type interactions" begin
    using BijectiveHilbert: Compact, index_type
    using BijectiveHilbert: encode_hilbert_zero, decode_hilbert_zero!
    using BijectiveHilbert: encode_hilbert, decode_hilbert!
    using Random

    # UInt128 index for sum(m) > 64; spot round-trips (no exhaustive loop).
    g = Compact(fill(9, 8))
    @test index_type(g) == UInt128
    Random.seed!(20260706)
    for _ in 1:64
        X = UInt16.(rand(0:511, 8))
        h = encode_hilbert_zero(g, X)
        @test h isa UInt128
        Y = zeros(UInt16, 8)
        decode_hilbert_zero!(g, Y, h)
        @test X == Y
    end

    # Int-h decode overload converts to the index type.
    gi = Compact(Int64, [3, 2, 4])
    Xi = zeros(Int, 3)
    decode_hilbert_zero!(gi, Xi, 5)                    # Int literal
    @test encode_hilbert_zero(gi, Xi) == index_type(gi)(5)

    # 1-based wrappers.
    g2 = Compact{UInt64, UInt32}([3, 2, 4])
    X1 = UInt32[6, 3, 12]                              # 1-based == 0-based [5,2,11]
    h1 = encode_hilbert(g2, X1)
    @test h1 >= 1
    Y1 = zeros(UInt32, 3)
    decode_hilbert!(g2, Y1, h1)
    @test X1 == Y1
end
