using TestItemRunner


@testitem "Compact constructor validation" begin
    using BijectiveHilbert: Compact

    # Valid construction
    c = Compact{UInt64, UInt32}([3, 2, 4])
    @test c.mmax == 4
    @test c.total_bits == 9
    @test length(c.m) == 3

    # Edge case: uniform dimensions
    c2 = Compact{UInt64, UInt32}([3, 3, 3])
    @test c2.mmax == 3
    @test c2.total_bits == 9

    # Edge case: single dimension
    c3 = Compact{UInt64, UInt32}([5])
    @test c3.mmax == 5
    @test c3.total_bits == 5

    # Edge case: all zeros (degenerate)
    c4 = Compact{UInt64, UInt32}([0, 0])
    @test c4.mmax == 0
    @test c4.total_bits == 0

    # Error: empty vector
    @test_throws ArgumentError Compact{UInt64, UInt32}(Int[])

    # Error: negative values
    @test_throws ArgumentError Compact{UInt64, UInt32}([-1, 2])

    # Error: too many bits for index type
    @test_throws ArgumentError Compact{UInt8, UInt32}([5, 5])  # 10 bits > 8

    # Error: too many bits for coordinate type
    @test_throws ArgumentError Compact{UInt64, UInt8}([10])  # 10 > 8
end


@testitem "Compact basic encode/decode" begin
    using BijectiveHilbert: Compact, encode_hilbert_zero, decode_hilbert_zero!

    # Uniform case
    c = Compact{UInt64, UInt32}([3, 3, 3])
    X = UInt32[5, 2, 7]
    h = encode_hilbert_zero(c, X)
    Y = zeros(UInt32, 3)
    decode_hilbert_zero!(c, Y, h)
    @test X == Y

    # Anisotropic case
    c2 = Compact{UInt64, UInt32}([3, 2, 4])
    X2 = UInt32[5, 2, 11]
    h2 = encode_hilbert_zero(c2, X2)
    Y2 = zeros(UInt32, 3)
    decode_hilbert_zero!(c2, Y2, h2)
    @test X2 == Y2

    # Edge case: single dimension
    c3 = Compact{UInt64, UInt32}([4])
    X3 = UInt32[9]
    h3 = encode_hilbert_zero(c3, X3)
    Y3 = zeros(UInt32, 1)
    decode_hilbert_zero!(c3, Y3, h3)
    @test X3 == Y3

    # Edge case: all zeros input
    c4 = Compact{UInt64, UInt32}([2, 2])
    X4 = UInt32[0, 0]
    h4 = encode_hilbert_zero(c4, X4)
    Y4 = zeros(UInt32, 2)
    decode_hilbert_zero!(c4, Y4, h4)
    @test X4 == Y4
end


@testitem "Compact is its own inverse (uniform)" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact
    for b in 1:4, n in 2:4
        c = Compact{UInt64, UInt32}(fill(b, n))
        @test HilbertTestSuite.check_own_inverse(c, b, n)
    end
end


@testitem "Compact is its own inverse (anisotropic)" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact

    # Test various anisotropic configurations
    test_cases = [
        [3, 2, 4],
        [2, 3],
        [1, 2, 3],
        [4, 2, 3, 1],
        [2, 2, 3],
        [1, 1, 1, 1],
    ]

    for ms in test_cases
        c = Compact{UInt64, UInt32}(ms)
        @test HilbertTestSuite.check_own_inverse(c, ms)
    end
end


@testitem "Compact is complete set (uniform)" setup=[HilbertTestSuite] begin
    using BijectiveHilbert: Compact
    for n in 2:3, b in 2:3
        c = Compact{UInt64, UInt32}(fill(b, n))
        @test HilbertTestSuite.check_complete_set(c, b, n)
    end
end


@testitem "Compact type stability" begin
    using BijectiveHilbert
    using BijectiveHilbert: Compact, encode_hilbert_zero, decode_hilbert_zero!
    using BijectiveHilbert: affine_apply, affine_apply_inv, child_entry, child_dir, embed_state

    # Test helper functions
    @test @inferred(affine_apply(UInt32(5), UInt32(3), 2, 4)) isa UInt32
    @test @inferred(affine_apply_inv(UInt32(5), UInt32(3), 2, 4)) isa UInt32
    @test @inferred(child_entry(UInt32(5), 4)) isa UInt32
    @test @inferred(child_dir(5, 4)) isa Int
    @test @inferred(child_dir(UInt32(5), 4)) isa Int

    # Test embed_state with vectors
    A_old = [1, 3]
    pos_new = [2, -1, 1]
    e_old, d_old = UInt32(0b01), 0
    result = @inferred(embed_state(A_old, 2, pos_new, e_old, d_old))
    @test result isa Tuple{UInt32, Int}

    # Test embed_state with views (as used in actual code)
    axes_level = [1 1; 3 2; 0 3]
    pos_level = [2 1; -1 2; 1 3]
    A_view = @view axes_level[1:2, 1]
    pos_view = @view pos_level[:, 2]
    result2 = @inferred(embed_state(A_view, 2, pos_view, UInt32(1), 0))
    @test result2 isa Tuple{UInt32, Int}

    # Test top-level encode/decode
    c = Compact{UInt64, UInt32}([3, 2, 4])
    X = UInt32[5, 2, 11]
    @test @inferred(encode_hilbert_zero(c, X)) isa UInt64

    Y = zeros(UInt32, 3)
    h = UInt64(42)
    @inferred decode_hilbert_zero!(c, Y, h)
end


@testitem "Compact 1-based API" begin
    using BijectiveHilbert: Compact, encode_hilbert, decode_hilbert!

    c = Compact{UInt64, UInt32}([3, 2, 4])

    # 1-based coordinates and index
    X = UInt32[6, 3, 12]  # 1-based: corresponds to 0-based [5, 2, 11]
    h = encode_hilbert(c, X)
    @test h >= 1  # 1-based index

    Y = zeros(UInt32, 3)
    decode_hilbert!(c, Y, h)
    @test X == Y
end


@testitem "Compact with different type combinations" begin
    using BijectiveHilbert: Compact, encode_hilbert_zero, decode_hilbert_zero!

    # UInt128 index for large grids
    c1 = Compact{UInt128, UInt64}([20, 20, 20])
    X1 = UInt64[12345, 67890, 11111]
    h1 = encode_hilbert_zero(c1, X1)
    Y1 = zeros(UInt64, 3)
    decode_hilbert_zero!(c1, Y1, h1)
    @test X1 == Y1

    # UInt32 index, UInt16 coordinates
    c2 = Compact{UInt32, UInt16}([4, 4])
    X2 = UInt16[10, 5]
    h2 = encode_hilbert_zero(c2, X2)
    Y2 = zeros(UInt16, 2)
    decode_hilbert_zero!(c2, Y2, h2)
    @test X2 == Y2
end
