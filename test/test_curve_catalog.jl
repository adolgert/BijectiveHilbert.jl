using TestItemRunner


@testitem "Curve catalog Gray code counts" begin
    using BijectiveHilbert: gray_count

    @test gray_count(3) == 1
    for k in 4:10
        @test gray_count(k) == 10
    end
    @test gray_count(2) == 0
    @test gray_count(11) == 0
end


@testitem "Curve catalog path counts" begin
    using BijectiveHilbert: gray_count, path_count

    # The single k=3 Gray code has 5 catalogued paths.
    @test path_count(3, 0) == 5

    # For k >= 4 each of the 10 Gray codes has 10 catalogued paths.
    for k in 4:10
        for gi in 0:(gray_count(k) - 1)
            @test path_count(k, gi) == 10
        end
    end

    # Out-of-range (k, gi) has no paths.
    @test path_count(2, 0) == 0
    @test path_count(3, 1) == 0
end


@testitem "Curve catalog gray_for spot check and range errors" begin
    using BijectiveHilbert: gray_for

    @test gray_for(3, 0) == UInt16[0, 4, 6, 2, 3, 7, 5, 1]

    @test_throws ArgumentError gray_for(2, 0)
    @test_throws ArgumentError gray_for(3, 1)
    @test_throws ArgumentError gray_for(11, 0)
end


@testitem "Curve catalog every Gray code is a valid cyclic Gray code" begin
    using BijectiveHilbert: gray_count, gray_for

    for k in 3:10
        for gi in 0:(gray_count(k) - 1)
            gray = gray_for(k, gi)
            n = 1 << k
            @test length(gray) == n
            @test gray[1] == 0
            @test length(Set(gray)) == n
            @test all(g -> g < n, gray)
            # Hamming-1 steps including the wraparound last -> first.
            for i in 1:n
                a = gray[i]
                b = gray[i == n ? 1 : i + 1]
                @test count_ones(a ⊻ b) == 1
            end
        end
    end
end


@testitem "Curve catalog path_child_entry shape and range errors" begin
    using BijectiveHilbert: gray_count, path_count, path_child_entry

    for k in 3:10
        n = 1 << k
        for gi in 0:(gray_count(k) - 1)
            for pi in 0:(path_count(k, gi) - 1)
                entry = path_child_entry(k, gi, pi)
                @test length(entry) == n
                @test all(e -> e < n, entry)
            end
        end
    end

    @test_throws ArgumentError path_child_entry(3, 0, 99)
    @test_throws ArgumentError path_child_entry(2, 0, 0)
end


@testitem "Curve catalog blob integrity" begin
    using BijectiveHilbert
    using SHA

    blob_path = joinpath(dirname(pathof(BijectiveHilbert)), "..", "data", "curve_paths.bin")
    bytes = read(blob_path)

    @test length(bytes) == BijectiveHilbert.CURVE_PATHS_BYTES
    @test bytes2hex(SHA.sha256(bytes)) == BijectiveHilbert.CURVE_PATHS_SHA256
end
