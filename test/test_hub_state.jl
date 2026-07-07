using TestItemRunner


@testitem "hub-state closed-form primitives" begin
    using BijectiveHilbert: brgc_rotated_gray, brgc_rotated_rank,
        glued_child_entry, glued_child_dir, validate_mismatch, axis_of_single_bit

    for k in 1:10
        n = 1 << k
        # brgc_rotated_gray / brgc_rotated_rank are inverse bijections.
        seen = falses(n)
        for w in 0:(n - 1)
            g = brgc_rotated_gray(w, k)
            @test 0 <= g < n
            @test brgc_rotated_rank(g, k) == w
            seen[Int(g) + 1] = true
        end
        @test all(seen)                          # permutation
        @test brgc_rotated_gray(0, k) == 0
        @test brgc_rotated_gray(n - 1, k) == 1

        # The closed-form glued child path produces a valid mismatch.
        gray = UInt64[brgc_rotated_gray(w, k) for w in 0:(n - 1)]
        ce = UInt64[glued_child_entry(w, k) for w in 0:(n - 1)]
        chi = UInt64[gray[i] ⊻ ce[i] for i in 1:n]
        @test validate_mismatch(k, gray, chi) === nothing

        # glued_child_dir agrees with the axis derived from chi.
        for w in 0:(n - 1)
            d = w + 1 < n ? axis_of_single_bit(chi[w + 1] ⊻ chi[w + 2]) :
                axis_of_single_bit(chi[w + 1])
            @test glued_child_dir(w, k) == d % k
        end
    end
end


@testitem "hub-state affine primitives are inverses" begin
    using BijectiveHilbert: glued_affine_apply, glued_affine_apply_inv
    using Random

    rng = MersenneTwister(1234)
    for _ in 1:5000
        k = rand(rng, 1:32)
        d = rand(rng, 0:(k - 1))
        mask = k == 64 ? ~UInt64(0) : (UInt64(1) << k) - UInt64(1)
        e = rand(rng, UInt64) & mask
        x = rand(rng, UInt64) & mask
        y = glued_affine_apply(x, e, d, k)
        @test y <= mask
        @test glued_affine_apply_inv(y, e, d, k) == x
    end
end


@testitem "build_hub_state on brgc_rotated matches closed forms" begin
    using BijectiveHilbert: build_hub_state, brgc_rotated_gray, brgc_rotated_rank,
        validate_mismatch

    for k in 1:10
        n = 1 << k
        gray = UInt64[brgc_rotated_gray(w, k) for w in 0:(n - 1)]
        t = build_hub_state(k, gray)

        @test t.k == k
        @test length(t.gray) == n
        @test length(t.gray_rank) == n
        @test length(t.child_entry) == n
        @test length(t.child_dir) == n

        # gray / gray_rank reproduce the closed forms.
        @test t.gray == gray
        for g in 0:(n - 1)
            @test t.gray_rank[Int(g) + 1] == brgc_rotated_rank(g, k)
        end

        # child_entry starts at the origin and every direction is reduced mod k.
        @test t.child_entry[1] == 0
        @test all(t.child_dir .< k)

        # The derived mismatch is valid.
        chi = UInt64[t.gray[i] ⊻ t.child_entry[i] for i in 1:n]
        @test validate_mismatch(k, t.gray, chi) === nothing
    end
end


@testitem "build_hub_state on a hand-made non-BRGC Gray code" begin
    using BijectiveHilbert: build_hub_state, validate_mismatch

    # A cyclic Gray code that is not the BRGC-rotated one.
    gray = UInt64[0, 4, 6, 2, 3, 7, 5, 1]
    t = build_hub_state(3, gray)
    @test t.gray == gray
    @test t.child_entry[1] == 0
    @test all(t.child_dir .< 3)
    chi = UInt64[t.gray[i] ⊻ t.child_entry[i] for i in 1:8]
    @test validate_mismatch(3, t.gray, chi) === nothing
end


@testitem "hub-state catalog: build_hub_state solves every Gray code" begin
    using BijectiveHilbert: build_hub_state, validate_mismatch, gray_count, gray_for

    for k in 3:10
        for gi in 0:(gray_count(k) - 1)
            gray = gray_for(k, gi)
            @test length(gray) == (1 << k)
            t = build_hub_state(k, gray)
            @test t.k == k
            @test all(t.child_dir .< k)
            chi = UInt64[UInt64(t.gray[i]) ⊻ t.child_entry[i] for i in 1:(1 << k)]
            @test validate_mismatch(k, t.gray, chi) === nothing
        end
    end
end


@testitem "hub-state catalog: build_from_child_entry solves every path" begin
    using BijectiveHilbert: build_from_child_entry, validate_mismatch,
        gray_count, gray_for, path_count, path_child_entry

    for k in 3:10
        for gi in 0:(gray_count(k) - 1)
            gray = gray_for(k, gi)
            for pi in 0:(path_count(k, gi) - 1)
                ce = path_child_entry(k, gi, pi)
                @test length(ce) == (1 << k)
                t = build_from_child_entry(k, gray, ce)
                @test t.k == k
                # build_from_child_entry reproduces the given child_entry.
                @test t.child_entry == UInt64.(ce)
                @test all(t.child_dir .< k)
                chi = UInt64[UInt64(t.gray[i]) ⊻ t.child_entry[i] for i in 1:(1 << k)]
                @test validate_mismatch(k, t.gray, chi) === nothing
            end
        end
    end
end


@testitem "hub-state validation rejects corrupted input" begin
    using BijectiveHilbert: validate_gray, validate_mismatch, build_from_child_entry,
        brgc_rotated_gray, glued_child_entry

    # validate_gray rejects a non-Gray step (differs in two bits).
    @test_throws ArgumentError validate_gray(3, UInt64[0, 3, 1, 2, 6, 7, 5, 4])
    # validate_gray rejects a non-permutation (repeated value).
    @test_throws ArgumentError validate_gray(3, UInt64[0, 1, 3, 2, 6, 7, 5, 5])
    # validate_gray rejects a wrong-length sequence.
    @test_throws ArgumentError validate_gray(3, UInt64[0, 1, 3, 2])

    gray = UInt64[brgc_rotated_gray(w, 3) for w in 0:7]
    ce = UInt64[glued_child_entry(w, 3) for w in 0:7]
    bad = copy(ce)
    bad[3] = bad[3] ⊻ UInt64(1)   # corrupt one entry

    badchi = UInt64[gray[i] ⊻ bad[i] for i in 1:8]
    @test_throws ArgumentError validate_mismatch(3, gray, badchi)
    @test_throws ArgumentError build_from_child_entry(3, gray, bad)

    # A child_entry whose entry corner is nonzero violates the entry constraint.
    badentry = copy(ce)
    badentry[1] = UInt64(1)
    @test_throws ArgumentError build_from_child_entry(3, gray, badentry)
end
