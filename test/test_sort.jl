using TestItemRunner
using LinearAlgebra

@testmodule LocalityMetrics begin
    using LinearAlgebra
    using Random
    using Statistics

    export mean_consecutive_distance, locality_ratio, bandwidth_sum, distance_rank_correlation

    """
        mean_consecutive_distance(pts)

    Compute the average Euclidean distance between consecutive points.
    Lower values indicate better spatial locality.
    """
    function mean_consecutive_distance(pts)
        n = length(pts)
        total = 0.0
        for i in 1:(n-1)
            total += norm(pts[i+1] - pts[i])
        end
        total / (n - 1)
    end

    """
        locality_ratio(pts_sorted, pts_random)

    Ratio of MCD for random ordering vs sorted ordering.
    Values >> 1 indicate good spatial locality (typically 2-10x).
    """
    function locality_ratio(pts_sorted, pts_random)
        mcd_sorted = mean_consecutive_distance(pts_sorted)
        mcd_random = mean_consecutive_distance(pts_random)
        mcd_random / mcd_sorted
    end

    """
        bandwidth_sum(pts, k=10)

    Sum of distances for points within index distance k.
    Lower values indicate better local clustering.
    """
    function bandwidth_sum(pts, k=10)
        n = length(pts)
        total = 0.0
        for i in 1:n
            for j in (i+1):min(i+k, n)
                total += norm(pts[i] - pts[j])
            end
        end
        total
    end

    """
        distance_rank_correlation(pts, sample_size=100)

    Pearson correlation between index distance and Euclidean distance.
    Positive values indicate locality preservation.
    """
    function distance_rank_correlation(pts, sample_size=100)
        n = length(pts)
        # Sample pairs to avoid O(n^2) computation
        rng = Random.MersenneTwister(42)
        index_dists = Float64[]
        euclidean_dists = Float64[]

        for _ in 1:sample_size
            i, j = rand(rng, 1:n, 2)
            push!(index_dists, abs(i - j))
            push!(euclidean_dists, norm(pts[i] - pts[j]))
        end

        # Pearson correlation
        cor(index_dists, euclidean_dists)
    end
end


@testitem "sort smoke" begin
    using BijectiveHilbert
    using StaticArrays
    using Random

    rng = Xoshiro(293482)
    n = 10
    pts_original = [SVector{2,Float64}(rand(rng, 2)...) for _ in 1:n]
    bare = copy(pts_original)
    hilbertsort!(bare)
    encoder = Simple2D(UInt128)
    enc = copy(pts_original)
    hilbertsort!(enc; encoder=encoder)
    @test enc == bare
    bits_per_axis = 52
    bits = copy(pts_original)
    hilbertsort!(bits; bits_per_axis=bits_per_axis)
    @test bits == bare
    together = copy(pts_original)
    hilbertsort!(together, encoder=encoder, bits_per_axis=bits_per_axis)
    @test together == bare
end


@testitem "2D Locality" setup=[LocalityMetrics] begin
    using BijectiveHilbert
    using StaticArrays
    using LinearAlgebra
    using Random

    # Generate random 2D points
    rng = Random.MersenneTwister(123)
    n = 1000
    pts_original = [SVector{2,Float64}(rand(rng, 2)...) for _ in 1:n]

    # Hilbert sorted
    pts_hilbert = copy(pts_original)
    BijectiveHilbert.hilbertsort!(pts_hilbert)

    # Random sorted
    pts_random = shuffle(rng, copy(pts_original))

    # Test 1: Mean Consecutive Distance
    mcd_hilbert = LocalityMetrics.mean_consecutive_distance(pts_hilbert)
    mcd_random = LocalityMetrics.mean_consecutive_distance(pts_random)

    @test mcd_hilbert < mcd_random

    # Test 2: Locality Ratio
    ratio = LocalityMetrics.locality_ratio(pts_hilbert, pts_random)
    @test ratio > 2.0  # Hilbert should be at least 2x better

    # Test 3: Bandwidth Sum
    bw_hilbert = LocalityMetrics.bandwidth_sum(pts_hilbert, 10)
    bw_random = LocalityMetrics.bandwidth_sum(pts_random, 10)
    @test bw_hilbert < bw_random

    # Test 4: Distance-Rank Correlation
    corr_hilbert = LocalityMetrics.distance_rank_correlation(pts_hilbert, 200)
    corr_random = LocalityMetrics.distance_rank_correlation(pts_random, 200)
    @test corr_hilbert > 0.3  # Should have positive correlation
    @test corr_hilbert > corr_random
end


@testitem "3D Locality" setup=[LocalityMetrics] begin
    using BijectiveHilbert
    using StaticArrays
    using LinearAlgebra
    using Random

    # Generate random 3D points
    rng = Random.MersenneTwister(456)
    n = 1000
    pts_original = [SVector{3,Float64}(rand(rng, 3)...) for _ in 1:n]

    # Hilbert sorted
    pts_hilbert = copy(pts_original)
    BijectiveHilbert.hilbertsort!(pts_hilbert)

    # Random sorted
    pts_random = shuffle(rng, copy(pts_original))

    # Test 1: Mean Consecutive Distance
    mcd_hilbert = LocalityMetrics.mean_consecutive_distance(pts_hilbert)
    mcd_random = LocalityMetrics.mean_consecutive_distance(pts_random)

    @test mcd_hilbert < mcd_random

    # Test 2: Locality Ratio
    ratio = LocalityMetrics.locality_ratio(pts_hilbert, pts_random)
    @test ratio > 2.0  # Hilbert should be at least 2x better

    # Test 3: Bandwidth Sum
    bw_hilbert = LocalityMetrics.bandwidth_sum(pts_hilbert, 10)
    bw_random = LocalityMetrics.bandwidth_sum(pts_random, 10)
    @test bw_hilbert < bw_random

    # Test 4: Distance-Rank Correlation
    corr_hilbert = LocalityMetrics.distance_rank_correlation(pts_hilbert, 200)
    corr_random = LocalityMetrics.distance_rank_correlation(pts_random, 200)
    @test corr_hilbert > 0.3  # Should have positive correlation
    @test corr_hilbert > corr_random
end


@testitem "Higher Dimensions" setup=[LocalityMetrics] begin
    using BijectiveHilbert
    using StaticArrays
    using LinearAlgebra
    using Random

    # Test that locality holds in higher dimensions
    for D in [4, 5, 6]
        rng = Random.MersenneTwister(100 + D)
        n = 500
        pts_original = [SVector{D,Float64}(rand(rng, D)...) for _ in 1:n]

        pts_hilbert = copy(pts_original)
        BijectiveHilbert.hilbertsort!(pts_hilbert)

        pts_random = shuffle(rng, copy(pts_original))

        ratio = LocalityMetrics.locality_ratio(pts_hilbert, pts_random)
        @test ratio > 1.5  # At least 1.5x better in higher dimensions
    end
end
