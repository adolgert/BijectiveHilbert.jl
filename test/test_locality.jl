using Test
using BijectiveHilbert
using StaticArrays
using LinearAlgebra
using Random

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

@testset "Spatial Locality Metrics" begin
    @testset "2D Locality" begin
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
        mcd_hilbert = mean_consecutive_distance(pts_hilbert)
        mcd_random = mean_consecutive_distance(pts_random)

        @test mcd_hilbert < mcd_random
        println("2D MCD - Hilbert: $(round(mcd_hilbert, digits=4)), Random: $(round(mcd_random, digits=4))")

        # Test 2: Locality Ratio
        ratio = locality_ratio(pts_hilbert, pts_random)
        @test ratio > 2.0  # Hilbert should be at least 2x better
        println("2D Locality Ratio: $(round(ratio, digits=2))x improvement")

        # Test 3: Bandwidth Sum
        bw_hilbert = bandwidth_sum(pts_hilbert, 10)
        bw_random = bandwidth_sum(pts_random, 10)
        @test bw_hilbert < bw_random
        println("2D Bandwidth(k=10) - Hilbert: $(round(bw_hilbert, digits=2)), Random: $(round(bw_random, digits=2))")

        # Test 4: Distance-Rank Correlation
        corr_hilbert = distance_rank_correlation(pts_hilbert, 200)
        corr_random = distance_rank_correlation(pts_random, 200)
        @test corr_hilbert > 0.3  # Should have positive correlation
        @test corr_hilbert > corr_random
        println("2D Correlation - Hilbert: $(round(corr_hilbert, digits=3)), Random: $(round(corr_random, digits=3))")
    end

    @testset "3D Locality" begin
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
        mcd_hilbert = mean_consecutive_distance(pts_hilbert)
        mcd_random = mean_consecutive_distance(pts_random)

        @test mcd_hilbert < mcd_random
        println("3D MCD - Hilbert: $(round(mcd_hilbert, digits=4)), Random: $(round(mcd_random, digits=4))")

        # Test 2: Locality Ratio
        ratio = locality_ratio(pts_hilbert, pts_random)
        @test ratio > 2.0  # Hilbert should be at least 2x better
        println("3D Locality Ratio: $(round(ratio, digits=2))x improvement")

        # Test 3: Bandwidth Sum
        bw_hilbert = bandwidth_sum(pts_hilbert, 10)
        bw_random = bandwidth_sum(pts_random, 10)
        @test bw_hilbert < bw_random
        println("3D Bandwidth(k=10) - Hilbert: $(round(bw_hilbert, digits=2)), Random: $(round(bw_random, digits=2))")

        # Test 4: Distance-Rank Correlation
        corr_hilbert = distance_rank_correlation(pts_hilbert, 200)
        corr_random = distance_rank_correlation(pts_random, 200)
        @test corr_hilbert > 0.3  # Should have positive correlation
        @test corr_hilbert > corr_random
        println("3D Correlation - Hilbert: $(round(corr_hilbert, digits=3)), Random: $(round(corr_random, digits=3))")
    end

    @testset "Higher Dimensions" begin
        # Test that locality holds in higher dimensions
        for D in [4, 5, 6]
            rng = Random.MersenneTwister(100 + D)
            n = 500
            pts_original = [SVector{D,Float64}(rand(rng, D)...) for _ in 1:n]

            pts_hilbert = copy(pts_original)
            BijectiveHilbert.hilbertsort!(pts_hilbert)

            pts_random = shuffle(rng, copy(pts_original))

            ratio = locality_ratio(pts_hilbert, pts_random)
            @test ratio > 1.5  # At least 1.5x better in higher dimensions
            println("$(D)D Locality Ratio: $(round(ratio, digits=2))x improvement")
        end
    end
end
