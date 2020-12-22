using Jute
import Random: MersenneTwister, rand!

@testcase "hilbert encode zero" begin
  zz = Set(Int64[])
  for x in 0:15
    for y in 0:15
      z = encode_hilbert_zero(x, y)
      # @show (x, y, z)
      @test z ∉ zz
      @test 0 ≤ z < 256
      push!(zz, z)
    end
  end
end

@testcase "gonna fail" begin
@test false
end

@testcase "hilbert decode zero" begin
  table2 = [
    15 12 11 10
    14 13  8  9
    1  2  7  6
    0  3  4  5
  ]
  table2 = table2[end:-1:1, :]

  table2ind = CartesianIndices(table2)
  for z in 0:15
    x, y = decode_hilbert_zero(z)
    found = table2ind[table2 .== z][1]
    @test found[2] - 1 == x
    @test found[1] - 1 == y
  end
end

@testcase "hilbert all values seen" begin
  xy = Set(Tuple{Int64, Int64}[])
  last = [-1, 0]
  for z in 0:255
    x, y = decode_hilbert_zero(z)
    @test (x, y) ∉ xy
    push!(xy, (x, y))
    # The next point in the grid must be beside the previous.
    # This is a unique feature of this version of the Hilbert curve.
    @test ((x - last[1])^2 == 1) || ((y - last[2])^2 == 1)
    last = [x, y]
  end
  for i in 0:15
    for j in 0:15
      @test (i, j) in xy
    end
  end
end


@testcase "hilbert is bijective" begin
  for z in 1:400
    x, y = decode_hilbert_zero(z)
    zz = encode_hilbert_zero(x, y)
    @test zz == z
  end
end


@testcase "hilbert order makes consecutive points near each other" begin
  point_to_point_distance(x) = sum(sqrt.(sum((x[:, 1:end-1] .- x[:, 2:end]).^2, dims = 1)))
  rng = MersenneTwister(984720987)
  points_in_space = zeros(2, 100)
  rand!(rng, points_in_space)
  points_reordered = points_in_space[:, hilbert_order(points_in_space, 50)]
  @test point_to_point_distance(points_reordered) < 0.25 * point_to_point_distance(points_in_space)
end