
@safetestset hilbert_decode_zero = "hilbert decode zero" begin
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


@safetestset simple2d_own_inverse = "simple2d is its own inverse" begin
  using BijectiveHilbert: Simple2D, check_own_inverse
  gg = Simple2D(Int)
  for b in 2:7
    @test check_own_inverse(gg, b, 2)
  end
end


@safetestset simple2d_complete_set = "simple2d is a complete set" begin
  using BijectiveHilbert: Simple2D, check_complete_set
  gg = Simple2D(Int)
  for b in 2:7
    @test check_complete_set(gg, b, 2)
  end
end
