@safetestset compare_with_cpp = "the c++ implementation is identical" begin
if isfile("test/libskilling.so")
function TransposetoAxes!(X::Vector{UInt}, b::Int, n::Int)
    ccall(
        (:TransposetoAxes, "test/libskilling"),
        Cvoid,
        (Ref{UInt}, Int, Int),
        X, b, n
        )
end

function AxestoTranspose!(X::Vector{UInt32}, b::Int32, n::Int32)
    ccall(
        (:AxestoTranspose, "test/libskilling"),
        Cvoid,
        (Ref{UInt32}, Int32, Int32),
        X, b, n
        )
end

using BijectiveHilbert: axes_to_transpose!, transpose_to_axes!
using Random
rng = MersenneTwister(29472349)
for i in 1:100
    n = Int32(rand(rng, 2:6))
    b = Int32(rand(rng, 2:8))

    X = convert(Vector{UInt32}, rand(rng, 0:(1<<b - 1), n))
    H = copy(X)
    axes_to_transpose!(H, b, n)
    XC = copy(X)
    AxestoTranspose!(XC, b, n)
    @test H == XC
end
end
end


# Implementation test
@safetestset interleave_single = "interleave single example" begin
using BijectiveHilbert
t = UInt8[0b1100, 0b0110, 0b0011]
b = 4
n = length(t)
h = BijectiveHilbert.interleave_transpose(t, b, n)
@test h == 0b001011110100
end


# Implementation test
@safetestset interleave_outerleave = "outerleave is opposite" begin
using BijectiveHilbert
n = 3
b = 4
X = zeros(UInt8, n)
for h in 0:(1<<b - 1)
    BijectiveHilbert.outerleave_transpose!(X, UInt64(h), b, n)
    h2 = BijectiveHilbert.interleave_transpose(X, b, n)
    @test h == h2
end
end


@safetestset inverse_of_itself = "GlobalGray is its own inverse" begin
using BijectiveHilbert
using Random
rng = MersenneTwister(29472349)
for i in 1:1000
    n = rand(rng, 2:6)
    b = rand(rng, 2:8)
    gg = GlobalGray(b, n)
    AT = axis_type(gg)
    TT = transpose_type(gg)
    X = convert(Vector{AT}, rand(rng, 0:(1<<b - 1), n))
    X0 = copy(X)
    h = encode_hilbert_zero!(gg, X)
    decode_hilbert_zero!(gg, X, h)
    @test X == X0
end
end


@safetestset hilbert_one_diff = "GlobalGray values next to each other" begin
  using BijectiveHilbert
  xy = Set(Tuple{Int64, Int64}[])
  last = [-1, 0]
  n = 3
  b = 5
  gg = GlobalGray(b, n)
  A = axis_type(gg)
  TT = transpose_type(gg)
  X = zeros(A, n)
  Y = copy(X)
  for h in 0:(1<<(n*b) - 1)
    decode_hilbert_zero!(gg, X, TT(h))
    tdiff = UInt64(0)
    for cmp_idx in 1:n
        if X[cmp_idx] > Y[cmp_idx]
            tdiff += X[cmp_idx] - Y[cmp_idx]
        else
            tdiff += Y[cmp_idx] - X[cmp_idx]
        end
    end
    @test tdiff == 1
    if tdiff > 1
        @show h, X, Y
        break
    end
    Y .= X
  end
end
