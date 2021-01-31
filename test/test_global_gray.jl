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
h = BijectiveHilbert.interleave_transpose(UInt64, t, b, n)
@test h == 0b100110011001
end


# Implementation test
@safetestset interleave_outerleave = "outerleave is opposite" begin
using BijectiveHilbert
n = 3
b = 4
X = zeros(UInt8, n)
for h in 0:(1<<(n*b) - 1)
    BijectiveHilbert.outerleave_transpose!(X, UInt64(h), b, n)
    h2 = BijectiveHilbert.interleave_transpose(UInt64, X, b, n)
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
    TT = index_type(gg)
    X = convert(Vector{AT}, rand(rng, 0:(1<<b - 1), n))
    X0 = copy(X)
    h = encode_hilbert_zero(gg, X)
    decode_hilbert_zero!(gg, X, h)
    @test X == X0
end
end


@safetestset hilbert_one_diff = "GlobalGray values next to each other" begin
using BijectiveHilbert: GlobalGray, check_complete_set
n = 3
b = 4
gg = GlobalGray(b, n)
@test check_complete_set(gg, b, n)
end


@safetestset simple2d_against_file = "simple2d agrees with C code" begin
function read_skill(fn)
    lines = readlines(fn)
    xyz = zeros(Int, 3, length(lines))
    hh = zeros(Int, 3, length(lines))
    idx = 0
    for ll in lines
        axmatch = r"^\((\d+), (\d+), (\d+)\)"
        hmatch = r"\) \((\d+), (\d+), (\d+)\)"
        if occursin(axmatch, ll)
            idx += 1
            mm = match(axmatch, ll)
            xyz[:, idx] = parse.(Int, mm.captures)
            hh[:, idx] = parse.(Int, match(hmatch, ll).captures)
        end
    end
    return xyz, hh
end

fn = "test/skill3_4_4_4.txt"
if !isfile(fn)
    fn = basename(fn)
end
if isfile(fn)
    xyz, hh = read_skill(fn)
end

using BijectiveHilbert
for check_idx in 1:size(xyz, 2)
    ax = convert(Vector{UInt8}, xyz[:, check_idx])
    h = convert(Vector{UInt8}, hh[:, check_idx])
    BijectiveHilbert.axes_to_transpose!(ax, 4, 3)
    @test ax == h
end

for check_idx in 1:size(xyz, 2)
    ax = convert(Vector{UInt8}, xyz[:, check_idx])
    h = convert(Vector{UInt8}, hh[:, check_idx])
    BijectiveHilbert.transpose_to_axes!(h, 4, 3)
    @test h == ax
end
end
