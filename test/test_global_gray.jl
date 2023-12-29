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
using BijectiveHilbert: GlobalGray
using ..HilbertTestSuite: check_complete_set
n = 3
b = 4
gg = GlobalGray(b, n)
@test check_complete_set(gg, b, n)
end


@safetestset globalgray_against_file = "globalgray agrees with C code" begin
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
end


@safetestset globalgray_type_interactions = "GlobalGray type interactions" begin
    using BijectiveHilbert
    using UnitTestDesign
    using Random
    rng = Random.MersenneTwister(9790323)
    for retrial in 1:5
        AxisTypes = shuffle(rng, [Int8, Int, UInt, Int128, UInt8, UInt128])
        IndexTypes = shuffle(rng, [Union{}, Int8, UInt8, Int, UInt, Int128, UInt128])
        Count= shuffle(rng, [0, 1])
        Dims = shuffle(rng, [2, 3, 4])
        Bits = shuffle(rng, [2, 3, 4, 5])
        test_set = all_pairs(
            AxisTypes, IndexTypes, Count, Dims, Bits;
        )
        for (A, I, C, D, B) in test_set
            if I == Union{}
                gg = GlobalGray(B, D)
                I = index_type(gg)
            else
                gg = GlobalGray(I, B, D)
            end
            if B * D > log2(typemax(I))
                continue
            end
            last = (one(I) << (B * D)) - one(I) + I(C)
            mid = one(I) << (B * D - 1)
            few = 5
            X = zeros(A, D)
            hlarr = vcat(C:min(mid, few), max(mid + 1, last - few):last)
            for hl in hlarr
                hli = I(hl)
                if C == 0
                    decode_hilbert_zero!(gg, X, hli)
                    hl2 = encode_hilbert_zero(gg, X)
                    @test hl2 == hli
                    @test typeof(hl2) == typeof(hli)
                else
                    decode_hilbert!(gg, X, hli)
                    hl2 = encode_hilbert(gg, X)
                    @test hl2 == hli
                    @test typeof(hl2) == typeof(hli)
                end
            end
        end
    end
end


@safetestset globalgray_own_inv = "GlobalGray is own inverse" begin
    using BijectiveHilbert: GlobalGray
    using ..HilbertTestSuite: check_own_inverse
    for n in 2:5
        for b in 2:4
            gg = GlobalGray(b, n)
            @test check_own_inverse(gg, b, n)
        end
    end
end


@safetestset globalgray_complete_set = "GlobalGray is complete set" begin
    using BijectiveHilbert: GlobalGray
    using ..HilbertTestSuite: check_complete_set
    for n in [2, 3, 5]
        for b in [2, 3]
            fc = GlobalGray(b, n)
            @test check_complete_set(fc, b, n)
        end
    end
    end
    