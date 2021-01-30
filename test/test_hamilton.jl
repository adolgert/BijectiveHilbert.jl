
@safetestset hi_g_definition = "g is the direction of gray code" begin
using BijectiveHilbert: brgc, hi_g
for i in 0:1000
    @test brgc(i) ⊻ brgc(i + 1) == 1<<hi_g(i)
    @test brgc(i) == brgc(i + 1) ⊻ 1<<hi_g(i)  # xor. It's how it works.
end
end


@safetestset binary_reflection_of_gc = "gray code is binary reflected" begin
using BijectiveHilbert: brgc
n = 10
for i in 0:(1<<n - 1)
    @test brgc(1<<n - 1 - i) == brgc(i) ⊻ 1<<(n-1)
end
end


@safetestset hi_g_matches_long_form = "hi_g is the same as its definition" begin
using BijectiveHilbert: hi_g, hi_g_orig
for i = 0:100
    @test hi_g(i) == hi_g_orig(i)
end
end


@safetestset hi_g_symmetry = "hi_g symmetry invariant" begin
using BijectiveHilbert: brgc, hi_g, bshow
# n = 3
# for i in 0:(1<<n)
#     println(
#         bshow(brgc(i)), " ",
#         bshow(brgc(i + 1)), " ",
#         bshow(brgc(i) ⊻ brgc(i + 1)), " ",
#         bshow(1<<hi_g(i)), " ",
#         bshow(i)
#     )
# end

# Test for the hi_g function.
n = 5
for i in 0:(1<<n - 2)
    # println(hi_g(i), " ", hi_g(1<<n - 2 - i))
    @test(hi_g(i) == hi_g(1<<n - 2 - i))
end
end


@safetestset efdg_match_table = "directions match figure 5 in hamilton report" begin
using BijectiveHilbert: hi_e, hi_f, hi_d, hi_g
n = 2
trials = [
    #i     e     f  d  g
    [0, 0b00, 0b01, 0, 0],
    [1, 0b00, 0b10, 1, 1],
    [2, 0b00, 0b10, 1, 0],
    [3, 0b11, 0b10, 0, 2]
]
for (i, e, f, d, g) in trials
    @test hi_e(i) == Int(e)
    @test hi_f(i, n) == Int(f)
    @test hi_d(i, n) == d
    if g != 42
        @test hi_g(i) == g
    end
end
end


# The 2d example in the paper has a simpler 1d equivalent,
# which is two segments in a row. d=0 for both. f is 1.
@safetestset efdg_match_1d = "directions match a one-dimensional equivalent" begin
using BijectiveHilbert: hi_e, hi_f, hi_d, hi_g
n = 1
trials = [
    #i  e  f  d  g
    [0, 0, 1, 0, 0],  # the first segment should match n=2
    [1, 0, 1, 0, 1]  # the second segment will be redirected for n=2
]
for (i, e, f, d, g) in trials
    @test hi_e(i) == e
    @test hi_f(i, n) == f
    @test hi_d(i, n) == d
    @test hi_g(i) == g
end
end


# A 3d example.
@safetestset efdg_match_3d = "directions match a three-dimensional equivalent" begin
using BijectiveHilbert: hi_e, hi_f, hi_d, hi_g
n = 3
trials = [
    #i      e      f  d  g
    [0, 0b000, 0b001, 0, 0], # segments 0-2 match n=2
    [1, 0b000, 0b010, 1, 1],
    [2, 0b000, 0b010, 1, 0],
    [3, 0b011, 0b111, 2, 2], # segment 3 gets redirected in z-direction
    [4, 0b011, 0b111, 2, 0], # segment 4 also in z direction
    [5, 0b110, 0b100, 1, 1], # then we repeat the 2d plane in reverse.
    [6, 0b110, 0b100, 1, 0],
    [7, 0b101, 0b100, 0, 3]
]
for (i, e, f, d, g) in trials
    @test hi_e(i) == Int(e)
    @test hi_f(i, n) == Int(f)
    @test hi_d(i, n) == d
    @test hi_g(i) == g
end
end


@safetestset hi_e_matches_float_version = "hi_e matches its floating point version" begin
using BijectiveHilbert: hi_e, hi_e_orig
# Show our implementation for hi_e matches the specification.
for i = 1:100
    @test(hi_e(i) == hi_e_orig(i))
end
end


@safetestset cubes_neighbor_along_g_coord = "cubes are neighbors along gth coordinate" begin
# From page 12, just below the figure. The definition of the directions.
using BijectiveHilbert: hi_e, hi_f, hi_g
for n in 2:5
    # Goes to 2^n-2, not 2^n-1 because there is no neighbor for the last point.
    for i in 0:(1<<n - 2)
        fside = hi_f(i, n) ⊻ (1<<hi_g(i))
        eside = hi_e(i + 1)
        @test fside == eside
        if fside != eside
            @show n, i, fside, eside
            break
        end
    end
end
end


@safetestset hi_d_symmetry_invariant = "hi_d symmetry invariant works corollary 2.7" begin
using BijectiveHilbert: hi_d, bshow
# n = 3
# for i = 0:(1<<n)
#     println(i, " ", bshow(hi_d(i, n)))
# end

# invariant for d on page 12. Corollary 2.7.
# Does not hold true for i=0 and i=1<<n - 1.
n = 5
for i = 0:(1<<n - 1)
    #println(hi_d(i, n), " ", hi_d((1<<n)-1-i, n))
    d = hi_d(i, n)
    dflip = hi_d((1<<n) - 1 - i, n)
    @test(d == dflip)
end
end


@safetestset hi_d_compare_g = "hi_d equals hi_g at 2n-1 lemma 2.8" begin
using BijectiveHilbert: hi_d, hi_g
for n = 1:5
    d = hi_d(1<<n - 1, n)
    g = hi_g(1<<n - 1)
    @test d % n == g % n
end
end


@safetestset hi_d_compares_with_g = "hi_d equals g for mod 2" begin
using BijectiveHilbert: hi_d, hi_g
for n = 2:5
    for i = 0:(1<<n - 1)
        d = hi_d(i, n)
        if i == 0
            @test d == 0
        elseif i % 2 == 0
            g = hi_g(i - 1)
            @test d == (g % n)
        else
            g = hi_g(i)
            @test d == (g % n)
        end
    end
end
end


@safetestset hi_edg_invariant = "hi e, d, g invariant is consistent" begin
using BijectiveHilbert: hi_d, hi_e, hi_g
# invariant for e, d, and g. pg. 11. Eq. 1.
# Contrary to the tech report, this is NOT true for 2^n-1.
n = 5
for i = 0:(1<<n - 2)
    @test(hi_e(i + 1) == hi_e(i) ⊻ (1 << hi_d(i, n)) ⊻ (1 << hi_g(i)))
end
end


@safetestset hi_e_reflection = "hi_e reflects to f" begin
using BijectiveHilbert: hi_d, hi_e, hi_f
n = 5
for i = 0:(1<<n - 1)
    e = hi_e(i)
    d = hi_d(i, n)
    f = hi_f(i, n)
    # Corollary 2.7, pg 12, states that e ⊻ d = f, but that's not true.
    # Page 11, first paragraph, says e(i) \xor f(i) = 2^d(i). That's this.
    @test e ⊻ (1<<d) == f
end
end


@safetestset rotation_invariant_zero = "rotation invariant on hi_e is zero" begin
using BijectiveHilbert: hi_T, hi_e, hi_d
# lemma 2.11, page 15
# Assert T_{e,d}(e) == 0
n = 3
for i = 0:(1<<n - 1)
    # println(hi_T(hi_e(i), hi_d(i, n), hi_e(i), n))
    @test(0 == hi_T(hi_e(i), hi_d(i, n), hi_e(i), n))
end
end


# "rotation invariant on hi_f is power of two" begin
# lemma 2.11, page 15
# Assert T_{e,d}(f) == 2^(n-1)
# Contrary to the tech report, this is NOT true.
# @test(0x1<<(n-1) == hi_T(hi_f(i, n), hi_d(i, n), hi_e(i), n))


# Top of page 16, stated invariant:
# (T_(e,d)(a) rotleft (d + 1)) ⊻ e == a
# Also NOT true


@safetestset ef_invariant = "ef relationship holds" begin
using BijectiveHilbert: hi_e, hi_f
# invariant bottom of pg. 11. Lemma 2.6.
# fails for i=0 and i=2^n
n = 5
for i = 0b1:(0b1<<(n - 1))
    e = hi_e(i)
    ff = hi_f((0b1<<n) - 0b1 - i, n) ⊻ (0b1<<(n-1))
    # println(e, " ", ff)
    @test(e == ff)
    @test typeof(e) == UInt8
    @test typeof(ff) == UInt8
end
end


@safetestset inverse_t_itself = "transform is its own inverse" begin
using BijectiveHilbert: hi_d, hi_e, hi_T, hi_T_inv
# Check that the inverse of T is an inverse.
n = 3
for i = 0x0:(0x1<<n - 0x1)
    for b = 0x0:(0x1<<n - 0x1)
        @test(typeof(i) == UInt8)
        d = UInt8(hi_d(i, n))
        e = UInt8(hi_e(i))
        a = hi_T(b, d, e, n)
        b1 = hi_T_inv(a, d, e, n)
        # println(join(string.(typeof.((i, b, a, b1))), " "))
        # println("b ", bitstring(b), " a ", bitstring(a), " b1 ", bitstring(b1))
        @test(b == b1)
    end
end
end


@safetestset get_ith_bit_trials = "trials show can get ith bit of vectors" begin
using BijectiveHilbert: ith_bit_of_indices, get_location

v = [0b01110100, 0b01101010, 0b01011001]
trials = [
    [0, "00000100"],
    [1, "00000010"],
    [2, "00000001"],
    [3, "00000110"],
    [4, "00000101"],
    [5, "00000011"],
    [6, "00000111"],
    [7, "00000000"]
]
for (idx, t0) in trials
    @test bitstring(ith_bit_of_indices(3, v, idx)) == t0
    # Look, Hamilton's version matches.
    @test bitstring(get_location(v, idx)) == t0
end
end


@safetestset hilbert_index_complete4d = "hilbert index is a complete set for 4d" begin
using BijectiveHilbert: SpaceGray, check_complete_set
b = 3
n = 4
gg = SpaceGray(b, n)
@test check_complete_set(gg, b, n)
end


# @safetestset hilbert_index_paper_tr = "paper and techreport agree" begin
# using BijectiveHilbert: hilbert_index, hilbert_index_paper
# They do NOT agree.
# m = 0x4  # One larger won't fit into a byte and must fail.
# for i in 0:(1<<m - 1)
#     for j in 0:(1<<m - 1)
#         h_tr = hilbert_index(0x2, m, UInt8[i, j])
#         h_p = hilbert_index_paper(0x2, m, UInt8[i, j])
#         @test typeof(h_tr) == typeof(h_p)
#         @test h_p == h_tr
#         if h_p != h_tr
#             break
#         end
#     end
# end
# end


@safetestset libhilbert_matches = "libhilbert output matches" begin
using BijectiveHilbert: SpaceGray, encode_hilbert_zero
n = 4
b = 5
gg = SpaceGray(b, n)
trials = [
    [10, 1, 1, 1, 1],
    [11, 1, 0, 1, 1],
    [12, 1, 0, 1, 0],
    [13, 1, 1, 1, 0],
    [14, 1, 1, 0, 0],
    [15, 1, 0, 0, 0],
    [16, 2, 0, 0, 0],
    [17, 2, 0, 1, 0],
    [18, 2, 0, 1, 1],
    [19, 2, 0, 0, 1]
]
for trial in trials
    v = convert(Vector{UInt8}, trial[2:end])
    result = Int(encode_hilbert_zero(gg, v))
    @test result == trial[1]
end
end


@safetestset correct_bits_set_in_indices = "indices bits match" begin
using BijectiveHilbert: set_indices_bits!
m = 3
i = 4
l = 0b110
p = zeros(Int, 3)
set_indices_bits!(p, l, m, i)
@test p == [0, 1<<i, 1<<i]
end


@safetestset paper_inv = "paper hilbert is own inverse" begin
using BijectiveHilbert: check_own_inverse, SpaceGray
for n in 2:5
    for b in 2:4
        gg = SpaceGray(b, n)
        @test check_own_inverse(gg, b, n)
    end
end
end


@safetestset brgc_mask_equal = "mask calculation matches paper" begin
using BijectiveHilbert: extract_mask, extract_mask_paper, hi_d
using Random
rng = MersenneTwister(974073242)
for trial in 1:1000
    n = rand(rng, 2:7)
    m = rand(rng, 2:5, n)
    i = rand(rng, 0:(n - 1))
    d = hi_d(i, n)
    a = extract_mask(m, n, d, i)
    b = extract_mask_paper(m, n, d, i)
    @test a == b
end
end


@safetestset libhilbert_hc_matches = "libhilbert compressed output matches" begin
using BijectiveHilbert: encode_hilbert_zero, Compact
n = 4
b = 5
trials = [
    [10, 1, 1, 1, 1],
    [11, 1, 0, 1, 1],
    [12, 1, 0, 1, 0],
    [13, 1, 1, 1, 0],
    [14, 1, 1, 0, 0],
    [15, 1, 0, 0, 0],
    [16, 2, 0, 0, 0],
    [17, 2, 0, 1, 0],
    [18, 2, 0, 1, 1],
    [19, 2, 0, 0, 1]
]
for trial in trials
    v = convert(Vector{UInt8}, trial[2:end])
    ms = [b, b, b, b]
    gg = Compact(ms)
    result = Int(encode_hilbert_zero(gg, v))
    @test result == trial[1]
end
end


@safetestset compact_index_is_its_inverse = "compact index is its inverse" begin
using BijectiveHilbert: Compact, check_own_inverse
using Random
rng = MersenneTwister(432479874)
for n = [5]
    for i in 1:1
        ms = rand(rng, 2:5, n)
        gg = Compact(ms)
        @test check_own_inverse(gg, ms, n)
    end
end
end


@safetestset compact_index_matches_examples = "compact index matches examples" begin
include(joinpath(dirname(@__FILE__), "check_hamilton.jl"))
using BijectiveHilbert: Compact, decode_hilbert_zero!, index_type, axis_type
fns = ["test/3_2_5_7.txt", "test/2_5_3.txt", "test/3_4_2_3.txt",
       "test/3_5_2_3.txt", "test/4_5_5_5.txt"]
for fn in fns
    if !isfile(fn)
        fn = basename(fn)
    end
    if isfile(fn)
        h, X, ms = hamilton_example(fn)
        n = size(X, 1)
        gg = Compact(ms)
        p = zeros(axis_type(gg), n)
        H = index_type(gg)
        for i in eachindex(h)
            decode_hilbert_zero!(gg, p, H(h[i]))
            @test convert(Vector{UInt8}, p) == X[:, i]
        end
    end
end
end


@safetestset compact_index_inv_matches_examples = "inverse compact index matches examples" begin
include(joinpath(dirname(@__FILE__), "check_hamilton.jl"))
using BijectiveHilbert: Compact, encode_hilbert_zero, axis_type, index_type
fns = ["test/3_2_5_7.txt", "test/2_5_3.txt", "test/3_4_2_3.txt",
       "test/3_5_2_3.txt", "test/4_5_5_5.txt"]
for fn in fns
    if !isfile(fn)
        fn = basename(fn)
    end
    if isfile(fn)
        h, X, ms = hamilton_example(fn)
        n = size(X, 1)
        gg = Compact(ms)
        H = index_type(gg)
        for i in eachindex(h)
            XX = convert(Vector{axis_type(gg)}, X[:, i])
            hc = encode_hilbert_zero(gg, XX)
            @test hc == H(h[i])
            if hc != H(h[i])
                @show n, ms, i, h[i], hc, X[:, i]
                break
            end
        end
    end
end
end


@safetestset compact_for_signed = "compact hilbert for signed integers" begin
using BijectiveHilbert
using Random
rng = MersenneTwister(9742439)
for i in 1:1000
    n = rand(rng, 3:5)
    ms = rand(rng, 2:7, n)
    gg = Compact(Int, ms)
    gg2 = Compact(ms)
    x = zeros(Int, n)
    A = axis_type(gg2)
    y = zeros(A, n)
    for i in 1:n
        x[i] = rand(rng, 0:(1<<ms[i] - 1))
        y[i] = A(x[i])
    end
    a = encode_hilbert_zero(gg, x)
    b = Int(encode_hilbert_zero(gg2, y))
    @test a == b
end
end