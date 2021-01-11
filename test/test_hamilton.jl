@safetestset gray_code_symmetry = "Gray code symmetry" begin
using BijectiveHilbert
# for i in 0:10
#     println(lpad(i, 3, " "), " ", lpad(string(BijectiveHilbert.brgc(i), base = 2), 8, "0"))
# end
# Here's a test of the gray code.
n = 5
for i in 0:(1 << n - 1)
    # println(BijectiveHilbert.brgc(1<<n - 1 - i), " ", BijectiveHilbert.brgc(i) ⊻ (1 << (n-1)))
    @test(BijectiveHilbert.brgc(1<<n - 1 - i) == BijectiveHilbert.brgc(i) ⊻ (1 << (n - 1)))
end

end


@safetestset brgc_own_inverse = "Gray code is own inverse" begin
using BijectiveHilbert
for i in 1:100
    @test(BijectiveHilbert.brgc_inv(BijectiveHilbert.brgc(i)) == i)
end
end


@safetestset hi_g_matches_long_form = "hi_g is the same as its definition" begin
using BijectiveHilbert: hi_g, hi_g_orig
for i = 1:100
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


@safetestset hi_e_matches_float_version = "hi_e matches its floating point version" begin
using BijectiveHilbert: hi_e, hi_e_orig
# Show our implementation for hi_e matches the specification.
for i = 1:100
    @test(hi_e(i) == hi_e_orig(i))
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
    @test d == g
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
            @test d == g
        else
            g = hi_g(i)
            @test d == g
        end
    end
end
end


@safetestset hi_edg_invariant = "hi e, d, g invariant is consistent" begin
using BijectiveHilbert: hi_d, hi_e, hi_g
# invariant for e, d, and g. pg. 11. Eq. 1.
n = 5
for i = 0:(1<<n - 1)
    @test(hi_e(i + 1) == hi_e(i) ⊻ (1 << hi_d(i, n)) ⊻ (1 << hi_g(i)))
end
end


@safetestset hi_e_reflection = "hi_e reflects to f" begin
using BijectiveHilbert: hi_e, hi_d, hi_f
n = 5
# XXX doesn't hold for 0
for i = 0:(1<<n - 1)
    e = hi_e(i)
    d = hi_d(i, n)
    f = hi_f(i, n)
    # Corollary 2.7, pg 12, states that e ⊻ d = f, but that's not true.
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


@safetestset rotation_invariant_2n = "rotation invariant on hi_f is power of two" begin
using BijectiveHilbert: hi_T, hi_e, hi_d, hi_f
# lemma 2.11, page 15
# Assert T_{e,d}(f) == 2^(n-1)
n = 0x5
for i = 0x0:(0x1<<n - 0x2)
    f = hi_f(i, n)
    d = hi_d(i, n)
    e = hi_e(i)
    v = hi_T(f, d, e, n)
    # println(bitstring(v), " ", v)
    @test(0x1<<(n-1) == hi_T(hi_f(i, n), hi_d(i, n), hi_e(i), n))
end
end


@safetestset bitrotate_invariant = "identity invariant for rotation" begin
using BijectiveHilbert: hi_d, hi_e, bitrotaten, hi_T
# Top of page 16, stated invariant:
# (T_(e,d)(a) rotleft (d + 1)) ⊻ e == a
n = 5
for i = 0x0:(0x1<<n - 0x1)
    d = hi_d(i, n)
    e = hi_e(i)
    a = bitrotaten(hi_T(i, d, e, n), d+0x1, n) ⊻ e
    @test typeof(a) == UInt8
    @test typeof(i) == UInt8
    @test(a == i)
end
end


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
using BijectiveHilbert: ith_bit_of_indices

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
end
end


@safetestset hilbert_index_complete2d = "hilbert index is a complete set for 2d" begin
using BijectiveHilbert: hilbert_index, hilbert_index_paper

m = 0x4  # One larger won't fit into a byte and must fail.
seen = Dict{UInt8,Tuple{UInt8,UInt8}}()
for i in 0:(1<<m - 1)
    for j in 0:(1<<m - 1)
        h = hilbert_index(0x2, m, UInt8[i, j])
        # println(bitstring(h))
        seen[h] = (i, j)
        # assert that hilbert indices are within the range.
        @test h >= 0
        @test h < 1<<(2m)
    end
end
# Assert that every unique value is seen.
@test(length(seen) == 1<<(2m))
x0, y0 = seen[0x0]
for hidx in 0x1:(0x1<<2m - 0x1)
    x1, y1 = seen[hidx]
    dx = (x1 > x0) ? x1 - x0 : x0 - x1  # because unsigned
    dy = (y1 > y0) ? y1 - y0 : y0 - y1
    # Assert that each x,y is one step away from previous value.
    @test(dx == 0x1 || dy == 0x1)
    x0, y0 = (x1, y1)
end
end


@safetestset hilbert_index_complete4d = "hilbert index is a complete set for 4d" begin
using BijectiveHilbert: hilbert_index, hilbert_index_paper
# Try a 4d hilbert curve.
dim_cnt = 4
m = 3
indices = Base.IteratorsMD.CartesianIndices(tuple(collect(1<<m for i in 1:dim_cnt)...))
seen2 = Dict{Int64, NTuple{dim_cnt,Int64}}()
for idx in indices
    h = hilbert_index(dim_cnt, m, Tuple(idx), 4)
    seen2[h] = Tuple(idx)
    @test h >= 0
    @test h < 1<<(dim_cnt * m)
end
@test(length(seen2) == 1<<(dim_cnt * m))
for hidx in 0:(1<<(dim_cnt*m) - 2)  # compare with next, so stop one early.
    differ = seen2[hidx] .!= seen2[hidx + 1]
    @test(sum(differ) == 1)
    choose = [x for x in 1:dim_cnt if differ[x]]
    a = collect(seen2[hidx])[choose]
    b = collect(seen2[hidx + 1])[choose]
    dx = (a > b) ? a - b : b - a
    if dx != 1
        @test dx == 1
        break
    end
end
end


@safetestset hilbert_index_paper_tr = "paper and techreport agree" begin
using BijectiveHilbert: hilbert_index, hilbert_index_paper

m = 0x4  # One larger won't fit into a byte and must fail.
for i in 0:(1<<m - 1)
    for j in 0:(1<<m - 1)
        h_tr = hilbert_index(0x2, m, UInt8[i, j])
        h_p = hilbert_index_paper(0x2, m, UInt8[i, j])
        @test typeof(h_tr) == typeof(h_p)
        @test h_p == h_tr
        if h_p != h_tr
            break
        end
    end
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
