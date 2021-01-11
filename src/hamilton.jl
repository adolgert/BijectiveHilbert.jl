# Compact Hilbert Indices by Chris Hamilton. Technical Report CS-2006-07.
# 6059 University Ave., Halifax, Nova Scotia, B3H 1W5, Canada.
#
# In this document, the notation is:
# ∧ = & for and
# ∨ = | for or
# ⊻ = for xor
# ▷ = >> for shift right
# ◁ = << for shift left
# ⌈x⌉ = ceil(x)
# ⌊x⌋ = floor(x)

# binary reflected gray code
brgc(i) = i ⊻ (i >> 1)

function brgc_inv(g::Integer)
    i = g
    m = log_base2(g)
    for j = 1:m
        i = i ⊻ (g >> j)
    end
    i
end


# hi_g(i) is the direction in which brgc changes.
# He uses a function g(i) = trailing_set_bits.
hi_g(i) = trailing_set_bits(i)
# This is the definition of g, to compare against for testing.
hi_g_orig(i) = floor(typeof(i), log2(brgc(i) ⊻ brgc(i + 1)))

# entry points. pg. 13 bottom.
# e = gc(2⌊(i-1)/2⌋) = gc(2*floor((i-1)/2)).
# domain is i=0 to i=2^n - 1.
function hi_e(i::Integer)
    i == zero(i) && return zero(i)
    brgc((i - one(i)) & ~one(i))
end

# This is verbatim from the article, for comparison.
function hi_e_orig(i::Base.BitInteger)
    i == zero(i) && return zero(i)
    brgc(2 * floor(typeof(i), (i - 1) / 2))
end


"""
page 12. directions
d(i) = 0 if i=0
    = g(i-1) mod n if i=0 mod 2
    = g(i) mod n if i=1 mod 2.
Domain is 0 <= i <= 2^n - 1.
"""
function hi_d(i, n)
    i == zero(i) && return zero(i)
    i & 1 == 0 && return mod1(hi_g(i - one(i)), n)
    mod1(hi_g(i), n)
end


"""
exit points
Corllary 2.7 on pg. 12 says:
    hi_f(i, n) = hi_e(i) ⊻ hi_d(i, n)
"""
hi_f(i, n) = hi_e(one(i)<<n - one(i) - i) ⊻ (one(i)<<(n-1))
# hi_f(i, n) = hi_e(i) ⊻ hi_d(i, n)
#hi_f(i, n) = hi_e(i) ⊻ (one(i)<<hi_d(i, n))

# bitrotate is a left shift, so negate d+1.
# T_{(e,d)}(b), so read right-to-left.
# The paper means to bit rotate over the n bits that are in use,
# not all n bits. This is not the usual bitrotate!
hi_T(b, d, e, n) = rotateright(b ⊻ e, d + one(d), n)

# The author's code differs with his paper. It doesn't add one.
# https://github.com/pdebuyl/libhilbert/blob/master/include/Hilbert/Algorithm.hpp
hi_T_inv(b, d, e, n) = rotateleft(b, d + one(d), n) ⊻ e

# Lemma 2.12, page 15.
# Is it -2 or -1?
# function hi_T_inv(b, d, e, n)
#     hi_T(b, n - d - one(b) - one(b), bitrotate(e, -(d + one(b)), n), n)
# end


"""
    ith_bit_of_indices(n, p, i)

Given `n` indices in `p`, take the `i`-th bit of each one
and place it so that the first vector's value is at the
0-th place, the second vector's value at the 1st place, and so-on.
`i` is zero-indexed.
"""
function ith_bit_of_indices(n, p, i)
    l = zero(eltype(p))
    for j = 1:n
        l |= (p[j] & (one(eltype(p))<<i)) >> (i - j + one(eltype(p)))
    end
    l
end


"""
    set_indices_bits(p, l, m, i)

Given bits in `l`, set each index in the array `p` according
to the bit in `l`. Set the i-th bit of each index in `p`.
"""
function set_indices_bits!(p, l, m, i)
    for j in 0:(m - 1)
        v = (l & (1<<j)) >>j
        p[j + 1] = p[j + 1] | (v<<i)
    end
end


"""
Return `m` bits from the `i`-th set of `m` bits
in the integer `h`. `i` is zero-based.
"""
function bitrange(h, n, i)
    v = zero(h)
    for j in 0:(n - 1)
        v = v | (h & (1<<(i*n + j)))
    end
    v
end


"""
    hilbert_index(n, m, p)

Hilbert index for an `n`-dimensional vector `p`, with each
component of extent less than 2^m. Algorithm 1 of Hamilton and
Rau-Chaplin.
"""
function hilbert_index_paper(n, m, p)
    h = zero(eltype(p))  # hilbert index
    e = zero(eltype(p))  # entry point
    d = zero(eltype(p))  # direction
    for i = (m - 1):-1:0  # i is an index. Can be any type.
        l = ith_bit_of_indices(n, p, i)
        t = rotateright(l ⊻ e, d, n)
        w = brgc_inv(t)
        h = (h << n) | w
        e = e ⊻ rotateleft(hi_e(w), d, n)
        d = mod1(d + hi_d(w, n) + one(d), n)  # n or m for hi_d?
    end
    h
end


# Make d an argument because it determines the initial direction and maybe
# it isn't initialized correctly. Maybe d and e need different initial values.
function hilbert_index(n, m, p, d = zero(eltype(p)))
    h = zero(eltype(p))  # hilbert index
    e = zero(eltype(p))  # entry point
    for i = (m - 1):-1:0  # i is an index. Can be any type.
        l = ith_bit_of_indices(n, p, i)
        # @show l
        l = hi_T(l, d, e, n)  # n or m?
        w = brgc_inv(l)
        h = (h << n) | w
        e = e ⊻ rotateleft(hi_e(w), d + one(d), n)
        d = mod1(d + hi_d(w, n) + one(d), n)  # n or m for hi_d?
    end
    h
end


function hilbert_index_inv(n, m, h)
    e = zero(h)  # entry point
    d = zero(h)  # direction
    p = zeros(typeof(h), m)
    for i = (m - 1):-1:0  # i is an index. Can be any type.
        w = bitrange(h, n, i)
        l = brgc(w)
        l = hi_T_inv(l, d, e, n)
        set_indices_bits!(p, l, m, i)
        e = e ⊻ rotateleft(hi_e(w), d + one(d), n)
        d = mod1(d + hi_d(w, n) + one(d), n)  # n or m for hi_d?
    end
    p
end
