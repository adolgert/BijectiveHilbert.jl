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


# hi_g(i) is the direction in which brgc changes.
# It is also tells us the axis along which the exit of one
# subcube meets the entrance of the next subcube.
hi_g(i) = trailing_set_bits(i)
# This is the definition of g, to compare against for testing.
hi_g_orig(i) = floor(typeof(i), log2(brgc(i) ⊻ brgc(i + 1)))


# e is the entry vertex of the ith sub-hypercube in a Gray code
# ordering of sub-hypercubes. If e is 0b011, then this is z-y-x
# within the subcube of 0 in the z, 1 in the y, 1 in the x.
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
f is the exit vertex of the ith sub-hypercube in a Gray code
ordering of the sub-hypercubes.
Corllary 2.7 on pg. 12 says:
    hi_f(i, n) = hi_e(i) ⊻ hi_d(i, n)
"""
hi_f(i, n) = hi_e(one(i)<<n - one(i) - i) ⊻ (one(i)<<(n-1))
# hi_f(i, n) = hi_e(i) ⊻ hi_d(i, n)
#hi_f(i, n) = hi_e(i) ⊻ (one(i)<<hi_d(i, n))


"""
page 12. directions
d(i) = 0 if i=0
    = g(i-1) mod n if i=0 mod 2
    = g(i) mod n if i=1 mod 2.
Domain is 0 <= i <= 2^n - 1.
"""
function hi_d(i, n)
    i == zero(i) && return zero(i)
    ((i & 1) == 0) && return hi_g(i - one(i)) % n
    hi_g(i) % n
end


# bitrotate is a left shift, so negate d+1.
# T_{(e,d)}(b), so read right-to-left.
# The paper means to bit rotate over the n bits that are in use,
# not all n bits. This is not the usual bitrotate!
# This had a +1 in the paper.
hi_T(b, d, e, n) = rotateright(b ⊻ e, d, n)

# The author's code differs with his paper. It doesn't add one.
# https://github.com/pdebuyl/libhilbert/blob/master/include/Hilbert/Algorithm.hpp
hi_T_inv(b, d, e, n) = rotateleft(b, d, n) ⊻ e

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


# hamilton version of getting the ith bit, zero-indexed i.
function get_location(p::Vector{T}, i) where {T}
    l = zero(T)
    for j = eachindex(p)
        if (p[j] & (one(T) << i)) != 0
            l |= (one(T) << (j - 1))
        end
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


function update1(l, t, w, n, e, d)
    e = l ⊻ (one(l) << d)
    d += one(d) + first_set_bit(t)
    while d >= n
        d -= n
    end
    if (w & one(w)) == zero(w)
        if d == zero(d)
            e ⊻= one(e) << (n - 1)
        else
            e ⊻= one(e) << (d - 1)
        end
    end
    e, d
end


function update2(l, t, w, n, e, d)
    e = l
    e ⊻= (one(e) << d)
    d += one(d) + first_set_bit(t)
    while d >= n
        d -= n
    end
    e, d
end


"""
    hilbert_index(n, m, p)

Hilbert index for an `n`-dimensional vector `p`, with each
component of extent less than 2^m. Algorithm 1 of Hamilton and
Rau-Chaplin.
"""
function hilbert_index_paper(n, m, p)
    T = UInt64
    h = zero(T)  # hilbert index
    e = zero(T)  # entry point
    d = one(T)  # direction
    nmask = fbvn1s(T, n)
    for i = (m - 1):-1:0  # i is an index. Can be any type.
        l = T(ith_bit_of_indices(n, p, i))
        t = hi_T(l, d, e, n)
        w = t
        if i < m - 1
            w ⊻= (one(w) << (n - 1))
        end
        # Concatenate to the index
        h |= (w & nmask) << (i * n)
        e, d = update2(l, t, w, n, e, d)
    end
    brgc_inv(h)
end


function hilbert_index_inv_paper!(n, m, h, p)
    T = UInt64
    e = zero(T)
    d = one(T)
    l = zero(T)
    p .= zero(eltype(p))
    nmask = fbvn1s(T, n)
    for i = (m - 1):-1:0
        w = (h >> (i * n)) & nmask
        t = brgc(w)
        l = hi_T_inv(t, d, e, n)
        set_indices_bits!(p, l, n, i)
        e, d = update1(l, t, w, n, e, d)
    end
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
        d = (d + hi_d(w, n) + one(d)) % n  # n or m for hi_d?
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
        d = (d + hi_d(w, n) + one(d)) % n  # n or m for hi_d?
    end
    p
end
