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

include("test/bitops.jl")

# binary reflected gray code
brgc(i) = i ⊻ (i >> 1)
for i in 0:10
    println(lpad(i, 3, " "), " ", lpad(string(brgc(i), base = 2), 8, "0"))
end

# Here's a test of the gray code.
n = 5
for i in 0:(1 << n - 1)
    println(brgc(1<<n - 1 - i), " ", brgc(i) ⊻ (1 << (n-1)))
    @assert(brgc(1<<n - 1 - i) == brgc(i) ⊻ (1 << (n - 1)))
end

function brgc_inv(g::Integer)
    i = g
    m = log_base2(g)
    for j = 1:m
        i = i ⊻ (g >> j)
    end
    i
end

for i in 1:100
    @assert(brgc_inv(brgc(i)) == i)
end

# hi_g(i) is the direction in which brgc changes.
# He uses a function g(i) = trailing_set_bits.
hi_g(i) = trailing_set_bits(i)
hi_g_orig(i) = floor(typeof(i), log2(brgc(i) ⊻ brgc(i + 1)))

for i = 1:100
    @assert hi_g(i) == hi_g_orig(i)
end

n = 3
for i in 0:(1<<n)
    println(
        bshow(brgc(i)), " ",
        bshow(brgc(i + 1)), " ",
        bshow(brgc(i) ⊻ brgc(i + 1)), " ",
        bshow(1<<hi_g(i)), " ",
        bshow(i)
    )
end

# Test for the hi_g function.
n = 5
for i in 0:(1<<n - 2)
    println(hi_g(i), " ", hi_g(1<<n - 2 - i))
    @assert(hi_g(i) == hi_g(1<<n - 2 - i))
end

# entry points. pg. 13 bottom.
# e = gc(2⌊(i-1)/2⌋) = gc(2*floor((i-1)/2)).
# domain is i=0 to i=2^n - 1.
function hi_e(i::Integer)
    i == zero(i) && return zero(i)
    brgc((i - one(i)) & ~one(i))
end

function hi_e_orig(i::Base.BitInteger)
    i == zero(i) && return zero(i)
    brgc(2 * floor(typeof(i), (i - 1) / 2))
end

# Show our implementation for hi_e matches the specification.
for i = 1:100
    @assert(hi_e(i) == hi_e_orig(i))
end


"""
page 12. directions
d(i) = 0 if i=0
    = g(i-1) mod n if i=0 mod 2
    = g(i) mod n if i=1 mod 2.
"""
function hi_d(i, n)
    i == zero(i) && return zero(i)
    i & 1 == 0 && return mod1(hi_g(i - one(i)), n)
    mod1(hi_g(i), n)
end

n = 3
for i = 0:(1<<n)
    println(i, " ", bshow(hi_d(i, n)))
end

# invariant for d on page 12
# Does not hold true for i=0 and i=1<<n - 1.
n = 5
for i = 1:(1<<n - 2)
    println(hi_d(i, n), " ", hi_d((1<<n)-1-i, n))
    @assert(hi_d(i, n) == hi_d((1<<n) - 1 - i, n))
end

# invariant for e, d, and g. pg. 11.
n = 5
for i = 0:(1<<n)
    @assert(hi_e(i + 1) == hi_e(i) ⊻ (1 << hi_d(i, n)) ⊻ (1 << hi_g(i)))
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
hi_T(b, d, e, n) = bitrotate(b ⊻ e, -d - one(b), n)

# From https://github.com/pdebuyl/libhilbert/blob/master/include/Hilbert/Algorithm.hpp
hi_T_inv(b, d, e, n) = bitrotate(b, d + one(b), n) ⊻ e

# Lemma 2.12, page 15.
# Is it -2 or -1?
# function hi_T_inv(b, d, e, n)
#     hi_T(b, n - d - one(b) - one(b), bitrotate(e, -(d + one(b)), n), n)
# end

# lemma 2.11, page 15
# Assert T_{e,d}(e) == 0
n = 3
for i = 0:(1<<n - 1)
    println(hi_T(hi_e(i), hi_d(i, n), hi_e(i), n))
    @assert(0 == hi_T(hi_e(i), hi_d(i, n), hi_e(i), n))
end

# lemma 2.11, page 15
# Assert T_{e,d}(f) == 2^(n-1)
n = 0x5
for i = 0x0:(0x1<<n - 0x2)
    f = hi_f(i, n)
    d = hi_d(i, n)
    e = hi_e(i)
    v = hi_T(f, d, e, n)
    println(bitstring(v), " ", v)
    @assert(0x1<<(n-1) == hi_T(hi_f(i, n), hi_d(i, n), hi_e(i), n))
end

# Top of page 16, stated invariant:
# (T_(e,d)(a) rotleft (d + 1)) ⊻ e == a
n = 5
for i = 0x0:(0x1<<n - 0x1)
    d = hi_d(i, n)
    e = hi_e(i)
    a = bitrotate(hi_T(i, d, e, n), d+0x1, n) ⊻ e
    @assert typeof(a) == UInt8
    @assert typeof(i) == UInt8
    @assert(a == i)
end

# invariant bottom of pg. 11. Lemma 2.6.
# fails for i=0 and i=2^n
n = 5
for i = 0b1:(0b1<<(n - 1))
    e = hi_e(i)
    ff = hi_f((0b1<<n) - 0b1 - i, n) ⊻ (0b1<<(n-1))
    println(e, " ", ff)
    @assert(e == ff)
    @assert typeof(e) == UInt8
    @assert typeof(ff) == UInt8
end

# Check that the inverse of T is an inverse.
n = 0x3
for i = 0x0:(0x1<<n - 0x1)
    for b = 0x0:(0x1<<n - 0x1)
        @assert(typeof(i) == UInt8)
        d = hi_d(i, n)
        @assert(typeof(d) == UInt8)
        e = hi_e(i)
        @assert(typeof(e) == UInt8)
        a = hi_T(b, d, e, n)
        b1 = hi_T_inv(a, d, e, n)
        println(join(string.(typeof.((i, b, a, b1))), " "))
        println("b ", bitstring(b), " a ", bitstring(a), " b1 ", bitstring(b1))
        @assert(b == b1)
    end
end

function vector_bits(n, p, i)
    l = zero(eltype(p))
    for j = 1:n
        l |= (p[j] & i) >> (i - j)
    end
    l
end


"""
    hilbert_index(n, m, p)

Hilbert index for an `n`-dimensional vector `p`, with each
component of extent less than 2^m.
"""
function hilbert_index(n, m, p)
    h = zero(eltype(p))  # hilbert index
    e = zero(eltype(p))  # entry point
    d = zero(eltype(p))  # direction
    for i = (m - 1):-1:0
        l = vector_bits(n, p, i)
        l = hi_T(l, d, e)
        w = brgc_inv(l)
        e = e ⊻ bitrotate(hi_e(w), d + 1)
        d += mod1(hi_d(w) + 1, n)
        h = min((h << n), w)
    end
    h
end
