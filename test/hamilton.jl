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

include("bitops.jl")

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

# Show our implementation for hi_e matches the specification.
for i = 1:20
    @assert((i - 1) & ~1 == 2 * floor(Int, (i - 1) / 2))
    println((i - 1) & ~1, " ", 2 * floor(Int, (i - 1) / 2))
end


# page 12. directions
function hi_d(i, n)
    i == zero(i) && return zero(i)
    i & 1 == 0 && return mod1(hi_g(i - 1), n)
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


# exit points
# hi_f(i, n) = hi_e(i) ⊻ hi_d(i, n)
hi_f(i, n) = hi_e(i) ⊻ (1<<hi_d(i, n))
# bitrotate is a left shift, so negate d+1.
# T_{(e,d)}(b), so read right-to-left.
# The paper means to bit rotate over the n bits that are in use,
# not all n bits.
hi_T(b, d, e) = bitrotate(b ⊻ e, -(d + 1))
hi_T_inv(b, d, e) = T(b, n - d - 1, bitrotate(e, -(d + 1)))

n = 3
for i = 0:(1<<n - 1)
    println(hi_T(hi_e(i), hi_d(i, n), hi_e(i)))
    @assert(0 == hi_T(hi_e(i), hi_d(i, n), hi_e(i)))
end

n = 3
for i = 0:(1<<n - 1)
    v = hi_T(hi_f(i, n), hi_d(i, n), hi_e(i))
    println(bitstring(v), " ", v)
    #@assert(0 == hi_T(hi_e(i), hi_d(i, n), hi_e(i)))
end


# invariant bottom of pg. 11
# XXX failing.
n = 3
for i = 0:(1<<n)
    println(hi_e(i), " ", hi_f((1<<n) - 1- i, n) ⊻ (1<<(n-1)))
    #@assert(hi_e(i) == hi_f(1<<n - 1 - i, n) ⊻ (1<<(n-1)))
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
