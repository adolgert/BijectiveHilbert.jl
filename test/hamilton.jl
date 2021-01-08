# Compact Hilbert Indices by Chris Hamilton. Technical Report CS-2006-07.
# 6059 University Ave., Halifax, Nova Scotia, B3H 1W5, Canada.

# binary reflected gray code
brgc(i) = i ⊻ (i >> 1)
for i in 1:10
    println(lpad(i, 3, " "), " ", lpad(string(brgc(i), base = 2), 8, "0"))
end


function log_base2(v::Integer)::Integer
    r = zero(v)
    while (v = v >> one(v)) > 0
      r += one(v)
    end
    r
end


function set_bits(v::Integer)
    c = zero(v)
    while v != 0
        c += v & 1
        v = v>>1
    end
    c
end

for i in 1:100
    m = log_base2(i) + 1
    println(i, " ", string(i, base = 2), " ", m, " ", ceil(Int, log(2, i)))
    @assert((1<<(m - 1)) & i != 0)
    @assert(1<<m & i == 0)
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

# He uses a function g(i) = trailing_set_bits.
hi_g(i) = set_bits(i) - 1

# entry points
# e = gc(2 * floor((i-1)/2))
function hi_e(i)
    i == zero(i) && return zero(i)
    brgc((i - one(i)) & ~one(i))
end
for i = 1:20
    println((i - 1) & ~1, " ", 2 * floor(Int, (i - 1) / 2))
end
# bitrotate(x, k_to_the_left)
# page 12. directions
function hi_d(i, n)
    i == zero(i) && return zero(i)
    i & 1 == 0 && return mod1(hi_g(i - 1), n)
    mod1(g(i), n)
end

# exit points
hi_f(i) = hi_e(i) ⊻ hi_d(i)
hi_T(b, d, e) = bitrotate(b ⊻ e, -(d + 1))
hi_T_inv(b, d, e) = T(b, n - d - 1, bitrotate(e, -(d + 1)))

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
