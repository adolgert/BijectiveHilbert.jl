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

function inv_brgc(g::Integer)
    i = g
    m = log_base2(g)
    for j = 1:m
        i = i ⊻ (g >> j)
    end
    i
end

for i in 1:100
    @assert(inv_brgc(brgc(i)) == i)
end

# He uses a function g(i) = trailing_set_bits.
g(i) = set_bits(i) - 1

# entry points
# e = gc(2 * floor((i-1)/2))
function e(i)
    i == zero(i) && return zero(i)
    brgc((i - one(i)) & ~one(i))
end
for i = 1:20
    println((i - 1) & ~1, " ", 2 * floor(Int, (i - 1) / 2))
end
# bitrotate(x, k_to_the_left)
# page 12. directions
function d(i, n)
    i == zero(i) && return zero(i)
    i & 1 == 0 && return mod1(g(i - 1), n)
    mod1(g(i), n)
end

# exit points
f(i) = e(i) ⊻ d(i)

# page 14. stopping
function hilbert_index(n, m, p)
    h = zero(eltype(p))
    e = zero(eltype(p))
    d = zero(eltype(p))
    for i = (m - 1):-1:0
        for j = 1:n
            l |= (p[j] & i) >> 
        end
    end
end
