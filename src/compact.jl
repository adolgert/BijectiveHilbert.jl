# Compact Hilbert curve for anisotropic grids (different resolution per axis)
# Based on C implementation in test/hilbert_affine.c

"""
    sort_axes_by_priority(m::Vector{Int})

Sort axes by (m[j], j) ascending. Returns a permutation where axes with
fewer bits come first, breaking ties by axis index.
"""
function sort_axes_by_priority(m::AbstractVector{<:Integer})
    n = length(m)
    perm = sortperm(collect(zip(m, 1:n)))
    return perm
end


"""
    Compact{T,B}(m::AbstractVector{<:Integer})
    Compact(T, m::AbstractVector{<:Integer})
    Compact(m::AbstractVector{<:Integer})

Compact Hilbert curve algorithm for anisotropic grids where each axis
can have a different number of bits.  If you don't specify type parameters
they are chosen for you depending on the size of `m`.

# Type Parameters
- `T` - Hilbert index type (e.g., UInt64, UInt128)
- `B` - Coordinate type (e.g., UInt32)

# Arguments
- `m` - Vector of bit counts, one per axis

# Example
```julia
c = Compact{UInt64, UInt32}([3, 2, 4])  # Dimension [2^3, 2^2, 2^4]
h = encode_hilbert_zero(c, UInt32[5, 2, 11])
```

This code starts with the tech report and paper by Chris Hamilton. That paper
and subsequent versions
make a Hilbert curve for unequal side lengths without a guarantee that consecutive
points are adjacent. This code will always produce a lattice-continuous Hilbert
curve.

Hamilton, Chris. "Compact hilbert indices." Dalhousie University, Faculty of Computer Science, Technical Report CS-2006-07 (2006).
"""
struct Compact{T, B} <: HilbertAlgorithm{T}
    m::Vector{Int}           # exponents, one per axis
    mmax::Int                # max(m), number of levels
    total_bits::Int          # sum(m), bits in Hilbert index
    k_level::Vector{Int}     # k_level[s] = count of active axes at level s
    axes_level::Matrix{Int}  # axes_level[j, s] = axis at position j, level s
    pos_level::Matrix{Int}   # pos_level[ax, s] = position of axis ax at level s (-1 if inactive)
end


function Compact{T,B}(m::AbstractVector{<:Integer}) where {T<:Unsigned, B<:Unsigned}
    # Validate inputs
    n = length(m)
    n > 0 || throw(ArgumentError("m must be non-empty"))
    all(x -> x >= 0, m) || throw(ArgumentError("m values must be non-negative"))
    mmax = maximum(m)
    total_bits = sum(m)
    total_bits <= 8 * sizeof(T) || throw(ArgumentError("total_bits exceeds index type capacity"))
    mmax <= 8 * sizeof(B) || throw(ArgumentError("max(m) exceeds coordinate type capacity"))

    # Handle edge case of all zeros
    if mmax == 0
        return Compact{T,B}(copy(m), 0, 0, Int[], zeros(Int, n, 0), zeros(Int, n, 0))
    end

    # Build active axes structure
    order = sort_axes_by_priority(m)

    k_level = zeros(Int, mmax)
    axes_level = zeros(Int, n, mmax)
    pos_level = fill(-1, n, mmax)

    for s in 1:mmax
        k = 0
        for i in 1:n
            ax = order[i]
            if m[ax] >= s
                k += 1
                axes_level[k, s] = ax
            end
        end
        k_level[s] = k
        for j in 1:k
            pos_level[axes_level[j, s], s] = j
        end
    end

    Compact{T,B}(copy(m), mmax, total_bits, k_level, axes_level, pos_level)
end


function Compact(m::AbstractVector{<:Integer})
    T = large_enough_unsigned(sum(m))       # index type from total bits
    B = large_enough_unsigned(maximum(m))   # coord type from max bits per axis
    Compact{T,B}(m)
end


function Compact(T, m::AbstractVector{<:Integer})
    B = large_enough_unsigned(maximum(m))   # coord type from max bits per axis
    Compact{T,B}(m)
end


axis_type(::Compact{T,B}) where {T,B} = B


# ============================================================================
# Helper functions
# ============================================================================

"""
    child_entry(w, k)

Hamilton entry vertex: the entry point of sub-hypercube w.
For w == 0, returns 0. Otherwise returns brgc((w - 1) & ~1) masked to k bits.
"""
function child_entry(w::B, k::Int)::B where {B<:Unsigned}
    w == zero(B) && return zero(B)
    brgc((w - one(B)) & ~one(B)) & fbvn1s(B, k)
end


"""
    child_dir(w, k)

Hamilton direction: the axis along which we exit sub-hypercube w.
Uses trailing_ones to compute efficiently.
"""
function child_dir(w::Integer, k::Int)::Int
    w == 0 && return 0
    isodd(w) ? trailing_ones(w) % k : trailing_ones(w - 1) % k
end


"""
    affine_apply(x, e, d, k)

Apply affine transformation S_{e,δ}(x) = rotl(x, d+1, k) ⊻ e
where δ = d + 1. Used in decode path.
"""
function affine_apply(x::B, e::B, d::Int, k::Int)::B where {B<:Unsigned}
    (rotateleft(x, (d + 1) % k, k) ⊻ e) & fbvn1s(B, k)
end


"""
    affine_apply_inv(y, e, d, k)

Inverse affine transformation S^{-1}(y) = rotr(y ⊻ e, d+1, k).
Used in encode path.
"""
function affine_apply_inv(y::B, e::B, d::Int, k::Int)::B where {B<:Unsigned}
    rotateright(y ⊻ e, (d + 1) % k, k) & fbvn1s(B, k)
end


"""
    embed_state(A_old, k_old, pos_new, e_old, d_old)

Map state (e, d) when transitioning from k_old active axes to more axes.
Scatters bits of e from old positions to new positions and maps the
direction axis through position lookup.

# Arguments
- `A_old` - view of axes at old level: axes_level[1:k_old, s]
- `k_old` - number of active axes at old level
- `pos_new` - view of positions at new level: pos_level[:, s-1]
- `e_old` - entry point at old level (k_old bits)
- `d_old` - direction at old level (0-based)

# Returns
- `(e_new, d_new)` - state mapped to new level
"""
function embed_state(A_old::AbstractVector{Int}, k_old::Int, pos_new::AbstractVector{Int},
                     e_old::B, d_old::Int)::Tuple{B,Int} where {B<:Unsigned}
    e_new = zero(B)
    for j in 1:k_old
        if (e_old >> (j - 1)) & one(B) != zero(B)
            new_pos = pos_new[A_old[j]]
            e_new |= one(B) << (new_pos - 1)
        end
    end
    # d is 0-based, so A_old[d_old + 1] gets the axis, then look up its new position
    dir_axis = A_old[d_old + 1]
    d_new = pos_new[dir_axis] - 1  # convert back to 0-based
    (e_new, d_new)
end


# ============================================================================
# Encode / Decode
# ============================================================================

"""
    encode_hilbert_zero(c::Compact, X)

Encode a point X (0-based coordinates) to a Hilbert index (0-based).
"""
function encode_hilbert_zero(c::Compact{T,B}, X::AbstractVector)::T where {T,B}
    (; mmax, k_level, axes_level, pos_level) = c
    mmax == 0 && return zero(T)

    e = zero(B)
    d = 0
    h = zero(T)

    for s in mmax:-1:1
        k = k_level[s]
        A = @view axes_level[1:k, s]
        mask = fbvn1s(B, k)

        e &= mask
        d = d % k

        # Gather bits from coordinates at level s
        plane = zero(B)
        for j in 1:k
            ax = A[j]
            plane |= ((B(X[ax]) >> (s - 1)) & one(B)) << (j - 1)
        end
        plane &= mask

        # Inverse affine transform, then Gray decode
        pre = affine_apply_inv(plane, e, d, k)
        w = brgc_inv(pre) & mask

        # Pack digit MSB-first
        h = (h << k) | T(w)

        # Update state
        entry = child_entry(w, k) & mask
        e = (e ⊻ rotateleft(entry, (d + 1) % k, k)) & mask
        d = (d + child_dir(w, k) + 1) % k

        # Embed state if k increases at next level
        if s > 1 && k_level[s - 1] > k
            pos_new = @view pos_level[:, s - 1]
            e, d = embed_state(A, k, pos_new, e, d)
        end
    end

    h
end


"""
    decode_hilbert_zero!(c::Compact, X, h)

Decode a Hilbert index h (0-based) to point X (0-based coordinates).
"""
function decode_hilbert_zero!(c::Compact{T,B}, X::AbstractVector, h::T) where {T,B}
    (; mmax, total_bits, k_level, axes_level, pos_level) = c

    fill!(X, zero(eltype(X)))
    mmax == 0 && return

    bit_pos = total_bits
    e = zero(B)
    d = 0

    for s in mmax:-1:1
        k = k_level[s]
        A = @view axes_level[1:k, s]
        mask = fbvn1s(B, k)

        e &= mask
        d = d % k

        # Extract k-bit digit
        bit_pos -= k
        w = B((h >> bit_pos) & T(mask))

        # Gray encode, then affine transform
        g = brgc(w) & mask
        plane = affine_apply(g, e, d, k)

        # Scatter bits to coordinates at level s
        for j in 1:k
            ax = A[j]
            X[ax] |= eltype(X)(((plane >> (j - 1)) & one(B)) << (s - 1))
        end

        # Update state (same as encode)
        entry = child_entry(w, k) & mask
        e = (e ⊻ rotateleft(entry, (d + 1) % k, k)) & mask
        d = (d + child_dir(w, k) + 1) % k

        # Embed state if k increases at next level
        if s > 1 && k_level[s - 1] > k
            pos_new = @view pos_level[:, s - 1]
            e, d = embed_state(A, k, pos_new, e, d)
        end
    end
end
