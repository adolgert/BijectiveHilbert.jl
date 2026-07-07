# Hub-state curve-table machinery.
#
# Port of the "hub state" algorithm from HilbertCurveCompact
# (src/hub_state.zig, src/gray_brgc.zig, src/affine.zig, src/bit_utils.zig),
# following the paper "Gluing the Seam of a Hilbert Curve".
#
# A Hilbert curve is a cyclic Gray code (the vertex ordering) glued to a
# "child path" (a per-vertex entry corner e and exit direction d). Given ANY
# cyclic Gray code, `build_hub_state` derives a valid child path in O(2^k) by
# walking the mismatch variable chi: a forward walk to the all-ones "hub", a
# reverse walk from the end, and a bridge between them. The child entry is then
# child_entry[w] = gray[w] xor chi[w].
#
# CONVENTION: this family rotates by `d` (the paper's convention), NOT `d+1`
# (Hamilton's). Do not confuse these with src/compact.jl's Hamilton-convention
# `affine_apply`/`child_entry`/`child_dir`.

# ---------------------------------------------------------------------------
# k-bit rotation and single-bit helpers (faithful port of bit_utils.zig).
#
# These reduce the rotation amount mod k and mask to k bits, so they accept
# k == 1 with a rotation of 1 (which src/bitops.jl `rotateleft` rejects via its
# `k < n` assertion). Local helpers are therefore used instead of bitops.jl.
# ---------------------------------------------------------------------------

"""
    rotl_k(x, r, k)

Rotate the low `k` bits of `x` left by `r`, wrapping. `r` is reduced mod `k`.
Faithful port of `bit_utils.rotl_k`. Returns a `UInt64` masked to `k` bits.
"""
function rotl_k(x::Integer, r::Integer, k::Integer)::UInt64
    kk = Int(k)
    @assert 0 < kk <= 64
    mask = kk == 64 ? ~UInt64(0) : (UInt64(1) << kk) - UInt64(1)
    rr = mod(Int(r), kk)
    xm = UInt64(x) & mask
    rr == 0 && return xm
    return ((xm << rr) | (xm >> (kk - rr))) & mask
end


"""
    rotr_k(x, r, k)

Rotate the low `k` bits of `x` right by `r`, wrapping. `r` is reduced mod `k`.
Faithful port of `bit_utils.rotr_k`. Returns a `UInt64` masked to `k` bits.
"""
function rotr_k(x::Integer, r::Integer, k::Integer)::UInt64
    kk = Int(k)
    @assert 0 < kk <= 64
    mask = kk == 64 ? ~UInt64(0) : (UInt64(1) << kk) - UInt64(1)
    rr = mod(Int(r), kk)
    xm = UInt64(x) & mask
    rr == 0 && return xm
    return ((xm >> rr) | (xm << (kk - rr))) & mask
end


"""
    axis_of_single_bit(x)

Return the 0-based index of the single set bit of `x`. Throws `ArgumentError`
if `x` is zero or has more than one bit set. Port of
`bit_utils.axis_of_single_bit`.
"""
function axis_of_single_bit(x::Integer)::Int
    xu = UInt64(x)
    (xu != 0 && (xu & (xu - one(UInt64))) == 0) ||
        throw(ArgumentError("value $x is not a single set bit"))
    return trailing_zeros(xu)
end


# ---------------------------------------------------------------------------
# Closed-form standard (BRGC-rotated) curve.
# ---------------------------------------------------------------------------

"""
    brgc_rotated_gray(w, k)

The Gray code used by the standard curve: `rotl_k(brgc(w), 1, k)`. Port of
`gray_brgc.brgc_rotated_g`.
"""
brgc_rotated_gray(w::Integer, k::Integer)::UInt64 = rotl_k(brgc(UInt64(w)), 1, k)


"""
    brgc_rotated_rank(g, k)

Inverse of [`brgc_rotated_gray`](@ref): `gray_decode(rotr_k(g, 1, k))`. Port of
`gray_brgc.brgc_rotated_r`.
"""
brgc_rotated_rank(g::Integer, k::Integer)::UInt64 = brgc_inv(rotr_k(UInt64(g), 1, k))


"""
    glued_child_entry(w, k)

Entry vertex of child `w` on the standard BRGC path (closed form, rotate-by-`d`
convention). Port of `gray_brgc.child_entry`.
"""
function glued_child_entry(w::Integer, k::Integer)::UInt64
    wu = UInt64(w)
    wu == 0 && return UInt64(0)
    return rotl_k(brgc((wu - one(UInt64)) & ~one(UInt64)), 1, k)
end


"""
    glued_child_dir(w, k)

Direction of child `w` on the standard BRGC path, reduced mod `k` (closed
form). Port of `gray_brgc.child_dir`. `trailing_zeros(~x)` reproduces the Zig
`trailing_ones` (`@ctz(~x)`).
"""
function glued_child_dir(w::Integer, k::Integer)::UInt8
    kk = Int(k)
    kk == 1 && return UInt8(0)
    w == 0 && return UInt8(1)
    wu = UInt64(w)
    t = (wu & one(UInt64)) != 0 ? trailing_zeros(~wu) : trailing_zeros(~(wu - one(UInt64)))
    d = t + 1
    d >= kk && (d -= kk)
    return UInt8(d)
end


# ---------------------------------------------------------------------------
# Affine state-machine primitives (k-bit, 0-based UInt64 values).
#
# NOTE: src/compact.jl defines `affine_apply`/`affine_apply_inv` with
# Hamilton's `d+1` convention. These hub-state versions use the paper's
# rotate-by-`d` convention, so they carry the `glued_` prefix (like
# `glued_child_entry`) to keep the two conventions from ever cross-dispatching.
# ---------------------------------------------------------------------------

"""
    glued_affine_apply(x, e, d, k)

Apply the affine map `rotl_k(x, d, k) xor e` (paper's rotate-by-`d`
convention). Port of `affine.CyclicAffine.apply`.
"""
glued_affine_apply(x::Integer, e::Integer, d::Integer, k::Integer)::UInt64 =
    rotl_k(x, d, k) ⊻ UInt64(e)


"""
    glued_affine_apply_inv(y, e, d, k)

Inverse of [`glued_affine_apply`](@ref): `rotr_k(y xor e, d, k)`. Port of
`affine.CyclicAffine.apply_inv`.
"""
glued_affine_apply_inv(y::Integer, e::Integer, d::Integer, k::Integer)::UInt64 =
    rotr_k(UInt64(y) ⊻ UInt64(e), d, k)


# ---------------------------------------------------------------------------
# Curve tables.
# ---------------------------------------------------------------------------

"""
    CurveTables

Tables describing one Hilbert curve on `k` axes.

- `gray[w+1]`       vertex (Gray code) for digit `w` (`w` 0-based)
- `gray_rank[g+1]`  the digit `w` whose vertex is `g`
- `child_entry[w+1]` entry corner `e` for child `w`
- `child_dir[w+1]`   exit direction `d` for child `w`, reduced mod `k`
"""
struct CurveTables
    k::Int
    gray::Vector{UInt64}
    gray_rank::Vector{UInt64}
    child_entry::Vector{UInt64}
    child_dir::Vector{UInt8}
end


# ---------------------------------------------------------------------------
# Validation. Both throw ArgumentError so a corrupted catalog can never yield
# a discontinuous curve.
# ---------------------------------------------------------------------------

"""
    validate_gray(k, gray)

Check that `gray` is a length-`2^k` sequence of distinct values in `[0, 2^k)`
whose consecutive entries differ in exactly one bit. Port of
`hub_state.validate_gray`.
"""
function validate_gray(k::Integer, gray::AbstractVector{<:Integer})
    kk = Int(k)
    n = 1 << kk
    length(gray) == n ||
        throw(ArgumentError("gray length $(length(gray)) != $n (k=$kk)"))
    seen = falses(n)
    for i in 1:n
        g = Int(gray[i])
        (0 <= g < n) ||
            throw(ArgumentError("gray value $g out of range [0,$n)"))
        seen[g + 1] && throw(ArgumentError("gray value $g is not unique"))
        seen[g + 1] = true
        if i < n
            diff = UInt64(gray[i]) ⊻ UInt64(gray[i + 1])
            count_ones(diff) == 1 ||
                throw(ArgumentError("gray adjacency at index $(i - 1) is not single-bit"))
        end
    end
    return nothing
end


"""
    validate_mismatch(k, gray, chi)

Check the mismatch sequence `chi` against `gray` exactly as
`hub_state.validate_mismatch`: `chi[0] == 0` (entry), `chi[last]` a single bit
(exit), each `gray` step's changed axis is set in the following `chi` (seam /
face constraint), and consecutive `chi` entries differ in exactly one bit.
"""
function validate_mismatch(k::Integer, gray::AbstractVector{<:Integer},
                           chi::AbstractVector{<:Integer})
    kk = Int(k)
    n = 1 << kk
    length(chi) == n || throw(ArgumentError("chi length $(length(chi)) != $n"))
    length(gray) == n || throw(ArgumentError("gray length $(length(gray)) != $n"))
    UInt64(chi[1]) == 0 || throw(ArgumentError("invalid entry: chi[0] != 0"))
    count_ones(UInt64(chi[n])) == 1 ||
        throw(ArgumentError("invalid exit: chi[last] is not a single bit"))
    for w in 0:(n - 2)
        diff = UInt64(gray[w + 1]) ⊻ UInt64(gray[w + 2])
        count_ones(diff) == 1 ||
            throw(ArgumentError("invalid gray adjacency at index $w"))
        d = axis_of_single_bit(diff)
        ((UInt64(chi[w + 2]) >> d) & 1) != 0 ||
            throw(ArgumentError("invalid seam at index $w: axis $d not set in chi"))
    end
    for w in 0:(n - 2)
        diff = UInt64(chi[w + 1]) ⊻ UInt64(chi[w + 2])
        count_ones(diff) == 1 ||
            throw(ArgumentError("invalid chi adjacency at index $w"))
    end
    return nothing
end


"""
    fill_from_chi!(gray, chi, child_entry, child_dir)

Derive `child_entry[w] = gray[w] xor chi[w]` and `child_dir` (the axis in which
consecutive `chi` entries differ, or the single bit of the final `chi`) from a
validated mismatch sequence. Port of `hub_state.fill_from_chi`.
"""
function fill_from_chi!(gray::AbstractVector{UInt64}, chi::AbstractVector{UInt64},
                        child_entry::AbstractVector{UInt64},
                        child_dir::AbstractVector{UInt8})
    n = length(gray)
    for i in 0:(n - 1)
        child_entry[i + 1] = gray[i + 1] ⊻ chi[i + 1]
        if i + 1 < n
            child_dir[i + 1] = UInt8(axis_of_single_bit(chi[i + 1] ⊻ chi[i + 2]))
        else
            child_dir[i + 1] = UInt8(axis_of_single_bit(chi[i + 1]))
        end
    end
    return nothing
end


# ---------------------------------------------------------------------------
# Table builders.
# ---------------------------------------------------------------------------

"""
    build_hub_state(k, gray)

Build [`CurveTables`](@ref) from any cyclic Gray code `gray` using the
three-phase hub-state walk (forward to the all-ones hub, reverse from the end,
bridge). Validates the input Gray code and the derived mismatch sequence,
throwing `ArgumentError` on failure. Port of `hub_state.build_hub_state`.
"""
function build_hub_state(k::Integer, gray_in::AbstractVector{<:Integer})::CurveTables
    kk = Int(k)
    validate_gray(kk, gray_in)

    n = 1 << kk
    gray = Vector{UInt64}(undef, n)
    for i in 1:n
        gray[i] = UInt64(gray_in[i])
    end

    gray_rank = Vector{UInt64}(undef, n)
    for i in 0:(n - 1)
        gray_rank[gray[i + 1] + 1] = UInt64(i)
    end

    # Adjacency axis d[w] = axis of the single changed bit between gray[w] and
    # gray[w+1] (0-based, length n-1; Julia dd[w+1] holds 0-based element w).
    dd = Vector{Int}(undef, n - 1)
    for i in 0:(n - 2)
        dd[i + 1] = axis_of_single_bit(gray[i + 1] ⊻ gray[i + 2])
    end

    chi = zeros(UInt64, n)
    first_hub_idx = kk
    second_hub_idx = n - kk

    # Phase 1: forward walk to the all-ones hub at index first_hub_idx.
    chi[1] = UInt64(0)
    for w in 0:(first_hub_idx - 1)
        d_w = dd[w + 1]
        d_bit = UInt64(1) << d_w
        current = chi[w + 1]
        if (current & d_bit) == 0
            chi[w + 2] = current | d_bit
        else
            for b in 0:(kk - 1)
                if b != d_w && ((current >> b) & 1) == 0
                    chi[w + 2] = current | (UInt64(1) << b)
                    break
                end
            end
        end
    end

    # Phase 3: reverse walk from the exit.
    chi[n] = UInt64(1) << dd[n - 1]   # dd[n-1] == 0-based d[n-2]
    for w in (n - 2):-1:second_hub_idx
        nxt = chi[w + 2]
        required_bit = w > second_hub_idx ? dd[w] : -1   # dd[w] == 0-based d[w-1]
        if required_bit >= 0 && ((nxt >> required_bit) & 1) == 0
            chi[w + 1] = nxt | (UInt64(1) << required_bit)
        else
            found = false
            for b in 0:(kk - 1)
                if ((nxt >> b) & 1) == 0
                    chi[w + 1] = nxt | (UInt64(1) << b)
                    found = true
                    break
                end
            end
            if !found
                chi[w + 1] = nxt ⊻ UInt64(1)
            end
        end
    end

    # Phase 2: bridge from the hub forward to meet the reverse walk.
    for w in first_hub_idx:(second_hub_idx - 1)
        d_w = dd[w + 1]
        d_bit = UInt64(1) << d_w
        current = chi[w + 1]
        if (current & d_bit) == 0
            chi[w + 2] = current | d_bit
        else
            if count_ones(current) == kk
                for b in 0:(kk - 1)
                    if b != d_w
                        chi[w + 2] = current ⊻ (UInt64(1) << b)
                        break
                    end
                end
            else
                for b in 0:(kk - 1)
                    if b != d_w && ((current >> b) & 1) == 0
                        chi[w + 2] = current | (UInt64(1) << b)
                        break
                    end
                end
            end
        end
    end

    validate_mismatch(kk, gray, chi)

    child_entry = Vector{UInt64}(undef, n)
    child_dir = Vector{UInt8}(undef, n)
    fill_from_chi!(gray, chi, child_entry, child_dir)

    return CurveTables(kk, gray, gray_rank, child_entry, child_dir)
end


"""
    build_from_child_entry(k, gray, child_entry)

Build [`CurveTables`](@ref) from a Gray code and an explicit `child_entry`
array (an embedded path from the catalog). The mismatch `chi[w] = gray[w] xor
child_entry[w]` is validated exactly like a hub-state solution, so a corrupted
catalog entry is rejected rather than producing a discontinuous curve.
`gray_rank` and `child_dir` are derived. Port of
`hub_state.build_from_child_entry`.
"""
function build_from_child_entry(k::Integer, gray_in::AbstractVector{<:Integer},
                                child_entry_in::AbstractVector{<:Integer})::CurveTables
    kk = Int(k)
    validate_gray(kk, gray_in)

    n = 1 << kk
    length(child_entry_in) == n ||
        throw(ArgumentError("child_entry length $(length(child_entry_in)) != $n"))

    gray = Vector{UInt64}(undef, n)
    for i in 1:n
        gray[i] = UInt64(gray_in[i])
    end

    gray_rank = Vector{UInt64}(undef, n)
    for i in 0:(n - 1)
        gray_rank[gray[i + 1] + 1] = UInt64(i)
    end

    chi = Vector{UInt64}(undef, n)
    for i in 1:n
        chi[i] = UInt64(gray_in[i]) ⊻ UInt64(child_entry_in[i])
    end

    validate_mismatch(kk, gray, chi)

    child_entry = Vector{UInt64}(undef, n)
    child_dir = Vector{UInt8}(undef, n)
    fill_from_chi!(gray, chi, child_entry, child_dir)

    return CurveTables(kk, gray, gray_rank, child_entry, child_dir)
end
