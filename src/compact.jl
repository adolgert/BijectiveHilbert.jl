# Compact Hilbert curve for anisotropic grids (different resolution per axis),
# built with the glued-seam construction.
#
# Port of the affine transducer from HilbertCurveCompact
# (src/hilbert_affine.zig `encode`/`decode`, src/domain.zig domain construction),
# following the paper "Gluing the Seam of a Hilbert Curve" (Andrew Dolgert, 2026,
# https://doi.org/10.1184/R1/32104066.v1).
#
# Builds on the already-merged hub-state machinery in src/hub_state.jl and the
# curve catalog in src/curve_catalog.jl. Axes are stable-sorted DESCENDING by
# bit count, so the active axes at every level are a prefix and no per-level
# state embedding is needed (contrast src/compact_hamilton.jl, which sorts
# ascending and uses `embed_state`). The affine state uses the paper's
# rotate-by-`d` convention (`glued_affine_apply`), NOT Hamilton's `d+1`
# convention.

"""
    build_glued_tables_for_k(k, gray, path, gray_index, path_index)

Resolve the [`CurveTables`](@ref) for a single level width `k`, or `nothing`
when the closed-form BRGC gray + standard path applies at that width. Mirrors
the per-k table resolution loop in HilbertCurveCompact's `domainCreate`
(src/domain.zig): random Gray codes come from the catalog with
`gray_index mod gray_count(k)` at widths `k >= 3`, catalogued child paths use
`path_index mod path_count(k, gi)`, and everything else derives a hub-state path
from the resolved Gray code.
"""
function build_glued_tables_for_k(k::Int, gray::Symbol, path::Symbol,
                                  gray_index::Int, path_index::Int)::Union{Nothing,CurveTables}
    nk = 1 << k
    from_catalog = false
    gi_eff = 0
    local gray_seq::Vector{UInt64}

    if gray == :random && k >= 3
        gcount = gray_count(k)
        if gcount > 0
            gi_eff = mod(gray_index, gcount)
            gray_seq = convert(Vector{UInt64}, gray_for(k, gi_eff))
            from_catalog = true
        end
    end
    if !from_catalog
        gray_seq = UInt64[brgc_rotated_gray(w, k) for w in 0:(nk - 1)]
    end

    if path == :catalog
        # No catalogued Gray code at this width -> closed-form BRGC collapse.
        from_catalog || return nothing
        pcount = path_count(k, gi_eff)
        if pcount > 0
            pi_eff = mod(path_index, pcount)
            ce = convert(Vector{UInt64}, path_child_entry(k, gi_eff, pi_eff))
            return build_from_child_entry(k, gray_seq, ce)
        end
        # No catalogued path at this width: derive one below ("make do").
    end

    return build_hub_state(k, gray_seq)
end


"""
    Compact{T,B}(m::AbstractVector{<:Integer}; kwargs...)
    Compact(T, m::AbstractVector{<:Integer}; kwargs...)
    Compact(m::AbstractVector{<:Integer}; kwargs...)

Glued-seam Hilbert curve for anisotropic grids where each axis can have a
different number of bits. Always produces a lattice-continuous Hilbert curve
(consecutive indices decode to points that differ in one axis by one). If you
don't specify the type parameters they are chosen for you from the size of `m`.

# Type Parameters
- `T` - Hilbert index type (e.g., `UInt64`, `UInt128`)
- `B` - Coordinate type (e.g., `UInt32`)

# Arguments
- `m` - Vector of bit counts, one per axis (axis `i` spans `[0, 2^m[i])`).

# Keyword arguments
Select which curve family is used. Indices are 0-based catalog ids.
- `gray::Symbol = :brgc` - Gray-code family, `:brgc` or `:random`.
- `path::Symbol` - child-path family, `:standard`, `:hub`, or `:catalog`.
  Defaults to `:standard` when `gray == :brgc`, otherwise `:hub`.
- `gray_index::Integer = 0` - which catalog Gray code (only for `gray = :random`).
- `path_index::Integer = 0` - which catalog child path (only for `path = :catalog`).

Only these `(gray, path)` combinations are legal:
`(:brgc, :standard)`, `(:brgc, :hub)`, `(:random, :hub)`, `(:random, :catalog)`.

The selection must exist at the widest level `k_max` (the number of axes active
at the top level): for `:random`, `gray_index < gray_count(k_max)` is required
(this also rejects `k_max > 10`, where the catalog is empty), and for
`:catalog`, `path_index < path_count(k_max, gray_index)`. At smaller level
widths `k`, `gray_index mod gray_count(k)` and `path_index mod path_count(k)`
are used, and where the catalog has no entry (`k <= 2`) the curve collapses to
the closed-form BRGC curve.

# Example
```julia
g = Compact([3, 2, 4])                                  # default BRGC curve
h = encode_hilbert_zero(g, UInt8[5, 2, 11])
# Four axes are active at the top level (k_max = 4), so gray_index may be 0-9
# and path_index 0-9. A three-axis domain (k_max = 3) has only gray_index = 0.
g2 = Compact([4, 3, 2, 2]; gray=:random, path=:catalog, gray_index=2, path_index=7)
```

# Relationship to [`CompactHamilton`](@ref)
`Compact` is the recommended anisotropic, lattice-continuous algorithm, and it
adds a choice of curve family (`gray`/`path`) that the legacy
[`CompactHamilton`](@ref) engine does not offer. `Compact` sorts axes by
descending bit count and uses the paper's rotate-by-`d` affine state, whereas
`CompactHamilton` uses Hamilton & Rau-Chaplin's ascending-sort, `d+1`
convention. The two produce identical indices on uniform grids, but their
anisotropic Hilbert indices generally differ; both are valid lattice-continuous
curves.

This code ports the affine transducer from the HilbertCurveCompact library and
follows:

 - Dolgert, Andrew (2026). Gluing the Seam of a Hilbert Curve. Carnegie Mellon
   University. Preprint. https://doi.org/10.1184/R1/32104066.v1
"""
struct Compact{T, B} <: HilbertAlgorithm{T}
    n::Int                                          # number of axes
    m_bits::Vector{Int}                             # bit counts, sorted non-increasing
    axis_perm::Vector{Int}                          # sorted position -> original 1-based axis
    k_levels::Vector{Int}                           # k_levels[s] = active axes at level s
    m_sum::Int                                      # sum(m), bits in Hilbert index
    max_m::Int                                      # max(m), number of levels
    tables::Vector{Union{Nothing, CurveTables}}     # tables[k]; nothing = closed-form BRGC
    gray::Symbol
    path::Symbol
    gray_index::Int
    path_index::Int
end


function Compact{T, B}(m::AbstractVector{<:Integer};
                       gray::Symbol = :brgc,
                       path::Symbol = (gray == :brgc ? :standard : :hub),
                       gray_index::Integer = 0,
                       path_index::Integer = 0) where {T <: Unsigned, B <: Unsigned}
    n = length(m)
    n > 0 || throw(ArgumentError("m must be non-empty"))
    all(x -> x >= 0, m) || throw(ArgumentError("m values must be non-negative"))
    n <= 64 || throw(ArgumentError("length(m) must be <= 64"))

    m_sum = sum(Int, m)
    max_m = maximum(Int, m)
    m_sum <= 8 * sizeof(T) || throw(ArgumentError("sum(m) exceeds index type capacity"))
    max_m <= 8 * sizeof(B) || throw(ArgumentError("max(m) exceeds coordinate type capacity"))

    # Legal (gray, path) combinations (mirrors domain.zig's legality rules).
    legal = (gray === :brgc && path === :standard) ||
            (gray === :brgc && path === :hub) ||
            (gray === :random && path === :hub) ||
            (gray === :random && path === :catalog)
    legal || throw(ArgumentError("illegal (gray, path) combination ($gray, $path)"))
    gray_index >= 0 || throw(ArgumentError("gray_index must be >= 0"))
    path_index >= 0 || throw(ArgumentError("path_index must be >= 0"))

    # Degenerate all-zero domain: single point, zero-width index.
    if max_m == 0
        return Compact{T, B}(n, fill(0, n), collect(1:n), Int[], 0, 0,
                             Union{Nothing, CurveTables}[], gray, path,
                             Int(gray_index), Int(path_index))
    end

    # Stable-sort axes descending by bit count (ties keep original axis order),
    # so active axes at every level are a prefix (no state embedding needed).
    order = collect(1:n)
    sort!(order; by = i -> -Int(m[i]), alg = Base.Sort.MergeSort)
    m_bits = Int[Int(m[order[j]]) for j in 1:n]
    axis_perm = order

    # k_levels[s] = number of axes with m >= s (a prefix count since sorted).
    k_levels = Vector{Int}(undef, max_m)
    for s in 1:max_m
        cnt = 0
        for j in 1:n
            m_bits[j] >= s && (cnt += 1)
        end
        k_levels[s] = cnt
    end

    k_max = k_levels[1]
    needs_tables = !(gray === :brgc && path === :standard)

    tables = Vector{Union{Nothing, CurveTables}}(nothing, k_max)
    if needs_tables
        # The selection must exist at the widest level width k_max; smaller
        # widths use index-mod-count fallback (handled per k below).
        if gray === :random && k_max >= 3
            gray_index < gray_count(k_max) ||
                throw(ArgumentError("gray_index $gray_index out of range at k_max=$k_max " *
                                    "(gray_count=$(gray_count(k_max)))"))
            if path === :catalog
                path_index < path_count(k_max, gray_index) ||
                    throw(ArgumentError("path_index $path_index out of range at k_max=$k_max, " *
                                        "gray_index=$gray_index (path_count=$(path_count(k_max, gray_index)))"))
            end
        end

        done = falses(k_max)
        for s in 1:max_m
            k = k_levels[s]
            (k == 0 || done[k]) && continue
            done[k] = true
            tables[k] = build_glued_tables_for_k(k, gray, path, Int(gray_index), Int(path_index))
        end
    end

    Compact{T, B}(n, m_bits, axis_perm, k_levels, m_sum, max_m, tables,
                  gray, path, Int(gray_index), Int(path_index))
end


function Compact(m::AbstractVector{<:Integer}; kwargs...)
    T = large_enough_unsigned(sum(m))       # index type from total bits
    B = large_enough_unsigned(maximum(m))   # coord type from max bits per axis
    Compact{T, B}(m; kwargs...)
end


function Compact(T, m::AbstractVector{<:Integer}; kwargs...)
    U = unsigned(T)
    B = large_enough_unsigned(maximum(m))   # coord type from max bits per axis
    Compact{U, B}(m; kwargs...)
end


axis_type(::Compact{T, B}) where {T, B} = B


function Base.show(io::IO, g::Compact{T, B}) where {T, B}
    m = zeros(Int, g.n)
    for j in 1:g.n
        m[g.axis_perm[j]] = g.m_bits[j]
    end
    print(io, "Compact{$T,$B}(", m, "; gray=:", g.gray, ", path=:", g.path,
          ", gray_index=", g.gray_index, ", path_index=", g.path_index, ")")
end


# ============================================================================
# Encode / Decode
# ============================================================================

"""
    encode_hilbert_zero(g::Compact, X)

Encode a point `X` (0-based coordinates) to a Hilbert index (0-based). Port of
`hilbert_affine.encode` (the computed, non-table path).
"""
function encode_hilbert_zero(g::Compact{T, B}, X::AbstractVector)::T where {T, B}
    max_m = g.max_m
    max_m == 0 && return zero(T)

    axis_perm = g.axis_perm
    k_levels = g.k_levels
    tables = g.tables

    st_e = UInt64(0)
    st_d = 0
    h = zero(T)
    for s in max_m:-1:1
        k = k_levels[s]
        tab = tables[k]

        # Gather bit-plane s: active axes at level s are the prefix j = 1:k
        # (axes are sorted descending by bit count), so plane bit (j-1) is
        # bit (s-1) of X[axis_perm[j]].
        plane = UInt64(0)
        for j in 1:k
            plane |= ((UInt64(X[axis_perm[j]]) >> (s - 1)) & one(UInt64)) << (j - 1)
        end

        pre = glued_affine_apply_inv(plane, st_e, st_d, k)
        w = tab === nothing ? brgc_rotated_rank(pre, k) : tab.gray_rank[pre + 1]

        h = (h << k) | T(w)
        s == 1 && break

        if tab === nothing
            entry = glued_child_entry(w, k)
            dir = Int(glued_child_dir(w, k))
        else
            entry = tab.child_entry[w + 1]
            dir = Int(tab.child_dir[w + 1])
        end
        st_e = glued_affine_apply(entry, st_e, st_d, k)
        st_d = (st_d + dir) % k
    end
    h
end


"""
    decode_hilbert_zero!(g::Compact, X, h)

Decode a Hilbert index `h` (0-based) into point `X` (0-based coordinates). Port
of `hilbert_affine.decode` (the computed, non-table path).
"""
function decode_hilbert_zero!(g::Compact{T, B}, X::AbstractVector, h::T) where {T <: Integer, B}
    fill!(X, zero(eltype(X)))
    max_m = g.max_m
    max_m == 0 && return

    axis_perm = g.axis_perm
    k_levels = g.k_levels
    tables = g.tables

    st_e = UInt64(0)
    st_d = 0
    bit_pos = g.m_sum
    for s in max_m:-1:1
        k = k_levels[s]
        tab = tables[k]

        bit_pos -= k
        maskT = (one(T) << k) - one(T)
        w = UInt64((h >> bit_pos) & maskT)

        gcode = tab === nothing ? brgc_rotated_gray(w, k) : tab.gray[w + 1]
        plane = glued_affine_apply(gcode, st_e, st_d, k)

        # Scatter bit-plane s: active axes at level s are the prefix j = 1:k
        # (zero-bit axes are never active, so they stay 0 from fill!).
        for j in 1:k
            ax = axis_perm[j]
            bitval = (plane >> (j - 1)) & one(UInt64)
            X[ax] |= eltype(X)(bitval << (s - 1))
        end
        s == 1 && break

        if tab === nothing
            entry = glued_child_entry(w, k)
            dir = Int(glued_child_dir(w, k))
        else
            entry = tab.child_entry[w + 1]
            dir = Int(tab.child_dir[w + 1])
        end
        st_e = glued_affine_apply(entry, st_e, st_d, k)
        st_d = (st_d + dir) % k
    end
    return
end


function decode_hilbert_zero!(g::Compact{T, B}, X::AbstractVector, h::Integer) where {T, B}
    decode_hilbert_zero!(g, X, T(h))
end
