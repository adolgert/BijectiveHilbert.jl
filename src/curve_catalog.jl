# Loader / accessor layer for the curve catalog. The data itself is in the
# generated file curve_catalog_data.jl (RANDOM_GRAY, PATH_INDEX) and in the
# committed blob data/curve_paths.bin. These accessors are internal (not
# exported); gi and pi are 0-based catalog ids, matching the Zig source.

"""
    gray_count(k::Integer)::Int

Number of catalog Gray codes for level width `k`. Returns 0 when no Gray codes
are catalogued for that `k`.
"""
function gray_count(k::Integer)::Int
    n = 0
    for g in RANDOM_GRAY
        g.k == k && (n += 1)
    end
    return n
end

"""
    gray_for(k::Integer, gi::Integer)::Vector{UInt16}

The catalog Gray code for level width `k` and 0-based id `gi`, as a cyclic Gray
code starting at 0. Throws `ArgumentError` when `(k, gi)` is not catalogued.
"""
function gray_for(k::Integer, gi::Integer)::Vector{UInt16}
    for g in RANDOM_GRAY
        if g.k == k && g.gi == gi
            return g.gray
        end
    end
    throw(ArgumentError("no catalog Gray code for k=$k, gi=$gi"))
end

"""
    path_count(k::Integer, gi::Integer)::Int

Number of catalogued child_entry paths for the Gray code `(k, gi)`. Returns 0
when `(k, gi)` has no catalogued paths.
"""
function path_count(k::Integer, gi::Integer)::Int
    n = 0
    for p in PATH_INDEX
        (p.k == k && p.gi == gi) && (n += 1)
    end
    return n
end

"""
    path_child_entry(k::Integer, gi::Integer, pi::Integer)::Vector{UInt16}

The child_entry values for the path `(k, gi, pi)`, where `pi` is the 0-based
path id. Reads the path blob data/curve_paths.bin on demand, decoding it as
little-endian `UInt16`. Throws `ArgumentError` when `(k, gi, pi)` is not
catalogued.
"""
function path_child_entry(k::Integer, gi::Integer, pi::Integer)::Vector{UInt16}
    ref = nothing
    for p in PATH_INDEX
        if p.k == k && p.gi == gi && p.pi == pi
            ref = p
            break
        end
    end
    ref === nothing && throw(ArgumentError("no catalog path for k=$k, gi=$gi, pi=$pi"))

    blob_path = joinpath(@__DIR__, "..", "data", "curve_paths.bin")
    bytes = read(blob_path)
    words = ltoh.(reinterpret(UInt16, bytes))
    off = Int(ref.off)
    len = Int(ref.len)
    return Vector{UInt16}(words[(off + 1):(off + len)])
end
