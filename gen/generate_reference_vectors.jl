# Generate committed cross-validation vectors for the GluedSeam Hilbert curve.
#
# This script builds the HilbertCurveCompact Zig shared library (the reference
# implementation), calls its C ABI via ccall, and writes test/glued_seam_vectors.txt
# with the reference Hilbert indices. It also cross-checks the Julia GluedSeam
# port against the Zig engine as it writes, so the emitted file is guaranteed to
# agree with both implementations at generation time.
#
# Usage:
#   julia --project=. gen/generate_reference_vectors.jl [/path/to/HilbertCurveCompact]
#
# Requirements:
#   - A `zig` compiler on PATH. The repo's `zig build` is broken under
#     zig 0.16-dev (the build.zig uses the pre-0.16 addCSourceFile API), so this
#     script builds the shared library DIRECTLY with `zig build-lib`, supplying
#     its own `build_options` module (curve_tier = .full) to enable the catalog
#     path table. It never edits the Zig repository.
#
# AXIS ORDER (the load-bearing detail):
#   Both sides take `m` and coordinates in USER axis order and internally
#   stable-sort axes DESCENDING by bit count (Julia: MergeSort by -m; Zig:
#   domain.zig stable_sort_indices_desc under HILB_AXIS_POLICY_SORTED). Zig's
#   encoder reads coords[axis_perm[j]] exactly as Julia reads X[axis_perm[j]],
#   so passing the identical `m` and identical (unpermuted) coordinates to both,
#   with the SORTED policy, feeds equivalent domains. No pre-permutation needed.
#
# The generated file is committed for reproducibility; regenerating it on the
# same inputs (fixed RNG seed MersenneTwister(20260706)) reproduces it exactly.

using Random
using Libdl
using BijectiveHilbert: GluedSeam, encode_hilbert_zero, decode_hilbert_zero!,
                       index_type, axis_type, gray_count, path_count

const DEFAULT_ZIG_REPO = "/Users/adolgert/dev/HilbertCurveCompact"
const ZIG_REPO = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_ZIG_REPO
const OUT_FILE = normpath(joinpath(@__DIR__, "..", "test", "glued_seam_vectors.txt"))
const SEED = 20260706

# ---------------------------------------------------------------------------
# Zig C ABI mirror (see include/hilbertcurve.h and src/domain.zig).
# ---------------------------------------------------------------------------

# extern struct hilb_domain_desc. Verified sizeof == 72 with the offsets:
#   struct_size@0 n@4 axis_bits@8 axis_sizes@16 axis_policy@24 axis_permutation@32
#   gray_family@40 path_family@44 seed@48 max_k_for_tables@56 flags@60
#   gray_index@64 path_index@66  (trailing pad -> 72)
struct DomainDesc
    struct_size::UInt32
    n::UInt32
    axis_bits::Ptr{UInt32}
    axis_sizes::Ptr{UInt64}
    axis_policy::UInt32           # HILB_AXIS_POLICY_SORTED = 0
    axis_permutation::Ptr{UInt32}
    gray_family::UInt32           # brgc_rotated=0, random=4
    path_family::UInt32           # brgc_standard=0, hub_state=3, table=5
    seed::UInt64
    max_k_for_tables::UInt32
    flags::UInt32
    gray_index::UInt16
    path_index::UInt16
end

# Julia (gray, path) symbols -> Zig enum values.
const GRAY_ENUM = Dict(:brgc => UInt32(0), :random => UInt32(4))
const PATH_ENUM = Dict(:standard => UInt32(0), :hub => UInt32(3), :catalog => UInt32(5))

# Resolved C ABI entry-point pointers (ccall on a Ptr{Cvoid} may use locals,
# unlike the (:sym, lib) form which forbids them).
struct ZigABI
    handle::Ptr{Cvoid}
    create::Ptr{Cvoid}
    free::Ptr{Cvoid}
    index_words::Ptr{Cvoid}
    encode::Ptr{Cvoid}
    decode::Ptr{Cvoid}
end

function load_zig_abi(libpath::AbstractString)
    h = dlopen(libpath)
    return ZigABI(h,
        dlsym(h, :hilb_domain_create),
        dlsym(h, :hilb_domain_free),
        dlsym(h, :hilb_domain_index_words),
        dlsym(h, :hilb_encode_point),
        dlsym(h, :hilb_decode_index))
end

function build_zig_lib(zig_repo::AbstractString)
    lib_src = joinpath(zig_repo, "src", "lib.zig")
    isfile(lib_src) || error("Zig source not found: $lib_src")
    builddir = mktempdir()
    bopts = joinpath(builddir, "build_options.zig")
    open(bopts, "w") do io
        # curve-tier = full embeds the child-path catalog needed for path=:catalog.
        println(io, "pub const CurveTier = enum { regular, full };")
        println(io, "pub const curve_tier: CurveTier = .full;")
    end
    libpath = joinpath(builddir, Sys.iswindows() ? "hilbertcurve.dll" :
                       (Sys.isapple() ? "libhilbertcurve.dylib" : "libhilbertcurve.so"))
    cmd = `zig build-lib -dynamic -lc -O ReleaseFast
           --dep build_options
           -Mroot=$lib_src
           -Mbuild_options=$bopts
           -femit-bin=$libpath`
    println("Building Zig reference library:\n  $cmd")
    run(cmd)
    isfile(libpath) || error("Zig library build produced no output at $libpath")
    println("Built: $libpath")
    return libpath
end

function zig_domain_create(abi::ZigABI, m::Vector{UInt32}, gray::Symbol, path::Symbol,
                           gi::Integer, pi::Integer)
    n = UInt32(length(m))
    dom = Ref{Ptr{Cvoid}}(C_NULL)
    st = GC.@preserve m begin
        desc = Ref(DomainDesc(UInt32(sizeof(DomainDesc)), n,
            pointer(m), Ptr{UInt64}(C_NULL), UInt32(0), Ptr{UInt32}(C_NULL),
            GRAY_ENUM[gray], PATH_ENUM[path], UInt64(0), UInt32(0), UInt32(0),
            UInt16(gi), UInt16(pi)))
        ccall(abi.create, Cint, (Ptr{DomainDesc}, Ptr{Ptr{Cvoid}}), desc, dom)
    end
    st == 0 || error("hilb_domain_create failed (status=$st) for m=$m gray=$gray path=$path gi=$gi pi=$pi")
    return dom[]
end

zig_index_words(abi::ZigABI, dom) = Int(ccall(abi.index_words, Csize_t, (Ptr{Cvoid},), dom))
zig_free(abi::ZigABI, dom) = ccall(abi.free, Cvoid, (Ptr{Cvoid},), dom)

function zig_encode(abi::ZigABI, dom, coords::Vector{UInt64}, words::Int)
    out = zeros(UInt64, words)
    st = ccall(abi.encode, Cint, (Ptr{Cvoid}, Ptr{UInt64}, Ptr{UInt64}), dom, coords, out)
    st == 0 || error("hilb_encode_point failed (status=$st)")
    return out
end

function zig_decode(abi::ZigABI, dom, index::Vector{UInt64}, n::Int)
    out = zeros(UInt64, n)
    st = ccall(abi.decode, Cint, (Ptr{Cvoid}, Ptr{UInt64}, Ptr{UInt64}), dom, index, out)
    st == 0 || error("hilb_decode_index failed (status=$st)")
    return out
end

# Assemble a little-endian u64 limb array into a single unsigned integer.
function limbs_to_uint(limbs::Vector{UInt64})
    v = UInt128(0)
    for i in eachindex(limbs)
        v |= UInt128(limbs[i]) << (64 * (i - 1))
    end
    return v
end

# ---------------------------------------------------------------------------
# Configuration matrix.
# ---------------------------------------------------------------------------

base_sels() = [(:brgc, :standard, 0, 0), (:brgc, :hub, 0, 0)]

# k_max = 3: only gray_index 0 exists; exercise every catalog path (0..4).
function k3_rich()
    s = Any[(:random, :hub, 0, 0)]
    for pi in 0:(path_count(3, 0) - 1)
        push!(s, (:random, :catalog, 0, pi))
    end
    return s
end

# k_max = 4: ten gray codes; exercise gray_index up to 9 and several paths.
function k4_rich()
    s = Any[]
    for gi in (0, 1, 9)
        push!(s, (:random, :hub, gi, 0))
    end
    for (gi, which) in ((0, :ends), (3, :mid), (9, :ends))
        pc = path_count(4, gi)
        pis = which === :ends ? unique([0, pc - 1]) : [pc ÷ 2]
        for pi in pis
            push!(s, (:random, :catalog, gi, pi))
        end
    end
    return s
end

# Each entry: (m, mode, npoints, selections). mode is :exhaustive or :random.
function config_matrix()
    return [
        # ---- Uniform, exhaustive ----
        ([2, 2],        :exhaustive, 0,  base_sels()),
        ([3, 3, 3],     :exhaustive, 0,  vcat(base_sels(), k3_rich())),
        ([2, 2, 2, 2],  :exhaustive, 0,  vcat(base_sels(), k4_rich())),
        # ---- Anisotropic, exhaustive (sum(m) <= 12) ----
        ([2, 5, 3],     :exhaustive, 0,  vcat(base_sels(),
            Any[(:random, :hub, 0, 0), (:random, :catalog, 0, 0), (:random, :catalog, 0, 2)])),
        ([3, 4, 2, 3],  :exhaustive, 0,  vcat(base_sels(),
            Any[(:random, :hub, 0, 0), (:random, :catalog, 9, path_count(4, 9) - 1)])),
        ([4, 3, 2],     :exhaustive, 0,  vcat(base_sels(), k3_rich())),
        ([1, 2, 3, 4],  :exhaustive, 0,  vcat(base_sels(),
            Any[(:random, :hub, 0, 0), (:random, :catalog, 0, 0),
                (:random, :catalog, 3, 5), (:random, :catalog, 9, 0)])),
        ([2, 0, 3],     :exhaustive, 0,  base_sels()),   # zero-bit axis (accepted by Zig)
        # ---- Wide index (sum = 80 > 64 bits), fixed-seed random ----
        ([20, 20, 20, 20], :random, 50, Any[(:brgc, :standard, 0, 0), (:random, :catalog, 0, 0)]),
    ]
end

# Enumerate points for a config. Exhaustive = odometer over the full box.
function enumerate_points(m::Vector{Int}, mode::Symbol, npoints::Int, rng)
    n = length(m)
    if mode === :exhaustive
        total = prod(1 << mi for mi in m)
        pts = Vector{Vector{UInt64}}(undef, total)
        p = zeros(UInt64, n)
        for idx in 1:total
            pts[idx] = copy(p)
            for i in 1:n
                p[i] += one(UInt64)
                p[i] < (one(UInt64) << m[i]) && break
                p[i] = zero(UInt64)
            end
        end
        return pts
    else
        pts = Vector{Vector{UInt64}}(undef, npoints)
        for j in 1:npoints
            pts[j] = UInt64[m[i] == 0 ? UInt64(0) : rand(rng, UInt64) & ((one(UInt64) << m[i]) - one(UInt64)) for i in 1:n]
        end
        return pts
    end
end

# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

function main()
    libpath = build_zig_lib(ZIG_REPO)
    abi = load_zig_abi(libpath)

    total_lines = 0
    total_mismatch = 0
    per_config = Tuple{String,Int}[]

    open(OUT_FILE, "w") do io
        println(io, "# GluedSeam cross-validation vectors.")
        println(io, "# Reference engine: HilbertCurveCompact Zig library (curve-tier=full),")
        println(io, "# generated by gen/generate_reference_vectors.jl (seed $SEED).")
        println(io, "#")
        println(io, "# Format: blocks are introduced by a header line")
        println(io, "#   # m=<m1,m2,...,mn> gray=<brgc|random> path=<standard|hub|catalog> gi=<int> pi=<int>")
        println(io, "# followed by one point per line: `x1 x2 ... xn h_hex` where the x_i are the")
        println(io, "# USER-ORDER 0-based coordinates and h_hex is the Hilbert index in lowercase hex")
        println(io, "# (no prefix). Indices wider than 64 bits are a single big hex value assembled")
        println(io, "# from the Zig little-endian u64 limb array (least-significant limb first).")
        println(io, "#")
        println(io, "# Axis order: m and coordinates are USER order on BOTH sides; each side stable-")
        println(io, "# sorts axes descending by bit count internally (Julia MergeSort by -m; Zig")
        println(io, "# HILB_AXIS_POLICY_SORTED). Zero-bit axes (e.g. m=2,0,3) are accepted by the Zig")
        println(io, "# ABI and are included.")
        println(io, "#")

        for (m, mode, npoints, sels) in config_matrix()
            n = length(m)
            m32 = UInt32.(m)
            mi = Int.(m)
            rng = MersenneTwister(SEED)          # reset per config -> deterministic points
            pts = enumerate_points(mi, mode, npoints, rng)

            for (gray, path, gi, pi) in sels
                g = GluedSeam(mi; gray = gray, path = path, gray_index = gi, path_index = pi)
                T = index_type(g)
                B = axis_type(g)
                dom = zig_domain_create(abi, m32, gray, path, gi, pi)
                words = zig_index_words(abi, dom)

                header = "# m=" * join(m, ",") * " gray=$gray path=$path gi=$gi pi=$pi"
                println(io, header)
                block_mismatch = 0
                for p in pts
                    zlimbs = zig_encode(abi, dom, p, words)
                    zh = limbs_to_uint(zlimbs)

                    # Cross-check Julia against Zig (encode equality).
                    X = B[B(p[i]) for i in 1:n]
                    jh = encode_hilbert_zero(g, X)
                    if UInt128(jh) != zh
                        block_mismatch += 1
                        total_mismatch += 1
                        if total_mismatch <= 10
                            println(stderr, "MISMATCH m=$m gray=$gray path=$path gi=$gi pi=$pi X=$(Int.(p))")
                            println(stderr, "  julia=$(UInt128(jh))  zig=$zh")
                        end
                    end
                    # Also confirm Julia decode inverts (the CI test asserts this too).
                    Y = zeros(B, n)
                    decode_hilbert_zero!(g, Y, T(jh))
                    if Y != X
                        println(stderr, "DECODE-INV FAIL m=$m sel=($gray,$path,$gi,$pi) X=$(Int.(p)) Y=$(Int.(Y))")
                    end

                    coordstr = join(Int.(p), " ")
                    println(io, coordstr, " ", string(zh; base = 16))
                    total_lines += 1
                end
                zig_free(abi, dom)
                push!(per_config, (chop(header; head = 2, tail = 0), length(pts)))
                if block_mismatch > 0
                    println(stderr, "  block had $block_mismatch mismatches: $header")
                end
            end
        end
    end

    println("\nWrote $OUT_FILE")
    println("Total vector lines: $total_lines")
    println("Blocks: $(length(per_config))")
    for (h, cnt) in per_config
        println("  $cnt  $h")
    end
    if total_mismatch == 0
        println("\nAll Julia GluedSeam indices match the Zig reference (0 mismatches).")
    else
        error("$total_mismatch Julia/Zig mismatches detected -- file is NOT consistent; investigate.")
    end
end

main()
