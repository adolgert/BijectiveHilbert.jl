# Manual live fuzz of the Compact Hilbert curve against the HilbertCurveCompact
# Zig reference engine. This is NOT part of the CI test suite (it needs a `zig`
# compiler and builds the reference shared library on the fly); it plays the same
# role as test/against_c.jl does for CompactHamilton. The committed, ccall-free CI check
# is test/test_compact_vectors.jl, which validates against vectors emitted by
# gen/generate_reference_vectors.jl.
#
# Run manually with:
#   julia --project=. test/against_zig.jl [num_configs] [/path/to/HilbertCurveCompact]
#
# It builds the Zig shared library directly with `zig build-lib` (the repo's
# `zig build` is broken under zig 0.16-dev), then draws random domains, random
# legal (gray, path, gray_index, path_index) selections, and random points, and
# checks that Julia encode == Zig encode and that decode inverts on both sides.
# It never edits the Zig repository.

using Random
using Libdl
using BijectiveHilbert: Compact, encode_hilbert_zero, decode_hilbert_zero!,
                        index_type, axis_type, gray_count, path_count

const DEFAULT_ZIG_REPO = "/Users/adolgert/dev/HilbertCurveCompact"
const NUM_CONFIGS = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
const ZIG_REPO = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_ZIG_REPO
const POINTS_PER_CONFIG = 40
const SEED = 20260706

# ---- Zig C ABI mirror (see include/hilbertcurve.h, src/domain.zig) ----

struct DomainDesc
    struct_size::UInt32
    n::UInt32
    axis_bits::Ptr{UInt32}
    axis_sizes::Ptr{UInt64}
    axis_policy::UInt32            # HILB_AXIS_POLICY_SORTED = 0
    axis_permutation::Ptr{UInt32}
    gray_family::UInt32           # brgc_rotated=0, random=4
    path_family::UInt32           # brgc_standard=0, hub_state=3, table=5
    seed::UInt64
    max_k_for_tables::UInt32
    flags::UInt32
    gray_index::UInt16
    path_index::UInt16
end

const GRAY_ENUM = Dict(:brgc => UInt32(0), :random => UInt32(4))
const PATH_ENUM = Dict(:standard => UInt32(0), :hub => UInt32(3), :catalog => UInt32(5))

struct ZigABI
    handle::Ptr{Cvoid}
    create::Ptr{Cvoid}
    free::Ptr{Cvoid}
    index_words::Ptr{Cvoid}
    encode::Ptr{Cvoid}
    decode::Ptr{Cvoid}
end

function build_zig_lib(zig_repo)
    lib_src = joinpath(zig_repo, "src", "lib.zig")
    isfile(lib_src) || error("Zig source not found: $lib_src")
    builddir = mktempdir()
    bopts = joinpath(builddir, "build_options.zig")
    open(bopts, "w") do io
        println(io, "pub const CurveTier = enum { regular, full };")
        println(io, "pub const curve_tier: CurveTier = .full;")
    end
    libpath = joinpath(builddir, Sys.iswindows() ? "hilbertcurve.dll" :
                       (Sys.isapple() ? "libhilbertcurve.dylib" : "libhilbertcurve.so"))
    cmd = `zig build-lib -dynamic -lc -O ReleaseFast
           --dep build_options -Mroot=$lib_src -Mbuild_options=$bopts
           -femit-bin=$libpath`
    println("Building Zig reference library:\n  $cmd")
    run(cmd)
    isfile(libpath) || error("no library produced at $libpath")
    return libpath
end

function load_zig_abi(libpath)
    h = dlopen(libpath)
    return ZigABI(h,
        dlsym(h, :hilb_domain_create), dlsym(h, :hilb_domain_free),
        dlsym(h, :hilb_domain_index_words),
        dlsym(h, :hilb_encode_point), dlsym(h, :hilb_decode_index))
end

function zig_domain_create(abi, m::Vector{UInt32}, gray, path, gi, pi)
    n = UInt32(length(m))
    dom = Ref{Ptr{Cvoid}}(C_NULL)
    st = GC.@preserve m begin
        desc = Ref(DomainDesc(UInt32(sizeof(DomainDesc)), n,
            pointer(m), Ptr{UInt64}(C_NULL), UInt32(0), Ptr{UInt32}(C_NULL),
            GRAY_ENUM[gray], PATH_ENUM[path], UInt64(0), UInt32(0), UInt32(0),
            UInt16(gi), UInt16(pi)))
        ccall(abi.create, Cint, (Ptr{DomainDesc}, Ptr{Ptr{Cvoid}}), desc, dom)
    end
    return (st, dom[])
end

zig_index_words(abi, dom) = Int(ccall(abi.index_words, Csize_t, (Ptr{Cvoid},), dom))
zig_free(abi, dom) = ccall(abi.free, Cvoid, (Ptr{Cvoid},), dom)

function zig_encode(abi, dom, coords::Vector{UInt64}, words)
    out = zeros(UInt64, words)
    st = ccall(abi.encode, Cint, (Ptr{Cvoid}, Ptr{UInt64}, Ptr{UInt64}), dom, coords, out)
    return (st, out)
end

function zig_decode(abi, dom, index::Vector{UInt64}, n)
    out = zeros(UInt64, n)
    st = ccall(abi.decode, Cint, (Ptr{Cvoid}, Ptr{UInt64}, Ptr{UInt64}), dom, index, out)
    return (st, out)
end

function limbs_to_uint(limbs)
    v = UInt128(0)
    for i in eachindex(limbs)
        v |= UInt128(limbs[i]) << (64 * (i - 1))
    end
    return v
end

# ---- Fuzz driver ----

# Draw a random legal selection for a domain whose active-axis count is k_max.
function random_selection(rng, k_max)
    if k_max >= 3 && rand(rng, Bool)
        gc = gray_count(k_max)
        if gc > 0
            gi = rand(rng, 0:(gc - 1))
            if rand(rng, Bool)
                pc = path_count(k_max, gi)
                if pc > 0
                    return (:random, :catalog, gi, rand(rng, 0:(pc - 1)))
                end
            end
            return (:random, :hub, gi, 0)
        end
    end
    return rand(rng, Bool) ? (:brgc, :standard, 0, 0) : (:brgc, :hub, 0, 0)
end

function run_fuzz()
    libpath = build_zig_lib(ZIG_REPO)
    abi = load_zig_abi(libpath)
    rng = MersenneTwister(SEED)

    total_points = 0
    mismatches = 0
    decode_fails = 0
    skipped = 0
    shown = 0

    for cfg in 1:NUM_CONFIGS
        n = rand(rng, 2:6)
        m = [rand(rng, 0:6) for _ in 1:n]
        sum(m) == 0 && (m[rand(rng, 1:n)] = rand(rng, 1:6))     # avoid all-zero domain
        sum(m) <= 128 || (skipped += 1; continue)               # keep within UInt128
        k_max = count(>(0), m)
        gray, path, gi, pi = random_selection(rng, k_max)

        g = Compact(m; gray = gray, path = path, gray_index = gi, path_index = pi)
        T = index_type(g); B = axis_type(g)
        st, dom = zig_domain_create(abi, UInt32.(m), gray, path, gi, pi)
        if st != 0
            skipped += 1
            continue
        end
        words = zig_index_words(abi, dom)

        for _ in 1:POINTS_PER_CONFIG
            p = UInt64[m[i] == 0 ? UInt64(0) :
                       rand(rng, UInt64) & ((one(UInt64) << m[i]) - one(UInt64)) for i in 1:n]
            est, zl = zig_encode(abi, dom, p, words)
            est == 0 || error("Zig encode status=$est for m=$m")
            zh = limbs_to_uint(zl)

            X = B[B(p[i]) for i in 1:n]
            jh = encode_hilbert_zero(g, X)
            if UInt128(jh) != zh
                mismatches += 1
                if shown < 15
                    shown += 1
                    println("ENCODE MISMATCH  m=$m gray=$gray path=$path gi=$gi pi=$pi")
                    println("  X=$(Int.(p))  julia=$(UInt128(jh))  zig=$zh")
                end
            end

            # Decode inversion on both sides.
            Y = zeros(B, n)
            decode_hilbert_zero!(g, Y, T(jh))
            dst, zc = zig_decode(abi, dom, zl, n)
            dst == 0 || error("Zig decode status=$dst")
            if Y != X || UInt64.(Y) != zc
                decode_fails += 1
                if shown < 15
                    shown += 1
                    println("DECODE MISMATCH  m=$m gray=$gray path=$path gi=$gi pi=$pi")
                    println("  X=$(Int.(p))  julia_dec=$(Int.(Y))  zig_dec=$(Int.(zc))")
                end
            end
            total_points += 1
        end
        zig_free(abi, dom)
    end

    println("\n" * "="^60)
    println("Compact vs Zig live fuzz")
    println("  configs requested:   $NUM_CONFIGS")
    println("  configs skipped:     $skipped (create rejected / out of range)")
    println("  points compared:     $total_points")
    println("  encode mismatches:   $mismatches")
    println("  decode mismatches:   $decode_fails")
    println("="^60)
    if mismatches == 0 && decode_fails == 0
        println("PASS: Julia Compact is bit-identical to the Zig reference.")
    else
        println("FAIL: divergence detected -- isolate the smallest failing case above.")
    end
end

run_fuzz()
