# Test Compact Julia implementation against C reference implementation
#
# This test is not part of the normal test suite because it requires:
# 1. A C compiler (gcc or clang)
# 2. Compiling test/hilbert_affine.c
#
# Run manually with:
#   julia --project test/against_c.jl

using Test
using Random
using BijectiveHilbert: Compact, encode_hilbert_zero, decode_hilbert_zero!

const TEST_DIR = @__DIR__
const C_SOURCE = joinpath(TEST_DIR, "hilbert_affine.c")
const C_LIB = joinpath(TEST_DIR, "libhilbert_affine.so")

function compile_c_library()
    if !isfile(C_SOURCE)
        error("C source file not found: $C_SOURCE")
    end

    # Try gcc first, then clang
    compiler = nothing
    for cc in ["gcc", "clang"]
        if success(`which $cc`)
            compiler = cc
            break
        end
    end

    if compiler === nothing
        error("No C compiler found (tried gcc, clang)")
    end

    # Compile shared library
    cmd = `$compiler -shared -fPIC -O2 -o $C_LIB $C_SOURCE`
    println("Compiling: $cmd")
    run(cmd)

    if !isfile(C_LIB)
        error("Failed to create shared library: $C_LIB")
    end

    println("Compiled successfully: $C_LIB")
end

# C function wrappers
function c_encode_64(point::Vector{UInt32}, m::Vector{Int32})::UInt64
    n = Int32(length(m))
    ccall((:hilbert_affine_encode_64, C_LIB), UInt64,
          (Ptr{UInt32}, Ptr{Int32}, Int32),
          point, m, n)
end

function c_decode_64!(point::Vector{UInt32}, h::UInt64, m::Vector{Int32})
    n = Int32(length(m))
    ccall((:hilbert_affine_decode_64, C_LIB), Cvoid,
          (UInt64, Ptr{Int32}, Int32, Ptr{UInt32}),
          h, m, n, point)
end

function c_encode_128(point::Vector{UInt32}, m::Vector{Int32})::UInt128
    n = Int32(length(m))
    h_lo = Ref{UInt64}(0)
    h_hi = Ref{UInt64}(0)
    ccall((:hilbert_affine_encode_128, C_LIB), Cvoid,
          (Ptr{UInt32}, Ptr{Int32}, Int32, Ptr{UInt64}, Ptr{UInt64}),
          point, m, n, h_lo, h_hi)
    return UInt128(h_hi[]) << 64 | UInt128(h_lo[])
end

function c_decode_128!(point::Vector{UInt32}, h::UInt128, m::Vector{Int32})
    n = Int32(length(m))
    h_lo = UInt64(h & 0xffffffffffffffff)
    h_hi = UInt64(h >> 64)
    ccall((:hilbert_affine_decode_128, C_LIB), Cvoid,
          (UInt64, UInt64, Ptr{Int32}, Int32, Ptr{UInt32}),
          h_lo, h_hi, m, n, point)
end

function test_encode_match(m::Vector{Int}, point::Vector{UInt32})
    m32 = Int32.(m)

    # Julia encode
    c = Compact{UInt64, UInt32}(m)
    h_julia = encode_hilbert_zero(c, point)

    # C encode
    h_c = c_encode_64(point, m32)

    if h_julia != h_c
        println("MISMATCH encode: m=$m, point=$point")
        println("  Julia: $h_julia")
        println("  C:     $h_c")
        return false
    end
    return true
end

function test_decode_match(m::Vector{Int}, h::UInt64)
    m32 = Int32.(m)
    n = length(m)

    # Julia decode
    c = Compact{UInt64, UInt32}(m)
    point_julia = zeros(UInt32, n)
    decode_hilbert_zero!(c, point_julia, h)

    # C decode
    point_c = zeros(UInt32, n)
    c_decode_64!(point_c, h, m32)

    if point_julia != point_c
        println("MISMATCH decode: m=$m, h=$h")
        println("  Julia: $point_julia")
        println("  C:     $point_c")
        return false
    end
    return true
end

function test_roundtrip_match(m::Vector{Int}, point::Vector{UInt32})
    m32 = Int32.(m)
    n = length(m)

    # Julia roundtrip
    c = Compact{UInt64, UInt32}(m)
    h_julia = encode_hilbert_zero(c, point)
    rt_julia = zeros(UInt32, n)
    decode_hilbert_zero!(c, rt_julia, h_julia)

    # C roundtrip
    h_c = c_encode_64(point, m32)
    rt_c = zeros(UInt32, n)
    c_decode_64!(rt_c, h_c, m32)

    ok = true
    if h_julia != h_c
        println("MISMATCH roundtrip encode: m=$m, point=$point")
        println("  Julia h: $h_julia")
        println("  C h:     $h_c")
        ok = false
    end
    if rt_julia != rt_c
        println("MISMATCH roundtrip decode: m=$m, point=$point")
        println("  Julia rt: $rt_julia")
        println("  C rt:     $rt_c")
        ok = false
    end
    if rt_julia != point
        println("MISMATCH roundtrip (not inverse): m=$m, point=$point")
        println("  Got: $rt_julia")
        ok = false
    end
    return ok
end

function exhaustive_test(m::Vector{Int})
    n = length(m)
    total_bits = sum(m)
    if total_bits > 20
        println("Skipping exhaustive test for m=$m (too many bits: $total_bits)")
        return true
    end

    c = Compact{UInt64, UInt32}(m)
    m32 = Int32.(m)

    # Test all possible points
    point = zeros(UInt32, n)
    num_points = prod(1 << mi for mi in m)

    for _ in 1:num_points
        # Test encode
        h_julia = encode_hilbert_zero(c, point)
        h_c = c_encode_64(point, m32)
        if h_julia != h_c
            println("MISMATCH exhaustive encode: m=$m, point=$point")
            println("  Julia: $h_julia, C: $h_c")
            return false
        end

        # Test decode
        rt_julia = zeros(UInt32, n)
        rt_c = zeros(UInt32, n)
        decode_hilbert_zero!(c, rt_julia, h_julia)
        c_decode_64!(rt_c, h_c, m32)
        if rt_julia != rt_c || rt_julia != point
            println("MISMATCH exhaustive decode: m=$m, h=$h_julia")
            println("  Julia: $rt_julia, C: $rt_c, expected: $point")
            return false
        end

        # Increment point (odometer style)
        for i in 1:n
            point[i] += one(UInt32)
            if point[i] < one(UInt32) << m[i]
                break
            end
            point[i] = zero(UInt32)
        end
    end

    return true
end

function run_tests()
    println("="^60)
    println("Testing Julia Compact vs C hilbert_affine")
    println("="^60)

    # Compile C library
    compile_c_library()

    @testset "Compact vs C" begin
        @testset "Uniform dimensions" begin
            for b in 1:4, n in 2:4
                m = fill(b, n)
                @test exhaustive_test(m)
            end
        end

        @testset "Anisotropic dimensions" begin
            test_cases = [
                [2, 3],
                [3, 2],
                [1, 2, 3],
                [3, 2, 1],
                [2, 3, 4],
                [4, 3, 2],
                [1, 2, 3, 4],
                [2, 2, 3],
                [3, 1, 2],
            ]
            for m in test_cases
                @test exhaustive_test(m)
            end
        end

        @testset "Edge cases" begin
            # Single dimension
            @test exhaustive_test([4])
            @test exhaustive_test([1])

            # All same small
            @test exhaustive_test([1, 1, 1])

            # One axis with 0 bits (degenerate but valid)
            # Note: This tests that axis with 0 bits is handled correctly
            m = [2, 0, 3]
            c = Compact{UInt64, UInt32}(m)
            m32 = Int32.(m)
            point = UInt32[2, 0, 5]
            h_julia = encode_hilbert_zero(c, point)
            h_c = c_encode_64(point, m32)
            @test h_julia == h_c
        end

        @testset "Larger random samples" begin
            # For larger bit counts, test random samples instead of exhaustive
            Random.seed!(42)

            test_configs = [
                [5, 5, 5],
                [6, 4, 5],
                [8, 8],
                [10, 5, 7],
            ]

            for m in test_configs
                c = Compact{UInt64, UInt32}(m)
                m32 = Int32.(m)
                n = length(m)

                # Test 1000 random points
                for _ in 1:1000
                    point = UInt32[rand(UInt32) & ((one(UInt32) << mi) - one(UInt32)) for mi in m]
                    @test test_roundtrip_match(m, point)
                end
            end
        end

        @testset "UInt128 index (large grids)" begin
            # Test cases that need more than 64 bits
            m = [20, 20, 20, 20]  # 80 bits total
            c = Compact{UInt128, UInt32}(m)
            m32 = Int32.(m)

            Random.seed!(123)

            for _ in 1:100
                point = UInt32[rand(UInt32) & ((one(UInt32) << mi) - one(UInt32)) for mi in m]

                # Julia
                h_julia = encode_hilbert_zero(c, point)
                rt_julia = zeros(UInt32, 4)
                decode_hilbert_zero!(c, rt_julia, h_julia)

                # C
                h_c = c_encode_128(point, m32)
                rt_c = zeros(UInt32, 4)
                c_decode_128!(rt_c, h_c, m32)

                @test h_julia == h_c
                @test rt_julia == rt_c
                @test rt_julia == point
            end
        end
    end

    println()
    println("="^60)
    println("All tests completed!")
    println("="^60)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_tests()

    # Clean up compiled library
    if isfile(C_LIB)
        rm(C_LIB)
        println("Cleaned up: $C_LIB")
    end
end
