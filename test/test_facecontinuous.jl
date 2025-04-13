using TestItemRunner


@testitem "FaceContinuous basic test" begin
# Test using the vector form of the Hilbert index so that we don't
# introduce another layer of error. When this works, delete it.
using BijectiveHilbert: FaceContinuous, H_encode!, H_decode!
# Works for UInt32 or UInt16.
A = UInt32
H = UInt32
fc1 = FaceContinuous(H, 9, 3)
fc2 = FaceContinuous(H, 8, 3)
X1 = zeros(A, 3)
X1[1] = A(0x8f)
X1[2] = A(0x8f)
X1[3] = A(0x8f)
hvec1 = zeros(A, 3)
hvec2 = zeros(A, 3)
H_encode!(fc1, X1, hvec1)
H_encode!(fc2, X1, hvec2)
@assert hvec1 == hvec2
X2 = zeros(A, 3)
H_decode!(fc1, hvec1, X1)
H_decode!(fc2, hvec2, X2)
@test X1 == X2
end


@testitem "FaceContinuous api test" begin
# Test using the vector form of the Hilbert index so that we don't
# introduce another layer of error. When this works, delete it.
using BijectiveHilbert: FaceContinuous, encode_hilbert_zero, decode_hilbert_zero!
# Works for UInt32 or UInt16.
A = UInt32
H = UInt32
fc = FaceContinuous(H, 9, 3)
X = zeros(A, 3)
X[1] = A(0x8f)
X[2] = A(0x85)
X[3] = A(0x43)
h = encode_hilbert_zero(fc, X)
X2 = zeros(A, 3)
decode_hilbert_zero!(fc, X2, h)
@test X == X2
X[1] = A(0x5)
X[2] = A(0x8)
X[3] = A(0x3)
h = encode_hilbert_zero(fc, X)
X3 = zeros(A, 3)
decode_hilbert_zero!(fc, X3, h)
@test X == X3
end


@testitem "FaceContinuous MArray" begin
    using StaticArrays
    bits = 4
    dimensions = 3
    gg = FaceContinuous(UInt32, bits, dimensions)
    X = @MArray [2, 3, 7]
    hilbert_idx = encode_hilbert(gg, X)
    fill!(X, 0)
    decode_hilbert!(gg, X, hilbert_idx)
    @test X[1] == 2
    @test X[2] == 3
    @test X[3] == 7
end


@testitem "FaceContinuous is its own inverse" setup=[HilbertTestSuite] begin
using BijectiveHilbert: FaceContinuous
for n in [2, 3, 5]
    for b in [2, 3]
        fc = FaceContinuous(b, n)
        @test HilbertTestSuite.check_own_inverse(fc, b, n)
    end
end
end


@testitem "FaceContinuous is complete set" setup=[HilbertTestSuite] begin
using BijectiveHilbert: FaceContinuous
for n in [2, 3, 5]
    for b in [2, 3]
        fc = FaceContinuous(b, n)
        @test HilbertTestSuite.check_complete_set(fc, b, n)
    end
end
end
