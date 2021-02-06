

@safetestset facecontinuous_smoke = "FaceContinuous basic test" begin
# Test using the vector form of the Hilbert index so that we don't
# introduce another layer of error. When this works, delete it.
using BijectiveHilbert: FaceContinuous, H_encode!, H_decode!, msb
# Works for UInt32 or UInt16.
A = UInt32
H = UInt32
fc = FaceContinuous(H, 9, 3)
X = zeros(A, 3)
X[1] = A(0x8f)
X[2] = A(0x8f)
X[3] = A(0x8f)
hvec = zeros(A, 3)
H_encode!(fc, X, hvec)
@show X, hvec, msb(hvec[1])
X .= zero(A)
H_decode!(fc, hvec, X)
@show X
end

@safetestset facecontinuous_api = "FaceContinuous api test" begin
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
@show X, h
X .= zero(A)
decode_hilbert_zero!(fc, X, h)
@show X
end


@safetestset facecontinuous_own_inverse = "FaceContinuous is its own inverse" begin
using BijectiveHilbert: FaceContinuous, check_own_inverse
for n in [2, 3, 5]
    for b in [2, 3]
        fc = FaceContinuous(b, n)
        @test check_own_inverse(fc, b, n)
    end
end
end


@safetestset facecontinuous_complete_set = "FaceContinuous is complete set" begin
using BijectiveHilbert: FaceContinuous, check_complete_set
for n in [2, 3, 5]
    for b in [2, 3]
        fc = FaceContinuous(b, n)
        @test check_complete_set(fc, b, n)
    end
end
end
