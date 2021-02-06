

@safetestset facecontinuous_smoke = "FaceContinuous basic test" begin
# Test using the vector form of the Hilbert index so that we don't
# introduce another layer of error. When this works, delete it.
using BijectiveHilbert: FaceContinuous, H_encode!, H_decode!
# Works for UInt32 or UInt16.
fc = FaceContinuous(UInt16, 16, 3)
X = zeros(UInt16, 3)
X[1] = UInt16(8)
X[2] = UInt16(5)
X[3] = UInt16(7)
hvec = zeros(UInt16, 3)
H_encode!(fc, X, hvec)
@show X, hvec
X .= zero(UInt16)
H_decode!(fc, hvec, X)
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
