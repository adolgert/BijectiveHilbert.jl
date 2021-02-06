

@safetestset facecontinuous_smoke = "FaceContinuous basic test" begin
# Test using the vector form of the Hilbert index so that we don't
# introduce another layer of error.
using BijectiveHilbert: FaceContinuous, H_encode!, H_decode!
fc = FaceContinuous(UInt32, 32, 2)
X = zeros(UInt32, 2)
X[1] = UInt32(2)
X[2] = UInt32(3)
hvec = zeros(UInt32, 2)
H_encode!(fc, X, hvec)
@show X, hvec
X .= zero(UInt32)
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
