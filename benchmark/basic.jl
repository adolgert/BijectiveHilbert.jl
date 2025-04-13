using BenchmarkTools
using BijectiveHilbert
using StaticArrays


println("Simple2d 2d encode")
fc = Simple2D(UInt64)
example1 = UInt8[1, 2]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2]
sut2 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example2)
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)


println("SpaceGray 2d encode")
fc = SpaceGray(UInt64, 8, 2)
example1 = UInt8[1, 2]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2]
sut2 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example2)
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)


println("SpaceGray 6d encode")
fc = SpaceGray(UInt64, 8, 6)
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2, 3, 4, 5, 6]
sut2 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example2)
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)


println("GlobalGray 6d encode")
gg = GlobalGray(UInt64, 8, 6)
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero!(gg, example1)
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2, 3, 4, 5, 6]
sut2 = () -> BijectiveHilbert.encode_hilbert_zero!(gg, example2)
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt8[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(gg, example3)
trial3 = @benchmark sut3()
println(trial3)


println("FaceContinuous 6d encode")

fc = FaceContinuous(UInt64, 8, 6)
example1 = UInt16[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2, 3, 4, 5, 6]
sut2 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example2)
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)
