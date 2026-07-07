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

println("Simple2d 2d decode")

xout = UInt8[0, 0]
sut4 = () -> decode_hilbert_zero!(fc, xout, UInt(372))
trial4 = @benchmark sut4()
println(trial4)

xout = @MArray [0, 0]
sut5 = () -> decode_hilbert_zero!(fc, xout, UInt(372))
trial5 = @benchmark sut5()
println(trial5)


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

println("SpaceGray 2d 1-based encode")
fc = SpaceGray(UInt64, 8, 2)
example1 = UInt8[1, 2]
sut1 = () -> BijectiveHilbert.encode_hilbert(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2]
sut2 = () -> BijectiveHilbert.encode_hilbert(fc, example2)
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2]
sut3 = () -> BijectiveHilbert.encode_hilbert(fc, example3)
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


println("SpaceGray 6d decode")
fc = SpaceGray(UInt64, 8, 6)
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2, 3, 4, 5, 6]
sut2 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example2, UInt64(137))
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example3, UInt64(137))
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

println("GlobalGray 6d decode")
gg = GlobalGray(UInt64, 8, 6)
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(gg, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2, 3, 4, 5, 6]
sut2 = () -> BijectiveHilbert.decode_hilbert_zero!(gg, example2, UInt64(137))
trial2 = @benchmark sut2()
println(trial2)


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


println("FaceContinuous 6d decode")

fc = FaceContinuous(UInt64, 8, 6)
example1 = UInt16[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

example2 = @MArray [1, 2, 3, 4, 5, 6]
sut2 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example2, UInt64(137))
trial2 = @benchmark sut2()
println(trial2)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example3, UInt64(137))
trial3 = @benchmark sut3()
println(trial3)


println("Compact 6d encode")
fc = Compact(UInt64, [8, 8, 8, 8, 8, 8])
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)

println("Compact 6d decode")
fc = Compact(UInt64, [8, 8, 8, 8, 8, 8])
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example3, UInt64(137))
trial3 = @benchmark sut3()
println(trial3)


println("GluedSeam 6d encode (closed-form BRGC)")
fc = GluedSeam(fill(8, 6))
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)

println("GluedSeam 6d decode (closed-form BRGC)")
fc = GluedSeam(fill(8, 6))
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example3, UInt64(137))
trial3 = @benchmark sut3()
println(trial3)


println("GluedSeam 6d encode (table-driven, random Gray + hub path)")
fc = GluedSeam(fill(8, 6); gray=:random, path=:hub)
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example3)
trial3 = @benchmark sut3()
println(trial3)

println("GluedSeam 6d decode (table-driven, random Gray + hub path)")
fc = GluedSeam(fill(8, 6); gray=:random, path=:hub)
example1 = UInt8[1, 2, 3, 4, 5, 6]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

example3 = UInt16[1, 2, 3, 4, 5, 6]
sut3 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example3, UInt64(137))
trial3 = @benchmark sut3()
println(trial3)


# Anisotropic head-to-head: Compact vs GluedSeam on axes of different sizes.
# (Uniform-m indices match exactly; anisotropic indices differ because
# GluedSeam sorts axes descending and Compact sorts ascending.)
println("Compact anisotropic [10,5,7] encode")
fc = Compact(UInt64, [10, 5, 7])
example1 = UInt16[3, 2, 5]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

println("GluedSeam anisotropic [10,5,7] encode")
fc = GluedSeam([10, 5, 7])
example1 = UInt16[3, 2, 5]
sut1 = () -> BijectiveHilbert.encode_hilbert_zero(fc, example1)
trial1 = @benchmark sut1()
println(trial1)

println("Compact anisotropic [10,5,7] decode")
fc = Compact(UInt64, [10, 5, 7])
example1 = UInt16[0, 0, 0]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)

println("GluedSeam anisotropic [10,5,7] decode")
fc = GluedSeam([10, 5, 7])
example1 = UInt16[0, 0, 0]
sut1 = () -> BijectiveHilbert.decode_hilbert_zero!(fc, example1, UInt64(137))
trial1 = @benchmark sut1()
println(trial1)
