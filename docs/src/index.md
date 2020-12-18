
# BijectiveHilbert

This package offers a function that turns two integer (i, j)-coordinates into a single integer z-coordinate, and vice versa.

```julia
julia> xy = zeros(Int, 8, 8)
julia> for y in 1:size(xy, 2)
           for x in 1:size(xy, 1)
               z = encode_hilbert(x, y)
               xy[x, y] = z
           end
       end
julia> @assert decode_hilbert(xy[5, 7]) == (5, 7)
julia> xy
8×8 Array{Int64,2}:
  1   2  15  16  17  20  21  22
  4   3  14  13  18  19  24  23
  5   8   9  12  31  30  25  26
  6   7  10  11  32  29  28  27
 59  58  55  54  33  36  37  38
 60  57  56  53  34  35  40  39
 61  62  51  52  47  46  41  42
 64  63  50  49  48  45  44  43
```
You use this function, called a [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve), if you have a two-dimensional array of data and want to improve memory locality when you access contiguous blocks of this data. It's used most often for geospatial work but is equally appropriate for dealing with large TIFF files. It helps memory locality because the order of the z-coordinate tends to pick values that are near each other in space. It belongs to the class of space-filling, self-avoiding, simple, and self-similar (FASS) curves, which includes Peano curves, and Morton z-curves.

This particular version of the function is easier to use than the traditional Hilbert curve because you don't have to tell the function how large the 2D array of coordinates is, or will eventually be. For a traditional Hilbert curve, the same value of z can refer to different values of (i, j). It comes from a paper by Chen, Wang, and Shi, titled, "A new algorithm for encoding and decoding the Hilbert order," in _Software Practice and Experience,_ 2007, vol. 37, pg 897-908.

Other space-filling, self-avoiding, simple, and self-similar curves:

* [HilbertSpaceFillingCurve](https://github.com/jonathanBieler/HilbertSpaceFillingCurve.jl) - An optimized version of the Hilbert curve.
* [Morton z-curve](https://github.com/JaneliaSciComp/Morton.jl) - Used for quadtrees and octrees.

We should make a group repo and include:

* The uneven-sided Hilbert curve: Hamilton, Chris H., and Andrew Rau-Chaplin. 2008. “Compact Hilbert Indices: Space-Filling Curves for Domains with Unequal Side Lengths.” Information Processing Letters 105 (5): 155–63.
* Peano curves
* Hilbert and Morton curves in arbitrary dimensions.

# Installation

```julia
pkg> add BijectiveHilbert
```
