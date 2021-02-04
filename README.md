# BijectiveHilbert

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://adolgert.github.io/BijectiveHilbert.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://adolgert.github.io/BijectiveHilbert.jl/dev)
[![Build Status](https://github.com/adolgert/BijectiveHilbert.jl/workflows/CI/badge.svg)](https://github.com/adolgert/BijectiveHilbert.jl/actions)
[![Coverage](https://codecov.io/gh/adolgert/BijectiveHilbert.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/adolgert/BijectiveHilbert.jl)

Data in two or more dimensions is stored linearly, so places nearby in the data can be far in memory or on disk. This package offers functions that sort multi-dimensional data to make storage and retrieval more efficient.

```julia
julia> using Pkg; Pkg.add("BijectiveHilbert")
julia> using BijectiveHilbert
julia> xy = zeros(Int, 8, 8)
julia> for y in 1:size(xy, 2)
           for x in 1:size(xy, 1)
               z = encode_hilbert(Simple2D(Int), [x, y])
               xy[x, y] = z
           end
       end
julia> X = zeros(Int, 2)
julia> decode_hilbert!(Simple2D(Int), X, xy[5, 7])
julia> X == [5, 7]
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
This function, called a [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve), is used most often for geospatial work or database implementation but is equally appropriate for dealing with large TIFF files. It helps memory locality because the order of the z-coordinate tends to pick values that are near each other in space. It belongs to the class of space-filling, self-avoiding, simple, and self-similar (FASS) curves, which includes Peano curves, and Morton z-curves.

Included are several variations of the Hilbert curve.

* [`Simple2D`](@ref), shown above, two-dimensional. Doesn't need to know axis dimensions.
* [`GlobalGray`](@ref), an n-dimensional curve where all axis dimensions must be equal.
* [`SpaceGray`](@ref), an n-dimensional curve with a different path. All axis dimensions must be equal.
* [`Compact`](@ref), an n-dimensional curve the permits each axis to be a different size.
