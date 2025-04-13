# BijectiveHilbert

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://adolgert.github.io/BijectiveHilbert.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://adolgert.github.io/BijectiveHilbert.jl/dev)
[![Build Status](https://github.com/adolgert/BijectiveHilbert.jl/workflows/CI/badge.svg)](https://github.com/adolgert/BijectiveHilbert.jl/actions)
[![Coverage](https://codecov.io/github/adolgert/BijectiveHilbert.jl/graph/badge.svg?token=gfF08n6KKg)](https://codecov.io/github/adolgert/BijectiveHilbert.jl)

This [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve) library encodes a multi-dimensional grid index into a single integer, such that nearby integers are nearby grid indices. See the [`Usage`](https://computingkitchen.com/BijectiveHilbert.jl/stable/usage/).

```julia
julia> using Pkg; Pkg.add("BijectiveHilbert")
julia> using BijectiveHilbert
julia> xy_coordinates = zeros(Int, 8, 8)
julia> for y in 1:size(xy_coordinates, 2)
           for x in 1:size(xy_coordinates, 1)
               z = encode_hilbert(Simple2D(Int), [x, y])
               xy_coordinates[x, y] = z
           end
       end
julia> xy_coordinates
8×8 Array{Int64,2}:
  1   2  15  16  17  20  21  22
  4   3  14  13  18  19  24  23
  5   8   9  12  31  30  25  26
  6   7  10  11  32  29  28  27
 59  58  55  54  33  36  37  38
 60  57  56  53  34  35  40  39
 61  62  51  52  47  46  41  42
 64  63  50  49  48  45  44  43
julia> X = zeros(Int, 2)
julia> decode_hilbert!(Simple2D(Int), X, xy_coordinates[5, 7])
julia> X == [5, 7]
```
This function, called a [Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve), is used most often for geospatial work or database implementation but is equally appropriate for dealing with large TIFF files. It belongs to the class of space-filling, self-avoiding, simple, and self-similar (FASS) curves, which includes Peano curves, and Morton z-curves.

Included are several variations of the Hilbert curve. They are type-stable and thoroughly tested, including bug fixes to amend the published algorithms.

* [`Simple2D`](https://computingkitchen.com/BijectiveHilbert.jl/stable/simple2d/) encodes two-dimensional Cartesian indices. Fastest and easiest to use.
* [`SpaceGray`](https://computingkitchen.com/BijectiveHilbert.jl/stable/compact/) encodes multi-dimensional Cartesian indices that all have the same power-of-two extent. The fastest multi-dimensional version.
* [`GlobalGray`](https://computingkitchen.com/BijectiveHilbert.jl/stable/globalgray/) like SpaceGray but a different pattern.
* [`FaceContinuous`](https://computingkitchen.com/BijectiveHilbert.jl/stable/facecontinuous/) like SpaceGray and GlobalGray, but another pattern again.


## Hilbert curves for computation

A video about using Hilbert curves:

[![Hilbert curves for computation](docs/src/hilbert_thumb.jpg)](https://youtu.be/MlfS7xo2L7w)
