```@meta
CurrentModule = BijectiveHilbert
```

# Usage

## Quick Start

The simplest case - 2D data:
```julia
using BijectiveHilbert
encoder = Simple2D(Int)
h = encode_hilbert(encoder, [x, y])
```

## Choose Your Encoder

### 2D data

[`Simple2D`](@ref) requires no setup - just specify the index type:
```julia
encoder = Simple2D(Int64)
```

### N-D, same-size axes

[`SpaceGray`](@ref) is the fastest N-D encoder. All axes must be the same power-of-two size:
```julia
# 4D space, each axis 0-7 (3 bits = 2^3 = 8 values)
encoder = SpaceGray(Int64, 3, 4)
```

[`GlobalGray`](@ref) and [`FaceContinuous`](@ref) use the same interface but produce different curve patterns:
```julia
encoder = GlobalGray(Int64, 3, 4)
encoder = FaceContinuous(Int64, 3, 4)
```

### N-D, different-size axes

[`Compact`](@ref) handles axes of different sizes (each must be a power of two):
```julia
# 3D space: axis 1 is 0-7 (3 bits), axis 2 is 0-3 (2 bits), axis 3 is 0-15 (4 bits)
encoder = Compact(Int64, [3, 2, 4])
```

## Encode and Decode

Convert a point to a Hilbert index and back:
```julia
point = [2, 1, 7, 3]
hilbert_index = encode_hilbert(encoder, point)

point_out = zeros(Int, 4)
decode_hilbert!(encoder, point_out, hilbert_index)
@assert point_out == point
```

For zero-based indexing, use the `_zero` variants:
```julia
point = [1, 0, 6, 2]
hilbert_index = encode_hilbert_zero(encoder, point)
decode_hilbert_zero!(encoder, point_out, hilbert_index)
```

## Sorting

Sort a `Vector{<:AbstractVector}` according to a Hilbert curve:
```julia
points = [rand(3) for _ in 1:1000]
hilbertsort!(points)
```
Methods for matrix arguments are also supported:
```julia
points_mat = permutedims(reduce(hcat, points))
hilbertsort!(points_mat)
```


## Performance Tips

- Use `StaticArrays.MVector` for the point vector when encoding/decoding in a tight loop
- The 1-based functions are thin wrappers around the 0-based implementations
