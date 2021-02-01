```@meta
CurrentModule = BijectiveHilbert
```

# Usage

All of the Hilbert algorithms have the same interface.

## Decide dimensions of the spatial axes

All of the algorithms, except [`Simple2D`](@ref), need to know ahead of time the extent of the coordinate system. If the Hilbert curve will be in a 32 x 32 space, then the dimension is 2, and it needs `log2(32)=5` bits of resolution in both directions. For [`GlobalGray`](@ref), that's `b=5`, `n=2`. If the sizes aren't powers of two or are uneven, then set the bits to cover the largest side, so (12 x 12 x 12 x 800) would be `b=10`, `n=4` because ``800 < 2^{10}``. There is one algorithm that deals better with uneven sides. The [`Compact`](@ref) algorithm can pack bits together so that the resulting Hilbert index takes less storage space. This is usually used for database storage. It could take (12 x 12 x 12 x 800) as `Compact([4, 4, 4, 10])` which results in a 24-bit integer that can be stored in a UInt32.

## Create an algorithm

For most Hilbert curve algorithms, you have to say, beforehand, how large the multidimensional coordinates are, in powers of two. For instance, a three-dimensional grid can have values from 1 to 16 in each dimension, so `n = 3` and `b = 4` because `16 = 2^b`.

```julia
using BijectiveHilbert
dimensions = 3
bits = 4
simple = Simple2D(Int)
gg = GlobalGray(UInt, bits, dimensions)
sg = SpaceGray(UInt, bits, dimensions)
compact = Compact(UInt, fill(bits, dimensions))
```

The first argument is a datatype for the Hilbert index. It should be large enough to hold all of the bits from the n-dimensional axes. If you don't specify one, it will use the smallest unsigned integer that can hold them.

Note that the [`Compact`](@ref) algorithm can have different sizes for each dimension, but they are all powers of two. It produces a Hilbert index that uses only as many bits as necessary.


## Encode and decode

You can encode from n-dimensions to the Hilbert index, or you can decode from a Hilbert index to n-dimensions.

```julia
for algorithm in [simple, gg, sg, compact]
    for k in 1:(1<<bits)
        for j in 1:(1<<bits)
            for i in 1:(1<<bits)
                X = [i, j, k]
                h = encode_hilbert(algorithm, X)
                X .= 0
                decode_hilbert!(algorithm, X, h)
            end
        end
    end
end
```

## Encode and decode with a zero-based value

The underlying algorithms use a zero-based axis and a zero-based Hilbert index. These are available, too.

```julia
for algorithm in [simple, gg, sg, compact]
    for k in 0:(1<<bits - 1)
        for j in 0:(1<<bits - 1)
            for i in 0:(1<<bits - 1)
                X = [i, j, k]
                h = encode_hilbert_zero(algorithm, X)
                X .= 0
                decode_hilbert_zero!(algorithm, X, h)
            end
        end
    end
end
```


# Index
```@index
```

```@autodocs
Modules = [BijectiveHilbert]
Private = false
```
