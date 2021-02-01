```@meta
CurrentModule = BijectiveHilbert
```

# Usage

All of the Hilbert algorithms have the same interface.

## Create an algorithm

For most Hilbert curve algorithms, you have to say, beforehand, how large the multidimensional coordinates are, in powers of two. For instance, a three-dimensional grid can have values from 1 to 32 in each dimension, so `n = 3` and `b = 5` because `32 = 2^b`.

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

Note that the `Compact` algorithm can have different sizes for each dimension, but they are all powers of two. It produces a Hilbert index that uses only as many bits as necessary.


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
    for k in 1:(1<<bits)
        for j in 1:(1<<bits)
            for i in 1:(1<<bits)
                X = [i, j, k] .- 1
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
