```@meta
CurrentModule = BijectiveHilbert
```

# Usage

## Create an encoder

```julia
    simple = Simple2D(HilbertType)
    sg = SpaceGray(HilbertType, bits, dimensions)
    gg = GlobalGray(HilbertType, bits, dimensions)
    sg = FaceContinuous(HilbertType, bits, dimensions)
```

  * `HilbertType` is an unsigned integer DataType large enough to hold the HilbertIndex. The returned HilbertIndex will have this type. Examples are `UInt32` and `UInt128`.

  * `bits` is the maximum bit width of a DataType that can hold the Cartesian index in each dimension. If `bits` is 7, that means all indices are between 1 and `128=2^7` or between 0 and `127=2^7-1`.

  * `dimension` is the integer number of Cartesian dimensions.

[`Simple2D`](@ref) doesn't need to know dimensions ahead of time.


## Encode and Decode

```julia
    hilbert_index = encode_hilbert(encoder, cartesian::AbstractVector)
    decode_hilbert!(encoder, cartesian, hilbert_index)
```
This converts from Cartesian indices in a vector to an integer Hilbert index. Then it converts back to Cartesian. decoding overwrites the cartesian array values. The input Cartesian vector may have any AbstractVector type. It can be faster if you use a `StaticArrays.MVector` of fixed size.

```julia
    hilbert_index = encode_hilbert_zero(encoder, cartesian::AbstractVector)
    decode_hilbert_zero!(encoder, cartesian, hilbert_index)
```
These are the same as above, but they use zero-based indices. The one-based are a layer on top of these zero-based implementations.


## Example bit calculations

If the Hilbert curve will be in a 32 x 32 space, then the dimension is 2, and it needs `log2(32)=5` bits of resolution in both directions. For [`GlobalGray`](@ref), that's `b=5`, `n=2`. If the sizes aren't powers of two or are uneven, then set the bits to cover the largest side, so (12 x 12 x 12 x 800) would be `b=10`, `n=4` because ``800 < 2^{10}``.


# Index
```@index
```

```@autodocs
Modules = [BijectiveHilbert]
Private = false
```
