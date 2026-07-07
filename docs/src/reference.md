```@meta
CurrentModule = BijectiveHilbert
```

# API Reference

## Algorithms

```@docs
Simple2D
SpaceGray
GlobalGray
FaceContinuous
Compact
```

## Deprecated

`CompactHamilton` is the pre-0.7 `Compact` algorithm (Hamilton & Rau-Chaplin
conventions). It is kept for one release cycle and will be removed in 0.8.0; new
code should use [`Compact`](@ref).

```@docs
CompactHamilton
```

## Functions

```@docs
encode_hilbert
decode_hilbert!
encode_hilbert_zero
decode_hilbert_zero!
hilbertsort
hilbertsort!
```
