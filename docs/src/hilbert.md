```@meta
CurrentModule = BijectiveHilbert
```

# Pseudo-Hilbert Curves for Computational Science

## Introduction

The Pseudo-Hilbert curve sorts regular arrays of points into a linear order that keeps nearby points in the array close to each other in the linear order. Scientific computing algorithms use this capability in a few ways.

1. An easy form of clustering.
2. Speed up queries for databases.
3. Allocate work for distributed computing.
4. Visualize correlation of large one-dimensional datasets.
5. Visualize change over time of two-dimensional datasets.
6. Apply a one-dimensional global optimization algorithm to a multi-dimensional dataset.
7. An R-tree algorithm that uses Hilbert curves to pick node locations.

These widely different applications all use the Hilbert curve, not for drawing, but to convert indices from multiple dimensions to a single dimension and to convert them back to multiple dimensions, as needed. For high performance computing, it's not a drawn curve but an algorithm.

## The algorithm

There are two functions. One accepts an array of natural numbers as input and returns a single Hilbert index as output. We call that encoding. The other function decodes the single Hilbert index into the original natural numbers. There are a lot of pairs of functions that can do this.

For instance, you may recall from your earliest math class that we use such a pair of functions to prove that the number of integers in the first quadrant of a two-dimensional graph are countably infinite. Starting from (0, 0), move to (0, 1) and then (1, 0). Then move from (2, 0) to (1, 1) to (0, 2). This weaves along diagonals (in boustrophedon order), sweeping out a one-dimensional covering of two dimensions. Each point in 2D has a corresponding, countable point in 1D. Nor is this the only example, by far.

Julia stores its matrices in column-major order. C stores its matrices in row-major order. Both are orderings of two-dimensional rectangles according to a linear storage order. So, too, is block storage of data. These are functions like the Hilbert curve because we can automatically translate from (i, j) to a single index in memory, if we know the dimensions of the matrix.

The Hilbert algorithm is equivalent to the use of column-major order with one addition. Two points near each other according to the Hilbert index will tend to be closer to each other if we measure the Euclidean distance between their two-dimensional indices. There isn't a guaranteed upper bound on how far apart two points can be, but experiments comparing nearness, called spatial locality, confirm that the Hilbert curve tends to arrange points closer to each other than does column-major storage, block storage, or other self-similar space-filling curves, such as the Peano curve.

That's why this algorithm is useful. Given a bunch of (i, j) indices, the Hilbert index helps to arrange them in a single vector such that nearby indices tend to be closer to each other. Using a Hilbert algorithm can be a little more complicated than that, only because tradition asks that you focus on counts of bits. Let's look at why this is an easy, if important, requirement.


## Counting bits

Like a column-major order for a matrix, a Hilbert algorithm needs to know, before doing a conversion, what extent the data might have. Most implementations of Hilbert curves assume that every dimension of the input coordinates will have the same size and will be a power of two. If your data is 25 x 38, that means you need a 64 x 64 grid on which to work. That means the i-axis requires 6 bits and the j-axis requires 6-bits. As a result, the Hilbert index will run from 0 to 4095, needing 12-bits in which to store its result.

If we don't use the full size of the allocated axes, are we not going to get as good spatial locality from the Hilbert curve? This is only approximate to start with. The way to use this algorithm is to create axes that are large enough and not worry about the unused bits of information.

The bit count matters for choosing types of function arguments. It's usually easy for the coordinates to hold their integer data, but higher dimensions can make the resulting Hilbert index too large for most data types. Eight dimensions of ten bits is 80 bits, which exceeds a typical Int64. That's fine in Julia, which has both Int128 and BigInt types. In addition, most of the math is defined for unsigned integers, but the algorithms in this library are tested to work up to the last sign-bit for signed integers, so they are OK to use.


## Properties of different Hilbert curves

There are multiple Hilbert curves. There is even a paper called, "How many three-dimensional Hilbert curves are there?" by Haverkort (2017). The different curves have different orientataions relative to their axes. They can be symmetric or asymmetric. They fill space in slightly different ways, (almost) all of which require domains that are powers of two in each direction. Haverkort also wrote a lovely paper describing these differences, in "Sixteen space-filling curves and traversals for d-dimensional cubes and simplices."

For use in scientific computing, it may be more important to judge different Hilbert curve algorithms on capability, speed, and relaxation of constraints.

* The [`Simple2D`](@ref) algorithm doesn't need to know how large the axes may be before you use it, but it only works in 2D.
* The [`GlobalGray`](@ref) algorithm is fast for an n-dimensional algorithm.
* The [`FaceContinuous`](@ref) algorithm is slower and is included because it has a different shape and is historically important as the first non-recursive n-dimensional algorithm.

In general, algorithms that are written explicitly for 2D are faster than n-dimensional equivalents.
