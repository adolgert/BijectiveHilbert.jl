# Implementations

See the [`Usage`](https://computingkitchen.com/BijectiveHilbert.jl/stable/usage/).

## The Simple2D Algorithm

If you want to sort some axes, then this Hilbert curve algorithm is the easiest to use. It doesn't need to know ahead of time how many bits it will need to generate a Hilbert index [^1]. As the length along each spatial axis grows, it creates gradually larger Hilbert indices to match it.

All of the algorithms use slightly different Hilbert curves. This one uses an asymmetric curve that shifts so that its endpoint is always an outside corder of each ``2^n`` x ``2^n`` tile. The next outer layer builds on the last.

## The SpaceGray Algorithm

The [`SpaceGray`](@ref) algorithm is the fastest way to generate multi-dimensional Hilbert curves. This library implements a classic space key algorithm that relies on Gray codes. It's not recursive. It follows four steps. Quoting the source paper [^2]:

1. Find the cell containing the point of interest.

2. Update the key (index) value appropriately.

3. Transform as necessary; and

4. Continue until sufficient precision has been obtained.

## Global Gray

See the [`Usage`](https://computingkitchen.com/BijectiveHilbert.jl/stable/usage/).

This is a very concise algorithm for Hilbert curve generation. It works in `n`-dimensions. It requires little code. It comes from a little paper [^3] behind a paywall, sadly.

Most algorithms for the Hilbert curve use Gray codes to generate the shape. He observed that, instead of using the space key algorithm, which dives to each level deeper and rotates the Gray code, the algorithm could use a global transformation of all values with a Gray code and then do a minor fix-up, afterwards, so untwist it. The resulting code is much simpler than earlier efforts.

For developers, note that this algorithm relies on encoding the Hilbert index in what, to me, was a surprising order. To understand the interleaving of the Hilbert index for this algorithm, start with a 2D value where higher bits are larger subscripts, ``(a_4a_3a_2a_1, b_4b_3b_2b_1)``. Skilling encodes this as ``a_4b_4a_3b_3a_2b_2a_1b_1``, which looks good on paper, but it means the first element of the vector has the higher bits.


## Face Continuous


This is the OG Hilbert code, where the first implementation of encoding came from a paper by Butz [^4] and, later, Lawder [^5] provided the decoding algorithm. It's called `FaceContinuous` because that was its main cited property in a review of Hilbert curves [^6].

For developers, there are two errors in the [code that Lawder corrected](http://www.dcs.bbk.ac.uk/~jkl/publications.html). The first is that there is a single-bit mask, called `mask`, that should be initialized from the number of levels, not from the size of the data type. This is true for both encoding and decoding. The second is that, during decoding, the first assignment, to the highest bit of the coordinates, assigns directly from P, the highest Hilbert index bits. It should assign from `A`, which is the binary-reflected Gray code of the highest bits. These problems wouldn't show up in testing unless the highest bits in the type were used, which is an understandable oversight.



## The Compact Algorithm

This algorithm can encode Hilbert indices for Cartesian domains that aren't square. You can make a Hilbert curve with dimensions `[32, 16, 8, 32]` or any combination of powers of two.

This code fixes Hamilton's original tech report [^7] by being more careful about how axes are embedded. Hamilton posted corrections to the article's code, [Hamilton's code](https://github.com/pdebuyl/libhilbert), but these seem to fail my unit tests, as well. I tried to work through this with Hamilton years ago, but only in 2025 did I figure this out.

[^1]: Chen, Ningtau; Wang, Nengchao; Shi, Baochang, "A new algorithm for encoding and decoding the Hilbert order," in Software---Practice and Experience, 2007, 37, 897-908.

[^2]: Hamilton, Chris H., and Andrew Rau-Chaplin. "Compact Hilbert indices for multi-dimensional data." First International Conference on Complex, Intelligent and Software Intensive Systems (CISIS'07). IEEE, 2007.

[^3]: Skilling, John. "Programming the Hilbert curve." AIP Conference Proceedings. Vol. 707. No. 1. American Institute of Physics, 2004.

[^4]: Butz, Arthur R. "Alternative algorithm for Hilbert's space-filling curve." IEEE Transactions on Computers 100.4 (1971): 424-426.

[^5]: Lawder, Jonathan K. "Calculation of mappings between one and n-dimensional values using the hilbert space-filling curve." School of Computer Science and Information Systems, Birkbeck College, University of London, London Research Report BBKCS-00-01 August (2000).

[^6]: Haverkort, Herman. "Sixteen space-filling curves and traversals for d-dimensional cubes and simplices." arXiv preprint arXiv:1711.04473 (2017).

[^7]: Hamilton, Chris. "Compact hilbert indices." Dalhousie University, Faculty of Computer Science, Technical Report CS-2006-07 (2006).
