# Global Gray

This is a very concise algorithm for Hilbert curve generation. It works in `n`-dimensions. It requires little code. It comes from a little paper [^1] behind a paywall, sadly.

Most algorithms for the Hilbert curve use Gray codes to generate the shape. He observed that, instead of using the space key algorithm, which dives to each level deeper and rotates the Gray code, the algorithm could use a global transformation of all values with a Gray code and then do a minor fix-up, afterwards, so untwist it. The resulting code is much simpler than earlier efforts.

For developers, note that this algorithm relies on encoding the Hilbert index in what, to me, was a surprising order. To understand the interleaving of the Hilbert index for this algorithm, start with a 2D value where higher bits are larger subscripts, ``(a_4a_3a_2a_1, b_4b_3b_2b_1)``. Skilling encodes this as ``a_4b_4a_3b_3a_2b_2a_1b_1``, which looks good on paper, but it means the first element of the vector has the higher bits.

[^1]: Skilling, John. "Programming the Hilbert curve." AIP Conference Proceedings. Vol. 707. No. 1. American Institute of Physics, 2004.
