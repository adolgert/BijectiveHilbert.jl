# Compact and SpaceGray

This algorithm is used for professional database implementation because it focuses on how to pack dimensions of different sizes into the smallest-possible Hilbert index. If the data is high-dimensional and one of the indices is larger than the others, then the usual Hilbert curve encodes all indices at the same size as the largest one. This means the Hilbert index needs to use a larger data type, which means it needs more storage. Hamilton and Rau-Chaplin used mathematical manipulation to drop unused bits from the Hilbert index while keeping its central feature, that nearby data remains nearby.

Both the [`Compact`](@ref) and [`SpaceGray`](@ref) algorithms use what's called the space key algorithm for Hilbert curves. It's not recursive. It follows four steps. Quoting their paper [^1]:

1. Find the cell containing the point of interest.

2. Update the key (index) value appropriately.

3. Transform as necessary; and

4. Continue until sufficient precision has been obtained.

This is a very typical way to generate Hilbert curves, and the algorithm, which I labeled SpaceGray, it their implementation of a classic space key algorithm that relies on Gray codes. They then layer the space key algorithm with a bit-packing algorithm that relies on Gray code rank. This is a mask, to exclude unused bits, and an ordering on the remaining bits that preserves the Hilbert structure.

As a note for developers, Hamilton's original tech report [^2] has errors that look, to me, like he developed the work for two dimensions and expanded it, incompletely, for n dimensions. It's impressively-detailed math that leads to a concise formulation. I wouldn't have figured out the problems, except that there is a copy of [Hamilton's code](https://github.com/pdebuyl/libhilbert), corrected, on Github.

[^1]: Hamilton, Chris H., and Andrew Rau-Chaplin. "Compact Hilbert indices for multi-dimensional data." First International Conference on Complex, Intelligent and Software Intensive Systems (CISIS'07). IEEE, 2007.

[^2]: Hamilton, Chris. "Compact hilbert indices." Dalhousie University, Faculty of Computer Science, Technical Report CS-2006-07 (2006).
