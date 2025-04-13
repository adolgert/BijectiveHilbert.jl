# SpaceGray and Compact

See the [`Usage`](https://computingkitchen.com/BijectiveHilbert.jl/stable/usage/).

## SpaceGray

The [`SpaceGray`](@ref) algorithm is the fastest way to generate multi-dimensional Hilbert curves. This library implements a classic space key algorithm that relies on Gray codes. It's not recursive. It follows four steps. Quoting the source paper [^1]:

1. Find the cell containing the point of interest.

2. Update the key (index) value appropriately.

3. Transform as necessary; and

4. Continue until sufficient precision has been obtained.


## Compact

The Compact algorithm was the star of this library, but every implementation I can find or make is broken, so I've removed it from the interface, but the papers and code are still there if someone knows better.

This algorithm can encode Hilbert indices for Cartesian domains that aren't square. This is very helpful for skewed coordinate systems.

As a note for developers, Hamilton's original tech report [^2] has errors that look, to me, like he developed the work for two dimensions and expanded it, incompletely, for n dimensions. It's impressively-detailed math that leads to a concise formulation. Hamilton posted corrections to the article's code, [Hamilton's code](https://github.com/pdebuyl/libhilbert), but these seem to fail my unit tests, as well.

[^1]: Hamilton, Chris H., and Andrew Rau-Chaplin. "Compact Hilbert indices for multi-dimensional data." First International Conference on Complex, Intelligent and Software Intensive Systems (CISIS'07). IEEE, 2007.

[^2]: Hamilton, Chris. "Compact hilbert indices." Dalhousie University, Faculty of Computer Science, Technical Report CS-2006-07 (2006).
