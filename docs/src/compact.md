# Compact

See the [`Usage`](https://computingkitchen.com/BijectiveHilbert.jl/stable/usage/).


## Compact

This algorithm can encode Hilbert indices for Cartesian domains that aren't square. You can make a Hilbert curve with dimensions `[32, 16, 8, 32]` or any combination of powers of two.

This code fixes Hamilton's original tech report [^1] by being more careful about how axes are embedded. Hamilton posted corrections to the article's code, [Hamilton's code](https://github.com/pdebuyl/libhilbert), but these seem to fail my unit tests, as well. I tried to work through this with Hamilton years ago, but only in 2025 did I figure this out.

[^1]: Hamilton, Chris. "Compact hilbert indices." Dalhousie University, Faculty of Computer Science, Technical Report CS-2006-07 (2006).
