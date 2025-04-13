# Simple2D

See the [`Usage`](https://computingkitchen.com/BijectiveHilbert.jl/stable/usage/).

If you want to sort some axes, then this Hilbert curve algorithm is the easiest to use. It doesn't need to know ahead of time how many bits it will need to generate a Hilbert index [^1]. As the length along each spatial axis grows, it creates gradually larger Hilbert indices to match it.

All of the algorithms use slightly different Hilbert curves. This one uses an asymmetric curve that shifts so that its endpoint is always an outside corder of each ``2^n`` x ``2^n`` tile. The next outer layer builds on the last.

[^1]: Chen, Ningtau; Wang, Nengchao; Shi, Baochang, "A new algorithm for encoding and decoding the Hilbert order," in Software---Practice and Experience, 2007, 37, 897-908.
