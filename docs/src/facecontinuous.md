# Face Continuous

This is the OG Hilbert code, where the first implementation of encoding came from a paper by Butz [^1] and, later, Lawder [^2] provided the decoding algorithm. It's called `FaceContinuous` because that was its main cited property in a review of Hilbert curves [^3].

For developers, there are two errors in the [code that Lawder corrected](http://www.dcs.bbk.ac.uk/~jkl/publications.html). The first is that there is a single-bit mask, called `mask`, that should be initialized from the number of levels, not from the size of the data type. This is true for both encoding and decoding. The second is that, during decoding, the first assignment, to the highest bit of the coordinates, assigns directly from P, the highest Hilbert index bits. It should assign from `A`, which is the binary-reflected Gray code of the highest bits. These problems wouldn't show up in testing unless the highest bits in the type were used, which is an understandable oversight.

[^1]: Butz, Arthur R. "Alternative algorithm for Hilbert's space-filling curve." IEEE Transactions on Computers 100.4 (1971): 424-426.

[^2]: Lawder, Jonathan K. "Calculation of mappings between one and n-dimensional values using the hilbert space-filling curve." School of Computer Science and Information Systems, Birkbeck College, University of London, London Research Report BBKCS-00-01 August (2000).

[^3]: Haverkort, Herman. "Sixteen space-filling curves and traversals for d-dimensional cubes and simplices." arXiv preprint arXiv:1711.04473 (2017).
