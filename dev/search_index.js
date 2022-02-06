var documenterSearchIndex = {"docs":
[{"location":"globalgray/#Global-Gray","page":"Global Gray","title":"Global Gray","text":"","category":"section"},{"location":"globalgray/","page":"Global Gray","title":"Global Gray","text":"This is a very concise algorithm for Hilbert curve generation. It works in n-dimensions. It requires little code. It comes from a little paper [1] behind a paywall, sadly.","category":"page"},{"location":"globalgray/","page":"Global Gray","title":"Global Gray","text":"Most algorithms for the Hilbert curve use Gray codes to generate the shape. He observed that, instead of using the space key algorithm, which dives to each level deeper and rotates the Gray code, the algorithm could use a global transformation of all values with a Gray code and then do a minor fix-up, afterwards, so untwist it. The resulting code is much simpler than earlier efforts.","category":"page"},{"location":"globalgray/","page":"Global Gray","title":"Global Gray","text":"For developers, note that this algorithm relies on encoding the Hilbert index in what, to me, was a surprising order. To understand the interleaving of the Hilbert index for this algorithm, start with a 2D value where higher bits are larger subscripts, (a_4a_3a_2a_1 b_4b_3b_2b_1). Skilling encodes this as a_4b_4a_3b_3a_2b_2a_1b_1, which looks good on paper, but it means the first element of the vector has the higher bits.","category":"page"},{"location":"globalgray/","page":"Global Gray","title":"Global Gray","text":"[1]: Skilling, John. \"Programming the Hilbert curve.\" AIP Conference Proceedings. Vol. 707. No. 1. American Institute of Physics, 2004.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"CurrentModule = BijectiveHilbert","category":"page"},{"location":"usage/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"All of the Hilbert algorithms have the same interface.","category":"page"},{"location":"usage/#Decide-dimensions-of-the-spatial-axes","page":"Usage","title":"Decide dimensions of the spatial axes","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"All of the algorithms, except Simple2D, need to know ahead of time the extent of the coordinate system. If the Hilbert curve will be in a 32 x 32 space, then the dimension is 2, and it needs log2(32)=5 bits of resolution in both directions. For GlobalGray, that's b=5, n=2. If the sizes aren't powers of two or are uneven, then set the bits to cover the largest side, so (12 x 12 x 12 x 800) would be b=10, n=4 because 800  2^10. There is one algorithm that deals better with uneven sides. The Compact algorithm can pack bits together so that the resulting Hilbert index takes less storage space. This is usually used for database storage. It could take (12 x 12 x 12 x 800) as Compact([4, 4, 4, 10]) which results in a 24-bit integer that can be stored in a UInt32.","category":"page"},{"location":"usage/#Create-an-algorithm","page":"Usage","title":"Create an algorithm","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"For most Hilbert curve algorithms, you have to say, beforehand, how large the multidimensional coordinates are, in powers of two. For instance, a three-dimensional grid can have values from 1 to 16 in each dimension, so n = 3 and b = 4 because 16 = 2^b.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"using BijectiveHilbert\ndimensions = 3\nbits = 4\nsimple = Simple2D(Int)\ngg = GlobalGray(UInt, bits, dimensions)\nsg = SpaceGray(UInt, bits, dimensions)\ncompact = Compact(UInt, fill(bits, dimensions))","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"The first argument is a datatype for the Hilbert index. It should be large enough to hold all of the bits from the n-dimensional axes. If you don't specify one, it will use the smallest unsigned integer that can hold them.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Note that the Compact algorithm can have different sizes for each dimension, but they are all powers of two. It produces a Hilbert index that uses only as many bits as necessary.","category":"page"},{"location":"usage/#Encode-and-decode","page":"Usage","title":"Encode and decode","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"You can encode from n-dimensions to the Hilbert index, or you can decode from a Hilbert index to n-dimensions.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"for algorithm in [simple, gg, sg, compact]\n    for k in 1:(1<<bits)\n        for j in 1:(1<<bits)\n            for i in 1:(1<<bits)\n                X = [i, j, k]\n                h = encode_hilbert(algorithm, X)\n                X .= 0\n                decode_hilbert!(algorithm, X, h)\n            end\n        end\n    end\nend","category":"page"},{"location":"usage/#Encode-and-decode-with-a-zero-based-value","page":"Usage","title":"Encode and decode with a zero-based value","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"The underlying algorithms use a zero-based axis and a zero-based Hilbert index. These are available, too.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"for algorithm in [simple, gg, sg, compact]\n    for k in 0:(1<<bits - 1)\n        for j in 0:(1<<bits - 1)\n            for i in 0:(1<<bits - 1)\n                X = [i, j, k]\n                h = encode_hilbert_zero(algorithm, X)\n                X .= 0\n                decode_hilbert_zero!(algorithm, X, h)\n            end\n        end\n    end\nend","category":"page"},{"location":"usage/#Index","page":"Usage","title":"Index","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Modules = [BijectiveHilbert]\nPrivate = false","category":"page"},{"location":"usage/#BijectiveHilbert.Bijective.Simple2D","page":"Usage","title":"BijectiveHilbert.Bijective.Simple2D","text":"Simple2D(::Type{T})\n\nThe type is a data type to hold the Hilbert index. It should have enough bits to hold all bits of integer axes to encode.\n\nThe variant algorithm used differs from the usual Hilbert code because it doesn't need to know the size of the whole grid before computing the code. It looks like a slightly-rotated version of the Hilbert curve, but it has the benefit that it is 1-1 between (x, y) and z, so you can translate back and forth.\n\nIf the size of the axis values or the size of the Hilbert index are too large to be stored in the data types of the axis vector or index, then you'll see an error from a failure of trunc, which tries to copy values into the smaller datatype. Solve this by using a larger datatype, so UInt64 instead of UInt32.\n\nIt comes from a paper:  N. Chen, N. Wang, B. Shi, A new algorithm for encoding and decoding the Hilbert order. Software—Practice and Experience 2007; 37(8): 897–908.\n\n\n\n\n\n","category":"type"},{"location":"usage/#BijectiveHilbert.Compact","page":"Usage","title":"BijectiveHilbert.Compact","text":"Compact(ms::Vector{Int})\nCompact(::Type{T}, ms::Vector{Int})\n\nThis algorithm is n-dimensional and permits dimensions to use different numbers of bits, specified in the ms vector. The type T is an optional data type for the Hilbert index. It should be greater than or equal to the sum of the bits.\n\nThis algorithm comes from three sources:\n\nA technical report, \"Compact Hilbert Indices\" by Chris Hamilton. Technical Report CS-2006-07. 6059 University Ave., Halifax, Nova Scotia, B3H 1W5, Canada. This report is informative but has many errors.\nA paper by Hamilton and Rau-Chaplin, \"Compact Hilbert Indices for Multi-Dimensional Data,\" 2007. Nice paper. Also wrong.\nThe libhilbert source code is a copy of Hamilton's work and has many corrections. This, ultimately, lead to the working code.\n\n\n\n\n\n","category":"type"},{"location":"usage/#BijectiveHilbert.FaceContinuous","page":"Usage","title":"BijectiveHilbert.FaceContinuous","text":"FaceContinuous(b::Int, n::Int)\nFaceContinuous(::Type{T}, b::Int, n::Int)\n\nFor n dimensions, use b bits of precision in this Hilbert curve. If you specify a type T, this will be used as the type of the Hilbert encoding. If not, the smallest unsigned integer that can hold n*b bits will be used for the Hilbert index data type.\n\nThis is the Butz algorithm, as presented by Lawder. Haverkort2017 says it is face continuous. The code is in lawder.c. The original paper had an error, and Lawder put a correction on his website. http://www.dcs.bbk.ac.uk/~jkl/publications.html\n\n\n\n\n\n","category":"type"},{"location":"usage/#BijectiveHilbert.GlobalGray","page":"Usage","title":"BijectiveHilbert.GlobalGray","text":"GlobalGray(b, n)\nGlobalGray(T, b, n)\n\nT is a data type for the Hilbert index. It can be signed or unsigned, as long as it has at least n * b bits. n is the number of dimensions, and b is the bits per dimension, so each axis value should be between 0 and 2^b - 1, inclusive, for the zero-based interface. They should be between 1 and 2^b, inclusive, for the one-based interface.\n\nThe GlobalGray algorithm is an n-dimensional Hilbert curve with a simplified implementation. It follows an article, \"Programming the Hilbert Curve,\" by John Skilling, 707 (2004), http://dx.doi.org/10.1063/1.1751381. I call it \"Global Gray\" because the insight of the article is that a single, global Gray code can be applied to all np bits of a Hilbert length.\n\n\n\n\n\n","category":"type"},{"location":"usage/#BijectiveHilbert.SpaceGray","page":"Usage","title":"BijectiveHilbert.SpaceGray","text":"SpaceGray(b, n)\nSpaceGray(::Type{T}, b, n)\n\nThis is an n-dimensional Hilbert curve where all n dimensions must have b bits in size. It was described in the same paper and examples as the Compact algorithm.\n\n\n\n\n\n","category":"type"},{"location":"usage/#BijectiveHilbert.axis_type-Tuple{FaceContinuous}","page":"Usage","title":"BijectiveHilbert.axis_type","text":"Used during testing to pick a type for the xyz coordinate.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.decode_hilbert!-Union{Tuple{T}, Tuple{A}, Tuple{BijectiveHilbert.HilbertAlgorithm{T}, Vector{A}, T}} where {A, T}","page":"Usage","title":"BijectiveHilbert.decode_hilbert!","text":"decode_hilbert!(ha::HilbertAlgorithm{T}, X::Vector{A})\n\nA 1-based Hilbert decode, from decode_hilbert_zero!. Both the Hilbert index and the axes start counting at 1 instead of 0.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.decode_hilbert_zero!-Union{Tuple{T}, Tuple{GlobalGray{T}, Vector, T}} where T","page":"Usage","title":"BijectiveHilbert.decode_hilbert_zero!","text":"decode_hilbert_zero!(ha::HilbertAlgorithm{T}}, X::Vector{A}, h::T)\n\nGiven a Hilbert index, h, computes an n-dimensional coordinate X. The type of the Hilbert index is large enought to contain the bits of all dimensions of the axis vector, X.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.encode_hilbert-Union{Tuple{T}, Tuple{A}, Tuple{BijectiveHilbert.HilbertAlgorithm{T}, Vector{A}}} where {A, T}","page":"Usage","title":"BijectiveHilbert.encode_hilbert","text":"encode_hilbert(ha::HilbertAlgorithm{T}, X::Vector{A})\n\nA 1-based Hilbert encoding. Both the Hilbert index and the axes start counting at 1 instead of 0.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.encode_hilbert_zero-Union{Tuple{T}, Tuple{GlobalGray{T}, Vector}} where T","page":"Usage","title":"BijectiveHilbert.encode_hilbert_zero","text":"encode_hilbert_zero(ha::HilbertAlgorithm{T}, X::Vector{A})\n\nTakes an n-dimensional vector X and returns a single integer of type T which orders X to improve spatial locality. The input X has multiple axes and the output is called a Hilbert index. This version is zero-based, so each axis counts from 0, and the smallest Hilbert index is 0.\n\n\n\n\n\n","category":"method"},{"location":"simple2d/#Simple2D","page":"Simple2D","title":"Simple2D","text":"","category":"section"},{"location":"simple2d/","page":"Simple2D","title":"Simple2D","text":"If you want to sort some axes, then this Hilbert curve algorithm is the easiest to use. It doesn't need to know ahead of time how many bits it will need to generate a Hilbert index [1]. As the length along each spatial axis grows, it creates gradually larger Hilbert indices to match it.","category":"page"},{"location":"simple2d/","page":"Simple2D","title":"Simple2D","text":"All of the algorithms use slightly different Hilbert curves. This one uses an asymmetric curve that shifts so that its endpoint is always an outside corder of each 2^n x 2^n tile. The next outer layer builds on the last.","category":"page"},{"location":"simple2d/","page":"Simple2D","title":"Simple2D","text":"[1]: Chen, Ningtau; Wang, Nengchao; Shi, Baochang, \"A new algorithm for encoding and decoding the Hilbert order,\" in Software–-Practice and Experience, 2007, 37, 897-908.","category":"page"},{"location":"facecontinuous/#Face-Continuous","page":"Face Continuous","title":"Face Continuous","text":"","category":"section"},{"location":"facecontinuous/","page":"Face Continuous","title":"Face Continuous","text":"This is the OG Hilbert code, where the first implementation of encoding came from a paper by Butz [1] and, later, Lawder [2] provided the decoding algorithm. It's called FaceContinuous because that was its main cited property in a review of Hilbert curves [3].","category":"page"},{"location":"facecontinuous/","page":"Face Continuous","title":"Face Continuous","text":"For developers, there are two errors in the code that Lawder corrected. The first is that there is a single-bit mask, called mask, that should be initialized from the number of levels, not from the size of the data type. This is true for both encoding and decoding. The second is that, during decoding, the first assignment, to the highest bit of the coordinates, assigns directly from P, the highest Hilbert index bits. It should assign from A, which is the binary-reflected Gray code of the highest bits. These problems wouldn't show up in testing unless the highest bits in the type were used, which is an understandable oversight.","category":"page"},{"location":"facecontinuous/","page":"Face Continuous","title":"Face Continuous","text":"[1]: Butz, Arthur R. \"Alternative algorithm for Hilbert's space-filling curve.\" IEEE Transactions on Computers 100.4 (1971): 424-426.","category":"page"},{"location":"facecontinuous/","page":"Face Continuous","title":"Face Continuous","text":"[2]: Lawder, Jonathan K. \"Calculation of mappings between one and n-dimensional values using the hilbert space-filling curve.\" School of Computer Science and Information Systems, Birkbeck College, University of London, London Research Report BBKCS-00-01 August (2000).","category":"page"},{"location":"facecontinuous/","page":"Face Continuous","title":"Face Continuous","text":"[3]: Haverkort, Herman. \"Sixteen space-filling curves and traversals for d-dimensional cubes and simplices.\" arXiv preprint arXiv:1711.04473 (2017).","category":"page"},{"location":"compact/#Compact-and-SpaceGray","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"","category":"section"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"This algorithm is used for professional database implementation because it focuses on how to pack dimensions of different sizes into the smallest-possible Hilbert index. If the data is high-dimensional and one of the indices is larger than the others, then the usual Hilbert curve encodes all indices at the same size as the largest one. This means the Hilbert index needs to use a larger data type, which means it needs more storage. Hamilton and Rau-Chaplin used mathematical manipulation to drop unused bits from the Hilbert index while keeping its central feature, that nearby data remains nearby.","category":"page"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"Both the Compact and SpaceGray algorithms use what's called the space key algorithm for Hilbert curves. It's not recursive. It follows four steps. Quoting their paper [1]:","category":"page"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"Find the cell containing the point of interest.\nUpdate the key (index) value appropriately.\nTransform as necessary; and\nContinue until sufficient precision has been obtained.","category":"page"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"This is a very typical way to generate Hilbert curves, and the algorithm, which I labeled SpaceGray, it their implementation of a classic space key algorithm that relies on Gray codes. They then layer the space key algorithm with a bit-packing algorithm that relies on Gray code rank. This is a mask, to exclude unused bits, and an ordering on the remaining bits that preserves the Hilbert structure.","category":"page"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"As a note for developers, Hamilton's original tech report [2] has errors that look, to me, like he developed the work for two dimensions and expanded it, incompletely, for n dimensions. It's impressively-detailed math that leads to a concise formulation. I wouldn't have figured out the problems, except that there is a copy of Hamilton's code, that he, himself, corrected.","category":"page"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"[1]: Hamilton, Chris H., and Andrew Rau-Chaplin. \"Compact Hilbert indices for multi-dimensional data.\" First International Conference on Complex, Intelligent and Software Intensive Systems (CISIS'07). IEEE, 2007.","category":"page"},{"location":"compact/","page":"Compact and SpaceGray","title":"Compact and SpaceGray","text":"[2]: Hamilton, Chris. \"Compact hilbert indices.\" Dalhousie University, Faculty of Computer Science, Technical Report CS-2006-07 (2006).","category":"page"},{"location":"hilbert/#Pseudo-Hilbert-Curves-for-Computational-Science","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"","category":"section"},{"location":"hilbert/#Introduction","page":"Pseudo-Hilbert Curves for Computational Science","title":"Introduction","text":"","category":"section"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"The Pseudo-Hilbert curve sorts regular arrays of points into a linear order that keeps nearby points in the array close to each other in the linear order. Scientific computing algorithms use this capability in a few ways.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"An easy form of clustering.\nSpeed up queries for databases.\nAllocate work for distributed computing.\nVisualize correlation of large one-dimensional datasets.\nVisualize change over time of two-dimensional datasets.\nApply a one-dimensional global optimization algorithm to a multi-dimensional dataset.\nAn R-tree algorithm that uses Hilbert curves to pick node locations.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"These widely different applications all use the Hilbert curve, not for drawing, but to convert indices from multiple dimensions to a single dimension and to convert them back to multiple dimensions, as needed. For high performance computing, it's not a drawn curve but an algorithm.","category":"page"},{"location":"hilbert/#The-algorithm","page":"Pseudo-Hilbert Curves for Computational Science","title":"The algorithm","text":"","category":"section"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"There are two functions. One accepts an array of natural numbers as input and returns a single Hilbert index as output. We call that encoding. The other function decodes the single Hilbert index into the original natural numbers. There are a lot of pairs of functions that can do this.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"For instance, you may recall from your earliest math class that we use such a pair of functions to prove that the number of integers in the first quadrant of a two-dimensional graph are countably infinite. Starting from (0, 0), move to (0, 1) and then (1, 0). Then move from (2, 0) to (1, 1) to (0, 2). This weaves along diagonals (in boustrophedon order), sweeping out a one-dimensional covering of two dimensions. Each point in 2D has a corresponding, countable point in 1D. Nor is this the only example, by far.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"Julia stores its matrices in column-major order. C stores its matrices in row-major order. Both are orderings of two-dimensional rectangles according to a linear storage order. So, too, is block storage of data. These are functions like the Hilbert curve because we can automatically translate from (i, j) to a single index in memory, if we know the dimensions of the matrix.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"The Hilbert algorithm is equivalent to the use of column-major order with one addition. Two points near each other according to the Hilbert index will tend to be closer to each other if we measure the Euclidean distance between their two-dimensional indices. There isn't a guaranteed upper bound on how far apart two points can be, but experiments comparing nearness, called spatial locality, confirm that the Hilbert curve tends to arrange points closer to each other than does column-major storage, block storage, or other self-similar space-filling curves, such as the Peano curve.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"That's why this algorithm is useful. Given a bunch of (i, j) indices, the Hilbert index helps to arrange them in a single vector such that nearby indices tend to be closer to each other. Using a Hilbert algorithm can be a little more complicated than that, only because tradition asks that you focus on counts of bits. Let's look at why this is an easy, if important, requirement.","category":"page"},{"location":"hilbert/#Counting-bits","page":"Pseudo-Hilbert Curves for Computational Science","title":"Counting bits","text":"","category":"section"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"Like a column-major order for a matrix, a Hilbert algorithm needs to know, before doing a conversion, what extent the data might have. Most implementations of Hilbert curves assume that every dimension of the input coordinates will have the same size and will be a power of two. If your data is 25 x 38, that means you need a 64 x 64 grid on which to work. That means the i-axis requires 6 bits and the j-axis requires 6-bits. As a result, the Hilbert index will run from 0 to 4095, needing 12-bits in which to store its result.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"If we don't use the full size of the allocated axes, are we not going to get as good spatial locality from the Hilbert curve? This is only approximate to start with. The way to use this algorithm is to create axes that are large enough and not worry about the unused bits of information.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"The bit count matters for choosing types of function arguments. It's usually easy for the coordinates to hold their integer data, but higher dimensions can make the resulting Hilbert index too large for most data types. Eight dimensions of ten bits is 80 bits, which exceeds a typical Int64. That's fine in Julia, which has both Int128 and BigInt types. In addition, most of the math is defined for unsigned integers, but the algorithms in this library are tested to work up to the last sign-bit for signed integers, so they are OK to use.","category":"page"},{"location":"hilbert/#Properties-of-different-Hilbert-curves","page":"Pseudo-Hilbert Curves for Computational Science","title":"Properties of different Hilbert curves","text":"","category":"section"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"There are multiple Hilbert curves. There is even a paper called, \"How many three-dimensional Hilbert curves are there?\" by Haverkort (2017). The different curves have different orientataions relative to their axes. They can be symmetric or asymmetric. They fill space in slightly different ways, (almost) all of which require domains that are powers of two in each direction. Haverkort also wrote a lovely paper describing these differences, in \"Sixteen space-filling curves and traversals for d-dimensional cubes and simplices.\"","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"For use in scientific computing, it may be more important to judge different Hilbert curve algorithms on capability, speed, and relaxation of constraints.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"The Compact algorithm creates a Hilbert curve where each axis can have a size that's a different power of two.\nThe Simple2D algorithm doesn't need to know how large the axes may be before you use it, but it only works in 2D.\nThe GlobalGray algorithm is fast for an n-dimensional algorithm.\nThe FaceContinuous algorithm is slower and is included because it has a different shape and is historically important as the first non-recursive n-dimensional algorithm.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"In general, algorithms that are written explicitly for 2D are faster than n-dimensional equivalents.","category":"page"},{"location":"hilbert/#Other-work-around-Hilbert-curves","page":"Pseudo-Hilbert Curves for Computational Science","title":"Other work around Hilbert curves","text":"","category":"section"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"This library provides two functions for each algorithm, one to encode and one to decode. There is considerable work around making more-efficient queries of data stored in Hilbert order. These algorithms take a single Hilbert index and return the Hilbert index in any coordinate direction. They don't need to convert to coordinates, find the neighbor in that space, and convert back to the Hilbert index. There is also work to make rectangular queries faster, in general. These find sequences of Hilbert indexes that fall within a rectangular region in coordinate space.","category":"page"},{"location":"hilbert/","page":"Pseudo-Hilbert Curves for Computational Science","title":"Pseudo-Hilbert Curves for Computational Science","text":"None of those are implemented here. They would be useful for advanced database work.","category":"page"},{"location":"#BijectiveHilbert","page":"Home","title":"BijectiveHilbert","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Data in two or more dimensions is stored linearly, so places nearby in the data can be far in memory or on disk. This package offers functions that make multi-dimensional data storage and retrieval more efficient by sorting them so that nearby data is nearby in memory.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.add(\"BijectiveHilbert\")\njulia> using BijectiveHilbert\njulia> xy = zeros(Int, 8, 8)\njulia> for y in 1:size(xy, 2)\n           for x in 1:size(xy, 1)\n               z = encode_hilbert(Simple2D(Int), [x, y])\n               xy[x, y] = z\n           end\n       end\njulia> X = zeros(Int, 2)\njulia> decode_hilbert!(Simple2D(Int), X, xy[5, 7])\njulia> X == [5, 7]\njulia> xy\n8×8 Array{Int64,2}:\n  1   2  15  16  17  20  21  22\n  4   3  14  13  18  19  24  23\n  5   8   9  12  31  30  25  26\n  6   7  10  11  32  29  28  27\n 59  58  55  54  33  36  37  38\n 60  57  56  53  34  35  40  39\n 61  62  51  52  47  46  41  42\n 64  63  50  49  48  45  44  43","category":"page"},{"location":"","page":"Home","title":"Home","text":"This function, called a Hilbert curve, is used most often for geospatial work or database implementation but is equally appropriate for dealing with large TIFF files. It belongs to the class of space-filling, self-avoiding, simple, and self-similar (FASS) curves, which includes Peano curves, and Morton z-curves.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Included are several variations of the Hilbert curve. They are type-stable and thoroughly tested, including bug fixes on three of the implementations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Simple2D, shown above, two-dimensional. Doesn't need to know axis dimensions, from Chen, Wang, and Shi.\nGlobalGray, an n-dimensional curve where all axis dimensions must be equal, from Skilling.\nSpaceGray, an n-dimensional curve with a different path. All axis dimensions must be equal, from Hamilton and Rau-Chaplin.\nCompact, an n-dimensional curve the permits each axis to be a different size, from Hamilton and Rau-Chaplin.\nFaceContinuous, an n-dimensional curve, the oldest version from Butz and Lawder.","category":"page"}]
}