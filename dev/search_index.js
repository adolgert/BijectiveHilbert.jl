var documenterSearchIndex = {"docs":
[{"location":"usage/","page":"Usage","title":"Usage","text":"CurrentModule = BijectiveHilbert","category":"page"},{"location":"usage/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Modules = [BijectiveHilbert]\nPrivate = false","category":"page"},{"location":"usage/#BijectiveHilbert.brgc-Tuple{Any}","page":"Usage","title":"BijectiveHilbert.brgc","text":"Binary-reflected Gray code.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.brgc_inv-Tuple{Integer}","page":"Usage","title":"BijectiveHilbert.brgc_inv","text":"This is how the paper describes the inverse of the binary-reflected gray code.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.brgc_inv-Tuple{UInt128}","page":"Usage","title":"BijectiveHilbert.brgc_inv","text":"Inverse of the binary-reflected Gray code. This takes you from the Gray code to its index.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.decode_hilbert!-Union{Tuple{T}, Tuple{A}, Tuple{BijectiveHilbert.HilbertAlgorithm{T},Array{A,1},T}} where T where A","page":"Usage","title":"BijectiveHilbert.decode_hilbert!","text":"decode_hilbert(x, y)\n\nA 1-based Hilbert decode, from decode_hilbert_zero.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.decode_hilbert-Tuple{Integer}","page":"Usage","title":"BijectiveHilbert.decode_hilbert","text":"decode_hilbert(x, y)\n\nA 1-based Hilbert decode, from decode_hilbert_zero.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.decode_hilbert_zero-Tuple{Any}","page":"Usage","title":"BijectiveHilbert.decode_hilbert_zero","text":"decode_hilbert_zero(z::Integer) -> (x, y)\n\nComputes the (x, y) from a Hilbert code.\n\nThis function is zero-based. 0 <= x < 2^n, 0 <= y < 2^n, and 0 <= z < 4^n.\n\nSee also: encode_hilbert_zero, decode_hilbert.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.encode_hilbert-Tuple{Integer,Integer}","page":"Usage","title":"BijectiveHilbert.encode_hilbert","text":"encode_hilbert(x, y)\n\nA 1-based Hilbert code, so x and y start at one, and the z that this returns also starts at one.\n\nSee also: encode_hilbert_zero.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.encode_hilbert-Union{Tuple{T}, Tuple{A}, Tuple{BijectiveHilbert.HilbertAlgorithm{T},Array{A,1}}} where T where A","page":"Usage","title":"BijectiveHilbert.encode_hilbert","text":"encode_hilbert(gg::HilbertAlgorithm{T}, X::Vector{A})\n\nA 1-based Hilbert encoding.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.encode_hilbert_zero-Tuple{Integer,Integer}","page":"Usage","title":"BijectiveHilbert.encode_hilbert_zero","text":"encode_hilbert_zero(x::Integer, y::Integer)\n\nComputes an integer Hilbert index for x and y using a variant algorithm.\n\nGiven two integer indices for a 2-dimensional plane, return a single index. This index is designed to increase locality for 1-dimensional access. It does this by keeping nearby points in 2 dimensions also nearby in 1 dimension.\n\nx and y need to be integers that have bit-shifting operations.\n\nThe variant algorithm used differs from the usual Hilbert code because it doesn't need to know the size of the whole grid before computing the code [1]. It looks like a slightly-rotated version of the Hilbert curve, but it has the benefit that it is 1-1 between (x, y) and z, so you can translate back and forth.\n\nThis function is zero-based. 0 <= x < 2^n, 0 <= y < 2^n, and the result is 0 <= z < 4^n.\n\nSee also: decode_hilbert_zero, encode_hilbert.\n\n[1]: N. Chen, N. Wang, B. Shi, A new algorithm for encoding and decoding the Hilbert order. Software—Practice and Experience 2007; 37(8): 897–908.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.hi_d-Tuple{Any,Any}","page":"Usage","title":"BijectiveHilbert.hi_d","text":"page 12. directions d(i) = 0 if i=0     = g(i-1) mod n if i=0 mod 2     = g(i) mod n if i=1 mod 2. Domain is 0 <= i <= 2^n - 1.\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.hi_f-Tuple{Any,Any}","page":"Usage","title":"BijectiveHilbert.hi_f","text":"f is the exit vertex of the ith sub-hypercube in a Gray code ordering of the sub-hypercubes. Corllary 2.7 on pg. 12 says:     hif(i, n) = hie(i) ⊻ hi_d(i, n)\n\n\n\n\n\n","category":"method"},{"location":"usage/#BijectiveHilbert.hilbert_order-Tuple{AbstractArray,Any}","page":"Usage","title":"BijectiveHilbert.hilbert_order","text":"hilbert_order(v::AbstractArray, subdivisions)\n\nCalculates the permutation of v that would order a 2D array by Hilbert curve, where v is an array of real numbers, dimension (2, N) and subdivisions is a real number that specifies how many boxes to place the v into, per side.\n\nThis tells you how to order an array of 2-dimensional numbers so that they have more memory locality in 1 dimension. Given an array of real numbers of dimension (2, n), subdivide them in each dimension by subdivisions, and assign each point a Hilbert code. Return the permutation that would sort the given array by that Hilbert code.\n\nSee also: encode_hilbert_zero.\n\nExample\n\nrng = MersenneTwister(984720987)\npoints_in_space = zeros(2, 100)\nrand!(rng, points_in_space)\npoints_reordered = points_in_space[:, hilbert_order(points_in_space, 50)]\n\n\n\n\n\n","category":"method"},{"location":"#BijectiveHilbert","page":"Home","title":"BijectiveHilbert","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package offers a function that turns two integer (i, j)-coordinates into a single integer z-coordinate, and vice versa.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> xy = zeros(Int, 8, 8)\njulia> for y in 1:size(xy, 2)\n           for x in 1:size(xy, 1)\n               z = encode_hilbert(x, y)\n               xy[x, y] = z\n           end\n       end\njulia> @assert decode_hilbert(xy[5, 7]) == (5, 7)\njulia> xy\n8×8 Array{Int64,2}:\n  1   2  15  16  17  20  21  22\n  4   3  14  13  18  19  24  23\n  5   8   9  12  31  30  25  26\n  6   7  10  11  32  29  28  27\n 59  58  55  54  33  36  37  38\n 60  57  56  53  34  35  40  39\n 61  62  51  52  47  46  41  42\n 64  63  50  49  48  45  44  43","category":"page"},{"location":"","page":"Home","title":"Home","text":"You use this function, called a Hilbert curve, if you have a two-dimensional array of data and want to improve memory locality when you access contiguous blocks of this data. It's used most often for geospatial work but is equally appropriate for dealing with large TIFF files. It helps memory locality because the order of the z-coordinate tends to pick values that are near each other in space. It belongs to the class of space-filling, self-avoiding, simple, and self-similar (FASS) curves, which includes Peano curves, and Morton z-curves.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This particular version of the function is easier to use than the traditional Hilbert curve because you don't have to tell the function how large the 2D array of coordinates is, or will eventually be. For a traditional Hilbert curve, the same value of z can refer to different values of (i, j). It comes from a paper by Chen, Wang, and Shi, titled, \"A new algorithm for encoding and decoding the Hilbert order,\" in Software Practice and Experience, 2007, vol. 37, pg 897-908.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Other space-filling, self-avoiding, simple, and self-similar curves:","category":"page"},{"location":"","page":"Home","title":"Home","text":"HilbertSpaceFillingCurve - An optimized version of the Hilbert curve.\nMorton z-curve - Used for quadtrees and octrees.","category":"page"},{"location":"","page":"Home","title":"Home","text":"We should make a group repo and include:","category":"page"},{"location":"","page":"Home","title":"Home","text":"The uneven-sided Hilbert curve: Hamilton, Chris H., and Andrew Rau-Chaplin. 2008. “Compact Hilbert Indices: Space-Filling Curves for Domains with Unequal Side Lengths.” Information Processing Letters 105 (5): 155–63.\nPeano curves\nHilbert and Morton curves in arbitrary dimensions.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add BijectiveHilbert","category":"page"}]
}
