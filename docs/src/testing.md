# Testing

## Hilbert Curve Properties

 - **Encode and decode are inverses.** - Show that if you decode and encode every Hilbert index in a domain you get the original point back.
 - **Domain is covered** - Test all of the points in the domain and ensure you get unique results.
 - **Consecutive Hilbert indices are neighboring points** - This is tested as a Manhattan index.

## Test Against Lean/C

 - **Lean proof of Butz algorithm** - This is used by Butz and Hamilton for the SpaceGray algorithm.
 - **Lean proof of the Compact algorithm** - This proves the Compact algorithm will work for all cases, not just all cases below 8 dimensions that fit into 32GB of memory.

 For these, I used Lean to create a callable library. Then I made a C program that matched every Hilbert index and point value. Lastly, that C program was used to match the Julia code.
 