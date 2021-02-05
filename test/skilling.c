#include <stdio.h>
//+++++++++++++++++++++++++++ PUBLIC-DOMAIN SOFTWARE ++++++++++++++++++++++++++
// Functions: TransposetoAxes AxestoTranspose
// Purpose: Transform in-place between Hilbert transpose and geometrical axes
// Example: b=5 bits for each of n=3 coordinates.
// 15-bit Hilbert integer = A B C D E F G H I J K L M N O is stored
// as its Transpose
// X[0] = A D G J M X[2]|
// X[1] = B E H K N <-------> | /X[1]
// X[2] = C F I L O axes |/
// high low 0------ X[0]
// Axes are stored conventially as b-bit integers.
// Author: John Skilling 20 Apr 2001 to 11 Oct 2003
//-----------------------------------------------------------------------------
typedef unsigned int coord_t;                  // char,short,int for up to 8,16,32 bits per word
void TransposetoAxes(coord_t *X, int b, int n) // position, #bits, dimension
{
    coord_t N = 2 << (b - 1), P, Q, t;
    int i;
    // Gray decode by H ^ (H/2)
    t = X[n - 1] >> 1;
    for (i = n - 1; i > 0; i--)
        X[i] ^= X[i - 1];
    X[0] ^= t;
    // Undo excess work
    for (Q = 2; Q != N; Q <<= 1)
    {
        P = Q - 1;
        for (i = n - 1; i >= 0; i--)
            if (X[i] & Q)
                X[0] ^= P; // invert
            else
            {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    } // exchange
}

void InverseUndo(coord_t *X, int b, int n)
{
    coord_t M = 1 << (b - 1), P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1)
    {
        P = Q - 1;
        for (i = 0; i < n; i++)
        {
            printf("Q=%u X[0]=%u X[%d]=%u\n", Q, X[0], i, X[i]);
            if (X[i] & Q)
            {
                X[0] ^= P; // invert
                printf("X[0]=%u P=%u\n", X[0], P);
            }
            else
            {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
                printf("X[0]=%u X[%d]=%u t=%u\n", X[0], i, X[i], t);
            }
        }
    } // exchange
}

void AxestoTranspose(coord_t *X, int b, int n) // position, #bits, dimension
{
    coord_t M = 1 << (b - 1), P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1)
    {
        P = Q - 1;
        for (i = 0; i < n; i++)
            if (X[i] & Q)
                X[0] ^= P; // invert
            else
            {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    } // exchange
    // Gray encode
    for (i = 1; i < n; i++)
        X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if (X[n - 1] & Q)
            t ^= Q - 1;
    for (i = 0; i < n; i++)
        X[i] ^= t;
}
