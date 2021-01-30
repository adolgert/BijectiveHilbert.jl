/*
"Calculation of Mappings Between One and n-dimensional Values
Using the Hilbert Space-filling Curve", by J K Lawder.
Technical Report JL1/00, Aug 15, 2000. School of Computer Science
and Information Systems, Birkbeck College, University of London, UK.

This code assumes the following:
The macro ORDER corresponds to the order of curve and is 32,
thus coordinates of points are 32 bit values.
A U_int should be a 32 bit unsigned integer.
The macro DIM corresponds to the number of dimensions in a
space.
The derived-key of a Point is stored in an Hcode which is an
array of U_int. The bottom bit of a derived-key is held in the
bottom bit of the hcode[0] element of an Hcode and the top bit
of a derived-key is held in the top bit of the hcode[DIM-1]
element of and Hcode.
g_mask is a global array of masks which helps simplyfy some
calculations - it has DIM elements. In each element, only one
bit is zeo valued - the top bit in element no. 0 and the bottom
bit in element no. (DIM - 1). eg.
#if DIM == 5 const U_int g_mask[] = {16, 8, 4, 2, 1}; #endif
#if DIM == 6 const U_int g_mask[] = {32, 16, 8, 4, 2, 1}; #endif
etc...
*/
#define DIM 3
typedef unsigned int U_int;
#define ORDER 32

typedef struct
{
    U_int hcode[DIM];
} Hcode;
typedef Hcode Point;

/*===========================================================*/
/* calc_P */
/*===========================================================*/
U_int calc_P(int i, Hcode H)
{
    int element;
    U_int P, temp1, temp2;
    element = i / ORDER;
    P = H.hcode[element];
    if (i % ORDER > ORDER - DIM)
    {
        temp1 = H.hcode[element + 1];
        P >>= i % ORDER;
        temp1 <<= ORDER - i % ORDER;
        P |= temp1;
    }
    else
        P >>= i % ORDER; /* P is a DIM bit hcode */
/* the & masks out spurious highbit values */
#if DIM < ORDER
    P &= (1 << DIM) - 1;
#endif
    return P;
}

/*===========================================================*/
/* calc_P2 */
/*===========================================================*/
U_int calc_P2(U_int S)
{
    int i;
    U_int P;
    P = S & g_mask[0];
    for (i = 1; i < DIM; i++)
        if (S & g_mask[i] ^ (P >> 1) & g_mask[i])
            P |= g_mask[i];
    return P;
}

/*===========================================================*/
/* calc_J */
/*===========================================================*/
U_int calc_J(U_int P)
{
    int i;
    U_int J;
    J = DIM;
    for (i = 1; i < DIM; i++)
        if ((P >> i & 1) == (P & 1))
            continue;
        else
            break;
    if (i != DIM)
        J -= i;
    return J;
}

/*===========================================================*/
/* calc_T */
/*===========================================================*/
U_int calc_T(U_int P)
{
    if (P < 3)
        return 0;
    if (P % 2)
        return (P - 1) ^ (P - 1) / 2;
    return (P - 2) ^ (P - 2) / 2;
}

/*===========================================================*/
/* calc_tS_tT */
/*===========================================================*/
U_int calc_tS_tT(U_int xJ, U_int val)
{
    U_int retval, temp1, temp2;
    retval = val;
    if (xJ % DIM != 0)
    {
        temp1 = val >> xJ % DIM;
        temp2 = val << DIM - xJ % DIM;
        retval = temp1 | temp2;
        retval &= ((U_int)1 << DIM) - 1;
    }
    return retval;
}

/*===========================================================*/
/* H_decode */
/*===========================================================*/
/* For mapping from one dimension to DIM dimensions */
Point H_decode(Hcode H)
{
    U_int mask = (U_int)1 << ORDER - 1,
          A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    Point pt = {0};
    int i = ORDER * DIM - DIM, j;
    P = calc_P(i, H);
    J = calc_J(P);
    xJ = J - 1;
    A = S = tS = P ^ P / 2;
    T = calc_T(P);
    tT = T;
    /*--- distrib bits to coords ---*/
    for (j = DIM - 1; P > 0; P >>= 1, j--)
        if (P & 1)
            pt.hcode[j] |= mask;
    for (i -= DIM, mask >>= 1; i >= 0; i -= DIM, mask >>= 1)
    {
        P = calc_P(i, H);
        S = P ^ P / 2;
        tS = calc_tS_tT(xJ, S);
        W ^= tT;
        A = W ^ tS;
        /*--- distrib bits to coords ---*/
        for (j = DIM - 1; A > 0; A >>= 1, j--)
            if (A & 1)
                pt.hcode[j] |= mask;
        if (i > 0)
        {
            T = calc_T(P);
            tT = calc_tS_tT(xJ, T);
            J = calc_J(P);
            xJ += J - 1;
        }
    }
    return pt;
}

/*===========================================================*/
/* H_encode */
/*===========================================================*/
/* For mapping from DIM dimensions to one dimension */
Hcode H_encode(Point pt)
{
    U_int mask = (U_int)1 << ORDER - 1, element,
          A, W = 0, S, tS, T, tT, J, P = 0, xJ;
    Hcode h = {0};
    int i = ORDER * DIM - DIM, j;
    for (j = A = 0; j < DIM; j++)
        if (pt.hcode[j] & mask)
            A |= g_mask[j];
    S = tS = A;
    P = calc_P2(S);
    /* add in DIM bits to hcode */
    element = i / ORDER;
    if (i % ORDER > ORDER - DIM)
    {
        h.hcode[element] |= P << i % ORDER;
        h.hcode[element + 1] |= P >> ORDER - i % ORDER;
    }
    else
        h.hcode[element] |= P << i - element * ORDER;
    J = calc_J(P);
    xJ = J - 1;
    T = calc_T(P);
    tT = T;
    for (i -= DIM, mask >>= 1; i >= 0; i -= DIM, mask >>= 1)
    {
        for (j = A = 0; j < DIM; j++)
            if (pt.hcode[j] & mask)
                A |= g_mask[j];
        W ^= tT;
        tS = A ^ W;
        S = calc_tS_tT(xJ, tS);
        P = calc_P2(S);
        /* add in DIM bits to hcode */
        element = i / ORDER;
        if (i % ORDER > ORDER - DIM)
        {
            h.hcode[element] |= P << i % ORDER;
            h.hcode[element + 1] |= P >> ORDER - i % ORDER;
        }
        else
            h.hcode[element] |= P << i - element * ORDER;
        if (i > 0)
        {
            T = calc_T(P);
            tT = calc_tS_tT(xJ, T);
            J = calc_J(P);
            xJ += J - 1;
        }
    }
    return h;
}
