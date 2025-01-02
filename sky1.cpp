/*==============================================================================
 * sky1.c - astronomical coordinate conversion routines, IAU 1980
 *
 * Author:  David Hoadley
 *
 * Description: (see sky1.h)
 *
 * Copyright (c) 2020, David Hoadley <vcrumble@westnet.com.au>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 *==============================================================================
 */
/*------------------------------------------------------------------------------
 * Notes:
 *      Character set: UTF-8. (Non ASCII characters appear in this file)
 *----------------------------------------------------------------------------*/

/* ANSI includes etc. */
#include "instead-of-math.h"
#include <stdlib.h>

/* Local and project includes */
#include "sky1.h"

#include "astron.h"
#include "general.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                        // For use by REQUIRE() - assertions.

/*      Convert from units of 0.1 milliarcsec to radians */
#define MILLIARCSECx10_TO_RAD   (PI / (180.0 * 3600.0 * 10000.0))

#define NUM_TERMS               106

/*      Coefficients of fundamental arguments */
typedef struct {
    int cd;         // coefficient of d
    int clp;        // coefficient of l'
    int cl;         // coefficient of l
    int cf;         // coefficient of f
    int com;        // coefficient of omega
} Fcoeffs;

/*      Coefficients of dPsi and dEps accumulations */
typedef struct {
    double ps1;
    double ps2;
    double ep1;
    double ep2;
} PEcoeffs;

/*
 * Prototypes for local functions (not called from other modules)
 */

/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */
/*      Periodic Terms for the nutation in longitude and obliquity, sorted
        roughly in descending order of PE coefficient magnitude.
        (Terms in each of the following two tables MUST be in the same order) */
LOCAL const Fcoeffs coeffs[NUM_TERMS] =
{ //   d  l' l  F  Om
    {  0, 0, 0, 0, 1 },         /*   1 */
    { -2, 0, 0, 2, 2 },         /*   9 */
    {  0, 0, 0, 2, 2 },         /*  31 */
    {  0, 0, 0, 0, 2 },         /*   2 */
    {  0, 1, 0, 0, 0 },         /*  10 */
    {  0, 0, 1, 0, 0 },         /*  32 */
    { -2, 1, 0, 2, 2 },         /*  11 */
    {  0, 0, 0, 2, 1 },         /*  33 */
    {  0, 0, 1, 2, 2 },         /*  34 */
    { -2, -1, 0, 2, 2 },        /*  12 */
    { -2, 0, 1, 0, 0 },         /*  35 */
    { -2, 0, 0, 2, 1 },         /*  13 */
    {  0, 0, -1, 2, 2 },        /*  36 */
    {  2, 0, 0, 0, 0 },         /*  37 */
    {  0, 0, 1, 0, 1 },         /*  38 */
    {  2, 0, -1, 2, 2 },        /*  40 */
    {  0, 0, -1, 0, 1 },        /*  39 */
    {  0, 0, 1, 2, 1 },         /*  41 */
    { -2, 0, 2, 0, 0 },         /*  14 */
    {  0, 0, -2, 2, 1 },        /*   3 */
    {  2, 0, 0, 2, 2 },         /*  42 */
    {  0, 0, 2, 2, 2 },         /*  45 */
    {  0, 0, 2, 0, 0 },         /*  43 */
    { -2, 0, 1, 2, 2 },         /*  44 */
    {  0, 0, 0, 2, 0 },         /*  46 */
    { -2, 0, 0, 2, 0 },         /*  15 */
    {  0, 0, -1, 2, 1 },        /*  47 */
    {  0, 2, 0, 0, 0 },         /*  16 */
    {  2, 0, -1, 0, 1 },        /*  48 */
    { -2, 2, 0, 2, 2 },         /*  18 */
    {  0, 1, 0, 0, 1 },         /*  17 */
    { -2, 0, 1, 0, 1 },         /*  49 */
    {  0, -1, 0, 0, 1 },        /*  19 */
    {  0, 0, 2, -2, 0 },        /*   4 */
    {  2, 0, -1, 2, 1 },        /*  50 */

    {  2, 0, 1, 2, 2 },         /*  54 */
    { -2, 1, 1, 0, 0 },         /*  51 */
    {  0, 1, 0, 2, 2 },         /*  52 */
    {  0, -1, 0, 2, 2 },        /*  53 */
    {  2, 0, 0, 2, 1 },         /*  58 */

    {  2, 0, -2, 0, 1 },        /*  20 */
    {  2, 0, 1, 0, 0 },         /*  55 */
    { -2, 0, 2, 2, 2 },         /*  56 */
    {  2, 0, 0, 0, 1 },         /*  57 */
    { -2, 0, 1, 2, 1 },         /*  59 */

    { -2, -1, 0, 2, 1 },        /*  21 */
    { -2, 0, 0, 0, 1 },         /*  60 */
    {  0, -1, 1, 0, 0 },        /*  61 */
    {  0, 0, 2, 2, 1 },         /*  62 */

    { -2, 0, 2, 0, 1 },         /*  22 */
    { -2, 1, 0, 2, 1 },         /*  23 */
    { -1, 0, 1, 0, 0 },         /*  24 */
    { -2, 1, 0, 0, 0 },         /*  63 */
    {  0, 0, 1, -2, 0 },        /*  64 */
    {  1, 0, 0, 0, 0 },         /*  65 */

    {  0, 0, -2, 2, 2 },        /*   5 */
    { -1, -1, 1, 0, 0 },        /*   6 */
    {  0, 1, 1, 0, 0 },         /*  66 */
    {  0, 0, 1, 2, 0 },         /*  67 */
    {  0, -1, 1, 2, 2 },        /*  68 */
    {  2, -1, -1, 2, 2 },       /*  69 */
    {  0, 0, 3, 2, 2 },         /*  71 */
    {  2, -1, 0, 2, 2 },        /*  72 */

    { -2, -2, 0, 2, 1 },        /*   7 */
    {  0, 0, -2, 0, 1 },        /*  70 */
    {  0, 1, 1, 2, 2 },         /*  73 */
    { -2, 0, -1, 2, 1 },        /*  74 */
    {  0, 0, 2, 0, 1 },         /*  75 */
    {  0, 0, 1, 0, 2 },         /*  76 */
    {  0, 0, 3, 0, 0, },        /*  77 */
    {  1, 0, 0, 2, 2, },        /*  78 */
    {  4, 0, -1, 2, 2 },        /*  82 */

    {  0, 0, 2, -2, 1 },        /*   8 */
    { -2, 1, 2, 0, 0 },         /*  25 */
    {  2, 0, 0, -2, 1 },        /*  26 */
    {  2, 1, 0, -2, 0 },        /*  27 */
    {  0, 1, 0, 0, 2 },         /*  28 */
    {  1, 0, -1, 0, 1 },        /*  29 */
    { -2, 1, 0, 2, 0 },         /*  30 */

    {  0, 0, -1, 0, 2 },        /*  79 */
    { -4, 0, 1, 0, 0 },         /*  80 */
    {  2, 0, -2, 2, 2 },        /*  81 */
    { -4, 0, 2, 0, 0 },         /*  83 */
    { -2, 1, 1, 2, 2 },         /*  84 */
    {  2, 0, 1, 2, 1 },         /*  85 */
    {  4, 0, -2, 2, 2 },        /*  86 */
    {  0, 0, -1, 4, 2 },        /*  87 */
    { -2, -1, 1, 0, 0 },        /*  88 */
    { -2, 0, 2, 2, 1 },         /*  89 */
    {  2, 0, 2, 2, 2, },        /*  90 */
    {  2, 0, 1, 0, 1 },         /*  91 */
    { -2, 0, 0, 4, 2 },         /*  92 */
    { -2, 0, 3, 2, 2 },         /*  93 */
    { -2, 0, 1, 2, 0 },         /*  94 */
    {  0, 1, 0, 2, 1 },         /*  95 */
    {  2, -1, -1, 0, 1 },       /*  96 */
    {  0, 0, 0, -2, 1 },        /*  97 */
    { -1, 0, 0, 2, 2 },         /*  98 */
    {  2, 1, 0, 0, 0 },         /*  99 */
    { -2, 0, 1, -2, 0 },        /* 100 */
    {  0, -1, 0, 2, 1 },        /* 101 */
    { -2, 1, 1, 0, 1 },         /* 102 */
    {  2, 0, 1, -2, 0 },        /* 103 */
    {  2, 0, 2, 0, 0 },         /* 104 */
    {  4, 0, 0, 2, 2 },         /* 105 */
    {  1, 1, 0, 0, 0 },         /* 106 */
};

LOCAL const PEcoeffs pec[NUM_TERMS] = {
    {-171996.0, -174.2, 92025.0, 8.9 }, /*   1 */
    {-13187.0, -1.6, 5736.0, -3.1 },    /*   9 */
    {-2274.0, -0.2, 977.0, -0.5 },      /*  31 */
    { 2062.0,  0.2, -895.0, 0.5 },      /*   2 */
    { 1426.0, -3.4, 54.0, -0.1 },       /*  10 */
    {  712.0,  0.1, -7.0,  0.0 },       /*  32 */
    { -517.0, 1.2, 224.0, -0.6 },       /*  11 */
    { -386.0, -0.4, 200.0, 0.0 },       /*  33 */
    { -301.0, 0.0, 129.0, -0.1 },       /*  34 */
    {  217.0, -0.5, -95.0, 0.3 },       /*  12 */
    { -158.0, 0.0, -1.0, 0.0 },         /*  35 */
    {  129.0, 0.1, -70.0, 0.0 },        /*  13 */
    {  123.0, 0.0, -53.0, 0.0 },        /*  36 */
    {  63.0, 0.0, -2.0, 0.0 },          /*  37 */
    {  63.0, 0.1, -33.0, 0.0 },         /*  38 */
    { -59.0, 0.0, 26.0, 0.0 },          /*  40 */
    { -58.0, -0.1, 32.0, 0.0 },         /*  39 */
    { -51.0, 0.0, 27.0, 0.0 },          /*  41 */
    {  48.0, 0.0, 1.0, 0.0 },           /*  14 */
    {  46.0, 0.0, -24.0, 0.0 },         /*   3 */
    { -38.0, 0.0, 16.0, 0.0 },          /*  42 */
    { -31.0, 0.0, 13.0, 0.0 },          /*  45 */
    {  29.0, 0.0, -1.0, 0.0 },          /*  43 */
    {  29.0, 0.0, -12.0, 0.0 },         /*  44 */
    {  26.0, 0.0, -1.0, 0.0 },          /*  46 */
    { -22.0, 0.0, 0.0, 0.0 },           /*  15 */
    {  21.0, 0.0, -10.0, 0.0 },         /*  47 */
    {  17.0, -0.1, 0.0, 0.0 },          /*  16 */
    {  16.0, 0.0, -8.0, 0.0 },          /*  48 */
    { -16.0, 0.1, 7.0, 0.0 },           /*  18 */
    { -15.0, 0.0, 9.0, 0.0 },           /*  17 */
    { -13.0, 0.0, 7.0, 0.0 },           /*  49 */
    { -12.0, 0.0, 6.0, 0.0 },           /*  19 */
    {  11.0, 0.0, 0.0, 0.0 },           /*   4 */
    { -10.0, 0.0, 5.0, 0.0 },           /*  50 */
                                                /*  35 terms to here */
    { -8.0, 0.0, 3.0, 0.0 },            /*  54 */
    { -7.0, 0.0, 0.0, 0.0 },            /*  51 */
    {  7.0, 0.0, -3.0, 0.0 },           /*  52 */
    { -7.0, 0.0, 3.0, 0.0 },            /*  53 */
    { -7.0, 0.0, 3.0, 0.0 },            /*  58 */
                                                /*  40 terms to here */
    { -6.0, 0.0, 3.0, 0.0 },            /*  20 */
    {  6.0, 0.0, 0.0, 0.0 },            /*  55 */
    {  6.0, 0.0, -3.0, 0.0 },           /*  56 */
    { -6.0, 0.0, 3.0, 0.0 },            /*  57 */
    {  6.0, 0.0, -3.0, 0.0 },           /*  59 */
                                                /*  45 terms to here */
    { -5.0, 0.0, 3.0, 0.0 },            /*  21 */
    { -5.0, 0.0, 3.0, 0.0 },            /*  60 */
    {  5.0, 0.0, 0.0, 0.0 },            /*  61 */
    { -5.0, 0.0, 3.0, 0.0 },            /*  62 */
                                                /*  49 terms to here */
    {  4.0, 0.0, -2.0, 0.0 },           /*  22 */
    {  4.0, 0.0, -2.0, 0.0 },           /*  23 */
    { -4.0, 0.0, 0.0, 0.0 },            /*  24 */
    { -4.0, 0.0, 0.0, 0.0 },            /*  63 */
    {  4.0, 0.0, 0.0, 0.0 },            /*  64 */
    { -4.0, 0.0, 0.0, 0.0 },            /*  65 */
                                                /*  55 terms to here */
    { -3.0, 0.0, 1.0, 0.0 },            /*   5 */
    { -3.0, 0.0, 0.0, 0.0 },            /*   6 */
    { -3.0, 0.0, 0.0, 0.0 },            /*  66 */
    {  3.0, 0.0, 0.0, 0.0 },            /*  67 */
    { -3.0, 0.0, 1.0, 0.0 },            /*  68 */
    { -3.0, 0.0, 1.0, 0.0 },            /*  69 */
    { -3.0, 0.0, 1.0, 0.0 },            /*  71 */
    { -3.0, 0.0, 1.0, 0.0 },            /*  72 */
                                                /*  63 terms to here */
    { -2.0, 0.0,  1.0, 0.0 },           /*   7 */
    { -2.0, 0.0,  1.0, 0.0 },           /*  70 */
    {  2.0, 0.0, -1.0, 0.0 },           /*  73 */
    { -2.0, 0.0,  1.0, 0.0 },           /*  74 */
    {  2.0, 0.0, -1.0, 0.0 },           /*  75 */
    { -2.0, 0.0,  1.0, 0.0 },           /*  76 */
    {  2.0, 0.0,  0.0, 0.0 },           /*  77 */
    {  2.0, 0.0, -1.0, 0.0 },           /*  78 */
    { -2.0, 0.0,  1.0, 0.0 },           /*  82 */
                                                /*  72 terms to here */
    {  1.0, 0.0,  0.0, 0.0 },           /*   8 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  25 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  26 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  27 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  28 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  29 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  30 */

    {  1.0, 0.0, -1.0, 0.0 },           /*  79 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  80 */
    {  1.0, 0.0, -1.0, 0.0 },           /*  81 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  83 */
    {  1.0, 0.0, -1.0, 0.0 },           /*  84 */
    { -1.0, 0.0,  1.0, 0.0 },           /*  85 */
    { -1.0, 0.0,  1.0, 0.0 },           /*  86 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  87 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  88 */
    {  1.0, 0.0, -1.0, 0.0 },           /*  89 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  90 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  91 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  92 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  93 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  94 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  95 */
    {  1.0, 0.0,  0.0, 0.0 },           /*  96 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  97 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  98 */
    { -1.0, 0.0,  0.0, 0.0 },           /*  99 */
    { -1.0, 0.0,  0.0, 0.0 },           /* 100 */
    { -1.0, 0.0,  0.0, 0.0 },           /* 101 */
    { -1.0, 0.0,  0.0, 0.0 },           /* 102 */
    { -1.0, 0.0,  0.0, 0.0 },           /* 103 */
    {  1.0, 0.0,  0.0, 0.0 },           /* 104 */
    { -1.0, 0.0,  0.0, 0.0 },           /* 105 */
    {  1.0, 0.0,  0.0, 0.0 }            /* 106 */
};

LOCAL const int nTerms[] = { NUM_TERMS, 72, 63, 49, 35 };

/*
 *==============================================================================
 *
 * Implementation
 *
 *==============================================================================
 *
 * Global functions callable by other modules
 *
 *------------------------------------------------------------------------------
 */
GLOBAL void sky1_frameBiasFK5(V3D_Matrix *biasM)
/*! Create the frame bias matrix from the IAU 2000 precession-nutation model.
    This matrix is used to convert coordinates from the Geocentric Celestial
    Reference System (GCRS) to FK5 reference system. (Although anything to do
    with the ICRS/GCRS is related to the IAU2000+ precession-nutation theory,
    this routine is included here to allow ICRS star catalogue positions to be
    converted to apparent coordinates using the IAU1980 precession-nutation
    theory, which is quite good enough for tracking. It is not good enough for
    astrometry, but it will get us to better than within an arcsecond of our
    desired apparent position).
 \param[out] biasM  Frame bias matrix \b B

 \par Reference:
    _Astronomical Almanac_ 2007, page B28
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    const double eta0_as    = -19.9 / 1000.0;
    const double xi0_as     =  9.1 / 1000.0;
    const double dAlpha0_as = -22.9 / 1000.0;
    V3D_Matrix   cM, bM, aM;

    REQUIRE_NOT_NULL(biasM);

    /* First intermediate result. [cM] = [R2(xi0)] x [R3(dAlpha0)] */
    v3d_createRotationMatrix(&aM, Zaxis, arcsecToRad(dAlpha0_as));
    v3d_createRotationMatrix(&bM, Yaxis, arcsecToRad(xi0_as));
    v3d_multMxM(&cM, &bM, &aM);

    /* Re-use aM for rotation by eta0, obtain
     *      [biasM] = [R1(eta0)] x [R2(xi0)] x [R3(dAlpha0)] */
    v3d_createRotationMatrix(&aM, Xaxis, arcsecToRad(-eta0_as));
    v3d_multMxM(biasM, &aM, &cM);
}



GLOBAL void sky1_precessionIAU1976(double t0_cy,
                                   double t1_cy,
                                   Sky1_Prec1976 *terms)
/*! This procedure calculates the equatorial precession parameters ζ, z, and θ
    which represent the rotation required to transform the FK5 equatorial
    reference system of Julian epoch t0 to that of Julian epoch t1, according to
    the IAU 1976 precession constants.
 \param[in]  t0_cy  Julian centuries since J2000.0 of initial epoch, TT
                    timescale
 \param[in]  t1_cy  Julian centuries since J2000.0 of final epoch, TT timescale
 \param[out] terms  The three precession angles ζ, z, θ.

 \par References:
    Lieske,J.H., _Astron. Astrophys_, Vol 73, pp282-284, 1979\n
    _Supplement to The Astronomical Almanac 1984_, HMNAO, ppS18-19

 \par When to call this function
    It is quite likely that you will not need to call this function directly.
    If you are tracking a celestial object, using star_getApparent() or
    star_getTopocentric() or calling planet_getApparent() or
    planet_getTopocentric(), they will call this routine for you.
    Likewise, if you are tracking the object using the skyfast module, the call
    to skyfast_init() will call either star_getApparent() or
    planet_getApparent(), and therefore call this routine for you.
 \par
    The values calculated by this routine change only slowly. So if you are
    calling it yourself, you can call it infrequently. Intervals of up to an
    hour between calls will not introduce much error. This function is also
    called by sky1_createNPmatrix().
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double t0Sq;
    double t, tSq, tCu;
    double temp_as;         // angle (arcseconds)

    REQUIRE_NOT_NULL(terms);

    /* Calculate interval between initial epoch T0 and final epoch T1 in Julian
       centuries */
    t = (t1_cy - t0_cy);

    /* Calculate squares and cubes of intervals */
    tSq = t * t;
    tCu = tSq * t;
    t0Sq = t0_cy * t0_cy;

    /* Calculate equatorial precession parameters */
    temp_as = (2306.2181 + 1.39656 * t0_cy - 0.000139 * t0Sq) * t
               + (0.30188 - 0.000344 * t0_cy) * tSq + 0.017998 * tCu;
    terms->zeta_rad = arcsecToRad(temp_as);

    temp_as = (2306.2181 + 1.39656 * t0_cy - 0.000139 * t0Sq) * t
               + (1.09468 + 0.000066 * t0_cy) * tSq + 0.018203 * tCu;
    terms->zed_rad = arcsecToRad(temp_as);

    temp_as = (2004.3109 - 0.85330 * t0_cy - 0.000217 * t0Sq) * t
               - (0.42665 + 0.000217 * t0_cy) * tSq - 0.041833 * tCu;
    terms->theta_rad = arcsecToRad(temp_as);
}



GLOBAL void sky1_createPrec1976Matrix(const Sky1_Prec1976 *terms,
                                      V3D_Matrix *precM)
/*! This routine calculates the precession matrix, based on angles ζ, z, and θ.
    \a precM is the combined orthogonal rotation matrix \b P, required for
    rigorous precession transformations using rectangular coordinates and matrix
    methods:
        V1 = P * V0
   It appears to be calculated from\n
        \b P = R3(-z) × R2(θ) × R3(-ζ)
 \param[in]  terms  The three precession angles ζ, z, θ, as returned by
                        sky1_precessionIAU1976()
 \param[out] precM  Precession rotation matrix

 \par References:
    Lieske,J.H., _Astron. Astrophys_, Vol 73, pp282-284, 1979\n
    _Supplement to The Astronomical Almanac 1984_, HMNAO, ppS18-19

 \par When to call this function
    Call this function after new values of ζ, z and θ have been calculated
    by routine sky1_precessionIAU1976(). So, again, it need only be called
    infrequently. But as mentioned for that function, it is quite likely that
    you will not need to call this function directly. Any of the functions
    sky1_createNPmatrix(), star_getApparent(), star_getTopocentric(),
    planet_getApparent() or planet_getTopocentric() will call this
    routine for you.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Matrix zedM;
    V3D_Matrix thetaM;
    V3D_Matrix zetaM;
    V3D_Matrix tempM;

    REQUIRE_NOT_NULL(terms);
    REQUIRE_NOT_NULL(precM);

    /* Calculate the precession matrix as multiplication of rotation matrices */
    v3d_createRotationMatrix(&zetaM,  Zaxis, -terms->zeta_rad);
    v3d_createRotationMatrix(&thetaM, Yaxis, terms->theta_rad);
    v3d_createRotationMatrix(&zedM,   Zaxis, -terms->zed_rad);
    v3d_multMxM(precM, &zedM, v3d_multMxM(&tempM, &thetaM, &zetaM));
}



GLOBAL void sky1_nutationIAU1980(double t_cy,
                                 int    precision,
                                 Sky1_Nut1980 *nut)
/*! Calculates the nutation in longitude and obliquity, according to the IAU
    1980 Nutation Theory. Calculates first the fundamental nutation arguments,
    and then a series of terms. There are 106 terms in the series, but fewer
    may be selected for faster execution.
 \param[in]  t_cy       Julian centuries since J2000.0, TT timescale
 \param[in]  precision  How much precision do you want?
                          Valid range [0, 4]. Values outside this range will be
                          clamped to the range.
                       - 0 = full precision, use full 106-term series
                       - 1 = ignore terms < 0.2 milliarcseconds. 72-term series
                       - 2 = ignore terms < 0.3 milliarcseconds. 63-term series
                       - 3 = ignore terms < 0.5 milliarcseconds. 49-term series
                       - 4 = ignore terms < 1.0 milliarcseconds. 35-term series
 \param[out] nut    field \a nut->dPsi_rad - Nutation in longitude Δψ (radian)\n
                    field \a nut->dEps_rad - Nutation in obliquity Δε (radian)

 \par References:
        Final report of the IAU Working Group on Nutation,
        chairman P.K.Seidelmann, 1980,
        published in _Celestial Mechanics_, Vol 27, pp79-106, 1982.\n
        Kaplan,G.H. _USNO circular no. 163_, pp A3-6, 1981.\n
        _Supplement to the Astronomical Almanac_ 1984.\n
        The form of the data tables was inspired by, and extended from,
        Reda, I. and Andreas, A. (2003), "Solar Position Algorithm
        for Solar Radiation Applications", National Renewable Energy Laboratory,
        NREL publication no. NREL/TP-560-34302.

 \note
  - With precision = 0, exact agreement with tabulated values in the
    _Astronomical Almanac_ 1987 pages B24 - B31 was found (to 3 decimal places,
    as tabulated in the Almanac). Likewise for the 1990 edition pages B24 - B31.
    But this IAU 1980 theory has been superseded by the more accurate IAU 2000
    theory. This is evident when comparing outputs with the tabulated values in
    the 2007 Almanac, pages B34 - B41 (to 4 decimal places). There are
    differences of several milliarcseconds.

  - With precision = 1, differences up to 1.0 milliarcsecond were seen compared
    to the precision = 0 values during testing.
  - With precision = 2, differences up to 1.2 milliarcseconds were seen
  - With precision = 3, differences up to 1.5 milliarcseconds were seen
  - With precision = 4, differences of up to 5 milliarcseconds were seen.

    For tracking purposes, these are still very small.
    The testing was not exhaustive, so take this as a rough guide only.

 \par When to call this function
    It is quite likely that you will not need to call this function directly.
    If you are tracking a celestial object, using star_getApparent() or
    star_getTopocentric() or calling planet_getApparent() or
    planet_getTopocentric(), they will call this routine for you.
    Likewise, if you are tracking the object using the skyfast module, the call
    to skyfast_init() will call either star_getApparent() or
    planet_getApparent(), and therefore call this routine for you.
 \par
    The values calculated by this routine change only slowly. So if you are
    calling it yourself, you can call it infrequently. Intervals of up to an
    hour between calls will not introduce much error. This function is also
    called by sky1_createNPmatrix().
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    const double Revs2Arcsec = 1296000.0;        // Revolutions to arc seconds

    /* Fundamental Nutation arguments at date */
    double l;   // l  - mean anomaly of the Moon (radian) (usually called M?)
    double lp;  // l' - mean anomaly of the Sun (Earth) (radian) (M'?)
    double om;  // Ω  - longitude of the ascending node of the Moon's mean orbit
                //      on the ecliptic, measured from the mean equinox of date
                //      (radian)
    double d;   // D  - mean elongation of the Moon from the Sun (radian)
    double f;   // F  - Mean longitude of the Moon (L) minus mean longitude of
                //      the Moon's node (radian), = L - Ω

    double psiSum_masx10;       // Nutation in longitude (units - 0.1 mas)
    double epsSum_masx10;       // Nutation in obliquity (units - 0.1 mas)
    int i;
    double a_rad;               // angle - summation of args (radian)

    REQUIRE_NOT_NULL(nut);

    if (precision < 0) { precision = 0; }
    if (precision >= ARRAY_SIZE(nTerms)) { precision = ARRAY_SIZE(nTerms) - 1; }

    // Calculate FUNDAMENTAL ARGUMENTS in the FK5 reference system

    // Mean elongation of the Moon from the Sun (called X0 in NREL SPA)
    d = arcsecToRad(1072261.307 + (1236.0 * Revs2Arcsec + 1105601.328
                                   + (-6.891 + 0.019 * t_cy) * t_cy) * t_cy);
    //+
    d = normalize(d, TWOPI);
    //-

    // Solar Mean Anomaly = Mean longitude of the Sun minus mean longitude of
    // the Sun's perigee (called X1 in NREL SPA)
    lp = arcsecToRad(1287099.804 + (99.0 * Revs2Arcsec + 1292581.224
                                    + (-0.577 - 0.012 * t_cy) * t_cy) * t_cy);
    //+
    lp = normalize(lp, TWOPI);
    //-

    // Lunar Mean Anomaly = Mean longitude of the Moon minus mean longitude of
    // the Moon 's perigee (called X2 in NREL SPA)
    l = arcsecToRad(485866.733 + (1325.0 * Revs2Arcsec + 715922.633
                                  + (31.31 + 0.064 * t_cy) * t_cy) * t_cy);
    //+
    l = normalize(l, TWOPI);
    //-

    // Mean longitude of the Moon minus mean longitude of the Moon's node
    // (called X3 in NREL SPA)
    f = arcsecToRad(335778.877 + (1342.0 * Revs2Arcsec + 295263.137
                                  + (-13.257 + 0.011 * t_cy) * t_cy) * t_cy);
    //+
    f = normalize(f, TWOPI);
    //-


    // Longitude of the mean ascending node of the lunar orbit on the
    // ecliptic, measured from the mean equinox of date (X4 in NREL SPA)
    om = arcsecToRad(450160.28 + (-5.0 * Revs2Arcsec - 482890.539
                                  + (7.455 + 0.008 * t_cy) * t_cy) * t_cy);
    //+
    //om = normalize(om, TWOPI) - TWOPI;
    //-


    // Multiply through the table of nutation co-efficients and add up all the
    // terms.
    psiSum_masx10 = 0.0;
    epsSum_masx10 = 0.0;
    for (i = nTerms[precision] - 1; i >= 0; i--) {
        a_rad = d * coeffs[i].cd + lp * coeffs[i].clp + l * coeffs[i].cl
                                  + f * coeffs[i].cf + om * coeffs[i].com;
        psiSum_masx10 += sin(a_rad) * (pec[i].ps1 + t_cy * pec[i].ps2);
        // No need to calculate cos if its coefficients are zero. And if ep1
        // is zero, so is ep2, so we don't need to test it.
        if (pec[i].ep1 != 0.0) {
            epsSum_masx10 += cos(a_rad) * (pec[i].ep1 + t_cy * pec[i].ep2);
        }
    }

    // Convert results to radians
    nut->dPsi_rad = psiSum_masx10 * MILLIARCSECx10_TO_RAD;
    nut->dEps_rad = epsSum_masx10 * MILLIARCSECx10_TO_RAD;
}



GLOBAL void sky1_epsilon1980(double t_cy, Sky1_Nut1980 *nut)
/*! Calculate the obliquity of the ecliptic and the equation of the equinoxes
 \param[in]     t_cy  Julian centuries since J2000.0, TT timescale
 \param[in,out] nut   [in]  field \a nut->dPsi_rad - Nutation in longitude Δψ,
                            as returned by function sky1_nutationIAU1980()
                            (radian)\n
                      [in]  field \a nut->dEps_rad - Nutation in obliquity Δε,
                            as returned by function sky1_nutationIAU1980()
                            (radian)\n
                      [out] field \a nut->eps0_rad - Mean obliquity of the
                            ecliptic ε0 (radian)\n
                      [out] field \a nut->eqEq_rad - Equation of the equinoxes
                              = Δψ * cos(ε0 + Δε) (radian)  Note: not seconds
 \par References
    _Astronomical Almanac_ 1987 page B18\n
    International Astronomical Union, Standards of Fundamental Astronomy
        software collection, release 2017-04-20. http://www.iausofa.org
        routine \c iauObl80()

 \par When to call this function
    Call this function after new values of Δψ, and Δε have been calculated
    by routine sky1_nutationIAU1980(). So, again, it need only be called
    infrequently. But as mentioned for that function, it is quite likely that
    you will not need to call this function directly. Any of the functions
    sky1_createNPmatrix(), star_getApparent(), star_getTopocentric(),
    planet_getApparent() or planet_getTopocentric() will call this
    routine for you.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double eps0_as;

    REQUIRE_NOT_NULL(nut);

    // Calculate Mean obliquity in radian
    eps0_as = 84381.448
              + (-46.8150 + (-0.00059 + 0.001813 * t_cy) * t_cy) * t_cy;
    nut->eps0_rad = arcsecToRad(eps0_as);

    // Calculate the Equation of the Equinoxes in radian (not in seconds, which
    // is more usually done)
    nut->eqEq_rad = nut->dPsi_rad * cos(nut->eps0_rad + nut->dEps_rad);
}



GLOBAL void sky1_createNut1980Matrix(const Sky1_Nut1980 *nut,
                                     V3D_Matrix *nutM)
/*! This routine calculates the Nutation matrix, using the nutation in
    longitude, the nutation in obliquity, and the mean obliquity of the equator.
    We calculate the Nutation matrix in a completely rigorous manner from the
    product of three rotation matrices:\n
        \b N = R1(-ε) × R3(-Δψ) × R1(ε0)
    where ε = ε0 + Δε   (true obliquity of the equator)
 \param[in]  nut   Angles Δψ, Δε and ε0, as returned by functions
                     sky1_nutationIAU1980() and sky1_epsilon1980()
 \param[out] nutM  Nutation rotation matrix  \b N

 \par References:
    _Astronomical Almanac_ 2007, page B31\n
     also
        Mueller, p75

 \par When to call this function
    Call it after you have called sky1_nutationIAU1980() and sky1_epsilon1980().
    As mentioned for the other functions above, this need only be done
    infrequently. But also as mentioned, it is quite likely that
    you will not need to call this function directly. Any of the functions
    sky1_createNPmatrix(), star_getApparent(), star_getTopocentric(),
    planet_getApparent() or planet_getTopocentric() will call this
    routine for you.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Matrix eps0M;
    V3D_Matrix psiM;
    V3D_Matrix epsM;
    V3D_Matrix tempM;

    REQUIRE_NOT_NULL(nut);
    REQUIRE_NOT_NULL(nutM);

    v3d_createRotationMatrix(&eps0M, Xaxis, nut->eps0_rad);
    v3d_createRotationMatrix(&psiM,  Zaxis, -nut->dPsi_rad);
    v3d_createRotationMatrix(&epsM, Xaxis, -(nut->eps0_rad + nut->dEps_rad));

    // Multiply the three matrices in the order epsM * psiM * eps0M
    v3d_multMxM(nutM, &epsM, v3d_multMxM(&tempM, &psiM, &eps0M));
}



GLOBAL void sky1_createNPmatrix(double t0_cy,
                                double t1_cy,
                                int    precision,
                                V3D_Matrix* npM)
/*! Call the various precession and nutation routines in this module to create
    a combined precession and nutation rotation matrix, suitable for converting
    a celestial object's coordinates from mean place at epoch \a t0_cy to
    apparent coordinates at epoch \a t1_cy.
 \param[in]  t0_cy      Julian centuries since J2000.0 of initial epoch,
                          TT timescale
 \param[in]  t1_cy      Julian centuries since J2000.0 of final epoch,
                          TT timescale
 \param[in]  precision  Precision specification to be passed to routine
                          sky1_nutationIAU1980(). See that routine for
                          explanation
 \param[out] npM        Combined precession and nutation matrix
                          \b NP = \b N × \b P

    This function calls routines sky1_precessionIAU1976(),
    sky1_createPrec1976Matrix(), sky1_nutationIAU1980(), sky1_epsilon1980()
    and sky1_createNut1980Matrix() to calculate the precession and nutation
    matrices.

 \par When to call this function
    Call this function if you are tracking non-solar system objects.
    The resulting matrix changes only very slowly, so this function needs to be
    called only once per hour.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Matrix nM, pM;
    Sky1_Prec1976 prec;
    Sky1_Nut1980  nut;

    REQUIRE_NOT_NULL(npM);

    sky1_precessionIAU1976(t0_cy, t1_cy, &prec);
    sky1_createPrec1976Matrix(&prec, &pM);
    sky1_nutationIAU1980(t1_cy, precision, &nut);
    sky1_epsilon1980(t1_cy, &nut);
    sky1_createNut1980Matrix(&nut, &nM);
    v3d_multMxM(npM, &nM, &pM);
}



GLOBAL double sky1_gmSiderealTimeIAU1982(double du)
/*! Calculate the Greenwich mean sidereal time.
 \returns    Greenwich Mean Sidereal Time (radian)
 \param[in]  du  Days since J2000.0, UT1 timescale

 \par Reference:
    _Astronomical Almanac_ 1987, page B6

 \par When to call this function
    If you are tracking a celestial object, you need not call this function
    directly. So long as you call sky1_appToTirs() every time around your
    control loop, it will call this function for you.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /*
     *  The expression given on page B6 of the Almanac, which is
     *     GMST = 24110.54841 + 8640184.812866 Tu + 0.093104 Tu^2 - 6.2e-6 Tu^3
     * is only defined for 0 h UT1 (midnight), meaning that it is only defined
     * for Du = 0.5, 1.5, 2.5 etc., (since Du is days since noon, 1st January
     * 2000). So we need to add a term which defines it for the rest of the day,
     * which means adding 2Pi radian per day. But this new term will add an
     *  offset of Pi at every midnight, so we need to take that out again.
     */

    // Time series coefficients for Greenwich Mean Sidereal Time
    const double B0 = secToRad(24110.54841) - PI;
    const double B1 = secToRad(8640184.812866);
    const double B2 = secToRad(0.093104);
    const double B3 = secToRad(-6.2e-6);

    double gmst_rad;        // Greenwich Mean Sidereal Time (radian)
    double tu;              // Julian centuries since J2000.0, UT1 timescale
    double dayfrac;         // fractional part of du

    tu = du / JUL_CENT;
    dayfrac = du - floor(du);

    gmst_rad = B0 + (tu * (B1 + (tu * (B2 + (tu * B3))))) + dayfrac * TWOPI;

    return normalize(gmst_rad, TWOPI);
}



GLOBAL void sky1_appToTirs(const V3D_Vector *appV,
                           double           j2kUT1_d,
                           double           eqEq_rad,
                           V3D_Vector *terInterV)
/*! Convert a position in geocentric apparent coordinates to geocentric
    coordinates in the Terrestrial Intermediate Reference System. This is the
    first stage of converting apparent coordinates to topocentric coordinates.
    The resulting vector depends upon the current rotational position of the
    Earth. (For the second stage, to obtain topocentric coordinates, call
    routine sky_siteTirsToTopo()).
 \param[in] appV     Position vector of apparent place
                     (unit vector in equatorial coordinates)
 \param[in] j2kUT1_d days since J2000.0, UT1 timescale, as returned by function
                     sky_updateTimes() in the \a j2kUT1_d field of the
                     Sky_Times struct.
 \param[in] eqEq_rad Equation of the equinoxes (radian), as returned by function
                     sky1_epsilon1980() in the \a eqEq_rad field of the
                     Sky1_Nut1980 struct.

 \param[out] terInterV  Position vector in Terrestrial Intermediate Ref System

 \par When to call this function
    When you have the position of a celestial object expressed in Apparent
    coordinates, use this function to convert it to Terrestrial coordinates at
    the relevant rotational position of the earth. (Effectively you are
    converting from Right Ascension and Declination to the negative of the
    Greenwich Hour Angle and Declination.)
    If you are running a control loop to enable continuous tracking
    of this object, you will need to call this function (once) every time around
    your control loop.
 \par
    Follow this function with a call to sky_siteTirsToTopo() to obtain the
    object's position in topocentric coordinates at the observing site.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double      gast_rad;   // Greenwich Apparent Sidereal Time (GAST)
    V3D_Matrix  earthRotM;  // rotation matrix for current GAST

    REQUIRE_NOT_NULL(appV);
    REQUIRE_NOT_NULL(terInterV);

    /* Get sidereal times */
    gast_rad = sky1_gmSiderealTimeIAU1982(j2kUT1_d) + eqEq_rad;

    /* Create earth rotation matrix from the current Apparent Sidereal Time */
    v3d_createRotationMatrix(&earthRotM, Zaxis, gast_rad);

    /* Convert apparent posn to posn in Terrestrial Intermediate Ref System */
    v3d_multMxV(terInterV, &earthRotM, appV);
}


/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules)
 *
 *------------------------------------------------------------------------------
 */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
