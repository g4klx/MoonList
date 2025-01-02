/*==============================================================================
 * sky0.c - astronomical coordinate conversion for NREL Sun Position Algorithm
 *
 * Author:  David Hoadley
 *
 * Description: (see sky0.h)
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

/* Local and project includes */
#include "sky0.h"

#include "astron.h"
#include "general.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       // For use by REQUIRE() - assertions.

/*      Convert from units of 0.1 milliarcsec to radians */
#define MILLIARCSECx10_TO_RAD  (PI / (180.0 * 3600.0 * 10000.0))

/*      Nutation constants from NREL SPA algorithm */
#define Y_COUNT 63
enum {TERM_X0, TERM_X1, TERM_X2, TERM_X3, TERM_X4, TERM_X_COUNT};
enum {TERM_PSI_A, TERM_PSI_B, TERM_EPS_C, TERM_EPS_D, TERM_PE_COUNT};
#define TERM_Y_COUNT TERM_X_COUNT


/*
 * Prototypes for local functions (not called from other modules)
 */

/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */
/*      Data tables from NREL SPA algorithm */
/*          Periodic Terms for the nutation in longitude and obliquity */
LOCAL const int Y_TERMS[Y_COUNT][TERM_Y_COUNT]=
{
    {0,0,0,0,1},
    {-2,0,0,2,2},
    {0,0,0,2,2},
    {0,0,0,0,2},
    {0,1,0,0,0},
    {0,0,1,0,0},
    {-2,1,0,2,2},
    {0,0,0,2,1},
    {0,0,1,2,2},
    {-2,-1,0,2,2},
    {-2,0,1,0,0},
    {-2,0,0,2,1},
    {0,0,-1,2,2},
    {2,0,0,0,0},
    {0,0,1,0,1},
    {2,0,-1,2,2},
    {0,0,-1,0,1},
    {0,0,1,2,1},
    {-2,0,2,0,0},
    {0,0,-2,2,1},
    {2,0,0,2,2},
    {0,0,2,2,2},
    {0,0,2,0,0},
    {-2,0,1,2,2},
    {0,0,0,2,0},
    {-2,0,0,2,0},
    {0,0,-1,2,1},
    {0,2,0,0,0},
    {2,0,-1,0,1},
    {-2,2,0,2,2},
    {0,1,0,0,1},
    {-2,0,1,0,1},
    {0,-1,0,0,1},
    {0,0,2,-2,0},
    {2,0,-1,2,1},
    {2,0,1,2,2},
    {0,1,0,2,2},
    {-2,1,1,0,0},
    {0,-1,0,2,2},
    {2,0,0,2,1},
    {2,0,1,0,0},
    {-2,0,2,2,2},
    {-2,0,1,2,1},
    {2,0,-2,0,1},
    {2,0,0,0,1},
    {0,-1,1,0,0},
    {-2,-1,0,2,1},
    {-2,0,0,0,1},
    {0,0,2,2,1},
    {-2,0,2,0,1},
    {-2,1,0,2,1},
    {0,0,1,-2,0},
    {-1,0,1,0,0},
    {-2,1,0,0,0},
    {1,0,0,0,0},
    {0,0,1,2,0},
    {0,0,-2,2,2},
    {-1,-1,1,0,0},
    {0,1,1,0,0},
    {0,-1,1,2,2},
    {2,-1,-1,2,2},
    {0,0,3,2,2},
    {2,-1,0,2,2},
};

LOCAL const double PE_TERMS[Y_COUNT][TERM_PE_COUNT]={
    {-171996,-174.2,92025,8.9},
    {-13187,-1.6,5736,-3.1},
    {-2274,-0.2,977,-0.5},
    {2062,0.2,-895,0.5},
    {1426,-3.4,54,-0.1},
    {712,0.1,-7,0},
    {-517,1.2,224,-0.6},
    {-386,-0.4,200,0},
    {-301,0,129,-0.1},
    {217,-0.5,-95,0.3},
    {-158,0,0,0},
    {129,0.1,-70,0},
    {123,0,-53,0},
    {63,0,0,0},
    {63,0.1,-33,0},
    {-59,0,26,0},
    {-58,-0.1,32,0},
    {-51,0,27,0},
    {48,0,0,0},
    {46,0,-24,0},
    {-38,0,16,0},
    {-31,0,13,0},
    {29,0,0,0},
    {29,0,-12,0},
    {26,0,0,0},
    {-22,0,0,0},
    {21,0,-10,0},
    {17,-0.1,0,0},
    {16,0,-8,0},
    {-16,0.1,7,0},
    {-15,0,9,0},
    {-13,0,7,0},
    {-12,0,6,0},
    {11,0,0,0},
    {-10,0,5,0},
    {-8,0,3,0},
    {7,0,-3,0},
    {-7,0,0,0},
    {-7,0,3,0},
    {-7,0,3,0},
    {6,0,0,0},
    {6,0,-3,0},
    {6,0,-3,0},
    {-6,0,3,0},
    {-6,0,3,0},
    {5,0,0,0},
    {-5,0,3,0},
    {-5,0,3,0},
    {-5,0,3,0},
    {4,0,0,0},
    {4,0,0,0},
    {4,0,0,0},
    {-4,0,0,0},
    {-4,0,0,0},
    {-4,0,0,0},
    {3,0,0,0},
    {-3,0,0,0},
    {-3,0,0,0},
    {-3,0,0,0},
    {-3,0,0,0},
    {-3,0,0,0},
    {-3,0,0,0},
    {-3,0,0,0},
};


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
GLOBAL void sky0_nutationSpa(double t_cy, Sky0_Nut1980 *nut)
/*! Calculates the nutation in longitude and obliquity, according to the
    algorithm set out in the NREL SPA document. This is a simplified version
    of the IAU 1980 Nutation Theory. Calculates first the fundamental nutation
    arguments, and then a series of terms. There are 106 terms in the full
    series, but this routine uses only the largest 63 of those terms.
 \param[in]  t_cy   centuries since J2000.0, TT timescale
 \param[out] nut    field \a nut->dPsi_rad - Nutation in longitude Δψ (radian)\n
                    field \a nut->dEps_rad - Nutation in obliquity Δε (radian)

 \par References:
    Reda, I. and Andreas, A. (2003), "Solar Position Algorithm
    for Solar Radiation Applications", National Renewable Energy Laboratory,
    NREL publication no. NREL/TP-560-34302, section 3.4.
    However, the theory on which it is based, and the full sequence of terms,
    can be found in:\n
        Final report of the IAU Working Group on Nutation,
        chairman P.K.Seidelmann, 1980,
        published in _Celestial Mechanics_, Vol 27, pp79-106, 1982.\n
   also
        Kaplan,G.H. _USNO circular no. 163_, pp A3-6, 1981.\n
        _Supplement to the Astronomical Almanac_ 1984.

 \par When to call this function
    It is quite likely that you will not need to call this function directly. It
    is used in the Solar Position Algorithm and the Moon Position Algorithm, so
    if you call sun_nrelApparent() or sun_nrelTopocentric(), or
    call moon_nrelApparent() or moon_nrelTopocentric(), those routines will call
    this routine for you. Likewise, if you are tracking the Sun or Moon using
    the skyfast module, the call to skyfast_init() will call either
    sun_nrelApparent() or moon_nrelApparent(), and therefore call this routine
    for you.
 \par
    The values calculated by this routine change only slowly. So if you are
    calling it yourself, you can call it infrequently. Intervals of up to an
    hour between calls will not introduce much error.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /* Fundamental Nutation arguments at date */
    double l;   // l  - mean anomaly of the Moon (radian) (usually called M?)
    double lp;  // l' - mean anomaly of the Sun (Earth) (radian) (M'?)
    double om;  /* Ω  - longitude of the ascending node of the Moon's mean orbit
                        on the ecliptic, measured from the mean equinox of date
                        (radian) */
    double d;   // D  - mean elongation of the Moon from the Sun (radian)
    double f;   // F  - Mean longitude of the Moon (L) minus mean longitude of
                //      the Moon's node (radian), = L - Ω

    double psiSum_masx10;       // Nutation in longitude (units - 0.1 mas)
    double epsSum_masx10;       // Nutation in obliquity (units - 0.1 mas)
    int i;
    double a_rad;               // angle - summation of args (radian)

    REQUIRE_NOT_NULL(nut);

    // Calculate FUNDAMENTAL ARGUMENTS in the FK5 reference system

    // Mean elongation of the moon from the sun (called X0 in NREL SPA)
    d = degToRad(297.85036 + t_cy * (445267.11148
                                     + t_cy * (-0.0019142
                                               + t_cy * (1.0/189474.0))));

    // Solar Mean Anomaly = Mean longitude of the sun minus mean longitude of
    // the sun's perigee (called X1 in NREL SPA)
    lp = degToRad(357.52772 + t_cy * (35999.05034
                                      + t_cy * (-0.0001603
                                                + t_cy * (-1.0/300000.0))));

    // Lunar Mean Anomaly = Mean longitude of the moon minus mean longitude of
    // the moon 's perigee (called X2 in NREL SPA)
    l = degToRad(134.96298 + t_cy * (477198.867398
                                     + t_cy * (0.0086972
                                               + t_cy * (1.0/56250.0))));

    // Mean longitude of the moon minus mean longitude of the moon's node
    // (called X3 in NREL SPA and (mistakenly, I think) called the Moon's
    // Argument of Latitude)
    f = degToRad(93.27191 + t_cy * (483202.017538
                                    + t_cy * (-0.0036825
                                              + t_cy * (1.0/327270.0))));

    // Longitude of the mean ascending node of the lunar orbit on the
    // ecliptic, measured from the mean equinox of date (X4 in NREL SPA)
    om = degToRad(125.04452 + t_cy * (-1934.136261
                                     + t_cy * (0.0020708
                                               + t_cy * (1.0/450000.0))));
#if 0
    printf("Nutation Args: d (X0) = %f°, lp (X1) = %f°, l (X2) = %f°\n"
           "               f (X3) = %f°, om (X4) = %f°\n",
           radToDeg(d), radToDeg(lp), radToDeg(l), radToDeg(f), radToDeg(om));
#endif

    // Multiply through the table of nutation co-efficients and add up all the
    // terms.
    psiSum_masx10 = 0.0;
    epsSum_masx10 = 0.0;
    for (i = 0; i < Y_COUNT; i++) {
        a_rad = d * Y_TERMS[i][0] + lp * Y_TERMS[i][1] + l * Y_TERMS[i][2]
                              + f * Y_TERMS[i][3] + om * Y_TERMS[i][4];
        psiSum_masx10 += sin(a_rad) * (PE_TERMS[i][TERM_PSI_A]
                                        + t_cy * PE_TERMS[i][TERM_PSI_B]);
        epsSum_masx10 += cos(a_rad) * (PE_TERMS[i][TERM_EPS_C]
                                        + t_cy * PE_TERMS[i][TERM_EPS_D]);
    }

    // Convert results to radians
    nut->dPsi_rad = psiSum_masx10 * MILLIARCSECx10_TO_RAD;
    nut->dEps_rad = epsSum_masx10 * MILLIARCSECx10_TO_RAD;
}



GLOBAL void sky0_epsilonSpa(double t_cy, Sky0_Nut1980 *nut)
/*! Calculate the obliquity of the ecliptic and the equation of the equinoxes
 \param[in]     t_cy    centuries since J2000.0, TT timescale
 \param[in,out] nut     [in]  field \a nut->dPsi_rad - Nutation in longitude Δψ,
                              as returned by function sky0_nutationSpa()
                              (radian)\n
                        [in]  field \a nut->dEps_rad - Nutation in obliquity Δε,
                              as returned by function sky0_nutationSpa()
                              (radian)\n
                        [out] field \a nut->eps0_rad - Mean obliquity of the
                              ecliptic ε0 (radian)\n
                        [out] field \a nut->eqEq_rad - Equation of the equinoxes
                              = Δψ * cos(ε0 + Δε) (radian) Note: not seconds

 \par Reference
    Reda, I. and Andreas, A. (2003), "Solar Position Algorithm
    for Solar Radiation Applications", National Renewable Energy Laboratory,
    NREL publication no. NREL/TP-560-34302, section 3.5.1

 \par When to call this function
    It is quite likely that you will not need to call this function directly. It
    is used in the Solar Position Algorithm and the Moon Position Algorithm, so
    if you call sun_nrelApparent() or sun_nrelTopocentric(), or
    call moon_nrelApparent() or moon_nrelTopocentric(), those routines will call
    this routine for you. Likewise, if you are tracking the Sun or Moon using
    the skyfast module, the call to skyfast_init() will call either
    sun_nrelApparent() or moon_nrelApparent(), and therefore call this routine
    for you.
 \par
    The values calculated by this routine change only slowly. So if you are
    calling it yourself, you can call it infrequently. Intervals of up to an
    hour between calls will not introduce much error.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double eps0_as;
    double u = t_cy/100.0;

    REQUIRE_NOT_NULL(nut);

    eps0_as = 84381.448
              + u * (-4680.93
                   + u * (-1.55
                        + u * (1999.25
                             + u * (-51.38
                                  + u * (-249.67
                                       + u * (-39.05
                                            + u * (7.12
                                                 + u * (27.87
                                                      + u * (5.79
                                                           + u * 2.45)))))))));
    nut->eps0_rad = arcsecToRad(eps0_as);
    // Calculate the Equation of the Equinoxes in radian (not in seconds, which
    // is more usually done)
   nut->eqEq_rad = nut->dPsi_rad * cos(nut->eps0_rad + nut->dEps_rad);
}



GLOBAL double sky0_gmSiderealTimeSpa(double du)
/*! Calculate the Greenwich mean sidereal time using the algorithm from the NREL
    SPA document. This is basically the IAU 1982 algorithm, but the various
    constants have been scaled in degrees.
 \returns    Greenwich Mean Sidereal Time (radian)
 \param[in]  du  days since J2000.0, UT1 timescale

 \par Reference
    Reda, I. and Andreas, A. (2003), "Solar Position Algorithm
    for Solar Radiation Applications", National Renewable Energy Laboratory,
    NREL publication no. NREL/TP-560-34302, section 3.8.1 & 3.8.2

 \par When to call this function
    If you are tracking a celestial object, you need not call this function
    directly. So long as you call sky0_appToTirs() every time around your
    control loop, it will call this function for you.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    static const double B0 = 280.46061837 * DEG2RAD;        // 67310.54841 s
    static const double B1 = 360.98564736629 * DEG2RAD;     // 86636.55536791 s
    static const double B2 = 0.000387933 * DEG2RAD;         // 0.093104 s
    static const double B3 = (-1.0 / 38710000.0) * DEG2RAD; // -6.2e-6 s
#if 0
    static const double B0d = 4.89496121274; // 3579;
    static const double B1d = 6.300388098984957;
    static const double B2d = 6.77070812713916E-06;
    static const double B3d = -4.50872966157151E-10;
#endif

    double gmst_rad;        // Greenwich Mean Sidereal Time (radian)
    double tu;              // Julian centuries since J2000.0, UT1 timescale

    tu = du / JUL_CENT;
    gmst_rad = B0 + B1 * du + (tu * tu * (B2 + (tu * B3)));

    return normalize(gmst_rad, TWOPI);
}



GLOBAL void sky0_appToTirs(const V3D_Vector *appV,
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
                     sky0_epsilonSpa() in the \a eqEq_rad field of the
                     Sky0_Nut1980 struct.

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
    gast_rad = sky0_gmSiderealTimeSpa(j2kUT1_d) + eqEq_rad;

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


