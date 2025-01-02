/*==============================================================================
 * sun.c - routines to calculate the Sun's position
 *
 * Author:  David Hoadley
 *
 * Description: (see sun.h)
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
 *      Character set: UTF-8. (Non-ASCII characters appear in this file)
 *----------------------------------------------------------------------------*/

/* ANSI includes etc. */
#include <float.h>
#include "instead-of-math.h"            /* for normalize() */

/* Local and project includes */
#include "sun.h"

#include "astron.h"
#include "general.h"
#include "vectors3d.h"

/*
 * Local #defines and typedefs 
 */
DEFINE_THIS_FILE;                       // For use by REQUIRE() - assertions.

/*      Constants from NREL SPA algorithm */
#define L_COUNT 6
#define B_COUNT 2
#define R_COUNT 5

typedef struct {
    double a;       // term A in (A * cos(B + C * t)). (called TERM_A in SPA)
    double cb;      // term B in (A * cos(B + C * t)). (called TERM_B in SPA)
    double cct;     // term C in (A * cos(B + C * t)). (called TERM_C in SPA)
} SunSeries;


/*
 * Prototypes for local functions (not called from other modules).
 */
LOCAL double sunLongitude(double t_ka);
LOCAL double sunLatitude(double t_ka);
LOCAL double sunDistance(double t_ka);
LOCAL double solarNoonApprox(double             noonGuess_d,
                             const Sky_DeltaTs  *deltas,
                             const Sky_SiteProp *site,
                             Sky_SiteHorizon *topo);
LOCAL double riseSetApprox(double             risesetGuess_d,
                           bool               getSunrise,
                           const Sky_DeltaTs  *deltas,
                           const Sky_SiteProp *site,
                           Sky_SiteHorizon *topo);

/*
 * Global variables accessible by other modules 
 */
/*      (none) */

/*
 * Local variables (not accessed by other modules)
 */
/*      (none) */


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
GLOBAL void sun_aaApparentApprox(double n,
                                 V3D_Vector *appV,
                                 double     *dist_au)
/*! This function calculates an approximate Sun position in apparent coordinates
    using the algorithm given in the _Astronomical Almanac_. The quoted
    accuracy is a precision of 0.01 degrees between 1950 and 2050.
 \param[in]  n        Days since J2000.0, UT1 timescale
 \param[out] appV     Position vector of Sun in apparent coordinates (unit
                        vector i.e. direction cosines)
 \param[out] dist_au  Geocentric distance of the Sun (Astronomical Units)
 
 \par When to call this routine
    Use this routine for a quick rough calculation of the Sun's position. If you
    need an accurate position, don't call this routine, call sun_nrelApparent()
    instead.

 \par Reference:
    The _Astronomical Almanac_, 2007, page C24
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double L_deg;           // Mean longitude, corrected for aberr (degrees)
    double g_deg;           // Mean anomaly (degrees)
    double lambda_deg;      // Ecliptic longitude (degrees)
    double epsilon_deg;     // Obliquity of ecliptic (degrees)

    REQUIRE_NOT_NULL(appV);
    REQUIRE_NOT_NULL(dist_au);

    L_deg = normalize(280.461 + (0.9856474 * n), 360.0);
    g_deg = normalize(357.529 + (0.9856003 * n), 360.0);
    lambda_deg = L_deg + 1.915 * sin(degToRad(g_deg))
                       + 0.02 * sin(degToRad(g_deg + g_deg));
    epsilon_deg = 23.439 - (0.0000004 * n);
    appV->a[0] = cos(degToRad(lambda_deg));
    appV->a[1] = cos(degToRad(epsilon_deg)) * sin(degToRad(lambda_deg));
    appV->a[2] = sin(degToRad(epsilon_deg)) * sin(degToRad(lambda_deg));
    
    *dist_au = 1.00014 - 0.01671 * cos(degToRad(g_deg))
                       - 0.00014 * cos(degToRad(g_deg + g_deg));
}



GLOBAL void sun_nrelApp2(double             t_cy,
                         const Sky0_Nut1980 *nut,
                         V3D_Vector *appV,
                         double     *dist_au)
/*! This function calculates the Sun position in apparent coordinates, using the
    NREL SPA algorithm (see reference). The quoted accuracy of this algorithm
    os 0.0003° (or about 1 arcsecond). It is much more computationally intensive
    than the approximate algorithm from the _Astronomical Almanac_ implemented
    by the routine sun_aaApparentApprox().
 \param[in]  t_cy     Julian centuries since J2000.0, TT timescale
 \param[in]  nut      Nutation terms and obliquity of the ecliptic, as returned
                      by functions sky0_nutationSpa() and sky0_epsilonSpa().

 \param[out] appV     Position vector of Sun in apparent coordinates (unit
                      vector i.e. direction cosines)
 \param[out] dist_au  Geocentric distance of the Sun (Astronomical Units)

 \par When to call this function
    In most cases, you will want to call sun_nrelApparent() or
    sun_nrelTopocentric() instead, and those functions call this function for
    you.

 \par Reference:
    Ibrahim Reda and Afshin Andreas,
    _Solar Position Algorithm for Solar Radiation Applications_
    National Renewable Energy Laboratory publication no. NREL/TP-560-34302,
    Revised January 2008
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double      t_ka;       // Millennia since J2000.0, TT timescale (J2KM)
    double      lambda_rad; // Sun longitude
    double      beta_rad;   // Sun latitude
    V3D_Matrix  epsM;       // Rotation by obliquity of ecliptic
    V3D_Vector  eclipV;     // Sun position in rectangular ecliptic coordinates

    REQUIRE_NOT_NULL(appV);
    REQUIRE_NOT_NULL(dist_au);

    t_ka = t_cy / 10.0;

    /* Calculate Sun longitude, latitude and distance from tables (steps 3.2 and
       3.3 of the algorithm in the SPA document). */
    *dist_au = sunDistance(t_ka);
    beta_rad = sunLatitude(t_ka);
    lambda_rad = sunLongitude(t_ka);

    /* Calculate the apparent longitude (in ecliptic coordinates), using the
       nutation in longitude and the light-time aberration correction. (The
       apparent longitude is the Sun's position 499 seconds earlier than the
       geometric position, because it took light that much time to reach Earth.)
       This is steps 3.6 and 3.7 of the SPA algorithm. */
    lambda_rad += nut->dPsi_rad - arcsecToRad(20.4898) / *dist_au;

    /* Convert to rectangular coordinates */
    v3d_polarToRect(&eclipV, lambda_rad, beta_rad);

    /* Rotate from ecliptic to equatorial coordinates */
    v3d_createRotationMatrix(&epsM, Xaxis, -(nut->eps0_rad + nut->dEps_rad));
    v3d_multMxV(appV, &epsM, &eclipV);
}



GLOBAL void sun_nrelApparent(double j2kTT_cy, Sky_TrueEquatorial *pos)
/*! Calculate the Sun's position as a unit vector and a distance, in apparent
    coordinates. It calls sun_nrelApp2() to obtain the Sun's position, after
    having called sky0_nutationSpa() to obtain the necessary nutation terms.
 \param[in]  j2kTT_cy   Julian centuries since J2000.0, TT timescale
 \param[out] pos        Timestamped structure containing position data and the
                        equation of the equinoxes

 \par When to call this function
 *  Because this function is computationally intensive, you may wish to limit
 *  your use of this function. 
 *      - if you want the Sun's position at multiple sites simultaneously at a
 *        single time, call this function, then follow it with a call to routine
 *        sky0_appToTirs(), and then make a separate call to
 *        sky_siteTirsToTopo() for each of one those sites.
 *      - if you want the Sun's position at one or more sites at closely spaced
 *        times (e.g. for tracking the Sun), pass this function to the
 *        skyfast_init() function.  skyfast_init() will
 *        call it for you several times, to fully calculate positions that will
 *        be saved and used later for interpolation by the skyfast_getApprox()
 *        function for tracking.
 * \par
 *  \em Alternatives:
 *      - If you want the Sun's position at a single site only at a single time,
 *        call sun_nrelTopocentric() instead, and it will call this function for
 *        you.
 *      - If you want the Sun's position at a single site for more than one time
 *        but the times are spaced more than a few hours apart, once again call
 *        sun_nrelTopocentric() instead.
 *      - If you want to track the sun, call skyfast_getApprox() instead.
 *        But this requires you to set up interpolation first with function
 *        skyfast_init() (and this function), as described above.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky0_Nut1980   nut;

    REQUIRE_NOT_NULL(pos);

    /* Calculate nutation (steps 3.4 of the algorithm in the SPA document) */
    sky0_nutationSpa(j2kTT_cy, &nut);
    
    /* Calculate the mean obliquity of the ecliptic (step 3.5) and the equation
       of the equinoxes */
    sky0_epsilonSpa(j2kTT_cy, &nut);
    pos->eqEq_rad = nut.eqEq_rad;

    /* Calculate sun apparent position */
    sun_nrelApp2(j2kTT_cy, &nut, &pos->appCirsV, &pos->distance_au);

    /* Now set the timestamp*/
    pos->timestamp_cy = j2kTT_cy;
}



GLOBAL void sun_nrelTopocentric(double             j2kUtc_d,
                                const Sky_DeltaTs  *deltas,
                                const Sky_SiteProp *site,
                                Sky_SiteHorizon *topo)
/*! Calls sun_nrelApparent() to calculate the Sun's position in apparent 
    coordinates using the NREL Sun Position Algorithm, and then converts this
    to topocentric horizon coordinates at the specified site.
 \param[in]  j2kUtc_d UTC time in "J2KD" form - i.e days since J2000.0
                      (= JD - 2 451 545.0)
 \param[in]  deltas   Delta T values, as set by the sky_initTime() (or 
                      sky_initTimeSimple() or sky_initTimeDetailed()) routines
 \param[in]  site     Properties of the observing site, particularly its
                      geometric location on the surface of the Earth, as set by
                      the sky_setSiteLocation() function (or sky_setSiteLoc2())
 \param[out] topo     Topocentric position, in both rectangular (unit vector)
                      form, and as Azimuth and Elevation (altitude) angles

 \par When to call this function
    Use this function if you are calculating the Sun topocentric position once,
    for a single site. But if you are going to be calculating it repeatedly, or
    for multiple sites, use of this function will cause you to perform a great
    many needless recalculations. Use skyfast_getApprox(), followed by
    sky0_appToTirs() and sky_siteTirsToTopo() instead.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_Times          atime;   // time, in various timescales
    Sky_TrueEquatorial pos;     // geocentric position of the Sun and distance
    V3D_Vector         terInterV; // unit vector in Terrestrial Intermed Ref Sys

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(site);
    REQUIRE_NOT_NULL(topo);

    sky_updateTimes(j2kUtc_d, deltas, &atime);

    sun_nrelApparent(atime.j2kTT_cy, &pos);

    /* Convert apparent position to topocentric Azimuth/Elevation coords */
    sky0_appToTirs(&pos.appCirsV, atime.j2kUT1_d, pos.eqEq_rad, &terInterV);
    sky_siteTirsToTopo(&terInterV, pos.distance_au, site, topo);
}



GLOBAL double sun_solarNoon(int                year,
                            int                month,
                            int                day,
                            const Sky_DeltaTs  *deltas,
                            const Sky_SiteProp *site,
                            Sky_SiteHorizon *topo)
/*! Routine to calculate the time of solar noon (Sun transit) for the day
    specified by \a year, \a month and \a day. This function uses the NREL SPA
    algorithm of sun_nrelTopocentric() to calculate the Sun's position.
 \returns              Solar noon for the day given in \a year, \a month and
                       \a day (returned as a J2KD date (= JD - 2 451 545.0),
                       UTC timescale). To view this as a local date and time,
                       add this value to \a site->timezone_d and pass the result
                       to function sky_j2kdToCalTime().
 \param[in] year, month, day
                       Date for which solar noon is desired
 \param[in]  deltas    Delta T values, as set by the sky_initTime() (or 
                       sky_initTimeSimple() or sky_initTimeDetailed()) routines
 \param[in]  site      Properties of the observing site, particularly its
                       geometric location on the surface of the Earth and its
                       time zone, as set by the sky_setSiteLocation() function
                       (or sky_setSiteLoc2())
 \param[out] topo      \b Optional. Topocentric position of the Sun at
                       transit, in both rectangular (unit vector) form, and as
                       Azimuth and Elevation (altitude) angles. If you are not
                       interested in these values, you can pass NULL to this
                       parameter.

 This routine uses an iterative approach. Two iterations is all that is needed
 to get a result within about 0.05 seconds.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_SiteHorizon topo1;      // Sun apparent position
    double          estimate_d; // estimate of the J2KD of solar noon

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(site);

    /* Make an initial guess of civil noon on the specified day */
    estimate_d = sky_calTimeToJ2kd(year, month, day,
                                   12, 0, 0.0, site->timeZone_d * 24.0);

    estimate_d = solarNoonApprox(estimate_d, deltas, site, &topo1);
    estimate_d = solarNoonApprox(estimate_d, deltas, site, &topo1);
 
    if (topo != NULL) {
        *topo = topo1;
    }
     return estimate_d;
}



GLOBAL double sun_riseSet(int                year,
                          int                month,
                          int                day,
                          bool               getSunrise,
                          const Sky_DeltaTs  *deltas,
                          const Sky_SiteProp *site,
                          Sky_SiteHorizon *topo)
/*! Routine to calculate the time of sunrise or sunset for the day specified by
    \a year, \a month and \a day. This function uses the NREL SPA algorithm of 
    sun_nrelTopocentric() to calculate the Sun's position.
 \returns                Sunrise (or sunset) time for the day given in \a year,
                         \a month and \a day (returned as a J2KD date 
                         (= JD - 2 451 545.0), UTC timescale).
                         To view this as a local date and time, add this
                         value to \a site->timezone_d and pass the result to
                         function sky_j2kdToCalTime().
 \param[in] year, month, day
                         Date for which sunrise or sunset time is desired
 \param[in]  getSunrise  If true, get sunrise time. If false, get sunset time
 \param[in]  deltas      Delta T values, as set by the sky_initTime() (or 
                         sky_initTimeSimple() or sky_initTimeDetailed())
                         routines
 \param[in]  site        Properties of the observing site, particularly its
                         geometric location on the surface of the Earth and its
                         time zone, as set by the sky_setSiteLocation() function
                         (or sky_setSiteLoc2()) and sky_setSiteTimeZone().
 \param[out] topo        \b Optional. Topocentric position of the Sun at rise or
                         set, in both rectangular (unit vector) form, and as
                         Azimuth and Elevation (altitude) angles. If you are not
                         interested in these values, you can pass NULL to this
                         parameter.

 This routine uses an iterative approach. It does two iterations.

 \note
   If the Sun does not rise (or set) on the specified day at the specified
   latitude, this function will return 0.0
 \note
 *  This routine assumes a standard refraction at the horizon of 34 arcminutes.
 *  Although this is typical of sunrise/sunset calculations, it is not
 *  necessarily the actual refraction that will occur at any given place or
 *  time. So the results should be considered approximate. Don't expect it to be
 *  be any better than the nearest minute to the time at which the sun appears
 *  or disappears.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_SiteHorizon topo1;      // Sun apparent position
    double          estimate_d; // estimate of the MJD of sunrise or sunset time
    Sky_SiteProp    st;         // copy of site details, without refraction

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(site);

   /* Sunrise and set calculations assume a standard refraction at the horizon
    * of 34 arcminutes, and a Sun semi-diameter of 16 arcminutes. So the 
    * calculation is based on the time at which the UNREFRACTED Sun position is
    * at -50 arcmin. So copy site information, and set refraction to zero for
    * that copy, which will be passed to riseSetApprox() */
    st = *site;
    st.refracPT = 0.0;

    /* Make an initial guess of 6 AM (or 6 PM) on the specified day */
    if (getSunrise) {
        estimate_d = sky_calTimeToJ2kd(year, month, day,
                                       6, 0, 0.0, site->timeZone_d * 24.0);
    } else {
        estimate_d = sky_calTimeToJ2kd(year, month, day,
                                       18, 0, 0.0, site->timeZone_d * 24.0);
    }

    estimate_d = riseSetApprox(estimate_d, getSunrise, deltas, &st, &topo1);
    if (estimate_d == 0.0) {
        /* Sun does not rise (or set) at the specified latitude on the requested
           day. Don't try another iteration. */
        return estimate_d;
    }
    estimate_d = riseSetApprox(estimate_d, getSunrise, deltas, &st, &topo1);

    if (topo != NULL) {
        *topo = topo1;
    }
    return estimate_d;  
}
                            

/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules).
 *
 *------------------------------------------------------------------------------
 */
LOCAL double sunLongitude(double t_ka)
/* This routine performs steps 3.2.1 to 3.2.4, 3.2.6, 3.3.1 and 3.3.2 of the
   algorithm outlined in the SPA document
 Returns - Sun geocentric longitude (radian) (geometric)
 Inputs
    t_ka - Julian ephemeris millennium, millennia since J2000.0, TT timescale

   Earth heliocentric longitude terms are stored in array L_TERMS.
   There are 6 intermediate terms to obtain (L0 -- L5), and each
   one is obtained by calculating
            L[i] = sum_over_rows_j( A[j] * cos(B[j] + C[j]*j_ka) )
   For L0, there are 64 rows j, for L1 there are 34 (etc - value
   is stored in lSubcount[i])
        Having obtained L0 -- L5, the longitude is obtained from
            (L0 + L1*j_ka + L2*j_ka^2 + ... + L5*j_ka^5) / 10^8
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /* Periodic terms for the Sun's ecliptic longitude */
    static const SunSeries l0[] = {
        {175347046.0, 0, 0},
        {  3341656.0, 4.6692568, 6283.07585},
        {    34894.0, 4.6261,   12566.1517 },
        {3497.0, 2.7441,  5753.3849},
        {3418.0, 2.8289,     3.5231},
        {3136.0, 3.6277, 77713.7715},
        {2676.0, 4.4181,  7860.4194},
        {2343.0, 6.1352,  3930.2097},
        {1324.0, 0.7425, 11506.7698},
        {1273.0, 2.0371,   529.691 },
        {1199.0, 1.1096,  1577.3435},
        {990, 5.233,  5884.927},
        {902, 2.045,    26.298},
        {857, 3.508,   398.149},
        {780, 1.179,  5223.694},
        {753, 2.533,  5507.553},
        {505, 4.583, 18849.228},
        {492, 4.205,   775.523},
        {357, 2.92,      0.067},
        {317, 5.849, 11790.629},
        {284, 1.899,   796.298},
        {271, 0.315, 10977.079},
        {243, 0.345,  5486.778},
        {206, 4.806,  2544.314},
        {205, 1.869,  5573.143},
        {202, 2.458,  6069.777},
        {156, 0.833,   213.299},
        {132, 3.411,  2942.463},
        {126, 1.083,    20.775},
        {115, 0.645,     0.98 },
        {103, 0.636,  4694.003},
        {102, 0.976, 15720.839},
        {102, 4.267,     7.114},
        { 99, 6.21,   2146.17},
        { 98, 0.68,    155.42},
        { 86, 5.98, 161000.69},
        { 85, 1.3,    6275.96},
        { 85, 3.67,  71430.7 },
        { 80, 1.81,  17260.15},
        { 79, 3.04,  12036.46},
        { 75, 1.76,   5088.63},
        { 74, 3.5,    3154.69},
        { 74, 4.68,    801.82},
        { 70, 0.83,   9437.76},
        { 62, 3.98,   8827.39},
        { 61, 1.82,   7084.9 },
        { 57, 2.78,   6286.6 },
        { 56, 4.39,  14143.5 },
        { 56, 3.47,   6279.55},
        { 52, 0.19,  12139.55},
        { 52, 1.33,   1748.02},
        { 51, 0.28,   5856.48},
        { 49, 0.49,   1194.45},
        { 41, 5.37,   8429.24},
        { 41, 2.4,   19651.05},
        { 39, 6.17,  10447.39},
        { 37, 6.04,  10213.29},
        { 37, 2.57,   1059.38},
        { 36, 1.71,   2352.87},
        { 36, 1.78,   6812.77},
        { 33, 0.59,  17789.85},
        { 30, 0.44,  83996.85},
        { 30, 2.74,   1349.87},
        { 25, 3.16,   4690.48}
    };
    static const SunSeries l1[] = {
        {628331966747.0, 0,     0},
        {206059.0, 2.678235, 6283.07585},
        {4303.0,   2.6351,  12566.1517},
        {425.0,    1.59,        3.523},
        {119.0,    5.796,      26.298},
        {109.0,    2.966,    1577.344},
        {93, 2.59, 18849.23},
        {72, 1.14,   529.69},
        {68, 1.87,   398.15},
        {67, 4.41,  5507.55},
        {59, 2.89,  5223.69},
        {56, 2.17,   155.42},
        {45, 0.4,    796.3 },
        {36, 0.47,   775.52},
        {29, 2.65,     7.11},
        {21, 5.34,     0.98},
        {19, 1.85,  5486.78},
        {19, 4.97,   213.3 },
        {17, 2.99,  6275.96},
        {16, 0.03,  2544.31},
        {16, 1.43,  2146.17},
        {15, 1.21, 10977.08},
        {12, 2.83,  1748.02},
        {12, 3.26,  5088.63},
        {12, 5.27,  1194.45},
        {12, 2.08,  4694   },
        {11, 0.77,   553.57},
        {10, 1.3,   6286.6 },
        {10, 4.24,  1349.87},
        { 9, 2.7,    242.73},
        { 9, 5.64,   951.72},
        { 8, 5.3,   2352.87},
        { 6, 2.65,  9437.76},
        { 6, 4.67,  4690.48}
    };
    static const SunSeries l2[] = {
        {52919.0, 0,         0},
        { 8720.0, 1.0721, 6283.0758},
        {  309.0, 0.867, 12566.152},
        {27, 0.05,    3.52 },
        {16, 5.19,    26.3 },
        {16, 3.68,   155.42},
        {10, 0.76, 18849.23},
        { 9, 2.06, 77713.77},
        { 7, 0.83,   775.52},
        { 5, 4.66,  1577.34},
        { 4, 1.03,     7.11},
        { 4, 3.44,  5573.14},
        { 3, 5.14,   796.3 },
        { 3, 6.05,  5507.55},
        { 3, 1.19,   242.73},
        { 3, 6.12,   529.69},
        { 3, 0.31,   398.15},
        { 3, 2.28,   553.57},
        { 2, 4.38,  5223.69},
        { 2, 3.75,     0.98}
    };
    static const SunSeries l3[] = {
        {289.0, 5.844, 6283.076},
        {35, 0,        0},
        {17, 5.49, 12566.15},
        { 3, 5.2,    155.42},
        { 1, 4.72,     3.52},
        { 1, 5.3,  18849.23},
        { 1, 5.97,   242.73}
    };
    static const SunSeries l4[] = {
        {114.0, 3.142, 0},
        {8, 4.13,  6283.08},
        {1, 3.84, 12566.15}
    };
    static const SunSeries l5[] = {
        {1, 3.14, 0}
    };
    /* Allow access to l0..l5 as if they form a single array lt[][] */
    static const SunSeries *lt[L_COUNT] = { l0, l1, l2, l3, l4, l5 };
    static const int lSubcount[L_COUNT] = { ARRAY_SIZE(l0),
                                            ARRAY_SIZE(l1),
                                            ARRAY_SIZE(l2),
                                            ARRAY_SIZE(l3),
                                            ARRAY_SIZE(l4),
                                            ARRAY_SIZE(l5)};

    double sum[L_COUNT];
    int    i, j;
    double earthSum;
    double tPower;

    /* Calculate the Earth heliocentric longitude (radian) */
    for (i = 0; i < L_COUNT; i++) {
        sum[i] = 0.0;
        for (j = 0; j < lSubcount[i]; j++) {
            sum[i] += lt[i][j].a * cos(lt[i][j].cb + lt[i][j].cct * t_ka);
        }
    }

    earthSum = 0.0;
    tPower = 1.0;
    for (i = 0; i < L_COUNT; i++) {
        /* replace  earthSum +=  sum[i] * pow(t_ka, i); with faster code */
        earthSum +=  sum[i] * tPower;
        tPower *= t_ka;
    }
    earthSum /= 1.0e8;
    
    /* Convert Earth heliocentric longitude to Sun geocentric longitude (radian)
       and force into the range 0 to TwoPi */
    earthSum += PI;
    
    return normalize(earthSum, TWOPI);
}



LOCAL double sunLatitude(double t_ka)
/* This routine performs steps 3.2.7 and 3.3.3 of the algorithm outlined
   in the SPA document
 Returns - Sun geocentric latitude (radian)
 Inputs
    t_ka - Julian ephemeris millennium, millennia since J2000.0, TT timescale

   Earth heliocentric latitude terms are stored in array B_TERMS.
   There are 2 intermediate terms to obtain (B0 -- B1), and each
   one is obtained by calculating
            B[i] = sum_over_rows_j( A[j] * cos(B[j] + C[j]*j_ka) )
   For B0, there are 5 rows j, for B1 there are 2 (value
   is stored in bSubcount(i))
        Having obtained B0 -- B1, the latitude is obtained from
            (B0 + L1*j_ka) / 10^8
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /* Periodic terms for the Sun's ecliptic latitude */
    static const SunSeries b0[] = {
        {280.0, 3.199, 84334.662},
        {102.0, 5.422, 5507.553},
        { 80.0, 3.88,  5223.69},
        { 44.0, 3.7,   2352.87},
        { 32.0, 4.0,   1577.34}
    };
    static const SunSeries b1[] = {
        {9.0, 3.9,  5507.55},
        {6.0, 1.73, 5223.69}
    };
    /* Allow access to b0..b1 as if they form a single array bt[][] */
    static const SunSeries *bt[B_COUNT] = { b0, b1 };
    static const int bSubcount[B_COUNT] = { ARRAY_SIZE(b0),
                                            ARRAY_SIZE(b1) };

    double sum[B_COUNT];
    int    i, j;
    double earthSum;

    /* Calculate the Earth heliocentric latitude (radian) */
    for (i = 0; i < B_COUNT; i++) {
        sum[i] = 0.0;
        for (j = 0; j < bSubcount[i]; j++) {
            sum[i] += bt[i][j].a * cos(bt[i][j].cb + bt[i][j].cct * t_ka);
        }
    }

    earthSum = sum[0] + sum[1] * t_ka;
    earthSum /= 1.0e8;

    /* Convert Earth heliocentric latitude to Sun geocentric latitude (radian)*/
    earthSum = -earthSum;
    
    return earthSum;
}



LOCAL double sunDistance(double t_ka)
/* This routine performs step 3.2.8 of the algorithm outlined in the SPA
   document.
 Returns - distance to the sun (astronomical units)
 Inputs
    t_ka - Julian ephemeris millennium, millennia since J2000.0, TT timescale

   Earth heliocentric distance terms are stored in array R_TERMS.
   There are 5 intermediate terms to obtain (R0 -- R4), and each
   one is obtained by calculating
            R[i] = sum_over_rows_j( A[j] * cos(B[j] + C[j]*j_ka) )
   For R0, there are 40 rows j, for R1 there are 10 (etc - value
   is stored in rSubcount(i))
        Having obtained R0 -- R4, the distance is obtained from
            (R0 + R1*j_ka + R2*j_ka^2 + ... + R4*j_ka^4) / 10^8
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /* Periodic terms for the Earth - Sun distance */
    static const SunSeries r0[] = {
        {100013989.0, 0.0, 0.0},
        {1670700.0, 3.0984635, 6283.07585},
        {13956.0,   3.05525,  12566.1517},
        {3084.0, 5.1985, 77713.7715},
        {1628.0, 1.1739, 5753.3849},
        {1576.0, 2.8469, 7860.4194},
        { 925.0, 5.453, 11506.77},
        { 542.0, 4.564,  3930.21},
        { 472.0, 3.661,  5884.927},
        { 346.0, 0.964,  5507.553},
        { 329.0, 5.9,    5223.694},
        { 307.0, 0.299,  5573.143},
        { 243.0, 4.273, 11790.629},
        { 212.0, 5.847,  1577.344},
        { 186.0, 5.022, 10977.079},
        { 175.0, 3.012, 18849.228},
        { 110.0, 5.055,  5486.778},
        {  98.0, 0.89,   6069.78},
        {  86.0, 5.69,  15720.84},
        {  86.0, 1.27, 161000.69},
        {  65.0, 0.27,  17260.15},
        {  63.0, 0.92,    529.69},
        {  57.0, 2.01,  83996.85},
        {  56.0, 5.24,  71430.7},
        {  49.0, 3.25,   2544.31},
        {  47.0, 2.58,    775.52},
        {  45.0, 5.54,   9437.76},
        {  43.0, 6.01,   6275.96},
        {  39.0, 5.36,   4694.0},
        {  38.0, 2.39,   8827.39},
        {  37.0, 0.83,  19651.05},
        {  37.0, 4.9,   12139.55},
        {  36.0, 1.67,  12036.46},
        {  35.0, 1.84,   2942.46},
        {  33.0, 0.24,   7084.9},
        {  32.0, 0.18,   5088.63},
        {  32.0, 1.78,    398.15},
        {  28.0, 1.21,   6286.6},
        {  28.0, 1.9,    6279.55},
        {  26.0, 4.59,  10447.39}
    };
    static const SunSeries r1[] = {
        {103019.0, 1.10749, 6283.07585},
        {  1721.0, 1.0644, 12566.1517},
        {702.0, 3.142,    0.0 },
        { 32.0, 1.02, 18849.23},
        { 31.0, 2.84,  5507.55},
        { 25.0, 1.32,  5223.69},
        { 18.0, 1.42,  1577.34},
        { 10.0, 5.91, 10977.08},
        {  9.0, 1.42,  6275.96},
        {  9.0, 0.27,  5486.78}
    };
    static const SunSeries r2[] = {
        {4359.0, 5.7846, 6283.0758},
        { 124.0, 5.579, 12566.152},
        {  12.0, 3.14,      0.0 },
        {   9.0, 3.63,  77713.77},
        {   6.0, 1.87,   5573.14},
        {   3.0, 5.47,  18849.23}
    };
    static const SunSeries r3[] = {
        {145.0, 4.273, 6283.076},
        {  7.0, 3.92, 12566.15 }
    };
    static const SunSeries r4[] = {
        {4.0, 2.56, 6283.08}
    };
    /* Allow access to r0..r4 as if they form a single array rt[][] */
    static const SunSeries *rt[R_COUNT] = { r0, r1, r2, r3, r4 };
    static const int rSubcount[R_COUNT] = { ARRAY_SIZE(r0),
                                            ARRAY_SIZE(r1),
                                            ARRAY_SIZE(r2),
                                            ARRAY_SIZE(r3),
                                            ARRAY_SIZE(r4)};

    double sum[R_COUNT];
    int    i, j;
    double earthSum;
    double tPower;

    for (i = 0; i < R_COUNT; i++) {
        sum[i] = 0.0;
        for (j = 0; j < rSubcount[i]; j++) {
            sum[i] += rt[i][j].a * cos(rt[i][j].cb + rt[i][j].cct * t_ka);
        }
    }
    
    earthSum = 0.0;
    tPower = 1.0;
    for (i = 0; i < R_COUNT; i++) {
        /* replace  earthSum +=  sum[i] * pow(t_ka, i); with faster code */
        earthSum +=  sum[i] * tPower;
        tPower *= t_ka;
    }
    earthSum /= 1.0e8;

    return earthSum;
}



LOCAL double solarNoonApprox(double             noonGuess_d,
                             const Sky_DeltaTs  *deltas,
                             const Sky_SiteProp *site,
                             Sky_SiteHorizon *topo)
/* Routine to calculate the time of solar noon for the day specified by
   noonGuess_d. The result returned is an approximate value, whose accuracy
   depends upon how close noonGuess_d is to true solar noon.
      If noonGuess_d is equal to civil noon in local time, the result can be
   expected to be within a few seconds of the true solar noon. Of course, if
   this routine is called again with the previous result used as the new guess,
   the new result will be more accurate. Generally, only one or two calls would
   ever be required.
 Returns         - Improved estimate of time of solar noon for the day & time
                   given in noonGuess_d, in the form J2KD (= JD - 2 451 545.0)
 Inputs
    noonGuess_d  - J2KD (= JD - 2 451 545.0) of the date at which solar noon is
                   desired, and the time of day of a guess of when solar noon
                   might be
    deltas       - delta times, to convert from UTC to TT
    site         - site properties
 Outputs
    topo         - Topocentric position of Sun at transit
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double          dec;        // ignored
    double          ha_rad;     // hour angle of the Sun (radian)
    
    sun_nrelTopocentric(noonGuess_d, deltas, site, topo);
    sky_siteAzElToHaDec(&topo->rectV, site, &ha_rad, &dec);

    /* HA gives time since celestial object passed meridian (so -ve HA gives
       time until object will pass meridian) in units of Sidereal time. But
       the sun moves against the fixed stars in that time, so, as an
       approximation, subtract HA from Solar time (UT1) instead. */
    return noonGuess_d - ha_rad / TWOPI;
}



LOCAL double riseSetApprox(double             risesetGuess_d,
                           bool               getSunrise,
                           const Sky_DeltaTs  *deltas,
                           const Sky_SiteProp *site,
                           Sky_SiteHorizon *topo)
/* Routine to calculate the time of Sun rise or set for the day specified by
   risesetGuess_d. The result returned is an approximate value, whose accuracy
   depends upon how close risesetGuess_d is to true Sun rise or set time.
   6 AM local time would be a good starting sunrise guess, and 6PM a good sunset
   guess.
   Of course, if this routine is called again with the previous result used as
   the new guess, the new result will be more accurate. Convergence is not
   quite as fast as the SolarNoonApprox() algorithm above.
 Returns            - Improved estimate of time of sunrise or sunset, (or zero
                      if the sun does not rise or does not set on this date)
 Inputs
    risesetGuess_d  - J2KD (= JD - 2 451 545.0) of date at which sunrise or
                      sunset time is desired, and time of day of a guess of when
                      sunrise or sunset might be
    getSunrise      - If true, get sunrise time. If false, get sunset time
    deltas          - delta times, to convert from UTC to TT
    site            - properties of the site. Note: it is assumed that the field
                      site->refracPT will be set to 0.0 before calling this
                      routine, in order to calculate an unrefracted position of
                      the Sun (as per note below).
 Output
    topo            - Topocentric position of Sun at rise or set

   Sunrise and set calculations assume a standard refraction at the horizon of
   34 arcminutes, and a sun semi-diameter of 16 arcminutes. So the calculation
   is based on the time at which the UNREFRACTED Sun position is at -50 arcmin.
 * 
 * TODO
 *  This routine uses the equation 
 *      cos(HA) = (sin(-50′) - sin(ϕA) sin(δ)) / (cos(ϕA) cos(δ))
 *  which means it will fail if astronomical latitude ϕA = π/2 (i.e. at the
 *  poles), and will presumably get less accurate the closer to the poles we
 *  get. Find a better expression.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double          ha1_rad;    // Hour angle of Sun at time riseSetGuess_d(rad)
    double          dec_rad;    // Declination of Sun (radian)
    double          ha2_rad;    // Hour angle of Sun at horizon (radian)
    double          cosHa2;     // Cos(Hour Angle) at horizon
    double          riseSetApprox_d;// Improved estimate of rise or set time

    sun_nrelTopocentric(risesetGuess_d, deltas, site, topo);
    sky_siteAzElToHaDec(&topo->rectV, site, &ha1_rad, &dec_rad);
    
    /* assuming Dec remains constant over the period, find out where dec circle
       intersects the Elevation = -50 arcminute circle */
    cosHa2 = (sin(degToRad(-50.0/60.0)) - sin(site->astLat_rad) * sin(dec_rad))
             / (cos(site->astLat_rad) * cos(dec_rad));
    /* If there is no intersection, the Sun either doesn't rise or doesn't set
       on this day at this latitude.  */
    if (cosHa2 > 1.0) {
        /* Sun does not rise */
        riseSetApprox_d = 0.0;
    } else if (cosHa2 < -1.0) {
        /* Sun does not set */
        riseSetApprox_d = 0.0;
    } else {
        if (getSunrise) {
            ha2_rad = -acos(cosHa2);
        } else {
            ha2_rad = acos(cosHa2);
        }
        /* Subtract difference in HA from Solar time (not sidereal time)
           i.e. this is the same approximation used in routine
           solarNoonApprox() above. */
        riseSetApprox_d = risesetGuess_d - (ha1_rad - ha2_rad) / TWOPI;
    }
    return riseSetApprox_d;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
