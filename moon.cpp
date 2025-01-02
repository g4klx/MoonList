/*==============================================================================
 * moon.c - routines to calculate the Moon's position
 *
 * Author:  David Hoadley
 *
 * Description: (see moon.h)
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
#include <float.h>
#include "instead-of-math.h"                /* for sincos() & normalize() */
#include <math.h>
#include <stdlib.h>

/* Local and project includes */
#include "moon.h"

#include "general.h"
///+
#include <stdio.h>
#include "skyio.h"
///-

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       /* For use by REQUIRE() - assertions. */

typedef struct {
    double lp_rad;  // L' - Mean Longitude of the Moon (radian)
    double d_rad;   // D  - Mean Elongation of the Moon from the Sun (radian)
    double m_rad;   // M  - Mean Anomaly of the Sun (radian)
    double mp_rad;  // M' - Mean Anomaly of the Moon (radian)
    double f_rad;   // F  - so-called Argument of Latitude of the Moon (radian)
    double e;       // E  - Eccentricity of the Earth's orbit around the Sun
} OrbTerms;

/* Note there is some inconsistency in the use of terms here, when compared
 * with the fundamental arguments of the 1980 nutation theory.
 * NREL Moon  Nutation1980    Term description
 *    L'         (L)           Mean Longitude of the Moon
 *    D           D            Mean Elongation of the Moon from the Sun
 *    M           l'           Mean Anomaly of the Sun
 *    M'          l            Mean Anomaly of the Moon
 *                Ω            Longitude of ascending node of Moon
 *    F           F            = L - Ω, Mean Longitude of Moon minus Longitude
 *                                    of ascending node.
 *
 * In the NREL SAMPA document, sections 3.2.5 and 3.2.4, the term F is referred
 * to as the "Argument of Latitude". I believe this is wrong. F is actually
 * the MEAN Longitude of Moon minus the Longitude of the Moon's ascending node
 * (F = L - Ω), or alternatively the MEAN anomaly plus the argument of periapsis
 * (F = M' + ω). But the "Argument of Latitude" is defined as the TRUE anomaly
 * plus the argument of periapsis - i.e. u = ν + ω.
 *
 * In the code that follows, I have kept this label, to make it easier for
 * someone reading the NREL SAMPA document to follow the code, but I have marked
 * it as "so-called" each time.
 */

typedef struct {
    double cd;
    int    cm;
    double cmp;
    double cf;
    double cl;
    double cr;
} LonDistTerms;


typedef struct {
    double cd;
    int    cm;
    double cmp;
    double cf;
    double cb;
} LatTerms;

/*
 * Prototypes for local functions (not called from other modules)
 */
LOCAL void moonOrbitals(double t_cy, OrbTerms *orb);
LOCAL void moonLongDist(const OrbTerms *orb, double *long_udeg, double *dist_m);
LOCAL double moonLatitude(const OrbTerms *orb);
LOCAL double riseSetApprox(double             risesetGuess_d,
                           bool               getMoonrise,
                           const Sky_DeltaTs  *deltas,
                           const Sky_SiteProp *site,
                           Sky_SiteHorizon    *topo);

/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */
/*      Constants found in the 2007 Astronomical Almanac, pages K6 & K7 */
LOCAL const double au_km = 1.49597871464e8;  // one Astronomical Unit (km)

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
GLOBAL void moon_nrelApp2(double             t_cy,
                          const Sky0_Nut1980 *nut,
                          V3D_Vector *appV,
                          double     *dist_au)
/*! Calculates the Moon's position in geocentric apparent coordinates, using the
    NREL Moon Position Algorithm.
 \param[in]  t_cy    Julian centuries since J2000.0, TT timescale
 \param[in]  nut     Nutation terms and obliquity of the ecliptic
 \param[out] appV    Position vector of Moon in apparent coordinates (unit
                       vector i.e. direction cosines)
 \param[out] dist_au Geocentric distance of the Moon (Astronomical Units)

 \par Reference
    Reda, I., "Solar Eclipse Monitoring for Solar Energy Applications Using the
    Solar and Moon Position Algorithms". National Renewable Energy Laboratory
    Technical Report NREL/TP-3B0-47681, March 2010
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    OrbTerms    orbt;               // Orbital terms for the Moon
    double      a1_rad, a2_rad, a3_rad;  // mysterious factors (radian)
    double      l_udeg, b_udeg;     // Moon long & lat terms (micro-degrees)
    double      r_m;                // Moon distance term (metres)
    double      dl_udeg, db_udeg;   // Δl & Δb (micro-degrees)
    double      lamdap_rad;         // λ' - Moon uncorrected longitude (radian)
    double      lamda_rad;          // λ  - Moon apparent longitude (radian)
    double      beta_rad;           // β  - Moon latitude (radian)
    double      mDelta_km;          // Δ  - Moon geocentric distance (km)
    V3D_Matrix  epsM;               // Rotation by obliquity of ecliptic
    V3D_Vector  eclV;               // Rectangular coordinates in ecliptic frame

    REQUIRE_NOT_NULL(appV);
    REQUIRE_NOT_NULL(dist_au);

    /* Steps 3.2.1 tp 3.2.5, and equation (15) of 3.2.6 */
    moonOrbitals(t_cy, &orbt);

    /* Steps 3.2.6 to 3.2.8 (the first!) */
    moonLongDist(&orbt, &l_udeg, &r_m);
    b_udeg = moonLatitude(&orbt);

    /* Steps 3.2.8 (the second!) 3.2.9 and 3.2.10 */
    a1_rad = degToRad(119.75 + 131.849 * t_cy);
    a2_rad = degToRad(53.09 + 479264.29 * t_cy);
    a3_rad = degToRad(313.45 + 481266.484 * t_cy);

    /* Steps 3.2.11 & 3.2.12 */
    dl_udeg = 3958.0 * sin(a1_rad) + 1962.0 * sin(orbt.lp_rad - orbt.f_rad)
              + 318.0 * sin(a2_rad);
    db_udeg = -2235.0 * sin(orbt.lp_rad) + 382.0 * sin(a3_rad)
              + 175.0 * sin(a1_rad - orbt.f_rad)
              + 175.0 * sin(a1_rad + orbt.f_rad)
              + 127.0 * sin(orbt.lp_rad - orbt.mp_rad)
              - 115.0 * sin(orbt.lp_rad + orbt.mp_rad);

    /* Steps 3.2.13 to 3.2.16 */
    lamdap_rad = normalize(orbt.lp_rad + degToRad((l_udeg + dl_udeg) / 1e6),
                           TWOPI);
    beta_rad = degToRad((b_udeg + db_udeg) / 1e6);
    mDelta_km = 385000.56 + r_m / 1000.0;
    /* Convert distance to the Moon from km to Astronomical Units, because all
       other celestial objects are specified that way, so our routines
       converting geocentric to topocentric place expect it to be that way. */
    *dist_au = mDelta_km / au_km;

    /* Step 3.6 */
    lamda_rad = lamdap_rad + nut->dPsi_rad;

    /* Convert to rectangular coordinates in ecliptic plane and rotate to
       equatorial plane (i.e to RA/Dec frame). */
    v3d_polarToRect(&eclV, lamda_rad, beta_rad);
    v3d_createRotationMatrix(&epsM, Xaxis, -(nut->eps0_rad + nut->dEps_rad));
    v3d_multMxV(appV, &epsM, &eclV);

    /* Now have vector of apparent coords. Could convert to RA and Dec using
       v3d_rectToPolar(&RA, &Dec, appV) but there isn't really any need.
       Had we done so we would have performed steps 3.8 and 3.9 of NREL SAMPA */
}



GLOBAL void moon_nrelApparent(double j2kTT_cy, Sky_TrueEquatorial *pos)
/*! Calculate the Moon's position as a unit vector and a distance, in apparent
    coordinates. It calls moon_nrelApp2() to obtain the Moon's position, after
    having called sky0_nutationSpa() to obtain the necessary nutation terms.
 \param[in]  j2kTT_cy   Julian centuries since J2000.0, TT timescale
 \param[out] pos        Timestamped structure containing position data and the
                        equation of the equinoxes

 \par When to call this function
    Because this function is computationally intensive, you may wish to limit
    your use of this function.
        - if you want the Moon's position at multiple sites simultaneously at a
          single time, call this function, then follow it with a call to routine
          sky0_appToTirs(), and then make a separate call to
          sky_siteTirsToTopo() for each of one those sites.
        - if you want the Moon's position at one or more sites at closely spaced
          times (e.g. for tracking the Moon), pass this function to the
          skyfast_init() function.  skyfast_init() will
          call it for you several times, to fully calculate positions that will
          be saved and used later for interpolation by the skyfast_getApprox()
          function for tracking.
 \par
    \em Alternatives:
        - If you want the Moon's position at a single site only at a single
          time, call moon_nrelTopocentric() instead, and it will call this
          function for you.
        - If you want the Moon's position at a single site for more than one
          time but the times are spaced more than an hour or so apart, once
          again call moon_nrelTopocentric() instead.
        - If you want to track the Moon, call skyfast_getApprox() instead.
          But this requires you to set up interpolation first with function
          skyfast_init() (and this function), as described above.
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

    /* Calculate Moon apparent position */
    moon_nrelApp2(j2kTT_cy, &nut, &pos->appCirsV, &pos->distance_au);

    /* Now set the timestamp*/
    pos->timestamp_cy = j2kTT_cy;
}



GLOBAL void moon_nrelTopocentric(double             j2kdUtc,
                                 const Sky_DeltaTs  *deltas,
                                 const Sky_SiteProp *site,
                                 Sky_SiteHorizon *topo)
/*! Calls moon_nrelApparent() to calculate the Moon's position in apparent
 *  coordinates using the NREL Moon Position Algorithm, and then converts this
 *  to topocentric horizon coordinates at the specified site.
 \param[in]  j2kdUtc  UTC time in "J2KD" form - i.e days since J2000.0
                      (= JD - 2 451 545.0)
 \param[in]  deltas   Delta T values, as set by the sky_initTime() (or
                      sky_initTimeSimple() or sky_initTimeDetailed()) routines
 \param[in]  site     Properties of the observing site, particularly its
                      geometric location on the surface of the Earth, as set by
                      the sky_setSiteLocation() function (or sky_setSiteLoc2())
 \param[out] topo     Topocentric position, in both rectangular (unit vector)
                      form, and as Azimuth and Elevation (altitude) angles
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_Times          atime;     // time, in various timescales
    Sky_TrueEquatorial pos;       // geocentric position of the Moon and dist.
    V3D_Vector         terInterV; // unit vector in Terrestrial Intermed Ref Sys

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(site);
    REQUIRE_NOT_NULL(topo);

    sky_updateTimes(j2kdUtc, deltas, &atime);

    /* Calculate Moon apparent position */
    moon_nrelApparent(atime.j2kTT_cy, &pos);

    /* Convert apparent position to topocentric Azimuth/Elevation coords */
    sky0_appToTirs(&pos.appCirsV, atime.j2kUT1_d, pos.eqEq_rad, &terInterV);
    sky_siteTirsToTopo(&terInterV, pos.distance_au, site, topo);
}



GLOBAL double moon_riseSet(int                   year,
                           int                   month,
                           int                   day,
                           bool                  getMoonrise,
                           const Sky_DeltaTs *deltas,
                           const Sky_SiteProp    *site,
                           Sky_SiteHorizon *topo)
/*! Routine to calculate the time of moonrise or moonset for the day specified
    by \a year, \a month and \a day. This function uses the NREL MPA algorithm
    of moon_nrelTopocentric() to calculate the Moon's position.
 \returns                Moonrise (or moonset) time for the day given in
                          \a year, \a month and \a day (returned as a J2KD date
                         (= JD - 2 451 545.0), UTC timescale).
                         To view this as a local date and time, add this
                         value to \a site->timezone_d and pass the result to
                         function sky_j2kdToCalTime().
 \param[in]  year, month, day
                         Date for which moonrise or moonset time is desired
 \param[in]  getMoonrise If true, get moonrise time. If false, get moonset time
 \param[in]  deltas      Delta T values, as set by the sky_initTime() (or
                         sky_initTimeSimple() or sky_initTimeDetailed())
                         routines
 \param[in]  site        Properties of the observing site, particularly its
                         geometric location on the surface of the Earth and its
                         time zone, as set by the sky_setSiteLocation() function
                         (or sky_setSiteLoc2()) and sky_setSiteTimeZone().
 \param[out] topo        \b Optional. Topocentric position of the Moon at rise
                         or set, in both rectangular (unit vector) form, and as
                         Azimuth and Elevation (altitude) angles. If you are not
                         interested in these values, you can pass NULL to this
                         parameter.

 This routine uses an iterative approach. Three iterations seems to be enough.

 \note
    If the Moon does not rise (or set) on the specified day, this routine will
    return the time just before midnight on the previous day, or just after
    midnight on the next day, when the rise (or set) actually occurs.
    If the routine encounters an error, it will return 0.0
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_SiteHorizon topo1;      // Moon apparent position
    double          estimate_d; // estimate of the MJD of Moon rise or set time
    Sky_SiteProp    st;         // copy of site details, without refraction
    int             i;          // iteration counter

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(site);

    /* Moonrise and set calculations assume a standard refraction at the horizon
    * of 34 arcminutes, and a Moon semi-diameter of 16 arcminutes. So the
    * calculation is based on the time at which the UNREFRACTED Moon position is
    * at -50 arcmin. So copy site information, and set refraction to zero for
    * that copy, which will be passed to riseSetApprox() */
    st = *site;
    st.refracPT = 0.0;

    /* Make an initial guess of civil noon on the specified day */
    estimate_d = sky_calTimeToJ2kd(year, month, day,
                                   12, 0, 0.0, site->timeZone_d * 24.0);

    /* Iterate to settle on a time. */
    for (i = 0; i < 3; i++) {
        estimate_d = riseSetApprox(estimate_d, getMoonrise, deltas,&st, &topo1);
        if (estimate_d == 0.0) {
            /* Error in calculations on the requested day.
               Don't try another iteration. */
            break;
        }
    }

    if (topo != NULL) {
        *topo = topo1;
    }
    return estimate_d;
}

/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules)
 *
 *------------------------------------------------------------------------------
 */
LOCAL void moonOrbitals(double t_cy, OrbTerms *orb)
/* Calculate the fundamental orbital parameters required to get the moon
   position. This is steps 3.2.1 to 3.2.5 of the NREL SAMPA document.
 Inputs
    t_cy    - julian centuries since J2000.0, TT timescale
 Outputs
    orb     - orbital terms for the moon, and eccentricity of Earth's orbit
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double temp_deg;

    /* Moon's Mean Longitude L' */
    temp_deg = 218.3164477
                + t_cy * (481267.88123421
                           + t_cy * (-0.0015786
                                     + t_cy * (1.0 / 538841.0
                                               + t_cy * (-1.0 / 65194000.0))));
    orb->lp_rad = normalize(degToRad(temp_deg), TWOPI);

    /* Mean Elongation of the Moon D */
    temp_deg = 297.8501921
                + t_cy * (445267.1114034
                          + t_cy * (-0.0018819
                                    + t_cy * (1.0 / 545868.0
                                              + t_cy * (-1.0 / 113065000.0))));
    orb->d_rad = normalize(degToRad(temp_deg), TWOPI);

    /* Sun's Mean Anomaly M */
    temp_deg = 357.5291092
                + t_cy * (35999.0502909
                          + t_cy * (-0.0001536
                                    + t_cy * (1.0 / 24490000.0)));
    orb->m_rad = normalize(degToRad(temp_deg), TWOPI);

    /* Moon's Mean Anomaly M' */
    temp_deg = 134.9633964
                + t_cy * (477198.8675055
                          + t_cy * (0.0087414
                                    + t_cy * (1.0 / 69699.0
                                              + t_cy * (-1.0 / 14712000.0))));
    orb->mp_rad = normalize(degToRad(temp_deg), TWOPI);

    /* Moon's so-called Argument of Latitude F */
    temp_deg = 93.2720950
                + t_cy * (483202.0175233
                          + t_cy * (-0.0036539
                                    + t_cy * (-1.0 / 3526000.0
                                              + t_cy * (1.0 / 863310000.0))));
    orb->f_rad = normalize(degToRad(temp_deg), TWOPI);

    /* Eccentricity of the Earth's orbit around the Sun */
    orb->e = 1.0 - 0.002516 * t_cy - 0.0000074 * t_cy * t_cy;
}



LOCAL void moonLongDist(const OrbTerms *orb, double *long_udeg, double *dist_m)
/* Accumulate longitude and distance terms from table lonD[]. This is steps
   3.2.6 and 3.2.7 of the NREL SAMPA document.
 Inputs
    orb         - orbital terms
 Outputs
    long_udeg   - accumulator of terms in longitude (micro-degrees)
    dist_m      - accumulator of terms in distance (metres)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    static const LonDistTerms lonD[] = {
        { 0, 0, 1, 0, 6288774,-20905355 },
        { 2, 0,-1, 0, 1274027,-3699111 },
        { 2, 0, 0, 0, 658314,-2955968 },
        { 0, 0, 2, 0, 213618,-569925 },
        { 0, 1, 0, 0,-185116, 48888 },
        { 0, 0, 0, 2,-114332,-3149 },
        { 2, 0,-2, 0, 58793, 246158 },
        { 2,-1,-1, 0, 57066,-152138 },
        { 2, 0, 1, 0, 53322,-170733 },
        { 2,-1, 0, 0, 45758,-204586 },
        { 0, 1,-1, 0,-40923,-129620 },
        { 1, 0, 0, 0,-34720, 108743 },
        { 0, 1, 1, 0,-30383, 104755 },
        { 2, 0, 0,-2, 15327, 10321 },
        { 0, 0, 1, 2,-12528, 0 },
        { 0, 0, 1,-2, 10980, 79661 },
        { 4, 0,-1, 0, 10675,-34782 },
        { 0, 0, 3, 0, 10034,-23210 },
        { 4, 0,-2, 0, 8548,-21636 },
        { 2, 1,-1, 0,-7888, 24208 },
        { 2, 1, 0, 0,-6766, 30824 },
        { 1, 0,-1, 0,-5163,-8379 },
        { 1, 1, 0, 0, 4987,-16675 },
        { 2,-1, 1, 0, 4036,-12831 },
        { 2, 0, 2, 0, 3994,-10445 },
        { 4, 0, 0, 0, 3861,-11650 },
        { 2, 0,-3, 0, 3665, 14403 },
        { 0, 1,-2, 0,-2689,-7003 },
        { 2, 0,-1, 2,-2602, 0 },
        { 2,-1,-2, 0, 2390, 10056 },
        { 1, 0, 1, 0,-2348, 6322 },
        { 2,-2, 0, 0, 2236,-9884 },
        { 0, 1, 2, 0,-2120, 5751 },
        { 0, 2, 0, 0,-2069, 0 },
        { 2,-2,-1, 0, 2048,-4950 },
        { 2, 0, 1,-2,-1773, 4130 },
        { 2, 0, 0, 2,-1595, 0 },
        { 4,-1,-1, 0, 1215,-3958 },
        { 0, 0, 2, 2,-1110, 0 },
        { 3, 0,-1, 0,-892, 3258 },
        { 2, 1, 1, 0,-810, 2616 },
        { 4,-1,-2, 0, 759,-1897 },
        { 0, 2,-1, 0,-713,-2117 },
        { 2, 2,-1, 0,-700, 2354 },
        { 2, 1,-2, 0, 691, 0 },
        { 2,-1, 0,-2, 596, 0 },
        { 4, 0, 1, 0, 549,-1423 },
        { 0, 0, 4, 0, 537,-1117 },
        { 4,-1, 0, 0, 520,-1571 },
        { 1, 0,-2, 0,-487,-1739 },
        { 2, 1, 0,-2,-399, 0 },
        { 0, 0, 2,-2,-381,-4421 },
        { 1, 1, 1, 0, 351, 0 },
        { 3, 0,-2, 0,-340, 0 },
        { 4, 0,-3, 0, 330, 0 },
        { 2,-1, 2, 0, 327, 0 },
        { 0, 2, 1, 0,-323, 1165 },
        { 1, 1,-1, 0, 299, 0 },
        { 2, 0, 3, 0, 294, 0 },
        { 2, 0,-1,-2, 0,   8752 },
    };
    static const int lonDistSize = ARRAY_SIZE(lonD);

    int i;
    double li;                  // lonD[i].cl, maybe scaled for eccentricity
    double ri;                  // lonD[i].cr, maybe scaled for eccentricity
    double a_rad, sinA, cosA;   // angle - summation of args (radian)


    *long_udeg = 0.0;
    *dist_m = 0.0;
    for (i = lonDistSize - 1; i >= 0; i--) {
        li = lonD[i].cl;
        ri = lonD[i].cr;
        if (lonD[i].cm != 0) {
            // Multiply terms by E or E^2 as appropriate
            li *= orb->e;
            ri *= orb->e;
            if (abs(lonD[i].cm) == 2) {
                li *= orb->e;
                ri *= orb->e;
            }
        }

        a_rad = lonD[i].cd * orb->d_rad + lonD[i].cm * orb->m_rad
                + lonD[i].cmp * orb->mp_rad + lonD[i].cf * orb->f_rad;
        sincos(a_rad, &sinA, &cosA);

        *long_udeg += li * sinA;
        *dist_m += ri * cosA;
    }
}



LOCAL double moonLatitude(const OrbTerms *orb)
/* Accumulate latitude terms from table lat[]. This is steps 3.2.8 (the first of
   the two sections marked 3.2.8!) of of the NREL SAMPA document.
 Returns    - accumulator of terms in latitude (micro-degrees)
 Input
    orb     - orbital terms
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    static const LatTerms lat[] = {
        { 0, 0, 0, 1, 5128122 },
        { 0, 0, 1, 1, 280602 },
        { 0, 0, 1,-1, 277693 },
        { 2, 0, 0,-1, 173237 },
        { 2, 0,-1, 1, 55413 },
        { 2, 0,-1,-1, 46271 },
        { 2, 0, 0, 1, 32573 },
        { 0, 0, 2, 1, 17198 },
        { 2, 0, 1,-1, 9266 },
        { 0, 0, 2,-1, 8822 },
        { 2,-1, 0,-1, 8216 },
        { 2, 0,-2,-1, 4324 },
        { 2, 0, 1, 1, 4200 },
        { 2, 1, 0,-1,-3359 },
        { 2,-1,-1, 1, 2463 },
        { 2,-1, 0, 1, 2211 },
        { 2,-1,-1,-1, 2065 },
        { 0, 1,-1,-1,-1870 },
        { 4, 0,-1,-1, 1828 },
        { 0, 1, 0, 1,-1794 },
        { 0, 0, 0, 3,-1749 },
        { 0, 1,-1, 1,-1565 },
        { 1, 0, 0, 1,-1491 },
        { 0, 1, 1, 1,-1475 },
        { 0, 1, 1,-1,-1410 },
        { 0, 1, 0,-1,-1344 },
        { 1, 0, 0,-1,-1335 },
        { 0, 0, 3, 1, 1107 },
        { 4, 0, 0,-1, 1021 },
        { 4, 0,-1, 1, 833 },
        { 0, 0, 1,-3, 777 },
        { 4, 0,-2, 1, 671 },
        { 2, 0, 0,-3, 607 },
        { 2, 0, 2,-1, 596 },
        { 2,-1, 1,-1, 491 },
        { 2, 0,-2, 1,-451 },
        { 0, 0, 3,-1, 439 },
        { 2, 0, 2, 1, 422 },
        { 2, 0,-3,-1, 421 },
        { 2, 1,-1, 1,-366 },
        { 2, 1, 0, 1,-351 },
        { 4, 0, 0, 1, 331 },
        { 2,-1, 1, 1, 315 },
        { 2,-2, 0,-1, 302 },
        { 0, 0, 1, 3,-283 },
        { 2, 1, 1,-1,-229 },
        { 1, 1, 0,-1, 223 },
        { 1, 1, 0, 1, 223 },
        { 0, 1,-2,-1,-220 },
        { 2, 1,-1,-1,-220 },
        { 1, 0, 1, 1,-185 },
        { 2,-1,-2,-1, 181 },
        { 0, 1, 2, 1,-177 },
        { 4, 0,-2,-1, 176 },
        { 4,-1,-1,-1, 166 },
        { 1, 0, 1,-1,-164 },
        { 4, 0, 1,-1, 132 },
        { 1, 0,-1,-1,-119 },
        { 4,-1, 0,-1, 115 },
        { 2,-2, 0, 1, 107 },
    };
    static const int latSize = ARRAY_SIZE(lat);

    int i;
    double bi;          // lat[i].cb, maybe scaled for eccentricity
    double lat_udeg;    // accumulator of terms in latitude (micro-degrees)

    lat_udeg = 0.0;
    for (i = latSize - 1; i >= 0; i--) {
        bi = lat[i].cb;
        if (lat[i].cm != 0) {
            // Multiply terms by E or E^2 as appropriate
            bi *= orb->e;
            if (abs(lat[i].cm) == 2) {
                bi *= orb->e;
            }
        }
        lat_udeg += bi * sin(lat[i].cd * orb->d_rad
                             + lat[i].cm * orb->m_rad
                             + lat[i].cmp * orb->mp_rad
                             + lat[i].cf * orb->f_rad);
    }

    return lat_udeg;
}



LOCAL double riseSetApprox(double             risesetGuess_d,
                           bool               getMoonrise,
                           const Sky_DeltaTs  *deltas,
                           const Sky_SiteProp *site,
                           Sky_SiteHorizon    *topo)
/* Routine to calculate the time of Moon rise or set for the day specified by
   MJDrisesetGuess. The result returned is an approximate value, whose accuracy
   depends upon how close MJDrisesetGuess is to true Moon rise or set time.

   Of course, if this routine is called again with the previous result used as
   the new guess, the new result will be more accurate.
 Returns            - Improved estimate of time of moonrise or moonset, (or zero
                      if the moon does not rise or does not set on this date)
 Inputs
    risesetGuess_d  - MJD of date at which moonrise or moonset time is desired,
                      and time of day of a guess of when moonrise or moonset
                      might be
    getMoonrise     - If true, get moonrise time. If false, get moonset time
    deltas          - delta times, to convert from UTC to TT
    site            - properties of the site. Note: it is assumed that the field
                      site->refracPT will be set to 0.0 before calling this
                      routine, in order to calculate an unrefracted position of
                      the Moon (as per note below).
 Output
    topo            - Topocentric position of Moon at approx. rise or set time

   Moonrise and set calculations assume a standard refraction at the horizon of
   34 arcminutes, and a moon semi-diameter of 16 arcminutes. So the calculation
   is based on the time at which the UNREFRACTED Moon position is at -50 arcmin.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double     ha1_rad;        // Hour angle of Moon at time riseSetGuess_d(rad)
    double     dec_rad;        // Declination of Moon (radian)
    double     ha2_rad;        // Hour angle of Moon at horizon (radian)
    double     cosHa2;         // Cos(Hour Angle) at horizon
    double     riseSetApprox_d;// Improved estimate of rise or set time

    moon_nrelTopocentric(risesetGuess_d, deltas, site, topo);
    sky_siteAzElToHaDec(&topo->rectV, site, &ha1_rad, &dec_rad);

    /* Assuming the Moon's declination remains constant over the period, find
     * out where dec circle intersects the Elevation = -50 arcminutes circle.
     * Of course the declination does not remain constant, so this estimate of
     * Hour Angle is approximate, but it will be a better estimate than
     * ha1_rad.  */
    cosHa2 = (sin(degToRad(-50.0/60.0)) - sin(site->astLat_rad) * sin(dec_rad))
             / (cos(site->astLat_rad) * cos(dec_rad));
    /* If there is no intersection, the Moon either doesn't rise or doesn't set
       on this day at this latitude. (In practice, does this ever happen?) */
    if (fabs(cosHa2) > 1.0) {
        riseSetApprox_d = 0.0;
    } else {
        if (getMoonrise) {
            ha2_rad = -acos(cosHa2);
        } else {
            ha2_rad = acos(cosHa2);
        }
        /* Subtract difference in HA from Solar time (not sidereal time). If
           the difference is large, there is a risk we will return the rise or
           set time for the day before or the day after the one requested. Make
           sure that doesn't happen first. The magic number of 1.03 is a little
           fudge factor that seems to help this routine converge a bit faster.*/
        if ((ha1_rad - ha2_rad) > PI)  { risesetGuess_d += 1.0; }
        if ((ha1_rad - ha2_rad) < -PI) { risesetGuess_d -= 1.0; }
        riseSetApprox_d = risesetGuess_d - (ha1_rad - ha2_rad) / TWOPI * 1.03;
    }
    return riseSetApprox_d;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
