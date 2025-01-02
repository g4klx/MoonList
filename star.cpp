/*==============================================================================
 * star.c - Astronomical conversion routines for stars and other objects
 *                 beyond the Solar System
 *
 * Author:  David Hoadley
 *
 * Description: (see star.h)
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
#include <ctype.h>
#include <errno.h>
#include "instead-of-math.h"            /* for sincos() & normalize() */
#include <limits.h>                     /* for INT_MAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local and project includes */
#include "star.h"

#include "sky1.h"
#include "astron.h"
#include "general.h"
#include "sky.h"
#include "skyio.h"
#include "vectors3d.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       /* For use by REQUIRE() - assertions. */

#define UNUSED_EPOCH_cy     (-99.9)     /* Set this when epoch is meaningless */
#define HALFPI_PLUS_SFA     (HALFPI + SFA)
#define TROP_YEAR           (TROP_CENT / 100.0) /* Length of Tropical year */
#define STRINGEMPTY         (-2)

/*
 * Prototypes for local functions (not called from other modules)
 */
LOCAL int extractObjectName(const char coordStr[],
                            char   nameStr[],
                            size_t nameStrSize,
                            const char **endPtr);
LOCAL int parseEpochStr(const char epochStr[],
                        Star_CoordSys *coordSys,
                        double        *epochT_cy,
                        const char    **endPtr);
LOCAL double myStrtod(const char str[], const char **endPtr, int *error);


#ifdef PREDEF_STANDARD_C_1999
/*          Compiler supports inline functions */
static inline char *saferStrncpy(char *dest, const char *src, size_t len)
/* The standard strncpy() function does not guarantee that the destination
   string will be NUL-terminated. This one does. So use it instead. */
{
    dest[0] = '\0';
    if (len > 0) {
        strncat(dest, src, len - 1);
    }
    return dest;
}

#else
/*          C89/C90 compiler only - no inline functions. Need macros instead */
/*      I swore I would never use the horrid comma operator. Here it is. */
  #define saferStrncpy(dest__, src__, len__)                                \
    ( (dest__)[0] = '\0',                                                   \
        ((len__) > 0) ? strncat((dest__), (src__), (len__) - 1) : (dest__) )
#endif


/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */
/*      Constants found in the 2007 Astronomical Almanac, pages K6 & K7 */
LOCAL const double lightTime_s = 499.0047863852; // time to travel 1 AU(seconds)
LOCAL const double au_km = 1.49597871464e8;      // one Astronomical Unit (km)

/*      Derived constants */
LOCAL const double invC_dpau = lightTime_s / 86400.0;// 1/c (light speed) (d/AU)
LOCAL const double auFactor = 86400.0 / au_km;       // Convert km/s to AU/d


LOCAL Star_CatalogPosn currentObject;   /* For use by star_getApparent() */

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
GLOBAL void star_setCurrentObject(const Star_CatalogPosn *c)
/*! Copies the catalogue information in \a c to internal storage for later use
    by star_getApparent()
 \param[in]   c  The catalogue information for the star, as set by one of the
                 three functions star_parseCoordString(), star_setCatalogPosn()
                 or star_setCatalogOldStyle().
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    currentObject = *c;
}



GLOBAL void star_catalogToApp(const Star_CatalogPosn *c,
                              double                 j2kTT_cy,
                              const Sky1_Nut1980     *nut,
                              V3D_Vector *appV,
                              double     *dist_au)
/*! Convert the catalogue coordinates for a star (or other object outside the
    Solar System) to apparent coordinates at time \a j2kTT_cy.
 \param[in]  c          Catalogue position and motion of object, as set by one
                        of the three functions star_parseCoordString(),
                        star_setCatalogPosn() or star_setCatalogOldStyle().
 \param[in]  j2kTT_cy   Julian centuries since J2000.0, TT timescale
 \param[in]  nut        Nutation angles and obliquity of the ecliptic at
                        time \a j2kTT_cy, as returned by functions
                        sky1_nutationIAU1980() and sky1_epsilon1980()
 \param[out] appV       Unit vector of geocentric apparent direction of
                        object, referred to the true equator and equinox at
                        time (\a j2kTT_cy)
 \param[out] dist_au    Distance to object (AU) as derived from
                        \a c->parallax_rad.

    This routine is simplified, in that:
        - It uses an approximate position for the Earth (calling star_earth())
        - It does not apply relativistic corrections for gravitational light
          deflection; this effect is ignored.
        - Neither does it apply a relativistic annual aberration correction - it
          uses a simple vector sum of the earth's velocity to calculate annual
          aberration (calling star_annAberr()).

 \par References:
    _Astronomical Almanac_ 1987 pages B39-40, or\n
    _Astronomical Almanac_ 2007 pages B66-67
 \note
    If the field \a c->cSys is equal to #FK4, this routine will return zero for
    the vector and distance. FK4 stellar positions are not supported (yet).
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky1_Prec1976   prec;       /* Precession angles */
    V3D_Matrix      cM, bM, aM; /* Used to obtain precession-nutation rotation*/
    V3D_Matrix      npM;        /* Precession/nutation combined matrix */
    V3D_Vector      ebV_au;     /* Barycentric position of Earth (AU) */
    V3D_Vector      ebdotV_aupd;/* Velocity of the Earth (AU/day) */
    V3D_Vector      sbV_au;     /* Barycentric position of Sun (AU) */
    double          tau;        /* Elapsed time for proper motion */
    V3D_Vector      qV;         /* Barycentric direction of object, (unit vector
                                   \b q in _Astronomical Almanac_) */
    V3D_Vector      mV_radpcy;  /* Space motion vector, (radian/Julian century)
                                   (vector \b m in _Astronomical Almanac_) */
    V3D_Vector      p0V;        /* Barycentric direction of celestial object */
    V3D_Vector      pV;         /* Geocentric geometric direction of object */
    V3D_Vector      p2V;        /* Geocentric "proper" direction of object */
    double          scale;

    REQUIRE_NOT_NULL(c);
    REQUIRE_NOT_NULL(appV);
    REQUIRE_NOT_NULL(dist_au);

    /* Obtain combined precession/nutation matrix from catalogue position to
       apparent coordinates. */
    switch (c->cSys) {
    case APPARENT:
    case INTERMEDIATE:
        /* Coordinates are already apparent. All we need to do is convert to
           a vector and return. Assume no diurnal parallax. */
        v3d_polarToRect(appV, c->ra_rad, c->dec_rad);
        *dist_au = 0.0;
        return;
        break;

    case FK4:
        /* Not supported yet. */
        appV->a[0] = appV->a[1] = appV->a[2] = 0.0;
        *dist_au = 0.0;
        return;
        break;

    case FK5:
        /* Create combined precession and nutation matrix */
        sky1_precessionIAU1976(c->eqnxT_cy, j2kTT_cy, &prec);
        sky1_createPrec1976Matrix(&prec, &bM);

        sky1_createNut1980Matrix(nut, &aM);
        v3d_multMxM(&npM, &aM, &bM);
        break;

    case ICRS:
        /* Create combined precession, nutation and frame bias matrix */
        sky1_frameBiasFK5(&aM);
        sky1_precessionIAU1976(c->eqnxT_cy, j2kTT_cy, &prec);
        sky1_createPrec1976Matrix(&prec, &bM);
        v3d_multMxM(&cM, &bM, &aM);

        sky1_createNut1980Matrix(nut, &aM);
        v3d_multMxM(&npM, &aM, &cM);
        break;

    default:
        break;
    }

    /* Obtain the position and velocity of the Earth, referred to the catalogue
       equator and equinox of our object of interest. */
    star_earth(j2kTT_cy, &npM, &ebV_au, &ebdotV_aupd, &sbV_au);

    /* Get catalogue position and space motion as vectors */
    star_catalogToVectors(c, &qV, &mV_radpcy);

    /* Calculate elapsed time from initial to final epoch in Julian centuries
       and apply space motion over that elapsed time */
    tau = j2kTT_cy - c->epochT_cy;
    p0V.a[0] = (qV.a[0] + tau * mV_radpcy.a[0]);
    p0V.a[1] = (qV.a[1] + tau * mV_radpcy.a[1]);
    p0V.a[2] = (qV.a[2] + tau * mV_radpcy.a[2]);

    /* Apply the correction for annual parallax to obtain the (geometric)
       geocentric position from the barycentric position. This is given
       rigorously by the vector sum of its barycentric position and the
       (negative of the) barycentric position of the earth.
            [objGeoV] = [objBarV] - annParallax_rad * [ebV_au] */
    pV.a[0] = (p0V.a[0] - c->parallax_rad * ebV_au.a[0]);
    pV.a[1] = (p0V.a[1] - c->parallax_rad * ebV_au.a[1]);
    pV.a[2] = (p0V.a[2] - c->parallax_rad * ebV_au.a[2]);
#if 0
    printVector("Eb", &ebV_au);
    printVector("Eb'", &ebdotV_aupd);
    printVector("Sb", &sbV_au);
    printMatrix("NP", &npM);
    printVector("q", &qV);
    printVector("m", &mV_radpcy);
    printf("tau = %.9f\n", tau);
    printVector("p0", &p0V);
    printVector("P", &pV);
#endif

    /*      Re-scale vector pV back to unity magnitude */
    scale = 1.0 / v3d_magV(&pV);
    pV.a[0] *= scale;
    pV.a[1] *= scale;
    pV.a[2] *= scale;

#ifdef RUN_RIDICULOUSLY_RIGOROUS_RELATIVISTIC_ROUTINE
    star_lightDeflection(&pV, &ebV_au, &ebdotV_aupd, &sbV_au, &p2V);
#else
    /* This version does not calculate vector p1 (light deflection) */

    /* Apply simple annual aberration correction to convert geometric position
       to "proper" position. */
    star_annAberr(&pV, &ebdotV_aupd, &p2V);
#endif
#if 0
    printVector("p", &pV);
    printVector("p2", &p2V);
#endif

    /* Calculate distance to object */
    if (c->parallax_rad <= 0.0) {
        *dist_au = 0.0;
    } else {
        *dist_au = RAD2ARCSEC / c->parallax_rad;
    }

    /* Convert to geocentric apparent coordinates by multiplying
       by matrices for precession and nutation */
    v3d_multMxV(appV, &npM, &p2V);
}



GLOBAL void star_getApparent(double j2kTT_cy, Sky_TrueEquatorial *pos)
/*! Calculate the position of the currently selected star (or other object
    outside the solar system) as a unit vector and a distance, in apparent
    coordinates. It calls star_catalogToApp() to obtain the star's position,
    after  having called sky1_nutationIAU1980() and sky1_epsilon1980() to obtain
    the necessary nutation terms.

    The star whose coordinates are obtained with this function is the star most
    recently specified with star_setCurrentObject()
 \param[in]  j2kTT_cy   Julian centuries since J2000.0, TT timescale
 \param[out] pos        Timestamped structure containing position data in
                        apparent coordinates and the equation of the equinoxes.

    This function is designed to be callable by the skyfast_init() and
    skyfast_backgroundUpdate() functions in a tracking application.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky1_Nut1980    nut;

    REQUIRE_NOT_NULL(pos);

    if (currentObject.cSys == INTERMEDIATE) {
        pos->eqEq_rad = 0.0;

    } else {
        /* Calculate nutation, the mean obliquity of the ecliptic and
           the equation of the equinoxes */
        sky1_nutationIAU1980(j2kTT_cy, 0, &nut);
        sky1_epsilon1980(j2kTT_cy, &nut);
        pos->eqEq_rad = nut.eqEq_rad;
    }

    star_catalogToApp(&currentObject, j2kTT_cy, &nut,
                      &pos->appCirsV, &pos->distance_au);

    /* Now set the timestamp*/
    pos->timestamp_cy = j2kTT_cy;
}



GLOBAL void star_getTopocentric(double             j2kUtc_d,
                                const Sky_DeltaTs  *deltas,
                                const Sky_SiteProp *site,
                                Sky_SiteHorizon *topo)
/*! Calls star_getApparent() to calculate the star's position in apparent
    coordinates, and then converts this to topocentric horizon coordinates at
    the specified site.

    The star whose coordinates are obtained with this function is the star most
    recently specified with star_setCurrentObject()
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
    Use this function if you are calculating the star topocentric position once,
    for a single site. But if you are going to be calculating it repeatedly, or
    for multiple sites, use of this function will cause you to perform a great
    many needless recalculations. Use skyfast_getApprox(), followed by
    sky1_appToTirs() and sky_siteTirsToTopo() instead.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_Times          atime;   // time, in various timescales
    Sky_TrueEquatorial pos;     // geocentric position of the Sun and distance
    V3D_Vector         terInterV; // unit vector in Terrestrial Intermed Ref Sys

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(site);
    REQUIRE_NOT_NULL(topo);

    sky_updateTimes(j2kUtc_d, deltas, &atime);

    star_getApparent(atime.j2kTT_cy, &pos);

    /* Convert apparent position to topocentric Azimuth/Elevation coords */
    sky1_appToTirs(&pos.appCirsV, atime.j2kUT1_d, pos.eqEq_rad, &terInterV);
    sky_siteTirsToTopo(&terInterV, pos.distance_au, site, topo);
}



GLOBAL void star_catalogToVectors(const Star_CatalogPosn *c,
                                  V3D_Vector *pV,
                                  V3D_Vector *vV_radpcy)
/*! This routine takes the mean equatorial polar coordinate of an object and its
    proper motions, annual parallax (distance) and radial velocity, and computes
    the rectangular position and velocity vectors with respect to the same
    coordinate system.
 \param[in]  c              Catalog position and motion of object

 \param[out] pV             Unit position vector (direction cosines)
 \param[out] vV_radpcy      Space velocity vector (radian/Julian century)

 \par Reference:
    _Astronomical Almanac_ 2007, pages B66 and B67
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double sinRa, cosRa, sinDec, cosDec;

    REQUIRE_NOT_NULL(c);
    REQUIRE_NOT_NULL(pV);
    REQUIRE_NOT_NULL(vV_radpcy);

    sincos(c->ra_rad, &sinRa, &cosRa);
    sincos(c->dec_rad, &sinDec, &cosDec);

    /* Form the position vector */
    pV->a[0] = cosRa * cosDec;
    pV->a[1] = sinRa * cosDec;
    pV->a[2] = sinDec;

    /* Form the velocity vector */
    /* vV = muRA * ∂([pV])/∂(RA) + muDec * ∂([pV])/∂(Dec)
     *                                           + parallax * radVel * [pV] */
    if (c->muRaInclCosDec) {
        // The cos(dec) term is already included in the figure for proper
        // motion in RA
        vV_radpcy->a[0] = -c->muRA_radpcy * sinRa
                         - c->muDec_radpcy * cosRa * sinDec
                         + c->parallax_rad * c->radVel_aupcy * pV->a[0];
        vV_radpcy->a[1] =  c->muRA_radpcy * cosRa
                         - c->muDec_radpcy * sinRa * sinDec
                         + c->parallax_rad * c->radVel_aupcy * pV->a[1];
    } else {
        // cos(dec) term was not included. Needs to be factored in here
        vV_radpcy->a[0] = -c->muRA_radpcy * pV->a[1]
                         - c->muDec_radpcy * cosRa * sinDec
                         + c->parallax_rad * c->radVel_aupcy * pV->a[0];
        vV_radpcy->a[1] =  c->muRA_radpcy * pV->a[0]
                         - c->muDec_radpcy * sinRa * sinDec
                         + c->parallax_rad * c->radVel_aupcy * pV->a[1];
    }

    vV_radpcy->a[2] =  c->muDec_radpcy * cosDec
                         + c->parallax_rad * c->radVel_aupcy * pV->a[2];
}



GLOBAL void star_earth(double           t1_cy,
                       const V3D_Matrix *npM,
                       V3D_Vector *ebV_au,
                       V3D_Vector *ebdotV_aupd,
                       V3D_Vector *sbV_au)
/*! Calculate the position and velocity of the Earth, for use in making annual
 *  parallax and aberration corrections.
 \param[in]  t1_cy        Epoch of time of observation (Julian centuries since
                          J2000.0, TT timescale)
 \param[in]  npM          Combined precession and nutation matrix \b NP that
                          will be used to convert the coordinates of the star
                          of interest from its catalogue mean place to
                          apparent coordinates. It will be used in this
                          routine to do the reverse - convert the apparent
                          coordinates of the Earth's position back into
                          mean coordinates at the catalogue epoch of the star.
 \param[out] ebV_au       Barycentric position of the Earth (AU) (vector \b Eb
                          in _Astronomical Almanac_)
 \param[out] ebdotV_aupd  Velocity of the Earth (AU/day) (vector \b Ėb
                          in _Astronomical Almanac_)
 \param[out] sbV_au       Barycentric position of the Sun (AU) (vector \b Sb
                          in _Astronomical Almanac_)

 *  In theory, this is the barycentric
 *  position, and the values are referred to the same system as is used for the
 *  catalogue position of the star (such as J2000, or
 *  the Barycentric Celestial Reference System/ICRS). This can be obtained from
 *  the IAU SOFA routine \c iauEpv00(), calculated in fabulously complicated
 *  detail. But the parallax and aberration corrections are small, which means
 *  that the accuracy required of this routine is not actually that high.
 *
 *  So to save on computing resources in a tracking application, this routine
 *  uses instead the simple algorithm in the _Astronomical Almanac_
 *  for the Sun's position, entitled "Low-precision formulas for the Sun's
 *  coordinates..." It differentiates the position vector to obtain the Earth's
 *  velocity. So there are a few approximations here:
 *  1. It calculates a heliocentric position, rather than a barycentric position
 *  2. The algorithm itself is approximate - accurate to 0.01° between 1950 and
 *     2050.
 *  3. It does not calculate the barycentric position of the Sun at all. So
 *     \a sbV_au will be returned as a zero vector.
 *
 *  But nonetheless, the effect on the final star position is very small. For
 *  example, the _Almanac_ gives an example of stellar reduction calculations
 *  for a fictitious star, (which appears to be very close to Alpha Centauri).
 *  This star has about as much parallax as one is ever likely to see.
 *
 *  In the 2007 _Almanac_, page B68, the apparent position of this star at
 *  2007 January 1, 0h TT , after performing a rigorous reduction, is
 *  -   RA = 14 40 03.4343 (hms), Dec = -60°51′37.770″
 *
 *  The position obtained by using the approximations of star_catalogToApp()
 *  (which calls this routine) and converts to apparent coordinates is
 *  -   RA = 14:40:03.4276 (hms), Dec = -60°51'37.782"
 *
 *  which is clearly good enough for tracking.

 *  In the 1987 _Almanac_, page B41, the apparent position of a (very slightly
 *  different) fictitious star at 1987 January 1, 0h TT is
 *  -   RA = 14 38 40.164 (hms), Dec = -60°46′44.82″
 *
 *  The position obtained using the routines here (as above) is
 *  -   RA = 14:38:40.157 (hms), Dec = -60°46'44.84"
 *
 * For a more detailed look at this, see \ref page-stellar-reduction-accuracy

 \par References:
    _Astronomical Almanac_ 1987 pages
        B39 (method), B41 (example) and C24 (Sun position) or\n
    _Astronomical Almanac_ 2007 pages
        B66 (method), B68 (example) and C24 (Sun position).
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    static const double g1_rad   = 0.9856003 * DEG2RAD;
    static const double l1_rad   = 0.9856474 * DEG2RAD;
    static const double lam1_rad = 1.915 * DEG2RAD;
    static const double lam2_rad = 0.020 * DEG2RAD;
    static const double r1       = 0.01671;
    static const double r2       = 0.00014;

    double L_rad;           // Mean longitude, corrected for aberration (radian)
    double g_rad;           // Mean anomaly (radian)
    double lambda_rad;      // Ecliptic longitude (radian)
    double epsilon_rad;     // Obliquity of ecliptic (radian)
    double sing, cosg;      // sine and cosine of g
    double sin2g, cos2g;    // sine and cosine of 2g
    double r_au;            // Sun-Earth distance (AU)
    double n;               // days since J2000.0
    double lambdadot;       // differentiation of lambda_rad
    double rdot;            // differentiation of r_au

    V3D_Vector posV;
    V3D_Vector velV;

    REQUIRE_NOT_NULL(npM);
    REQUIRE_NOT_NULL(ebV_au);
    REQUIRE_NOT_NULL(ebdotV_aupd);
    REQUIRE_NOT_NULL(sbV_au);

    n = t1_cy * JUL_CENT;
    // The following lines are taken from the Astronomical Almanac 2007 p C24
    // "Low precision formulas for the Sun's coordinates and the equat. of time"
    L_rad = normalize(degToRad(280.461) + (l1_rad * n), TWOPI);
    g_rad = normalize(degToRad(357.529) + (g1_rad * n), TWOPI);

    sincos(g_rad, &sing, &cosg);
    sincos(g_rad + g_rad, &sin2g, &cos2g);

    lambda_rad = L_rad + lam1_rad * sing + lam2_rad * sin2g;
    epsilon_rad = degToRad(23.439 - (0.0000004 * n));
    r_au = 1.00014 - r1 * cosg - r2 * cos2g;

    // Earth position vector is negative of vector to sun, which is why each of
    // the following terms is negated
    posV.a[0] = -r_au * cos(lambda_rad);
    posV.a[1] = -r_au * cos(epsilon_rad) * sin(lambda_rad);
    posV.a[2] = -r_au * sin(epsilon_rad) * sin(lambda_rad);
    // Rotate back to mean coordinates - using inverse of NP matrix
    v3d_multMtransxV(ebV_au, npM, &posV);

    // from the above data, calculate a velocity vector
    lambdadot = l1_rad + lam1_rad * g1_rad * cosg
                 + 2.0 * lam2_rad * g1_rad * cos2g;
    rdot = -r1 * g1_rad * sing - 2.0 * r2 * sin2g;

    velV.a[0] = r_au * sin(lambda_rad) * lambdadot - cos(lambda_rad) * rdot;
    velV.a[1] = -cos(epsilon_rad) * (r_au * cos(lambda_rad) * lambdadot
                                     + sin(lambda_rad) * rdot);
    velV.a[2] = -sin(epsilon_rad) * (r_au * cos(lambda_rad) * lambdadot
                                     + sin(lambda_rad) * rdot);
    // Rotate back to mean coordinates - using inverse of NP matrix
    v3d_multMtransxV(ebdotV_aupd, npM, &velV);

    // Barycentric Sun position.
    // Of course with this simple algorithm we are not actually calculating this
    sbV_au->a[0] = 0.0;
    sbV_au->a[1] = 0.0;
    sbV_au->a[2] = 0.0;
}



GLOBAL void star_annAberr(const V3D_Vector *p1V,
                          const V3D_Vector *earthVelV_aupd,
                          V3D_Vector *p2V)
/*! Apply the annual aberration correction, to convert an object's coordinates
    from geocentric geometric direction to "proper" direction (still geocentric)
 \param[in]  p1V            Unit vector pointing to object's geometric direction
                            as viewed from the centre of the earth.
 \param[in]  earthVelV_aupd Earth velocity vector, referred to the J2000 mean
                            equator and equinox, as returned by routine
                            star_earth() (AU/day - not a unit vector)
 \param[out] p2V            Unit vector pointing to object's proper
                            direction
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Vector aberrV;

    REQUIRE_NOT_NULL(earthVelV_aupd);
    REQUIRE_NOT_NULL(p2V);

    aberrV.a[0] = earthVelV_aupd->a[0] * invC_dpau;
    aberrV.a[1] = earthVelV_aupd->a[1] * invC_dpau;
    aberrV.a[2] = earthVelV_aupd->a[2] * invC_dpau;

    *p2V = *p1V;      // copy vector
    v3d_addToUVfast(p2V, &aberrV);
}



GLOBAL char *star_equinoxToStr(const Star_CatalogPosn *coordBlock,
                               char   equinoxStr[],
                               size_t eqnxStrSize)
/*! Write the coordinate system and equinox out in a standard form.
 \param[in]  coordBlock  Coordinate block describing the celestial object
 \param[out] equinoxStr  String to contain the coordinate system and equinox
 \param[in]  eqnxStrSize Size of equinoxStr (i.e. maximum number of characters
                         available for the output

 *  The output form depends upon the value of \a coordBlock->cSys:
 *  |\a coordBlock->cSys| Output                |
 *  |-------------------|-----------------------|
 *  | #APPARENT         | Apparent              |
 *  | #INTERMEDIATE     | Intermediate          |
 *  | #FK4              | Bnnnn.n (n is numeric)|
 *  | #FK5              | Jnnnn.n (n is numeric)|
 *  | #ICRS             | ICRS                  |
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int epoch_ax10;         /* Integer form of year, multiplied by 10 */

    REQUIRE_NOT_NULL(coordBlock);
    REQUIRE_NOT_NULL(equinoxStr);
    REQUIRE(eqnxStrSize < INT_MAX); /* Was a -ve value passed to this param? */

    switch (coordBlock->cSys) {
    case APPARENT:
        saferStrncpy(equinoxStr, "Apparent", eqnxStrSize);
        break;

    case INTERMEDIATE:
        saferStrncpy(equinoxStr, "Intermediate", eqnxStrSize);
        break;

    case FK4:
        /* Write out year to 1 decimal place without using the %f formatting
           specifier to snprintf(). (Many embedded systems don't support it) */
        epoch_ax10 = (int)( ((coordBlock->eqnxT_cy * JUL_CENT
                              + MJD_J2000 - MJD_B1900) * 10.0 / TROP_YEAR)
                           + 19000.5);
        snprintf(equinoxStr, eqnxStrSize, "B%05d ", epoch_ax10);
        if (eqnxStrSize > 7) { equinoxStr[6] = equinoxStr[5]; }
        if (eqnxStrSize > 6) { equinoxStr[5] = '.'; }
        break;

    case FK5:
        /* As above, don't use %f */
        epoch_ax10 = (int)(coordBlock->eqnxT_cy * 1000.0 + 20000.5);
        snprintf(equinoxStr, eqnxStrSize, "J%05d ", epoch_ax10);
        if (eqnxStrSize > 7) { equinoxStr[6] = equinoxStr[5]; }
        if (eqnxStrSize > 6) { equinoxStr[5] = '.'; }
        break;

    case ICRS:
        saferStrncpy(equinoxStr, "ICRS", eqnxStrSize);
        break;

    default:
        equinoxStr[0] = '\0';
        break;
    }
    return equinoxStr;
}



GLOBAL char *star_epochToStr(const Star_CatalogPosn *coordBlock,
                             char   epochStr[],
                             size_t epochStrSize)
/*! Write the epoch used for proper motion calculations out in a standard form.
 \param[in]  coordBlock   Coordinate block describing the celestial object
 \param[out] epochStr     String to contain the coordinate system and epoch
 \param[in]  epochStrSize Size of epochStr (i.e. maximum number of characters
                          available for the output

 *  The output form depends upon the value of \a coordBlock->cSys:
 *  |\a coordBlock->cSys| Output                              |
 *  |-------------------|-------------------------------------|
 *  | #APPARENT         | (empty string - epoch is irrelevant)|
 *  | #INTERMEDIATE     | (empty string - epoch is irrelevant)|
 *  | #FK4              | Bnnnn.n (n is numeric)              |
 *  | #FK5              | Jnnnn.n (n is numeric)              |
 *  | #ICRS             | Jnnnn.n (n is numeric)              |
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double  epochT_cy;
    int epoch_ax10;         /* Integer form of year, multiplied by 10 */

    REQUIRE_NOT_NULL(coordBlock);
    REQUIRE_NOT_NULL(epochStr);
    REQUIRE(epochStrSize < INT_MAX); /* Was a -ve value passed to this param? */

    if (coordBlock->epochSpecified) {
        epochT_cy = coordBlock->epochT_cy;
    } else {
        epochT_cy = coordBlock->eqnxT_cy;
    }
    switch (coordBlock->cSys) {
    case APPARENT:
    case INTERMEDIATE:
        epochStr[0] = '\0';
        break;

    case FK4:
        /* Write out year to 1 decimal place without using the %f formatting
           specifier to snprintf(). (Many embedded systems don't support it) */
        epoch_ax10 = (int)(( (epochT_cy * JUL_CENT + MJD_J2000 - MJD_B1900)
                            * 10.0 / TROP_YEAR) + 19000.5);
        snprintf(epochStr, epochStrSize, "B%05d ", epoch_ax10);
        if (epochStrSize > 7) { epochStr[6] = epochStr[5]; }
        if (epochStrSize > 6) { epochStr[5] = '.'; }
        break;

    case FK5:
    case ICRS:
        /* As above, don't use %f */
        epoch_ax10 = (int)(epochT_cy * 1000.0 + 20000.5);
        snprintf(epochStr, epochStrSize, "J%05d ", epoch_ax10);
        if (epochStrSize > 7) { epochStr[6] = epochStr[5]; }
        if (epochStrSize > 6) { epochStr[5] = '.'; }
        break;

    default:
        epochStr[0] = '\0';
        break;
    }
    return epochStr;
}



GLOBAL int star_setCatalogPosn(const char objectName[],
                               double ra_h,
                               double dec_deg,
                               Star_CoordSys coordSys,
                               double equinoxYear,
                               double epochYear,
                               double muRaSplat_maspa,
                               double muDec_maspa,
                               double annParallax_as,
                               double radVel_kmps,
                               Star_CatalogPosn *coordBlock)
/*! This is one of three alternative routines for loading the information for a
    celestial object into a block of form Star_CatalogPosn. (The other two
    routines are star_setCatalogOldStyle() and star_parseCoordString().)
 \returns                       An error code containing one of the values in
                                Star_CoordErrors
 \param[in]  objectName         Name of the object (optional)
 \param[in]  ra_h               Right ascension (hours)
 \param[in]  dec_deg            Declination (degrees)
 \param[in]  coordSys           Coordinate system for RA and Dec
 \param[in]  equinoxYear        Year of epoch of the mean equinox (if coordSys
                                is #FK4 or #FK5)
 \param[in]  epochYear          Year of epoch for proper motion (if coordSys is
                                #ICRS, #FK4 or #FK5)
 \param[in]  muRaSplat_maspa    μα* - proper motion in RA, including the cos(δ)
                                factor (milliarcseconds/year)
 \param[in]  muDec_maspa        μδ - proper motion in declination
                                (milliarcseconds/year)
 \param[in]  annParallax_as     π - annual parallax (arcseconds)
 \param[in]  radVel_kmps        ν - radial velocity (km/s), positive away from
                                the Earth
 \param[out] coordBlock         The block containing the data from all the input
                                parameters, basically scaled into radians

    The difference between this routine and star_setCatalogOldStyle() is in the
    way the proper motions in RA and Dec are specified. In this routine, μα* is
    specified in milliarcseconds per year, after scaling to absolute
    displacement on the sky by multiplying by cos(δ). μδ is also specified in
    milliarcseconds per year.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(coordBlock);

    saferStrncpy(coordBlock->objectName,
                 objectName,
                 sizeof(coordBlock->objectName));
    coordBlock->ra_rad = hrsToRad(ra_h);
    coordBlock->dec_rad = degToRad(dec_deg);

    coordBlock->cSys = coordSys;
    switch (coordBlock->cSys) {
    case APPARENT:
    case INTERMEDIATE:
        coordBlock->eqnxT_cy = UNUSED_EPOCH_cy;
        coordBlock->epochT_cy = UNUSED_EPOCH_cy;
        break;
    case FK4:
        coordBlock->eqnxT_cy = ((equinoxYear - 1900.0) * TROP_YEAR
                                + (MJD_B1900 - MJD_J2000)) / JUL_CENT;
        coordBlock->epochT_cy = ((epochYear - 1900.0) * TROP_YEAR
                                 + (MJD_B1900 - MJD_J2000)) / JUL_CENT;
        break;
    case FK5:
        coordBlock->eqnxT_cy = (equinoxYear - 2000.0) / 100.0;
        coordBlock->epochT_cy = (epochYear - 2000.0) / 100.0;
        break;
    case ICRS:
        coordBlock->eqnxT_cy = 0.0;
        coordBlock->epochT_cy = (epochYear - 2000.0) / 100.0;
        break;
    default:
        break;
    }

    coordBlock->muRA_radpcy = arcsecToRad(muRaSplat_maspa * 0.1);
    coordBlock->muDec_radpcy = arcsecToRad(muDec_maspa * 0.1);
    coordBlock->parallax_rad = arcsecToRad(annParallax_as);
    coordBlock->radVel_aupcy = radVel_kmps * auFactor * JUL_CENT;
    coordBlock->muRaInclCosDec = true;
    coordBlock->epochSpecified = (fabs(coordBlock->eqnxT_cy
                                       - coordBlock->epochT_cy) > SFA);

    /* A bit of basic sanity checking */
    if ((coordBlock->ra_rad < 0.0) || (coordBlock->ra_rad >= TWOPI)) {
        return STAR_RARANERR;
    }
    if (fabs(coordBlock->dec_rad) > HALFPI_PLUS_SFA) {
        return STAR_DECRANERR;
    }
    if (   ((coordBlock->cSys == FK4) || (coordBlock->cSys == FK5))
        && ((equinoxYear < 0.0) || (equinoxYear >  9999.9)))
    {
        return STAR_IVEPOCH;
    }
    if (   (   (coordBlock->cSys == ICRS) || (coordBlock->cSys == FK4)
            || (coordBlock->cSys == FK5))
        && ((equinoxYear < 0.0) || (equinoxYear >  9999.9)))
    {
        return STAR_IVEPOCH;
    }
    return STAR_NORMAL;
}



GLOBAL int star_setCatalogOldStyle(const char objectName[],
                                   double ra_h,
                                   double dec_deg,
                                   Star_CoordSys coordSys,
                                   double equinoxYear,
                                   double epochYear,
                                   double muRa_spcy,
                                   double muDec_aspcy,
                                   double annParallax_as,
                                   double radVel_kmps,
                                   Star_CatalogPosn *coordBlock)
/*! This is one of three alternative routines for loading the information for a
    celestial object into a block of form Star_CatalogPosn. (The other two
    routines are star_setCatalogPosn() and star_parseCoordString().)
 \returns                       An error code containing one of the values in
                                Star_CoordErrors
 \param[in]  objectName         Name of the object (optional)
 \param[in]  ra_h               Right ascension (hours)
 \param[in]  dec_deg            Declination (degrees)
 \param[in]  coordSys           Coordinate system for RA and Dec
 \param[in]  equinoxYear        Year of epoch of the mean equinox (if coordSys
                                is #FK4 or #FK5)
 \param[in]  epochYear          Year of epoch for proper motion (if coordSys is
                                #ICRS, #FK4 or #FK5)
 \param[in]  muRa_spcy          μα - proper motion in RA, \b not including the
                                cos(δ) factor (seconds/century)
 \param[in]  muDec_aspcy        μδ - proper motion in declination
                                (arcseconds/century)
 \param[in]  annParallax_as     π - annual parallax (arcseconds)
 \param[in]  radVel_kmps        ν - radial velocity (km/s), positive away from
                                the Earth
 \param[out] coordBlock         The block containing the data from all the input
                                parameters, basically scaled into radians

    The difference between this routine and star_setCatalogPosn() is in the way
    the proper motions in RA and Dec are specified. In this routine, μα is
    specified in seconds (of time) per century, and is not scaled to absolute
    displacement on the sky. μδ is specified in arcseconds per century.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(coordBlock);

    saferStrncpy(coordBlock->objectName,
                 objectName,
                 sizeof(coordBlock->objectName));
    coordBlock->ra_rad = hrsToRad(ra_h);
    coordBlock->dec_rad = degToRad(dec_deg);

    coordBlock->cSys = coordSys;
    switch (coordBlock->cSys) {
    case APPARENT:
    case INTERMEDIATE:
        coordBlock->eqnxT_cy = UNUSED_EPOCH_cy;
        coordBlock->epochT_cy = UNUSED_EPOCH_cy;
        break;
    case FK4:
        coordBlock->eqnxT_cy = ((equinoxYear - 1900.0) * TROP_YEAR
                                + (MJD_B1900 - MJD_J2000)) / JUL_CENT;
        coordBlock->epochT_cy = ((epochYear - 1900.0) * TROP_YEAR
                                 + (MJD_B1900 - MJD_J2000)) / JUL_CENT;
        break;
    case FK5:
        coordBlock->eqnxT_cy = (equinoxYear - 2000.0) / 100.0;
        coordBlock->epochT_cy = (epochYear - 2000.0) / 100.0;
        break;
    case ICRS:
        coordBlock->eqnxT_cy = 0.0;
        coordBlock->epochT_cy = (epochYear - 2000.0) / 100.0;
        break;
    default:
        break;
    }

    coordBlock->muRA_radpcy  = secToRad(muRa_spcy);
    coordBlock->muDec_radpcy = arcsecToRad(muDec_aspcy);
    coordBlock->parallax_rad = arcsecToRad(annParallax_as);
    coordBlock->radVel_aupcy = radVel_kmps * auFactor * JUL_CENT;
    coordBlock->muRaInclCosDec = false;
    coordBlock->epochSpecified = (fabs(coordBlock->eqnxT_cy
                                       - coordBlock->epochT_cy) > SFA);

    /* A bit of basic sanity checking */
    if ((coordBlock->ra_rad < 0.0) || (coordBlock->ra_rad >= TWOPI)) {
        return STAR_RARANERR;
    }
    if (fabs(coordBlock->dec_rad) > HALFPI_PLUS_SFA) {
        return STAR_DECRANERR;
    }
    if (   ((coordBlock->cSys == FK4) || (coordBlock->cSys == FK5))
        && ((equinoxYear < 0.0) || (equinoxYear >  9999.9)))
    {
        return STAR_IVEPOCH;
    }
    if (   (   (coordBlock->cSys == ICRS) || (coordBlock->cSys == FK4)
            || (coordBlock->cSys == FK5))
        && ((equinoxYear < 0.0) || (equinoxYear >  9999.9)))
    {
        return STAR_IVEPOCH;
    }
    return STAR_NORMAL;
}



GLOBAL int star_parseCoordString(const char inStr[],
                                 Star_CatalogPosn *coordBlock)
/*! This is one of three alternative routines for loading the information for a
    celestial object into a block of form Star_CatalogPosn. (The other two
    routines are star_setCatalogPosn() and star_setCatalogOldStyle().)
 \returns               An error code containing one of the values in
                        Star_CoordErrors
 \param[in]  inStr      The input string.
 \param[out] coordBlock The block containing the data from the input
                        string, basically scaled into radians

    This routine decodes an input string which is of the form
 \verbatim
    [name], RA, DEC[, Coords,[Epoch][, Mu_ra*, Mu_dec[, Pi[, Vr]]]]
 \endverbatim
    where the items shown inside square brackets are optional. This is basically
    a comma separated value (CSV) format. The text fields
    are as follows:
    - \b name   - Object name. If the name is longer than the field
                  Star_CatalogPosn.objectName, the extra characters will be
                  lost.
    - \b RA     - Right Ascension. This has one to three numeric subfields,
                  depending upon where the decimal point (if any) is placed. So
                  it may be in the form of
                + decimal hours (e.g. 12.582) (one subfield only)
                + hours and decimal minutes (e.g. 12:34.933 or 12 34.933) (two
                  subfields)
                + hours, minutes and seconds (e.g. 12:34:56 or 12 34 56, or
                  12:34:56.78 or 12:34:56:78)
    - \b DEC    - Declination. This also has one to three numeric subfields, as
                  above. it may be in the form of
                + decimal degrees (e.g. ±21.625 or ±21.625°)
                + degrees and decimal arcminutes (e.g. ±21 37.5  or ±21°37.5′)
                + degrees, arcminutes and arcseconds (e.g. ±21 37 30  or
                  ±21°37′30″, or ±21 37 30.0  or ±21°37′30.0″)

    In both of the above cases, the subfields may be separated by spaces,
    colons, degree, minute and seconds symbols, single and double quotes,
    or anything really. The parser is not fussy, and does not check.
    - \b Coords - the coordinate system reference for the coordinates. Valid
                  values for this field are
                + \b Apparent - referred to true equinox and equator of date
                + \b Jnnnn.n - FK5 position referred to mean equator and equinox
                  of the specified Julian epoch (e.g. J2000)
                + \b Bnnnn.n - FK4 position referred to mean equator and equinox
                  of the specified Besselian epoch (e.g. B1950)
                + \b nnnn.n - If not specified with a B or a J, values here
                  before 1984.0 will be assumed to be FK4 (Besselian) and those
                  after to be FK5 (Julian)
                + \b Intermediate - Celestial Intermediate Reference System.
                  This is referred to the true equator of date and the Celestial
                  Intermediate Origin (CIO).
                + \b ICRS - International Celestial Reference System
                .
                If this field is not specified, it will be assumed to be
                Apparent, and all following fields will be ignored.
    - \b Epoch  - Initial epoch for proper motion. If this field is omitted, the
                  epoch is assumed to be the same as the coordinate system
                  specification that precedes it. If that coordinate system is
                  ICRS, then the epoch is assumed to be J2000.0 if is not
                  specified. If the coordinate system is Apparent or
                  Intermediate, the Epoch is ignored.
    - \b Mu_ra* - proper motion in right ascension (milliarcseconds/year),
                  including the cos(DEC) factor.
                  This field will be ignored if Coords is set to Apparent or
                  Intermediate
    - \b Mu_dec - proper motion in declination (milliarcseconds/year).
                  This field will be ignored if Coords is set to Apparent or
                  Intermediate
    - \b Pi     - annual parallax (arcseconds)
                  This field will be ignored if Coords is set to Apparent or
                  Intermediate
    - \b Vr     - radial velocity (km/s), positive for objects moving away.
                  This field will be ignored if Coords is set to Apparent or
                  Intermediate

    The string format is fairly flexible. For example, fields may be separated
    with spaces rather than commas, so long as any omitted fields are at the end
    of the string. But if any fields not at the end are omitted, commas must be
    used to indicate the empty fields.
 \note
    The facilities provided by the C language for text handling are rudimentary,
    and as a result, this routine is slightly fragile. In trying to be flexible
    but not too complex, the following limitation applies: Do not ever have
    white space immediately preceding a comma. If you do, the routine may go
    wrong. So, in particular this means that omitted fields (not at the end of
    the string) must be denoted by two successive commas with no space in
    between.
 \verbatim
    Dummy, 12:34:56.789, +89°56′43.210″, J2000.0,, 23.455, 12.766
    Barnard's star, 17 57 48.500 +04 41 36.111 ICRS,, -802.803, 10362.542, 0.5474506, -110.353
 \endverbatim
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int         decodeError;
    double      tempVal;            /* temporary value, decoded from string */
    const char  *startPtr;
    const char  *endPtr;
    bool        muRAsupplied = false;
    bool        muDecSupplied = false;
    bool        piSupplied = false;
    Star_CoordSys dummyCoordSys;

    REQUIRE_NOT_NULL(inStr);
    REQUIRE_NOT_NULL(coordBlock);

    /* Start by assuming that we will be given a minimal string - just RA and
       Dec in Apparent coordinates. */
    memset(coordBlock, 0, sizeof(*coordBlock));
    coordBlock->eqnxT_cy = UNUSED_EPOCH_cy;
    coordBlock->epochT_cy = UNUSED_EPOCH_cy;

    startPtr = inStr;

    /* Extract the object name */
    decodeError = extractObjectName(startPtr,
                                    coordBlock->objectName,
                                    sizeof(coordBlock->objectName),
                                    &endPtr);
    if (decodeError != STAR_NORMAL) {
        return decodeError;
    }

    if (*endPtr == ',') { endPtr++; }
    startPtr = endPtr;

    /* Decode the right ascension */
    tempVal = skyio_sxStrToAng(startPtr, &endPtr, &decodeError);
    if (decodeError == 0) {
        /* RA is specified in hours. Convert to radians */
        coordBlock->ra_rad = hrsToRad(tempVal);
    } else if (decodeError == INVALID_ANGLE) {
        return STAR_ERRINRA;
    } else {
        return STAR_NORA;
    }
    if ((coordBlock->ra_rad < 0.0) || (coordBlock->ra_rad >= TWOPI)) {
        return STAR_RARANERR;
    }

    /* Decode the declination */
    startPtr = endPtr;
    tempVal = skyio_sxStrToAng(startPtr, &endPtr, &decodeError);
    if (decodeError == 0) {
        /* Dec is specified in degrees. Convert to radians */
        coordBlock->dec_rad = degToRad(tempVal);
    } else if (decodeError == INVALID_ANGLE) {
        return STAR_ERRINDEC;
    } else {
        return STAR_NODEC;
    }
    if (fabs(coordBlock->dec_rad) > HALFPI_PLUS_SFA) {
        return STAR_DECRANERR;
    }

    /* Decode the equinox/coordinate system specification */
    startPtr = endPtr;
    decodeError = parseEpochStr(startPtr,
                                &coordBlock->cSys,
                                &coordBlock->eqnxT_cy,
                                &endPtr);
    if (decodeError == STRINGEMPTY) {
        coordBlock->cSys = APPARENT;
        coordBlock->eqnxT_cy = UNUSED_EPOCH_cy;
        coordBlock->epochT_cy = UNUSED_EPOCH_cy;

    } else if (decodeError != 0) {
        return decodeError;
    } else {
        ;
    }
    if ((coordBlock->cSys == APPARENT) || (coordBlock->cSys == INTERMEDIATE)) {
        return STAR_NORMAL;
    }

    /* Decode the epoch specification (for proper motion) if it is present */
    if (*endPtr == ',') { endPtr++; }
    startPtr = endPtr;
    /* Is there a left parenthesis at this point? If so, look for a right
       parenthesis and decode the contents in between as an epoch specification
     */
    decodeError = parseEpochStr(startPtr,
                                &dummyCoordSys,
                                &coordBlock->epochT_cy,
                                &endPtr);
    if (decodeError == STRINGEMPTY) {
        coordBlock->epochSpecified = false;
        coordBlock->epochT_cy = coordBlock->eqnxT_cy;

    } else if (decodeError != 0) {
        return decodeError;
    } else {
        coordBlock->epochSpecified = true;
    }

    /* Are there any fields left? If so, we expect at least two items i.e.
       proper motion in RA and Dec. */
    if (*endPtr == ',') { endPtr++; }
    startPtr = endPtr;
    if (startPtr[0] != '\0') {
        tempVal = myStrtod(startPtr, &endPtr, &decodeError);
        if (decodeError == STRINGEMPTY) {
            ;
        } else if (decodeError != 0) {
            return STAR_ERRINMURA;
        } else {
            /* Convert from milliarcseconds per year to radians per century */
            coordBlock->muRA_radpcy = arcsecToRad(tempVal * 0.1);
            coordBlock->muRaInclCosDec = true;
            muRAsupplied = true;
        }

        if (*endPtr == ',') { endPtr++; }
        startPtr = endPtr;
        tempVal = myStrtod(startPtr, &endPtr, &decodeError);
        if (decodeError == STRINGEMPTY) {
            ;
        } else if (decodeError != 0) {
            return STAR_ERRINMUDEC;
        } else {
            /* Convert from milliarcseconds per year to radians per century */
            coordBlock->muDec_radpcy = arcsecToRad(tempVal * 0.1);
            muDecSupplied = true;
        }
        if (muRAsupplied && !muDecSupplied) {
            return STAR_NOMUDEC;
        }
        if (muDecSupplied && !muRAsupplied) {
            return STAR_NOMURA;
        }

    } else {
        return STAR_NORMAL;
    }

    /* Anything left? It will be annual parallax */
    if (*endPtr == ',') { endPtr++; }
    startPtr = endPtr;
    if (startPtr[0] != '\0') {
        tempVal = myStrtod(startPtr, &endPtr, &decodeError);
        if (decodeError == STRINGEMPTY) {
            ;
        } else if (decodeError != 0) {
            return STAR_ERRINPARAL;
        } else {
            /* Pi is specified in arcseconds. Convert to radians */
            coordBlock->parallax_rad = arcsecToRad(tempVal);
            piSupplied = true;
        }
        if (piSupplied && !muRAsupplied && !muDecSupplied) {
            return STAR_PINEEDSPM;
        }

    } else {
        return STAR_NORMAL;
    }

    /* Anything left? This is radial velocity */
    if (*endPtr == ',') { endPtr++; }
    startPtr = endPtr;
    if (startPtr[0] != '\0') {
        tempVal = myStrtod(startPtr, &endPtr, &decodeError);
        if (decodeError == STRINGEMPTY) {
            return STAR_NORMAL;
        } else if (decodeError != 0) {
            return STAR_ERRINVR;
        } else {
            /* Vr is specified in kilometres/second. Convert to AU per century*/
            coordBlock->radVel_aupcy = tempVal * auFactor * JUL_CENT;
        }
        if (muRAsupplied && muDecSupplied && piSupplied) {
            return STAR_NORMAL;
        } else {
            return STAR_VRNEEDSPM;
        }

    } else {
        return STAR_NORMAL;
    }
}



#if 0
GLOBAL void star_lightDeflection(const V3D_Vector *pV,
                                 const V3D_Vector *earthBV_au,
                                 const V3D_Vector *earthVelV_aupd,
                                 const V3D_Vector *sunBV_au,
                                 V3D_Vector *p2V)
/*! Applies the full relativistic corrections to the star's position for the
    effects of (a) light deflection (by the Sun's gravity) and (b) annual
    aberration. For tracking a star, this routine is basically unnecessary.
    The magnitude of light deflection when viewing a star at night is for all
    practical purposes negligible, and the very simple non-relativistic
    correction for annual aberration performed by the alternative routine
    star_annAberr() is within about 10 milliarcseconds (i.e. quite good enough).
    But here it is for completeness.
 \param[in]  pV             Unit vector pointing to object's geometric direction
                            as viewed from the centre of the earth.
 \param[in]  earthBV_au     Barycentric position of the Earth, referred to the
                            J2000 mean equator and equinox, as returned by
                            routine star_earth()
 \param[in]  earthVelV_aupd Earth velocity vector, referred to the J2000 mean
                            equator and equinox, as returned by routine
                            star_earth() (AU/day - not a unit vector)
 \param[in]  sunBV_au       Barycentric position of the Sun, referred to the
                            J2000 mean equator and equinox, as returned by
                            routine star_earth()
 \param[out] p2V            Unit vector pointing to object's proper
                            direction
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Vector      p1V;        /* pV, corrected for light deflection */
    V3D_Vector      dpV;
    V3D_Vector      e0V_au, eV; /* Heliocentric Earth direction */
    double          eMag_au;    /* Sun-Earth distance, |e0V_au| */
    double          ep;         /* dot product of eV and pV */
    V3D_Vector      vV;         /* Velocity vector, as fraction of light speed*/
    double          vMag;       /* Magnitude of velocity |vV| */
    double          p1dotV;     /* dot product of p1V & vV */
    double          invbeta, p1vbeta;
    double          scale;

    REQUIRE_NOT_NULL(pV);
    REQUIRE_NOT_NULL(earthBV_au);
    REQUIRE_NOT_NULL(earthVelV_aupd);
    REQUIRE_NOT_NULL(sunBV_au);
    REQUIRE_NOT_NULL(p2V);

    v3d_subtractV(&e0V_au, earthBV_au, sunBV_au);
    eMag_au = v3d_magV(&e0V_au);
    scale = 1.0 / eMag_au;
    eV.a[0] = e0V_au.a[0] * scale;
    eV.a[1] = e0V_au.a[1] * scale;
    eV.a[2] = e0V_au.a[2] * scale;

    ep = v3d_dotProductV(pV, &eV);

    dpV.a[0] = (2.0 * 9.87e-9 / eMag_au) * (eV.a[0] - ep * pV->a[0]) / (1 + ep);
    dpV.a[1] = (2.0 * 9.87e-9 / eMag_au) * (eV.a[1] - ep * pV->a[1]) / (1 + ep);
    dpV.a[2] = (2.0 * 9.87e-9 / eMag_au) * (eV.a[2] - ep * pV->a[2]) / (1 + ep);

    p1V.a[0] = pV->a[0] + dpV.a[0];
    p1V.a[1] = pV->a[1] + dpV.a[1];
    p1V.a[2] = pV->a[2] + dpV.a[2];
#if 1
    printVector("E", &e0V_au);
    printf("|E| = %.9f\n", eMag_au);
    printVector("e", &eV);
    printf("e.p = %.9f\n", ep);
    printVector("dp", &dpV);
    printVector("p1", &p1V);
#endif

    vV.a[0] = invC_dpau * earthVelV_aupd->a[0];
    vV.a[1] = invC_dpau * earthVelV_aupd->a[1];
    vV.a[2] = invC_dpau * earthVelV_aupd->a[2];
    vMag = v3d_magV(&vV);
    invbeta = sqrt(1 - vMag * vMag);
    p1dotV = v3d_dotProductV(&p1V, &vV);
    p1vbeta = 1.0 + p1dotV / (1.0 + invbeta);
#if 1
    printVector("V", &vV);
    printf("|V| = %.9f\n", vMag);
    printf("invbeta = %.9f\n", invbeta);
    printf("p1.V = %.9f\n", p1dotV);
    printf("1 + (p1.V)/(1 + invbeta) = %.9f\n", p1vbeta);
#endif
    p2V->a[0] = (invbeta * p1V.a[0] + p1vbeta * vV.a[0]) / (1.0 + p1dotV);
    p2V->a[1] = (invbeta * p1V.a[1] + p1vbeta * vV.a[1]) / (1.0 + p1dotV);
    p2V->a[2] = (invbeta * p1V.a[2] + p1vbeta * vV.a[2]) / (1.0 + p1dotV);
}
#endif


/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules)
 *
 *------------------------------------------------------------------------------
 */

LOCAL int extractObjectName(const char coordStr[],
                            char   nameStr[],
                            size_t nameStrSize,
                            const char **endPtr)
/*! Extract the celestial object's name from the coordinate string. If the name
    is present, it will be the first field in the string. It must be separated
    from the fields that follow by a comma. It may be enclosed in double quote
    characters, in which case a comma must immediately follow the closing double
    quote character. The double quotes will be removed from the output string.
    [in]  coordStr      The input coordinate string
    [out] nameStr       The object name
    [in]  nameStrSize   The maximum number chars that can be written to nameStr
    [out] endPtr        Points to first char in coordStr after the object name
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    size_t      n;
    const char  *sc, *ec;       /* Point to start char, end char of name */

    sc = coordStr;
    if (coordStr[0] == '"') {
        sc++;
        ec = strchr(sc, '"');
        *endPtr = ec;
        /* Move endPtr to char immediately after that second double-quote */
        if (*endPtr != NULL) { (*endPtr)++; }
    } else {
        ec = strchr(coordStr, ',');
        *endPtr = ec;
    }
    if (ec == NULL) {
        return STAR_ERRINOBJ;
    }

    n = (size_t)(ec - sc) + 1;  /* Include room for the '\0' at the end */
    if (n > nameStrSize) {
        n = nameStrSize;
    }
    saferStrncpy(nameStr, sc, n);

    return STAR_NORMAL;
}



LOCAL int parseEpochStr(const char epochStr[],
                        Star_CoordSys *coordSys,
                        double        *epochT_cy,
                        const char    **endPtr)
/*! Decodes a character string representing an Epoch into its component data;
    the coordinate system specification (APPARENT, FK4 or FK5 etc.), and the
    year. The string \a epochStr must be without spaces or tabs and be like
    one of the following format examples:
    - Besselian:            B1950  B1950.0  B1967.4  1975  1976.5
    - Julian:               J2000  J2000.0  J1985.1  1986  1987.5
    - Apparent Place:       APPARENT    (or abbreviation of at least 3 letters)
    - Intermediate Place:   INTERMEDIATE (or abbreviation of at least 3 letters)
    - ICRS Place:           ICRS (or abbreviation of at least 3 letters)
 \returns              An error code. STAR_NORMAL if OK, STAR_IVEPOCH if not, or
                       STRINGEMPTY if \a epochStr was empty
 \param[in]  epochStr  The input string containing the epoch or equinox
                       specification
 \param[out] coordSys  Coordinate system that has been implied by the equinox
                       specification
 \param[out] epochT_cy Time of epoch (Julian centuries since J2000.0)
 \param[out] endPtr    A pointer to the end of that part of \a epochStr that was
                       read to obtain the results. This may be pointing to some
                       white space between the results just decoded and further
                       data, or it may be the end of the string. A NULL
                       pointer may be passed to this parameter if you do not
                       need this value.

    If the year is specified without a preceding 'B' or 'J' character, it will
    be interpreted as a Besselian epoch if it is less than 1984.0, and as a
    Julian epoch if it is from 1984.0 onwards.

    If a Julian epoch is specified, \a coordSys will be set to #FK5.
    If a Besselian epoch is specified, \a coordSys will be set to #FK4.
    If APPARENT or INTERMEDIATE is specified, the epochT_cy value is meaningless
    so it is set to -99.9.
    If ICRS is specified, epochT_cy will be set to 0.0 (i.e. J2000.0)

 \note
    There are two independent reasons for specifying an epoch string. One is to
    identify the coordinate system in use and the equinox that applies to the
    celestial object's right ascension and declination, if that coordinate
    system is FK4 or FK5. The other is to specify the time zero for a celestial
    object's proper motion. In this case the value of \a coordSys returned by
    this function is irrelevant, and should be ignored.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    static const char *apparent = "APPARENT";
    static const char *intermediate = "INTERMEDIATE";
    static const char *icrs = "ICRS";
    double            tempEpoch;    /* Read into this variable and check range
                                       before inserting into epochT_cy */
    int               sError = 0;   /* Error status of myStrtod() calls */
    char              upEpStr[16];  /* Upper case version of epochStr[] */
    size_t            i, j;

    /* Trim leading space and convert epochStr to upper case before doing
       comparisons */
    i = 0;
    j = 0;
    while (isspace(epochStr[i])) {
        i++;
    }
    while(   (j < sizeof(upEpStr) - 1) && (epochStr[i] != '\0')
          && !isspace(epochStr[i]) && (epochStr[i] != ','))
    {
        upEpStr[j] = (char)toupper(epochStr[i]);
        i++;
        j++;
    }
    upEpStr[j] = '\0';
    if (endPtr != NULL) {
        *endPtr = &epochStr[i];
    }

    /* Was supplied string completely empty or blank? */
    if (upEpStr[0] == '\0') {
        *epochT_cy = UNUSED_EPOCH_cy;
        return STRINGEMPTY;
    }

        /* Does the string specify "apparent" or "intermediate" or "ICRS" with
       at least 3 characters? */
    if ((j >= 3) && (strstr(apparent, upEpStr) == apparent)) {
        *coordSys = APPARENT;
        *epochT_cy = UNUSED_EPOCH_cy;

    } else if ((j >= 3) && (strstr(intermediate, upEpStr) == intermediate)) {
        *coordSys = INTERMEDIATE;
        *epochT_cy = UNUSED_EPOCH_cy;

    } else if ((j >= 3) && (strstr(icrs, upEpStr) == icrs)) {
        *coordSys = ICRS;
        *epochT_cy = 0.0;

    } else {
        /* Must be a numeric epoch spec. Work out if it is Besselian or Julian*/
        if (upEpStr[0] == 'B') {
            *coordSys = FK4;
            tempEpoch = myStrtod(&epochStr[i - j + 1], NULL, &sError);

        } else if (upEpStr[0] == 'J') {
            *coordSys = FK5;
            tempEpoch = myStrtod(&epochStr[i - j + 1], NULL, &sError);

        } else if ((epochStr[i - j] >= '0') && (epochStr[i - j] <= '9')) {
            tempEpoch = myStrtod(&epochStr[0], NULL, &sError);
            if (tempEpoch < 1984.0) {
                *coordSys = FK4;
            } else {
                *coordSys = FK5;
            }

        } else {
            sError = STRINGEMPTY;   // String must be invalid.
        }

        /* Was a sensible epoch specified? */
        if ((sError != 0) || (tempEpoch < 0.0) || (tempEpoch >  9999.9)) {
            return STAR_IVEPOCH;

        } else {
            if (*coordSys == FK4) {
                *epochT_cy = ((tempEpoch - 1900.0) * TROP_YEAR
                              + (MJD_B1900 - MJD_J2000)) / JUL_CENT;
            } else if (*coordSys == FK5) {
                *epochT_cy = (tempEpoch - 2000.0) / 100.0;
            } else {
                return STAR_IVEPOCH;
            }
        }
    }
    return STAR_NORMAL;
}



LOCAL double myStrtod(const char str[], const char **endPtr, int *error)
/* Convert string to double-precision float.
Returns    - the result of the conversion as a floating point number
 Inputs
    str     - the input text string, to be passed to strtod()
 Output (optional)
    endPtr  - pointer to the first character after the last character that was
              read to obtain the number.
              A NULL pointer can be passed to this parameter if you don't care
              about this.
    error   - conversion status. One of
        0           : success
        STRINGEMPTY : No characters could be converted into a number. There was
                      no valid data in the string
        ERANGE      : the strtod() function detected overflow or underflow and
                      set errno to ERANGE. The returned result will be plus or
                      minus HUGE_VAL or zero.
              A NULL pointer can be passed to this parameter if you don't care
              about the conversion status.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    char        *endPtr2;
    double      num;
    int         status;
    int         savedErrno;

    savedErrno = errno;
    errno = 0;

    num = strtod(str, &endPtr2);
    if (endPtr2 == str) {
        status = STRINGEMPTY;   /* No conversion done; nothing valid found */
    } else {
        status = errno;         /* Might be ERANGE */
    }

    errno = savedErrno;
    if (endPtr != NULL) { *endPtr = endPtr2; }
    if (error != NULL) { *error = status; }
    return num;
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */



/*! \page page-stellar-reduction-accuracy Stellar reduction accuracy
 *
 *  Several tests were performed to compare the calculations of the routines in
 *  star.h / star.c with those given in examples in section B of the
 *  _Astronomical Almanac_ (various years).
 *
 *  The 1987 _Astronomical Almanac_ uses a fictitious star for its example
 *  calculations of stellar reduction. The catalog position of the star (given
 *  on page B40 of the _Almanac_) is
 *
 *      14 39 36.087, -60°50′07.14″, J2000,, -49.486 s/cy, +69.60″/cy, 0.752″, -22.2 km/s
 *
 *  The table below shows the results of converting this position to apparent
 *  coordinates (at the date shown in the table heading), using various
 *  approximations, and compared with the results from the _Almanac_ itself.
 *
 * |Apparent coords at 1987-01-01 0h TT |RA             |Dec             |
 * |------------------------------------|---------------|----------------|
 * |Almanac page B41                    | 14:38:40.164  | -60°46′44.82″  |
 * |Confirmation (see 1. below)         | 14:38:40.1641 | -60°46′44.822″ |
 * |No light deflection (see 2.)        | 14:38:40.1652 | -60°46′44.820″ |
 * |Simpler nutation (see 3.)           | 14:38:40.1639 | -60°46′44.820″ |
 * |IAU 2000 prec&nut                   | 14:38:40.1632 | -60°46′44.829″ |
 * |Approx Earth (see 4.)               | 14:38:40.1564 | -60°46′44.843″ |
 * |All of the above (see 5.)           | 14:38:40.1572 | -60°46′44.841″ |
 *
 *  The simplest calculation (see 5. below) gives a position error in RA of
 *  0.0069 s and in Dec of 0.019″, giving a total angular position error on the
 *  sky of 0.054″.
 *
 *  The same star appears on page B40 of the 1990 _Astronomical Almanac_.
 *  Likewise, the table below compares the calculations here with the results
 *  from the _Almanac_ itself.
 *
 * |Apparent coords at 1990-01-01 0h TT |RA             |Dec             |
 * |------------------------------------|---------------|----------------|
 * |Almanac page B41                    | 14:38:54.112  | -60°47′32.78″  |
 * |Confirmation (see 1. below)         | 14:38:54.1120 | -60°47′32.780″ |
 * |No light deflection (see 2.)        | 14:38:54.1131 | -60°47′32.778″ |
 * |Simpler nutation (see 3.)           | 14:38:54.1119 | -60°47′32.779″ |
 * |IAU2000 prec&nut                    | 14:39:54.1109 | -60°47′32.784″ |
 * |Approx Earth (see 4.)               | 14:38:54.1046 | -60°47′32.790″ |
 * |All of the above (see 5.)           | 14:38:54.1056 | -60°47′32.789″ |
 *
 *  The simple calculation (see 5. below) gives a position error in RA of
 *  0.0064 s and in Dec of 0.010″, giving a total angular position error on the
 *  sky of 0.048″.
 *
 *  The 2007 _Astronomical Almanac_ introduced a very slightly different
 *  fictitious star for its example calculations on page B68. This one has a
 *  catalog position of
 *
 *      14 39 36.4958, -60 50 02.309 ICRS, J2000.0, -3678.06 mas/yr, 482.87 mas/yr, 0.742″, -21.6 km/s
 *
 *  Again, the table below compares the calculations here with the results
 *  from the _Almanac_ itself.
 *
 * |Apparent coords at 2007-01-01 0h TT |RA             |Dec             |
 * |------------------------------------|---------------|----------------|
 * |Almanac page B68                    | 14:40:03.4343 | -60°51′37.770″ |
 * |Nutation IAU2000B (see 6. below)    | 14:40:03.4342 | -60°51′37.770″ |
 * |No light deflection (see 7.)        | 14:40:03.4353 | -60°51′37.769″ |
 * |Approx Earth (see 8.)               | 14:40:03.4265 | -60°51′37.784″ |
 * |All of the above (see 9.)           | 14:40:03.4276 | -60°51′37.782″ |
 * |Nutation IAU1980 with precision 0   | 14:40:03.4400 | -60°51′37.784″ |
 * |All above and Nutation 1980         | 14:40:03.4334 | -60°51′37.796″ |
 * |All above and simpler Nutation 1980 | 14:40:03.4336 | -60°51′37.796″ |
 *
 *  The simple calculation (see 9. below) gives a position error in RA of
 *  0.0067 s and in Dec of 0.012″, giving a total angular position error on the
 *  sky of 0.051″.  Or if we use 1980 nutation, 0.0009s and 0.026″,
 *
 *  These three results show very small position errors.
 *
 *  Test conditions
 *      1. Confirmation: sky1_nutationIAU1980() called with \a precision set to
 *         0, vectors \b Eb, \b Ėb, and \b Sb taken from the relevant _Almanac_,
 *         page B39, and star_lightDeflection() used to calculate aberration
 *         etc.
 *      2. Ignore light deflection. That is, conditions as for 1, but the much
 *         simpler routine star_annAberr() called instead of
 *         star_lightDeflection().
 *      3. Simpler nutation. Conditions as for 1, but sky1_nutationIAU1980()
 *         called with \a precision set to 4.
 *      4. Approx Earth. Conditions as for 1, but the vectors \b Eb, \b Ėb, and
 *         \b Sb are calculated by the very approximate routine star_earth()
 *      5. All simplifications: sky1_nutationIAU1980() called with \a precision
 *         set to 4, vectors \b Eb, \b Ėb, and \b Sb obtained approximately from
 *         star_earth(), and star_annAberr() used to calculate annual aberration
 *      6. Approximate confirmation: astc2_nutationIAU2000B() called to
 *         calculate nutation, vectors \b Eb, \b Ėb, and \b Sb taken from the
 *         2007 _Almanac_ page B68, and star_lightDeflection() used to calculate
 *         aberration etc.
 *      7. Ignore light deflection. That is, conditions as for 6, but the much
 *         simpler routine star_annAberr() called instead of
 *         star_lightDeflection().
 *      8. Approx Earth. Conditions as for 6, but the vectors \b Eb, \b Ėb, and
 *         \b Sb are calculated by the very approximate routine star_earth()
 *      9. All simplifications: astc2_nutationIAU2000B() called, vectors \b Eb,
 *         \b Ėb, and \b Sb obtained approximately from star_earth(), and
 *         star_annAberr() used to calculate annual aberration
 *
 *  */


