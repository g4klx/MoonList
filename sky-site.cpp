/*==============================================================================
 * sky-site.c - astronomical routines related to an observing site
 *
 * Author:  David Hoadley
 *
 * Description: (see the "Site routines" sections of sky.h)
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
 *      Things you might want to edit: definition of macro SPA_COMPARISONS
 *----------------------------------------------------------------------------*/


/* ANSI includes etc. */
#include "instead-of-math.h"                /* for sincos() */
#include <math.h>

/* Local and project includes */
#include "sky.h"
#include "general.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       /* For use by REQUIRE() - assertions */

/*      Un-comment the following if you are running comparison tests vs the
        NREL SPA code for Sun positions. It changes the constants describing
        the Earth to match those used in the SPA itself (instead of using the
        (better) constants from the Astronomical Almanac). Also it disables
        the correction for diurnal aberration (because SPA does not calculate
        this correction). */
/*--- #define SPA_COMPARISONS ---*/

/*
 * Prototypes for local functions (not called from other modules).
 */
LOCAL void createAzElBaseM(Sky_SiteProp *site);

/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */
#ifndef SPA_COMPARISONS
/*      Constants found in the 2007 Astronomical Almanac, pages K6 & K7 */
LOCAL const double c_kmps = 299792.458;      // speed of light (km/s)
LOCAL const double f = 0.0033528197;         // flattening factor of the geoid
LOCAL const double ae_km = 6378.1366;        // equatorial radius of geoid (km)
LOCAL const double omega_radps = 7.292115e-5;// mean earth rotation rate(rad/s)
LOCAL const double au_km = 1.49597871464e8;  // one Astronomical Unit (km)
#else
/*      Constants found in the NREL SPA document */
LOCAL const double f = 1.0 - 0.99664719;     // = 0.0033528100
LOCAL const double ae_km = 6378.140;
LOCAL const double au_km = 1.4960039e8;      // implied by sol parallax=8.794"
/*      Constants not found in NREL, so copied from above */
LOCAL const double c_kmps = 299792.458;      // speed of light (km/s)
LOCAL const double omega_radps = 7.292115e-5;// mean earth rotation rate(rad/s)
#endif
/*      Derived constants */
LOCAL const double esq = f + f - f * f;      // square of eccentricity of geoid

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
GLOBAL void sky_setSiteLocation(double latitude_deg,
                                double longitude_deg,
                                double height_m,
                                Sky_SiteProp *site)
/*! Initialise the \a site structure by calculating those site-related values
    that do not change with time.
 \param[in]  latitude_deg   Latitude of site (ϕ) (degrees)
 \param[in]  longitude_deg  Longitude of site (degrees)
 \param[in]  height_m       Height above ellipsoid (e.g. GPS height) (metres)

 \param[out] site           All fields initialised

    This is a simplified version of sky_setSiteLoc2(), which is good enough
    for most purposes.

    This function calculates the quantities that will be used to convert
    geocentric coordinates to topocentric ones. This involves calculating the
    position of the site relative to the centre of the earth, and calculating
    constants that will be used for diurnal parallax and aberration corrections.

    This is complicated by the fact that the earth is elliptical rather than
    perfectly spherical. See section K of The Astronomical Almanac (pages K11
    & K12 in the 2007 edition) for details of the geometry. This function
    effectively calculates the Geocentric latitude (ϕ′) from the Geodetic
    latitude (ϕ), and then sin(ϕ - ϕ′) and cos(ϕ - ϕ′). But in practice, the
    sin and cos terms can be obtained without explicitly calculating ϕ′ itself,
    and that is done here.

    The \a site->refracPT field is initialised assuming that the site pressure
    is 1010 hPa and the temperature is 10 °C. If other values are required,
    follow this call by a call to sky_setSiteTempPress().

    The \a site->timezone_d field is initialised to 0, assuming that the site
    time zone is UTC. If another value is required, follow this call by a call
    to sky_setSiteTimeZone().

 \par When to call this function
    At program initialisation time only
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /* Use same values for astronomical and geodetic lat/long */
    sky_setSiteLoc2(latitude_deg,
                         longitude_deg,
                         latitude_deg,
                         longitude_deg,
                         height_m,
                         site);
}



GLOBAL void sky_setSiteLoc2(double astLat_deg,
                            double astLong_deg,
                            double geodLat_deg,
                            double geodLong_deg,
                            double height_m,
                            Sky_SiteProp *site)
/*! Alternative initialisation function to sky_setSiteLocation(). It initialises
 *  the \a site structure, but it supports a distinction between Astronomical
 *  latitude & longitude and Geodetic latitude & longitude.
 \param[in]  astLat_deg     Astronomical Latitude of site (ϕA) (degrees)
 \param[in]  astLong_deg    Astronomical Longitude of site (degrees)
 \param[in]  geodLat_deg    Geodetic Latitude of site (GPS latitude) (ϕ)
                            (degrees)
 \param[in]  geodLong_deg   Geodetic Longitude of site (GPS longitude) (degrees)
 \param[in]  height_m       Height above ellipsoid (e.g. GPS height) (metres)

 \param[out] site           All fields initialised

    This function calculates the same quantities as described for
    sky_setSiteLocation().
    The only difference is that the \a site.azElBaseM rotation matrix is based
    on Astronomical latitude and latitude (instead of Geodetic).

    Astronomical latitude and longitude are coordinates as determined by the
    local gravity vector and astronomical observation.
    Geodetic latitude and longitude are coordinates as measured by some
    geodetic system - e.g. a map reference, or GPS coordinates.

 \note Distinguishing between Astronomical lat & lon and Geodetic lat & lon
    is only required for very precise work. It is most likely that you won't
    actually know the Astronomical coords, and it won't matter. Use function
    sky_setSiteLocation() instead - it sets the Astronomical coords to the same
    values as the  Geodetic or GPS coords.

    The \a site->refracPT field is initialised assuming that the site pressure
    is 1010 hPa and the temperature is 10 °C. If other values are required,
    follow this call by a call to sky_setSiteTempPress().

    The \a site->timezone_d field is initialised to 0, assuming that the site
    time zone is UTC. If another value is required, follow this call by a call
    to sky_setSiteTimeZone().

 \par When to call this function
    At program initialisation time only, if you need to distinguish between
    Astronomical latitude/longitude and Geodetic latitude/longitude. Otherwise
    call sky_setSiteLocation() instead.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double     geodLat_rad;         // Geodetic Latitude (radian)
    double     geodLong_rad;        // Geodetic Longitude (radian)
    double     sinPhi, cosPhi;      // sin and cos of Geodetic Latitude
    double     invC;                // = 1/C, (for C see Astron. Almanac p. K11)
    double     aeC_km;              // = ae * C (km)
    double     rhoCosPhiP_km;       // = ae * rho * cos(phi') (km)
    double     Kc;
    V3D_Matrix mat1, mat2;

    REQUIRE_NOT_NULL(site);

    /* 1. Calculate radian values of some of the site constants */
    site->astLat_rad = degToRad(astLat_deg);
    site->astLong_rad = degToRad(astLong_deg);

    geodLat_rad = degToRad(geodLat_deg);
    geodLong_rad = degToRad(geodLong_deg);

    /* 2. Create rotation matrix from astronomical latitude and longitude */
    createAzElBaseM(site);

    /* 3. Calculate the geocentric radius of the observatory and the two
       geocentric parallax quantities
       ae*rho*sin(Phi-Phi') and -ae*rho*cos(Phi-Phi')
        (The geodetic latitude is phi, the geocentric latitude is phi'.
         The angle (phi - phi') is known as the "angle of the vertical".) */
    sincos(geodLat_rad, &sinPhi, &cosPhi);
    invC = sqrt(1.0 - esq * (sinPhi * sinPhi));
    aeC_km = ae_km / invC;
    Kc = sqrt(1.0 - esq * (2.0 - esq) * (sinPhi * sinPhi));
    site->geocRadius_km = aeC_km * Kc + height_m / 1000.0;

    //  Parallax corrections are offsets in x (north) and z (zenith) directions
    site->rhoSin_au = aeC_km * esq * sinPhi * cosPhi / au_km;      // x (north)
    site->rhoCos_au = -(ae_km * invC + height_m / 1000.0) / au_km; // z (zenith)

    /* 4. Calculate the diurnal aberration correction */
    rhoCosPhiP_km = (aeC_km + height_m / 1000.0) * cosPhi;
    site->diurnalAberr = omega_radps * rhoCosPhiP_km / c_kmps;     // y (east)

    /* 5. Might want to convert from Az-El frame to HA-Dec frame. Create
       rotation matrix to do this. Matrix must first rotate around Y axis by
       Pi/2 - Latitude. Second rotation is around Z' axis by Pi. */
    v3d_createRotationMatrix(&mat1, Yaxis, HALFPI - site->astLat_rad);
    v3d_createRotationMatrix(&mat2, Zaxis, PI);
    v3d_multMxM(&site->haDecM, &mat2, &mat1);

    /* 6. Other initialisations, in case sky_setSiteTempPressure() or
       sky_setSiteTimeZone() don't get called. */
    site->refracPT = 1.0;               // Default to 10 °C and 1010 hPa
    site->timeZone_d = 0.0;             // Default to UTC.

    site->azElM = &site->azElBaseM;     // Assume no polar motion correction
}



GLOBAL void sky_setSiteTempPress(double temperature_degC,
                                 double pressure_hPa,
                                 Sky_SiteProp *site)
/*! Set refraction coefficients based on atmospheric temperature and pressure
    at the site.
 \param[in]  temperature_degC - average annual air temperature at the site (°C)
 \param[in]  pressure_hPa     - average annual air pressure at the site
                                (hPa = mbar)
 \param[out] site             - field \a refracPT, a combined refraction
                                coefficient for pressure & temperature

    This coefficient will be used in the simple refraction algorithm that is
    called from routine sky_siteTirsToTopo().

 \par When to call this function
    1. At program initialisation time, after calling sky_setSiteLocation() (or
        sky_setSiteLoc2())
    2. If you have updated values of pressure or temperature
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double p, t;

    REQUIRE_NOT_NULL(site);
    REQUIRE(temperature_degC > -100.0);         // what were you thinking?

    t = 283.0 / (273.0 + temperature_degC);
    p = pressure_hPa / 1010.0;
    site->refracPT = p * t;
}



GLOBAL void sky_setSiteTimeZone(double timeZone_h,
                                Sky_SiteProp *site)
/*! Set the time zone offset for this site.
 \param[in]  timeZone_h  Zonal correction (hours). Time zones east of
                           Greenwich are positive (e.g. Australian Eastern
                           Standard Time is +10.0) and those west of Greenwich
                           are negative (e.g. US Pacific Daylight Time is -7.0)

 \param[out] site        field \a timeZone_d, time zone scaled to fractions of
                           a day

 \par When to call this function
    1. At program initialisation time, after calling sky_setSiteLocation() (or
        sky_setSiteLoc2())
    2. If daylight saving begins or ends during the time the program is running
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(site);

    site->timeZone_d = timeZone_h / 24.0;
}



GLOBAL void sky_setupSiteSurface(double azimuth_deg,
                                 double slope_deg,
                                 Sky_SiteHorizon *surface)
/*! Set the orientation and slope of a surface (such as a solar panel) for which
    you want to calculate the incidence angle of incoming radiation.
 \param[in]  azimuth_deg  The azimuth angle of the normal to the surface
                          (degrees). This is measured clockwise from North, as
                          for all other azimuth angles.
 \param[in]  slope_deg    The slope of the surface measured from the horizontal
                          (degrees). Alternatively, this is the zenith distance
                          of the normal to the surface.
 \param[out] surface      The coordinates of the normal to the surface, in both
                          polar and rectangular forms. The rectangular form will
                          be passed later to sky_siteIncidence_rad().

 \par When to call this function
    At program initialisation time only, if you have a surface of interest.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(surface);

    surface->azimuth_rad = degToRad(azimuth_deg);
    surface->elevation_rad = degToRad(90.0 - slope_deg);
    v3d_polarToRect(&surface->rectV,
                    surface->azimuth_rad,
                    surface->elevation_rad);
}



GLOBAL void sky_adjustSiteForPolarMotion(const Sky_PolarMot *polar,
                                         Sky_SiteProp *site)
/*! Modify our azEl rotation matrix to incorporate a polar motion rotation
    matrix
 \param[in]     polar   Polar motion parameters, as set by function
                        sky_setPolarMotion()
 \param[in,out] site    Fields \a azElM and possibly \a azElPolM are updated

 \par When to call this function
    Polar motion is such a tiny effect that you may well decide not to bother
    with it. So you only need to call this routine if:
        1. you are bothering about polar motion at all, and
        2. changes were made to the \a polar structure by a recent call to
            function sky_setPolarMotion(). This routine is only ever to be
            called after that routine has run. This is expected to be very
            infrequently - once per day is more than enough.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(polar);
    REQUIRE_NOT_NULL(site);

    if (polar->correctionInUse) {
        v3d_multMxM(&site->azElPolM, &site->azElBaseM, &polar->corrM);
        site->azElM = &site->azElPolM;  // Use combined matrix from now on
    } else {
        site->azElM = &site->azElBaseM; // Use matrix without polar motion
    }
}



GLOBAL void sky_siteTirsToTopo(const V3D_Vector   *terInterV,
                               double             dist_au,
                               const Sky_SiteProp *site,
                               Sky_SiteHorizon *topo)
/*! Transform a coordinate vector from the Terrestrial Intermediate Reference
    System to topocentric Az/El coordinates for the observing site whose
    properties are described in parameter "site". Note: there is a mixture of
    vectors in right-handed and left-handed coordinate systems here
 \param[in]  terInterV  Position vector in Terrestrial Intermediate Reference
                        System. (This will have been obtained by calling either
                        routine astsite_appToTirs() or routine
                        astsite_intermedToTopo().)
 \param[in]  dist_au    Geocentric Distance to object (astronomical units).
                        Note: for far distant objects outside of the solar
                        system, you can supply 0.0 for this value. It will be
                        treated as infinity.
 \param[in]  site       Block of data describing the observing site, as
                        initialised by one of the functions
                        sky_setSiteLocation() or sky_setSiteLoc2().
 \param[out] topo       Topocentric position, both as a vector in horizon
                        coordinates, and as azimuth (radian) and elevation
                        (radian).

 \par When to call this function
    When you have the position of a celestial object expressed in the
    Terrestrial Intermediate Reference System, use this function to convert it
    to topocentric coordinates (the coordinates as seen from a particular
    observing location on the Earth's surface). If you are running a control
    loop to enable continuous tracking of this object, you will need to call
    this function every time around your control loop.
 \par
    If you wish to calculate the position of the object as seen from multiple
    sites simultaneously, you need to call this function once for each of those
    sites, passing the relevant \a site data block to each call.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double dEl_rad;         // Change in elevation from refraction (radian)
    double w;
    double tanZd;           // Tangent of zenith distance

    REQUIRE_NOT_NULL(terInterV);
    REQUIRE_NOT_NULL(site);
    REQUIRE_NOT_NULL(topo);

    /* Rotate from the equatorial system to the horizon coordinate system. The
       resultant position vector, referred to North, East and Zenith, is in a
       left-handed system, like the Azimuth & Elevation angle system to which it
       will eventually be converted. */
    v3d_multMxV(&topo->rectV, site->azElM, terInterV);
    /* Note, at this point topoV is not actually topocentric yet - it is still
       a geocentric point-of-view, even though the coordinate system is now
       specific to our particular site. */

    /* Transform from Geocentric to Topocentric coordinates - i.e. correct for
       Diurnal Aberration and Geocentric Parallax. This may be done by
       calculating a combined correction as a vector addition.
       (The rhoSin and rhoCos terms here are in a geodetic coordinate system
       [i.e. referred to the geodetic vertical] whereas the vector topoV is
       referred to the astronomical vertical. However the error introduced by
       simply adding the correction vector [as though the two vectors were in
       the same coordinate system] is minuscule, since the "deflection of the
       vertical" is so small - almost certainly < 20 arcseconds.) */
#ifndef SPA_COMPARISONS
    topo->rectV.a[1] += site->diurnalAberr;
#else
#warning "Correction for diurnal aberration is not being applied"
#endif
    if (dist_au > 0.0) {
        topo->rectV.a[0] += site->rhoSin_au / dist_au;
        topo->rectV.a[2] += site->rhoCos_au / dist_au;
    }
    // else
    //      We treat 0.0 (or -ve values) as meaning "infinitely far away".
    //      Objects that far away have no parallax, so we need do nothing here

    v3d_rectToPolar(&topo->azimuth_rad, &topo->elevation_rad, &topo->rectV);

#if 0
    /* Correct for atmospheric refraction. Use the simpler NREL SPA calculation
       instead of the more detailed atmospheric model used by Stromlo.
       Unfortunately, this formula is expressed in degrees, so we have to
       convert back and forth. */
    {
        double e0_deg;      // Elevation (not corrected for refraction, degrees)

        e0_deg = radToDeg(topo->elevation_rad);
        if (e0_deg > -2.0) {
            dEl_rad = degToRad(1.02 / (60.0 * tan(degToRad(e0_deg + 10.3
                                                           / (e0_deg + 5.11)))))
                      * site->refracPT;
        } else {
            dEl_rad = 0.0;
        }
    }
#elif 0
    /* Correct for atmospheric refraction using a radian version of the above */
    if (topo->elevation_rad > (-2.0 * DEG2RAD) {
        dEl_rad = 0.000296706 / tan(topo->elevation_rad + 0.00313756
                                           / (topo->elevation_rad + 0.0891863))
                  * site->refracPT;
    } else {
        dEl_rad = 0.0;
    }
#else
    /* Correct for atmospheric refraction using two simpler formulae */
    w = sqrt(topo->rectV.a[0] * topo->rectV.a[0]
             + topo->rectV.a[1] * topo->rectV.a[1]);
    if (topo->rectV.a[2] >= 0.268 * w) {
        /* Elevation is greater than 15° */
        tanZd = w / topo->rectV.a[2];
        dEl_rad = (2.8253e-4 * tanZd - 3.9948e-7 * tanZd * tanZd * tanZd)
                  * site->refracPT;

    } else if (topo->elevation_rad > (-2.0 * DEG2RAD)) {
        /* Low elevation, use this approx formula instead */
        dEl_rad =  (8.3323e-3 + 3.1786e-2 * topo->elevation_rad
                    + 2.0746e-2 * topo->elevation_rad * topo->elevation_rad)
                 / (1 + 20.995 * topo->elevation_rad
                    + 160.31 * topo->elevation_rad * topo->elevation_rad)
                 * site->refracPT;

    } else {
        dEl_rad = 0.0;
    }
#endif
    topo->elevation_rad += dEl_rad;

    /* Convert back to rectangular coords */
    v3d_polarToRect(&topo->rectV, topo->azimuth_rad, topo->elevation_rad);
}



GLOBAL void sky_siteAzElToHaDec(const V3D_Vector   *topoV,
                                const Sky_SiteProp *site,
                                double *hourAngle_rad,
                                double *declination_rad)
/*! Take a topocentric position vector in Azimuth/Elevation frame and use it to
    calculate the observed Hour Angle and Declination.
 \param[in]  topoV            Topocentric position vector in Azimuth/Elevation
                              frame (as returned by the sky_siteTirsToTopo()
                              function)
 \param[in]  site             Field \a haDecM - rotation matrix to HA/Dec
                              coordinates
 \param[out] hourAngle_rad    Hour angle (radian - note: not in hours)
 \param[out] declination_rad  Declination (radian)

    The hour angle and declination returned by this function will be those of
    the refracted position of the object.

 \par When to call this function
    After each call to sky_siteTirsToTopo(), but only if you want hour angle and
    declination.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Vector haDecV;  // topocentric posn vector in Hour Angle/Decl. frame
                        // like topoV, this vector is a left-handed set

    REQUIRE_NOT_NULL(topoV);
    REQUIRE_NOT_NULL(site);
    REQUIRE_NOT_NULL(hourAngle_rad);
    REQUIRE_NOT_NULL(declination_rad);

    v3d_multMxV(&haDecV, &site->haDecM, topoV);
    v3d_rectToPolar(hourAngle_rad, declination_rad, &haDecV);
}



GLOBAL double sky_siteIncidence_rad(const V3D_Vector *topoV,
                                    const V3D_Vector *surfaceV)
/*! Calculate the incidence angle of rays from the celestial object described
    by \a topoV (such as the sun) falling onto a surface described by
    \a surfaceV (such as a solar panel).
 \returns    Incidence angle of the radiation (radian)
 \param[in]  topoV      Unit vector pointing to the celestial object, as
                        returned by the sky_siteTirsToTopo() routine (in the
                        \a topo->rectV field). Must be of unity magnitude,
                        otherwise an assertion failure occurs.
 \param[in]  surfaceV   Unit vector of normal to the surface, as returned by
                        routine sky_setupSiteSurface() (in the \a surface->rectV
                        field). Must be of unity magnitude, otherwise an
                        assertion failure occurs.

 \par When to call this function
    After each call to sky_siteTirsToTopo(), but only if you have a surface that
    you want the incidence angle for.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    // This function assumes that both input vectors are unit vectors.
    REQUIRE(fabs(v3d_magVSq(topoV) - 1.0) < SFA);
    REQUIRE(fabs(v3d_magVSq(surfaceV) - 1.0) < SFA);

    // Dot product of two unit vectors gives cosine of angle between them
    return acos(v3d_dotProductV(topoV, surfaceV));
}

/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules).
 *
 *------------------------------------------------------------------------------
 */
LOCAL void createAzElBaseM(Sky_SiteProp *site)
/* Calculates the transformation matrix to convert a coordinate vector
   expressed in a terrestrial equatorial coordinate frame to one expressed in
   the horizon coordinate frame.

   The terrestrial equatorial frame has the geocentre as its origin, the
        x axis points to 0° Lat and 0° Long (called the terrestrial intermediate
               origin, or TIO)
        y axis points to 0° Lat and 90° East Longitude
        z axis points to the (north) celestial intermediate pole (CIP)
   If converted to polar coordinates, the angles are the negative of Greenwich
   Hour angle (GHA) and Declination.

   The horizon frame typically is centered on the observing site, and the vector
   is a Left-handed set, unlike the vector above. The
        x axis points to the North horizon,
        y axis points to the East horizon,
        z axis points to the zenith
   If converted to polar coordinates, the angles are Azimuth and Elevation.

   Note: this routine creates a matrix for rotation only; the resulting vector
   is in horizon coordinates, but it is still a geocentric position. Correcting
   for displacement to the earth's surface (diurnal parallax and aberration)
   needs to be carried out separately. (But that is a relatively simple
   operation in horizon coordinates.)

   To convert from one frame to the other, we must
    1. Rotate around Z axis by astronomical longitude
    2. Rotate around new Y' axis by -ve of astronomical latitude (phiA)
    3. swap the resulting Z'' and X'' axes to obtain the X' and Z' axes
       described above.
   Thus          azElBaseM = J * R2(-PhiA) * R3(longitude)
   where J = 0 0 1
             0 1 0
             1 0 0
 Inputs
    site->astLong_rad    - astronomical longitude of site (radian)
    site->astLat_rad     - astronomical latitude of site (radian)
 Output
    site->azElBaseM      - 3x3 rotation matrix, -GHA/Dec to Az/El
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double sinLon, cosLon;
    double sinLat, cosLat;

    sincos(site->astLong_rad, &sinLon, &cosLon);
    sincos(site->astLat_rad, &sinLat, &cosLat);
    site->azElBaseM.a[0][0] = -sinLat * cosLon;
    site->azElBaseM.a[0][1] = -sinLat * sinLon;
    site->azElBaseM.a[0][2] = cosLat;

    site->azElBaseM.a[1][0] = -sinLon;
    site->azElBaseM.a[1][1] = cosLon;
    site->azElBaseM.a[1][2] = 0.0;

    site->azElBaseM.a[2][0] = cosLat * cosLon;
    site->azElBaseM.a[2][1] = cosLat * sinLon;
    site->azElBaseM.a[2][2] = sinLat;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*! \page page-about About this code, and how to use it
 *
 *  This code has been written with the aim of providing astronomical routines
 *  for tracking the Sun, Moon, planets and stars using a small processor. It
 *  started as a project to implement accurate tracking of the Sun for solar
 *  heliostat applications. What was required was code that would execute more
 *  quickly than the traditional Solar Position Algorithm (see reference), but
 *  without the compromises in accuracy that many simpler algorithms suffer
 *  from.
 *
 *  Some of the principles of that implementation, in particular the use of
 *  rectangular coordinates for converting from geocentric apparent coordinates
 *  to site-specific topocentric coordinates, and the use of interpolation
 *  between apparent position vectors, are just as applicable to the tracking of
 *  any other celestial object. So gradually code to support tracking the Moon,
 *  planets and stars has been added.
 *
 *  There are a couple of use cases for this code. One could be a field of
 *  heliostats which need accurate solar tracking but are each controlled by a
 *  small cheap processor.
 *  Another could be a home project to automate a home telescope.
 *
 *  This code is intended to handle all the intricacies of accurate celestial
 *  calculations, leaving you free to concentrate on the code to actually drive
 *  your device to the desired position.
 *  The code is designed to provide accurate conversion of astronomical
 *  coordinate systems, such as is used on quite large telescopes. Tracking
 *  accuracy of better than one arcsecond can be expected, (although the
 *  positions of the Moon and other planets are not calculated to that level of
 *  accuracy).
 *      - \subpage page-design-choices
 *      - \subpage page-how-to-use
 *      - \subpage page-licensing
 *
 * \par Reference:
 *          Reda, Ibrahim and Andreas, Afshin.
 *          _Solar Position Algorithm for Solar Radiation Applications._
 *          National Renewable Energy Laboratory, publication no.
 *          NREL/TP-560-34302, June 2003, revised January 2008
 */

/*! \page page-conventions Conventions used in this code
 *      - \subpage page-c-style
 *      - \subpage page-var-suffixes
 *      - \subpage page-sec-arcsec
 *      - \subpage page-sign-conventions
 *      - \subpage page-interval-notation
 *  */

/*! \page page-time About time
 *      - \subpage page-timescales
 *      - \subpage page-time-variables
 *  */

/*! \page page-misc Miscellaneous
 *      - \subpage page-vernal-equinox
 *      - \subpage page-why-struct-array
 *  */

/*! \page page-design-choices Design choices
 *  To make the code suitable for a small processor, the following choices have
 *  been made.
 *      - The code can be compiled with a C89/C90 compliant compiler. (Some
 *        embedded systems still use compilers of this vintage.) But if a C99
 *        or later compiler is being used, the code will use a few inline
 *        functions in place of function-like macros.
 *      - Does not use any heap storage (dynamic memory). That is, there are no
 *        \c malloc() calls.
 *      - Does not use any variable-length arrays. Thus stack usage is
 *        predictable.
 *      - \c printf() calls are basically limited to debugging code.
 *      - Very limited use of \c snprintf(). Does not use the \%f specifier in
 *        any \c snprintf().
 *      - Uses an interpolation process to dramatically reduce the calling of
 *        trigonometric functions during tracking, with a very small loss of
 *        accuracy (see \ref page-interpolation).
 */

/*! \page page-how-to-use How to use this code
 *  The code has been provided as a set of separate files for you to include in
 *  your project. It has not been configured as a monolithic library, because
 *  there are likely to be parts that you don't want to include. But
 *  importantly, it is assumed that you will want to make some choices, and you
 *  will need to do that by making some simple edits to the source code of
 *  various modules. Here are some guidelines.
 *
 *  ###Required files
 *      - For all applications:
 *          + The basic sky routines - sky.h & sky-time.c & sky-site.c
 *          + Low level support routines and definitions -
 *            general.h, astron.h, instead-of-math.h, vectors3d.h & vectors3d.c.
 *      - For solar position calculations: the "all applications" set above,
 *        plus sky0.h & sky0.c, and sun.h & sun.c
 *      - For Sun tracking: the above, plus skyfast.h & skyfast.c
 *      - For Moon position calculations: the "all applications" set above, plus
 *        sky0.h & sky0.c, and moon.h & moon.c
 *      - For Moon tracking: the above, plus skyfast.h & skyfast.c
 *      - For planet position calculations: the "all applications" set above,
 *        plus sky1.h & sky1.c, and planet.h & planet.c
 *      - For tracking a planet: the above, plus skyfast.h & skyfast.c
 *      - For a star: the "all applications" set above, plus
 *        sky1.h & sky1.c, sky2.h & sky2.c, star.h & star.c, and
 *        skyio.h & skyio.c
 *      - For tracking a star: the above, plus skyfast.h & skyfast.c
 *
 *  You may also want to use skyio.h and skyio.c to format the output angles for
 *  display for the Sun, Moon or planets.
 *
 *  ###Editing the code
 *  You are likely to want to make some basic edits to the code, as listed
 *  below. But before you do, make sure your editor can handle files in UTF-8
 *  format. MacOS and Linux systems certainly will do this without you needing
 *  to think about it. Only Windows systems still using old 8-bit code pages
 *  like the Windows-1252 page are likely to have any issues. (You should move
 *  away from them to UTF-8 in any case. UTF-8 can handle any character; those
 *  old code pages don't.)
 *      - \subpage page-general-h
 *      - \subpage page-sky-h
 *      - \subpage page-skyfast-c
 */

/*! \page page-general-h Edits you may want to make to general.h
 *
 *  ####Integer types with a C89/C90 compiler.
 *  If you are compiling your code using a C89/C90 compiler, you will need to
 *  check that the following lines are correct for your target processor. (If
 *  you have a C99 or later compiler, you can skip this section.)
 \verbatim
    typedef signed short int    int16_t;
    typedef unsigned short int  uint16_t;
    typedef int                 int32_t;
    typedef unsigned int        uint32_t;
 \endverbatim
 *  If, for example, your processor defines an int as being 16 bits wide, you
 *  will get some strange warning messages. One possibility I have seen is
 *  "general.h:105:1: warning: division by zero [-Wdiv-by-zero]".
 *  Another possibility is
 *  "general.h:105:1: error: expression is not an integer constant expression".
 *  The important thing to do is to look at the line number where the error
 *  was detected, where you will see a line like
 \verbatim
    static_assert(sizeof(int16_t) == 2,  "typedef of int16_t is wrong");
 \endverbatim
 *  To fix the problem, don't change the line containing the static_assert(),
 *  go back and change the line(s) containing the typedef(s), so that these four
 *  types (\c int16_t, \c uint16_t, \c int32_t & \c uint32_t) are correctly
 *  defined for your processor.
 *
 *  ####Assertion checking
 *  Many routines in this collection use assertions to check their parameters.
 *  In particular, parameters which are passed by address are checked to make
 *  sure that a NULL pointer has not been passed to them. You can control this
 *  checking by defining three macros within this file. They are:
 *      - ASSERTION_LEVEL. Use this to set how much checking is done. While you
 *        are still writing and testing your own code which uses this code, you
 *        are strongly recommended to keep this defined at a non-zero value.
 *      - ENABLE_NULL_CHECKING. So long as ASSERTION_LEVEL is not set to zero,
 *        this macro enables the checking of parameters which are passed by
 *        reference, to make sure they are not NULL. (There is the occasional
 *        parameter which is allowed to be NULL - such a parameter is of course
 *        not subject to this check.) Once again, while you
 *        are still writing and testing your own code which uses this code, you
 *        are strongly recommended to keep this macro defined.
 *      - USE_STANDARD_ASSERT. The code is designed to be run on an embedded
 *        processor, where the traditional assert() behaviour -- write a message
 *        and exit the program -- does not make any sense. On an embedded system
 *        you need to write your own assertion handling routine onAssert__(),
 *        which handles the information in a way that helps you debug the
 *        problem. But if you are running this code on a processor which has an
 *        operating system (say, like a Raspberry Pi runs Linux), you don't need
 *        to go to this extra trouble. You can use the traditional assert()
 *        behaviour quite happily. Define this macro to enable this behaviour.
 */


/*! \page page-licensing Open Source Licensing
 *
 * All code here except for the routines planet_getHeliocentric() and
 * sky2_nutationIAU2000B() are covered by the MIT license, as follows:
 *
 * \copyright
 * \parblock
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
 * \endparblock
 *
**  =====================
 *
 *  The routines planet_getHeliocentric() and sky2_nutationIAU2000B() are
 *  covered by the SOFA license, as follows.
 *  See the description of each of those functions, which contains the necessary
 *  text to comply with section 3 of the license below:
 *
 * \copyright
 * \verbatim
**
**  Copyright (C) 2017
**  Standards Of Fundamental Astronomy Board
**  of the International Astronomical Union.
**
**  =====================
**  SOFA Software License
**  =====================
**
**  NOTICE TO USER:
**
**  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
**  CONDITIONS WHICH APPLY TO ITS USE.
**
**  1. The Software is owned by the IAU SOFA Board ("SOFA").
**
**  2. Permission is granted to anyone to use the SOFA software for any
**     purpose, including commercial applications, free of charge and
**     without payment of royalties, subject to the conditions and
**     restrictions listed below.
**
**  3. You (the user) may copy and distribute SOFA source code to others,
**     and use and adapt its code and algorithms in your own software,
**     on a world-wide, royalty-free basis.  That portion of your
**     distribution that does not consist of intact and unchanged copies
**     of SOFA source code files is a "derived work" that must comply
**     with the following requirements:
**
**     a) Your work shall be marked or carry a statement that it
**        (i) uses routines and computations derived by you from
**        software provided by SOFA under license to you; and
**        (ii) does not itself constitute software provided by and/or
**        endorsed by SOFA.
**
**     b) The source code of your derived work must contain descriptions
**        of how the derived work is based upon, contains and/or differs
**        from the original SOFA software.
**
**     c) The names of all routines in your derived work shall not
**        include the prefix "iau" or "sofa" or trivial modifications
**        thereof such as changes of case.
**
**     d) The origin of the SOFA components of your derived work must
**        not be misrepresented;  you must not claim that you wrote the
**        original software, nor file a patent application for SOFA
**        software or algorithms embedded in the SOFA software.
**
**     e) These requirements must be reproduced intact in any source
**        distribution and shall apply to anyone to whom you have
**        granted a further right to modify the source code of your
**        derived work.
**
**     Note that, as originally distributed, the SOFA software is
**     intended to be a definitive implementation of the IAU standards,
**     and consequently third-party modifications are discouraged.  All
**     variations, no matter how minor, must be explicitly marked as
**     such, as explained above.
**
**  4. You shall not cause the SOFA software to be brought into
**     disrepute, either by misuse, or use for inappropriate tasks, or
**     by inappropriate modification.
**
**  5. The SOFA software is provided "as is" and SOFA makes no warranty
**     as to its use or performance.   SOFA does not and cannot warrant
**     the performance or results which the user may obtain by using the
**     SOFA software.  SOFA makes no warranties, express or implied, as
**     to non-infringement of third party rights, merchantability, or
**     fitness for any particular purpose.  In no event will SOFA be
**     liable to the user for any consequential, incidental, or special
**     damages, including any lost profits or lost savings, even if a
**     SOFA representative has been advised of such damages, or for any
**     claim by any third party.
**
**  6. The provision of any version of the SOFA software under the terms
**     and conditions specified herein does not imply that future
**     versions will also be made available under the same terms and
**     conditions.
*
**  In any published work or commercial product which uses the SOFA
**  software directly, acknowledgement (see www.iausofa.org) is
**  appreciated.
**
**  Correspondence concerning SOFA software should be addressed as
**  follows:
**
**      By email:  sofa@ukho.gov.uk
**      By post:   IAU SOFA Center
**                 HM Nautical Almanac Office
**                 UK Hydrographic Office
**                 Admiralty Way, Taunton
**                 Somerset, TA1 2DN
**                 United Kingdom
**
 * \endverbatim
 */

/*! \page page-c-style C coding style
 *
 *  Code here largely conforms to [this style guide](../c-style-guide.pdf)
 */

/*! \page page-var-suffixes Suffixes on variable and parameter names
 *
 *  In accordance with the naming convention spelt out in the
 *  [style guide](../c-style-guide.pdf),
 *  many variables and constants have a suffix indicating their
 *  units. Mostly these are SI units, for example
 *      -  _d (days), _h (hours) _s (seconds) _ns (nanoseconds)
 *      -  _m (metres), _km (kilometres), _kmps (kilometres/second)
 *      -  _rad (radian), _deg (degrees), _radps (radian/second)
 *      -  _degC (degrees Celsius)
 *
 *  But we also have units for astronomical quantities as recommended by the
 *  International Astronomical Union (IAU) style manual
 *  (https://www.iau.org/publications/proceedings_rules/units/), for example
 *      -  _as (arcseconds) _mas (milliarcseconds), (see \ref page-sec-arcsec)
 *      -  _au (Astronomical Unit (of distance))
 *      -  _a (Julian year)
 *
 *  This last unit implies that we could use "ka" for millennia, and we do. But
 *  we also have quantities measured in Julian centuries. The IAU guidelines
 *  specifically rule out using "ha" for this (quite right: "ha" means
 *  hectares), and they oppose "cy", but they make no recommendation of an
 *  alternative. But the _Astronomical Almanac_ does use "cy" for a Julian
 *  century, so I will use suffix _cy.
 *      -  _ka (Julian millennium), _cy (Julian century)
 *
 *  Another convention that applies here: Objects whose name ends in a capital V
 *  are 3-dimensional vectors (often unit position vectors, or direction
 *  cosines) but not exclusively so. Objects whose name ends in a capital M are
 *  matrices, typically 3x3 rotation matrices. Occasionally, this is combined
 *  with a suffix. So a name like \c earthV_au means a vector whose components
 *  are scaled in Astronomical Units.
 */

/*! \page page-sec-arcsec Seconds versus arcseconds
 *  The fact that time in hours is subdivided into minutes and seconds, and
 *  angles in degrees are also subdivided into minutes and seconds, normally
 *  does not cause confusion. But astronomers measure angles on the celestial
 *  sphere using both units of time (for Right Ascension and Hour Angles) and
 *  units of degrees (for many other angles). So if we talk about an angle of
 *  n minutes, or n seconds, do we mean a fraction of an hour, or a degree?
 *
 *  For this reason, astronomers refer to the fractions of a degree as
 *  arcminutes, and reserve the term "minutes" for fractions of an hour. And
 *  fractions of an arcminute are called arcseconds rather than seconds. This
 *  program follows that convention, particularly in the use of suffixes on
 *  variable names (see page \ref page-var-suffixes). Variables in units of
 *  seconds of time have the suffix _s, those in units of arcseconds have the
 *  suffix _as.
 *
 *  In general, within this program most angles are stored in units of radians,
 *  and only converted to degrees-arcminutes-arcseconds, or
 *  hours-minutes-seconds when being written out in text form.
 *  The symbol ′ is used for arcminutes only, never for minutes of time.
 *  The symbol ″ is used for arcseconds only, never for seconds of time.
 */

/*! \page page-interval-notation Interval notation
 *  Frequently throughout the documentation and comments, valid ranges of
 *  function parameters are expressed in interval notation. This is a pair of
 *  numbers representing the two ends of the range, separated by a comma. The
 *  the pair of numbers is surrounded by brackets. Square brackets
 *  indicate an inclusive end-of-range, and round brackets (parentheses)
 *  indicate an exclusive end-of range.
 *
 *  Examples:
 *  - [0,10] means all values from 0 to 10 inclusive
 *  - [0.0, 360.0) means all values from 0.0 up to but not including 360.0
 *      itself.
 *  - (-π, +π] means all values from -π, but not including -π itself, up to and
 *      including +π.
 */

/*! \page page-sign-conventions Sign and direction conventions
 *  The following conventions apply throughout this software.
 *      - Terrestrial coordinates (latitude and longitude).
 *          - Polar form\n
 *            Longitudes east of Greenwich are positive, those west are
 *            negative.
 *            Latitudes north of the equator are positive, those south are
 *            negative. When the Earth is viewed from above the north pole
 *            looking towards the centre of the earth, longitude increases in an
 *            anticlockwise direction.
 *          - Rectangular form\n
 *            The position is indicated by a vector, whose components are as
 *            follows:
 *              - the x axis points from the centre of the earth towards the
 *                equator at 0° longitude,
 *              - the y axis points to 90° longitude, and
 *              -  the z axis points to the north pole.
 *              .
 *            This is a right-handed set.
 *
 *      - Equatorial coordinates (Right Ascension and declination).
 *          - Polar form\n
 *            The equator is the projection of the Earth's equator out into
 *            space, and the north celestial pole is the projection of the north
 *            pole out into space. When viewed from north celestial pole looking
 *            down towards the origin, Right Ascension increases in an
 *            anticlockwise direction (i.e. eastwards, like longitude), and
 *            north declinations are positive, south are negative.
 *          - Rectangular form\n
 *              - x points to the March Equinox (0° RA),
 *                  (\ref page-vernal-equinox "unfortunately" known as the
 *                  vernal equinox).
 *              - y points to 90° RA, and
 *              - z points to the north celestial pole.
 *              .
 *            (This is also a right-handed set.)
 *
 *      - Ecliptic coordinates. (Ecliptic latitude and longitude).\n
 *        These are exactly analogous to terrestrial coordinates, in that
 *        longitude increases anticlockwise, and the rectangular coordinate
 *        vector is a right-handed set.
 *
 *      - Horizon coordinates (azimuth and elevation (or altitude)).
 *          - Polar form\n
 *            Azimuth increases in a clockwise direction from True North (just
 *            like a compass bearing).
 *          - Rectangular form\n
 *              - x points to the True North horizon,
 *              - y  points the East horizon, and
 *              - z points to the zenith.
 *              .
 *            This is a left-handed set. (Compared to a right-handed set, the x
 *            and y axes are swapped.)
 *
 *      - (Hour angle and declination).
 *          - Polar form\n
 *            Hour angles are measured positive westward. This means it
 *            increases clockwise, unlike Right Ascension.
 *          - Rectangular form\n
 *              - x points to the intersection of the local meridian and the
 *                celestial equator,
 *              - y points to HA = 90°, and
 *              - z points to the north celestial pole.
 *              .
 *            Like Horizon coordinates, this is a left-handed set.
 */

/*! \page page-vernal-equinox Vernal equinox
 *  The equinox that occurs in March is normally referred to as the vernal
 *  equinox. This is unfortunate, because of the meaning of the word vernal:
 * \par
 *      \b ver′nal \em adj. of or occurring in spring.\n
 *      -- The Collins Concise Dictionary, 4th edn. 1999
 *
 *  Obviously it never occurred to the person who coined this term for the
 *  equinox that the Southern Hemisphere exists. Because \em everybody in the
 *  Southern Hemisphere knows that the spring equinox occurs in September, not
 *  March.
 */

