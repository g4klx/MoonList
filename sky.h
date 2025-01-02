#ifndef SKY_H
#define SKY_H
/*============================================================================*/
/*!\file    
 * \brief sky.h - structures and routines for astronomical observing & tracking
 *
 * \author  David Hoadley
 *
 * \details
 *      This collection is in two parts:
 *          - time routines. Routines for handling the different timescales used
 *            in astronomy, and for calculating the rotational orientation of
 *            the Earth. The timescales involved are TT, UT1 and UTC
 *            (see \ref page-timescales). Also included here is a routine to set
 *            polar motion parameters.
 *          - site routines.  Data that is specific to an observing site on the
 *            Earth, and routines to convert astronomical coordinates to
 *            site-specific (i.e. topocentric) coordinates at the site. More
 *            than one site may be supported simultaneously, if required.
 * 
 *      The routines are designed to provide an efficient implementation of
 *      the necessary calculations. When combined with the routines in the
 *      skyfast.h and skyfast.c module, they enable accurate tracking with
 *      even a small processor.
 *
 *      Other modules required:
 *          - sun.h (and sun.c) if you want to track the Sun, using either the
 *            NREL Solar Position Algorithm, or a simplified approximate one.
 *          - moon.h (and moon.c) if you want to track the Moon, using the
 *            NREL Moon Position Algorithm
 *          - (future) a module for tracking other planets
 *          - (future) a module for tracking stars
 *          - skyio.h (and skyio.c) for input-output routines: conversions of
 *            angles to and from strings in sexagesimal (e.g. degrees, minutes
 *            and seconds) form; writing out a day number as a date & time.
 *          - vectors3d.h (and vectors3d.c) which implements the rectangular
 *            matrix and vector operations used by these algorithms,
 *
 *      You will need to separately include one of the following modules:
 *          - sky0.h (and sky0.c) for nutation routines and sidereal time
 *            routines suitable for tracking the Sun or the Moon using the
 *            NREL Solar Position Algorithm and Moon Position Algorithm
 *          - (future) sky1.h (and sky1.c) for IAU 1980 precession, nutation and
 *            sidereal time routines, suitable for tracking stars
 *          - (future) sky2.h (and sky2.c) for IAU 2000 precession, nutation and
 *            sidereal time routines, suitable for tracking stars
 *
 *==============================================================================
 */
/*
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
 */
/*------------------------------------------------------------------------------
 * Notes:
 *      Character set: UTF-8. (Non ASCII characters appear in this file)
 *      Things you might want to edit: definition of macro INCLUDE_MJD_ROUTINES
 *                                     definition of macro POSIX_SYSTEM
 *----------------------------------------------------------------------------*/

#include <time.h>

#include "general.h"
#include "vectors3d.h"

/*
 * Global #defines and typedefs
 */

/*! Struct used for holding an object's coordinates in equatorial apparent or
 *  Intermediate coordinates.
 *  Apparent coordinates are those referred to the true equator and equinox of
 *  the time indicated in field #timestamp_cy. 
 *  Intermediate coordinates are in the Celestial Intermediate Reference System
 *  (CIRS), referred to the true equator of time #timestamp_cy and to the
 *  Celestial Intermediate Origin (CIO) instead of the equinox.
 *  If the object is in Apparent coordinates, the Equation of the Equinoxes
 *  (#eqEq_rad) field is required as part of converting to topocentric
 *  coordinates. If the object is in CIRS coordinates, field #eqEq_rad can be
 *  ignored. 
 *
 *  This structure is returned by the sun_nrelApparent(), moon_nrelApparent(),
 *  (future)planet_getApparent() and (future)star_getApparent() functions, and
 *  importantly for tracking, the skyfast_getApprox() function.*/
typedef struct {
    double      timestamp_cy;   /*!< Time applying to the other figures in
                                     this struct (centuries since J2000.0, TT
                                     timescale) */
    V3D_Vector  appCirsV;       /*!< Direction of object in apparent or CIRS 
                                     coordinates (effectively a unit vector). */
    double      distance_au;    /*!< Distance to object (Astronomical Units) or
                                     0.0 for far distant objects (that is, those
                                     with negligible diurnal parallax) */
    double      eqEq_rad;       /*!< Equation of the Equinoxes (radian). */
} Sky_TrueEquatorial;

/*!     Coordinates of a celestial object in the horizon frame, in both
        rectangular and polar forms.
        - The rectangular coordinate vector has the orientation
           +  x points to North horizon,
           +  y points to East horizon,
           +  z points to the zenith.
          (This is a left-handed set.)
        - In the polar form, the azimuth is measured clockwise from North
          (i.e. east = +90°).
        - This combined rectangular & polar form is used because it happens to
          be most convenient to calculate both forms at once. */
typedef struct {
    V3D_Vector rectV;         /*!< unit vector in horizon coordinates */
    double     azimuth_rad;   /*!< azimuth (radian) */
    double     elevation_rad; /*!< elevation (or altitude) (radian) */
} Sky_SiteHorizon;


#ifdef __cplusplus
extern "C" {
#endif


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *                      The TIME routines and structs
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*
 * Global #defines and typedefs
 */

/*      Un-comment the following line if you want routines that convert calendar
        dates and times to and from Modified Julian Date (MJD) format (where 
        MJD = JD - 2 400 000.5).  */
/*--- #define INCLUDE_MJD_ROUTINES ---*/

/*      Un-comment the following line if you want to use either of the routines
        that use the struct timespec data type (available on POSIX-compliant
        systems). These are sky_unixTimespecToJ2kd() & sky_unixTimespecToMjd()*/
/*--- #define POSIX_SYSTEM ---*/

/*!     This structure contains relatively constant data, and is set up by one
        of the three functions sky_initTime(), sky_initTimeSimple() or
        sky_initTimeDetailed(). The data which can vary
        is not expected to vary any more frequently than once per day or even
        less. Do not modify any of the fields in this structure directly; use
        the routines in this file to make all changes. In general, you won't
        need to access any of the individual fields here. */
typedef struct {
    double   deltaUT_d;     //!< UT1 - UTC, scaled to days
    double   deltaT_d;      //!< TT - UT1, scaled to days
    double   deltaTT_d;     //!< TT - UTC, scaled to days
} Sky_DeltaTs;

/*!     This structure contains the continuously varying time (and earth
        rotation) data, in various forms that we will find useful. 

        Do not modify any of these fields directly - use
        the sky_updateTimes() function to update them.
 
        But you will almost certainly want to read the values of any of the
        individual fields, and/or pass them to functions. In particular, the
        field #j2kTT_cy is passed to many routines.

        Notation used in section B of the _Astronomical Almanac_ 2007 is shown
        in [square brackets]. */
typedef struct {
    double     mjdUtc;    //!< Modified Julian Date (= JD - 2 400 000.5), UTC
    double     j2kUT1_d;  //!< days since J2000.0, UT1 timescale           [Du]
    double     j2kTT_d;   //!< days since J2000.0, TT timescale             [D]
    double     j2kTT_cy;  //!< Julian centuries since J2000.0, TT timescale [T]
    double     era_rad;   //!< Earth Rotation Angle (radian)                [θ]
} Sky_Times;

/*!     This structure contains polar motion parameters and a rotation
        matrix. Do not modify any of these fields directly - use the
        sky_setPolarMotion() function to do that. In general, you won't
        need to access any of the individual fields here. */
typedef struct {
    bool       correctionInUse; // if false, polar motion is being ignored
    double     xPolar_as;       // polar motion in x (arcseconds)
    double     yPolar_as;       // polar motion in y (arcseconds)
    V3D_Matrix corrM;           // polar motion rotation matrix
} Sky_PolarMot;

/*
 * Global functions available to be called by other modules
 */
/*      1. Initialise this module by choosing one of the following 3 fns */
void sky_initTime(int deltaAT_s, double deltaUT_s, Sky_DeltaTs *d);
void sky_initTimeSimple(Sky_DeltaTs *d);
void sky_initTimeDetailed(double mjdUtc,
                          double usnoMjdBase,
                          double usnoCoeffC11,
                          double usnoCoeffC12,
                          int    deltaAT_s,
                          Sky_DeltaTs *d);

/*      2. Get the time in "J2KD" form from one of the following functions. */
double sky_calTimeToJ2kd(int year, int month, int day,
                         int hour, int minute, double second,
                         double tz_h);
double sky_unixTimeToJ2kd(time_t unixTime);
#ifdef POSIX_SYSTEM
double sky_unixTimespecToJ2kd(struct timespec uTs);
#endif

#ifdef INCLUDE_MJD_ROUTINES
/*      2b. Alternatively get time in MJD form from one of the following fns */
double sky_calTimeToMjd(int year, int month, int day,
                        int hour, int minute, double second,
                        double tz_h);
double sky_unixTimeToMjd(time_t unixTime);
#ifdef POSIX_SYSTEM
double sky_unixTimespecToMjd(struct timespec uTs);
#endif
#endif


/*      3. Call the following to update the other astronomical times */
void sky_updateTimes(double            j2kUtc_d,
                     const Sky_DeltaTs *d,
                     Sky_Times *t);

#ifdef INCLUDE_MJD_ROUTINES
/*      3b. Update the other astronomical times */
void sky_updateTimesFromMjd(double            mjdUtc,
                            const Sky_DeltaTs *d,
                            Sky_Times *t);
#endif


/*      4. Call routines (from other modules, not this one) to obtain the
        position of your chosen star, planet or the sun in apparent coordinates
        (i.e. w.r.t true equator and equinox at the time of interest). As part
        of that process, you will have obtained the equation of the equinoxes.*/

/*      5. Call one of the following to update the rotational position of the
        Earth and convert your apparent position in to coordinates in the
        Terrestrial Intermediate Reference system. Choose one of
            sky0_appToTirs()   - if you are tracking the Sun or Moon using
                                  the NREL SPA or SAMPA functions.
            astc1_appToTirs()   - if you are using the IAU 1980 nutation
                                  routines etc.
            astc2_appToTirs()   - if you are using the IAU 2000 routines
*/

/*      6. Call sky_siteTirsToTopo() for each observing site or solar panel
        location to convert the TIRS coordinates to topocentric coordinates at
        the site or sites. */
 

/*      Repeat 2, 3, 4, 5 and 6 in a loop at whatever rate suits you */

/*      6. From time to time, (no more than daily) call the following to update
        polar motion parameters, if you need to take polar motion into account*/
void sky_setPolarMotion(double xPolar_as,
                        double yPolar_as,
                        double t_cy,
                        Sky_PolarMot *polar);
/*      6b. If you have called the routine above, you must also call the routine
        sky_adjustSiteForPolarMotion() for every observing site you are
        calculating positions for */

/*      Other routines */
void sky_j2kdToCalTime(double j2k_d,
                       int    *year,
                       int    *month,
                       int    *day,
                       int    *hour,
                       int    *minute,
                       double *second);

#ifdef INCLUDE_MJD_ROUTINES
void sky_mjdToCalTime(double mjd,
                      int    *year,
                      int    *month,
                      int    *day,
                      int    *hour,
                      int    *minute,
                      double *second);
#endif


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *                      The SITE routines and structs
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*
 * Global #defines and typedefs
 */
/*!     Site properties. Declare one object of the following type for each site
        that you want to calculate sky positions for (typically one site). Do
        not modify any of the fields in this structure directly; use the
        routines in this file to make all changes. In general, you won't need
        to access any of the individual fields here (except possibly
        #timeZone_d) */
typedef struct {
    double     astLat_rad;    //!< Astronomical latitude of site (ϕA) (radian)
    double     astLong_rad;   //!< Astronomical longitude of site (radian)
    double     geocRadius_km; //!< Approx geocentric radius of site (≈ ae*ρ)(km)
    double     rhoSin_au;     //!< ae*ρ*sin(ϕ - ϕ′) geocentre-to-site x (AU)
    double     rhoCos_au;     //!< -ae*ρ*cos(ϕ - ϕ′) geocentre-to-site z (AU)
    double     diurnalAberr;  //!< Diurnal aberration: caused by earth rotation
    double     refracPT;      //!< Refraction correction: pressure & temperature
    double     timeZone_d;    //!< time zone offset from UTC (fraction of a day)
    V3D_Matrix *azElM;        //!< points to either azElPolM or azElBaseM
    V3D_Matrix azElPolM;      //!< rotation matrix from TIRS to Az/El coords
    V3D_Matrix azElBaseM;     //!< as above, but excluding polar motion correctn
    V3D_Matrix haDecM;        //!< rotation matrix from Az/El to HA/Dec coords
} Sky_SiteProp;


/*
 * Global functions available to be called by other modules
 */
/*      To initialise the Sky_SiteProp structure, call one of the following two
        sky_setSiteLocation() functions for each site, and then call 
        sky_setSiteTempPressure() and sky_setSiteTimeZone() as required */
void sky_setSiteLocation(double latitude_deg,
                         double longitude_deg,
                         double height_m,
                         Sky_SiteProp *site);
void sky_setSiteLoc2(double astLat_deg,
                     double astLong_deg,
                     double geodLat_deg,
                     double geodLong_deg,
                     double height_m,
                     Sky_SiteProp *site);
void sky_setSiteTempPress(double temperature_degC,
                          double pressure_hPa,
                          Sky_SiteProp *site);
void sky_setSiteTimeZone(double timeZone_h,
                         Sky_SiteProp *site);

/*      Call the following if you have a surface for which you want to calculate
        the incidence angle of the rays from the celestial object being tracked.
        Typically this surface would be a solar panel, and you want the angle
        at which the Sun's rays are striking the surface. */
void sky_setupSiteSurface(double azimuth_deg,
                          double slope_deg,
                          Sky_SiteHorizon *surface);

/*      Call the following if you have called sky_setPolarMotion() to change
        polar motion parameters. */
void sky_adjustSiteForPolarMotion(const Sky_PolarMot *polar,
                                  Sky_SiteProp *site);


/*      You will need to call one or more of the following functions each time
        around your main loop, to convert to site-specific coords. (Typically
        you will call sky_siteTirsToTopo() to convert from terrestrial
        intermediate coordinates to topocentric coordinates.) */
void sky_siteTirsToTopo(const V3D_Vector   *terInterV,
                        double             dist_au,
                        const Sky_SiteProp *site,
                        Sky_SiteHorizon *topo);
void sky_siteAzElToHaDec(const V3D_Vector   *topoV,
                         const Sky_SiteProp *site,
                         double *hourAngle_rad,
                         double *declination_rad);
double sky_siteIncidence_rad(const V3D_Vector *topoV,
                             const V3D_Vector *surfaceV);


/*
 * Global variables accessible by other modules
 */

#ifdef __cplusplus
}
#endif

/*! \page page-sky-h Edits you may want to make to sky.h
 *
 *  Define the macro INCLUDE_MJD_ROUTINES if you have use for any of the
 *  routines sky_calTimeToMjd(), sky_unixTimeToMjd(), sky_unixTimespecToMjd(),
 *  sky_updateTimesFromMjd(), or sky_mjdToCalTime(). None of these routines are
 *  essential for tracking celestial objects.
 * 
 *  Define the macro POSIX_SYSTEM if you are using a system that supports the
 *  POSIX standard (a common Unix standard) and you have a use for either of the
 *  routines sky_unixTimespecToJ2kd() or sky_unixTimespecToMjd()
 */

#endif /* SKY_H */

