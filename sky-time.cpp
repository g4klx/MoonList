/*==============================================================================
 * sky-time.c - astronomical time routines
 *
 * Author:  David Hoadley
 *
 * Description: (see the "Time routines" sections of sky.h)
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
#include <math.h>
#include <time.h>

/* Local and project includes */
#include "sky.h"

#include "astron.h"
#include "general.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       // For use by REQUIRE() - assertions.

#ifdef INCLUDE_MJD_ROUTINES
#define MJD_1JAN1970    40587
#endif
#if 0
#define MJD_1JAN1904    16480
#define MJD_0JAN1900    15019
#define MJD_30DEC1899   15018
#define MJD_BASE        2400000.5   /* 17-Nov-1858, as a Julian Date */
#endif

#define START_GREGORIAN 2299160.0   /* Start Gregorian calendar 15-Oct-1582 */
#define TIME_T_J2000    946728000   /* J2000.0 in Unix/C time_t format (s) */

/*
 * Prototypes for local functions (not called from other modules).
 */
LOCAL double getDeltaUT1_s(double mjdUtc,
                           double usnoMjdBase,
                           double usnoCoeffC11,
                           double usnoCoeffC12);
LOCAL double calendarToJ2kd(int year, int month, int day);

/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */

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
GLOBAL void sky_initTime(int deltaAT_s, double deltaUT_s, Sky_DeltaTs *d)
/*! This is one of three alternative routines for setting up the various delta
    time values for use in ongoing calculations. (The other two are
    sky_initTimeDetailed() and sky_initTimeSimple()). For an explanation of
    these delta times, see \ref page-timescales.
 \param[in]  deltaAT_s  (= TAI - UTC). Cumulative number of leap seconds
                         (seconds)
 \param[in]  deltaUT_s  (= UT1 - UTC) (seconds). Valid range [-0.9,+0.9]
 \param[out] d          Fields initialised as follows:
                        - \a d->deltaUT_d initialised to \a deltaUT_s / 86400
                        - \a d->deltaT_d  initialised to
                            (\a deltaAT_s + 32.184 - \a deltaUT_s) / 86400
                        - \a d->deltaTT_d initialised to
                            (\a deltaAT_s + 32.184) / 86400

 \par When to call this function
    Call at program startup time.
    Use this routine if you have opted for option 3A, 3B or 3D in the note at
    the bottom of this file (see \ref page-timescales).
    (Use sky_initTimeDetailed() if you have opted for option 3C.)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    const double deltaTAI_s = 32.184;   // TT - TAI (seconds)
    double deltaTT_s;                   // TT - UTC (seconds)
    double deltaT_s;                    // TT - UT1 (seconds)

    REQUIRE_NOT_NULL(d);

    deltaTT_s = deltaTAI_s + deltaAT_s;
    d->deltaTT_d = deltaTT_s / 86400.0;

    deltaT_s = deltaTT_s - deltaUT_s;
    d->deltaT_d = deltaT_s / 86400.0;
    d->deltaUT_d = deltaUT_s / 86400.0;
}



GLOBAL void sky_initTimeDetailed(double mjdUtc,
                                 double usnoMjdBase,
                                 double usnoCoeffC11,
                                 double usnoCoeffC12,
                                 int    deltaAT_s,
                                 Sky_DeltaTs *d)
/*! This is one of three alternative routines for setting up the various delta
    time values for use in ongoing calculations. (The other two are
    sky_initTime() and sky_initTimeSimple()). This routine uses a prediction
    formula to calculate the value of delta_UT, rather than having it specified
    directly. For an explanation of these delta times, see \ref page-timescales.
 \param[in]  mjdUtc    Modified Julian Date (= JD - 2 400 000.5), UTC timescale
 \param[in]  usnoMjdBase,
             usnoCoeffC11,
             usnoCoeffC12
                       US Naval Observatory prediction formula values. These
                       three values are obtained from "Bulletin A" of the
                       International Earth Rotation Service (IERS), updated at
                       the IERS (& USNO) website from time to time
                       (see \ref page-timescales)
 \param[in]  deltaAT_s (= TAI - UTC). Cumulative number of leap seconds
                        (seconds)
 \param[out] d         Fields initialised as follows:
                        - \a d->deltaUT_d initialised using USNO prediction
                          formula
                        - \a d->deltaT_d  initialised to
                            (\a deltaAT_s + 32.184) / 86400 - \a d->deltaUT_d
                        - \a d->deltaTT_d initialised to
                            (\a deltaAT_s + 32.184) / 86400

 \par When to call this function
    Call this function at program startup time.
    Use this routine if you have opted for option 3C in the note at the bottom
    of this file (see \ref page-timescales). (Use sky_initTime() instead, if you
    have opted for option 3A, 3B or 3D.) If your program runs uninterrupted for
    many days, you may wish to call this function once per day with an updated
    value for \a mjdUtc.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double deltaUT_s;                   // UT1 - UTC (seconds)

    deltaUT_s = getDeltaUT1_s(mjdUtc, usnoMjdBase, usnoCoeffC11, usnoCoeffC12);
    sky_initTime(deltaAT_s, deltaUT_s, d);
}



GLOBAL void sky_initTimeSimple(Sky_DeltaTs *d)
/*! This is one of three alternative routines for setting up the various delta
    time values for use in ongoing calculations. (The other two are
    sky_initTime() and sky_initTimeDetailed()). This function initialises the
    delta times with simnple default values. For an explanation of
    these delta times, see \ref page-timescales.
 \param[out] d    Fields initialised as follows:
                        - \a d->deltaUT_d initialised to 0.0
                        - \a d->deltaT_d  initialised to (37 + 32.184) / 86400
                        - \a d->deltaTT_d initialised to (37 + 32.184) / 86400

 \par When to call this function
    Call this function at program startup time.
    Use this routine if you don't need the high accuracy of sub-second time. But
    if you do need it, call sky_initTime() or sky_initTimeDetailed() instead.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    /* Use the 2019 value of deltaAT_s. This value will slowly go out of date,
       but in practice, accuracy will hardly be affected.
       Assume deltaUT = 0.0 */
    sky_initTime(37, 0.0, d);
}



GLOBAL void sky_updateTimes(double            j2kUtc_d,
                            const Sky_DeltaTs *d,
                            Sky_Times *t)
/*! Convert the given "J2KD" in the UTC timescale to the other timescales, and
    pre-calculate some other quantities
 \param[in]  j2kUtc_d  Date in "J2KD" form - i.e. the number of days elapsed
                       since 2000 Jan 1.5 (= JD - 2 451 545.0), UTC timescale,
                       as returned by sky_calTimeToJ2kd(), sky_unixTimeToJ2kd()
                       or sky_unixTimespecToJ2kd().
 \param[in]  d         The various delta T values as set by one of
                       sky_initTime(), sky_initTimeDetailed() or
                       sky_initTimeSimple()
 \param[out] t         All fields are updated

 \par Reference
    _Astronomical Almanac_ 2007, page B9 (for Earth Rotation Angle)

 \par When to call this function
    Call this routine before calling any routine that calculates a celestial
    position.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(d);
    REQUIRE_NOT_NULL(t);

    t->mjdUtc = j2kUtc_d + MJD_J2000;

    t->j2kUT1_d = j2kUtc_d + d->deltaUT_d;
    t->j2kTT_d  = j2kUtc_d + d->deltaTT_d;
    t->j2kTT_cy = t->j2kTT_d / JUL_CENT;
    t->era_rad = (0.7790572732640 + 1.00273781191135488 * t->j2kUT1_d) * TWOPI;
}



GLOBAL double sky_calTimeToJ2kd(int year, int month, int day,
                                int hour, int minute, double second,
                                double tz_h)
/*! Return the number of days (and fraction of a day) since noon 2000 Jan 1
 *  (UTC) of the given calendar date and time. This is one of three alternative
 *  routines that return a J2KD - the other two are sky_unixTimeToJ2kd() and
 *  sky_unixTimespecToJ2kd().
 \returns days since Julian date 2 451 545.0, UTC timescale
 \param[in] year, month, day calendar date
 \param[in] hour             valid range [0, 23]
 \param[in] minute           valid range [0, 59]
 \param[in] second           valid range [0.0, 60.0)
 \param[in] tz_h             time zone offset (hours), positive for zones
                             east of Greenwich (e.g. Australian time zones
                             AEST = +10.0, ACST = +9.5, AEDT = +11.0;
                             UTC = 0.0;
                             USA time zones are negative)
 \par When to call this routine
    Call this (or one of the alternative routines) before each call to
    sky_updateTimes().
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double jd2k;
    double timeOfDay_d;

    jd2k = calendarToJ2kd(year, month, day);

    timeOfDay_d = (second + 60.0 * (minute + 60.0 * (hour - tz_h))) / 86400.0;
    jd2k += timeOfDay_d;
    return jd2k;
}



GLOBAL double sky_unixTimeToJ2kd(time_t unixTime)
/*! Convert a time in Unix system time format to days since 2000 Jan 1, noon
    UTC. This is one of three alternative
    routines that return a J2KD - the other two are sky_calTimeToJ2kd() and
    sky_unixTimespecToJ2kd().
 \returns  days since Julian Date 2 451 545.0, UTC timescale
 \param[in] unixTime - time in C standard \c time_t (or unix) format
                        ((sort of) seconds since 1-Jan-1970 UTC)

 \note This routine has no better resolution than one second. Use routine
   sky_unixTimespecToMjd() or sky_calTimeToMjd() instead, if you need
   sub-second resolution.
 \par When to call this routine
    Call this (or one of the alternative routines) before each call to
    sky_updateTimes().
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    return (double)(unixTime - TIME_T_J2000) / 86400.0;
}



#ifdef POSIX_SYSTEM
GLOBAL double sky_unixTimespecToJ2kd(struct timespec uTs)
/*! Convert a time in Unix timespec format to days since 2000 Jan 1, noon UTC.
    This is one of three alternative routines that return a J2KD - the other two
    are sky_calTimeToJ2kd() and sky_unixTimeToJ2kd().
 \returns  days since Julian Date 2 451 545.0, UTC timescale
 \param[in] uTs - time in "timespec" format, a combination of (sort of) seconds
                  since 1-Jan-1970 UTC, and nanoseconds of the current second.
                  This is the form returned by the POSIX \c clock_gettime()
                  function.
 \note
    The macro POSIX_SYSTEM must be defined at compile time to use this function.
 \par When to call this routine
    Call this (or one of the alternative routines) before each call to
    sky_updateTimes().
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    return ((double)(uTs.tv_sec - TIME_T_J2000) + (double)uTs.tv_nsec / 1e9)
            / 86400.0;
}
#endif



GLOBAL void sky_j2kdToCalTime(double j2k_d,
                              int *year,
                              int *month,
                              int *day,
                              int *hour,
                              int *minute,
                              double *second)
/*! This procedure converts the integral part of the date in "J2KD" form to a
    calendar date and the fractional part to hours, minutes and seconds. If the
    returned date is from 1582-10-15 or later, it is a Gregorian calendar date.
    If the returned date is 1582-10-04 or earlier, it is a Julian calendar date.
 \param[in] j2k_d   Days since J2000.0 (= Julian Date - 2 451 545.0)
                    Valid range: j2k_d >= -2447065, otherwise incorrect results

 \param[out] year, month, day      calendar date
 \param[out] hour, minute, second  time of day

 \par Reference
    The algorithm is based on a method of D A Hatcher, _Quarterly Journal of the
    Royal Astronomical Society_ 1984, Vol 25, pp 53-55. It is valid for dates
    after JD = 4480 (7 April 4701 BC, Julian calendar).

 \note To obtain calendar date and time for any timezone other than UTC, add
    the time zone offset (in units of fraction of a day) to \a mjd before
    calling this routine. The \a timezone_d field of the site properties struct
    (type Sky_SiteProp) provides this value.
 \par
    If \a second turns out to be within half a millisecond of the next round
    minute, this routine rounds time upwards.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int32_t j;
    int32_t n4;
    int32_t n10;
    double  timeOfDay;
    double  j2kdate;

    j2k_d += 0.5;       // Change from elapsed since noon to elapsed since 00:00
    j2kdate = floor(j2k_d);
    timeOfDay = (j2k_d - j2kdate) * 24.0;
    *hour = (int)timeOfDay;
    timeOfDay = (timeOfDay - *hour) * 60.0;
    *minute = (int)timeOfDay;
    *second = (timeOfDay - *minute) * 60.0;

    /* Round up if within half a millisecond of a round minute */
    if ((int)(*second + 0.0005) == 60) {
        *second = 0.0;
        (*minute)++;
        if (*minute == 60) {
            *minute = 0;
            (*hour)++;
            if (*hour == 24) {
                *hour = 0;
                j2kdate++;
            }
        }
    }


    j = (int32_t)j2kdate + 2451545;

    if (j <= START_GREGORIAN) {
        n4 = 4 * j;
    } else {
        n4 = 4 * (j + ((((4 * j - 17918) / 146097) * 3 + 2) >> 2) - 37);
    }

    n10 = 10 * (((n4 - 237) % 1461) >> 2) + 5;

    *year = (int)(n4 / 1461 - 4712);
    *month = (n10 / 306 + 2) % 12 + 1;
    *day = (n10 % 306) / 10 + 1;
}



#ifdef INCLUDE_MJD_ROUTINES
GLOBAL double sky_calTimeToMjd(int year, int month, int day,
                               int hour, int minute, double second,
                               double tz_h)
/*! Return the Modified Julian Date (= Julian Date - 2 400 000.5) of a calendar
    date and time.
 \returns the Modified Julian Date, UTC timescale
 \param[in] year, month, day calendar date
 \param[in] hour             valid range [0, 23]
 \param[in] minute           valid range [0, 59]
 \param[in] second           valid range [0.0, 60.0)
 \param[in] tz_h             time zone offset (hours), positive for zones
                             east of Greenwich (e.g. Australian time zones
                             AEST = +10.0, ACST = +9.5, AEDT = +11.0;
                             UTC = 0.0;
                             USA time zones are negative)
 \note
    The macro INCLUDE_MJD_ROUTINES must be defined at compile time to use this
    function.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double mjd;
    double timeOfDay_d;

    mjd = calendarToJ2kd(year,month, day) + MJD_J2000;

    timeOfDay_d = (second + 60.0 * (minute + 60.0 * (hour - tz_h))) / 86400.0;
    mjd += timeOfDay_d;
    return mjd;
}



GLOBAL double sky_unixTimeToMjd(time_t unixTime)
/*! Convert a time in Unix system time format to Modified Julian Date
 \returns  the Modified Julian Date (= Julian Date - 2 400 000.5), UTC timescale
 \param[in] unixTime - time in C standard \c time_t (or unix) format
                        ((sort of) seconds since 1-Jan-1970 UTC)

 \note This routine has no better resolution than one second. Use routine
   sky_unixTimespecToMjd() or sky_calTimeToMjd() instead, if you need
   sub-second resolution
 \note
    The macro INCLUDE_MJD_ROUTINES must be defined at compile time to use this
    function.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    return (double)unixTime / 86400.0 + MJD_1JAN1970;
}



#ifdef POSIX_SYSTEM
GLOBAL double sky_unixTimespecToMjd(struct timespec uTs)
/*! Convert a time in Unix timespec format to Modified Julian Date
 \returns  the Modified Julian Date (= Julian Date - 2 400 000.5), UTC timescale
 \param[in] uTs - time in "timespec" format, a combination of (sort of) seconds
                  since 1-Jan-1970 UTC, and nanoseconds of the current second.
                  This is the form returned by the POSIX \c clock_gettime()
                  function.
 \note
    The macros INCLUDE_MJD_ROUTINES and POSIX_SYSTEM must be defined at compile
    time to use this function.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    return ((double)uTs.tv_sec + (double)uTs.tv_nsec / 1e9) / 86400.0
            + MJD_1JAN1970;
}
#endif



GLOBAL void sky_updateTimesFromMjd(double mjdUtc,
                                const Sky_DeltaTs *d,
                                Sky_Times *t)
/*! Convert the given MJD in the UTC timescale to the other timescales, and
    pre-calculate some other quantities
 Inputs
 \param[in]  mjdUtc  Modified Julian Date (= JD - 2 400 000.5), UTC timescale
 \param[in]  d       The various delta T values as set by one of
                     sky_initTime(), sky_initTimeDetailed() or
                     sky_initTimeSimple()
 \param[out] t       All fields except \a gmst_rad and \a gast_rad are updated

 \par Reference
    _Astronomical Almanac_ 2007, page B9 (for Earth Rotation Angle)
 \note
    The macro INCLUDE_MJD_ROUTINES must be defined at compile time to use this
    function.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double mjdUT1;  // Modified Julian Date (= JD - 2 400 000.5), UT1 timescale
    double mjdTT;   // MJD as above, TT timescale

    REQUIRE_NOT_NULL(d);
    REQUIRE_NOT_NULL(t);

    t->mjdUtc = mjdUtc;
    mjdUT1    = mjdUtc + d->deltaUT_d;
    mjdTT     = mjdUtc + d->deltaTT_d;

    t->j2kUT1_d = mjdUT1 - MJD_J2000;
    t->j2kTT_d  = mjdTT - MJD_J2000;
    t->j2kTT_cy = t->j2kTT_d / JUL_CENT;
    t->era_rad = (0.7790572732640 + 1.00273781191135488 * t->j2kUT1_d) * TWOPI;
}



GLOBAL void sky_mjdToCalTime(double mjd,
                             int *year,
                             int *month,
                             int *day,
                             int *hour,
                             int *minute,
                             double *second)
/*! This procedure converts the integral part of the Modified Julian Date to a
    calendar date and the fractional part to hours, minutes and seconds. If the
    returned date is from 1582-10-15 or later, it is a Gregorian calendar date.
    If the returned date is 1582-10-04 or earlier, it is a Julian calendar date.
 \param[in] mjd     - Modified Julian Date (= Julian Date - 2 400 000.5)
                      Valid range: mjd > -2395520, otherwise incorrect results

 \param[out] year, month, day     - calendar date
 \param[out] hour, minute, second - time of day

 \par Reference
    The algorithm is based on a method of D A Hatcher, _Quarterly Journal of the
    Royal Astronomical Society_ 1984, Vol 25, pp 53-55. It is valid for dates
    after MJD = -2395520 (1st March 4701 BC).

 \note To obtain calendar date and time for any timezone other than UTC, add
    the time zone offset (in units of fraction of a day) to \a mjd before
    calling this routine. The \a timezone_d field of the site properties struct
    (type Sky_SiteProp) provides this value.

 \note
    The macro INCLUDE_MJD_ROUTINES must be defined at compile time to use this
    function.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int32_t j;
    int32_t n4;
    int32_t n10;
    double  timeOfDay;
    double  mjdate;

    mjdate = floor(mjd);
    timeOfDay = (mjd - mjdate) * 24.0;
    *hour = (int)timeOfDay;
    timeOfDay = (timeOfDay - *hour) * 60.0;
    *minute = (int)timeOfDay;
    *second = (timeOfDay - *minute) * 60.0;

    j = (int32_t)mjdate + 2400001;

    if (j <= START_GREGORIAN) {
        n4 = 4 * j;
    } else {
        n4 = 4 * (j + (((4 * j - 17918) / 146097) * 3 + 2) / 4 - 37);
    }

    n10 = 10 * (((n4 - 237) % 1461) / 4) + 5;

    *year = (int)(n4 / 1461 - 4712);
    *month = (n10 / 306 + 2) % 12 + 1;
    *day = (n10 % 306) / 10 + 1;
}
#endif



GLOBAL void sky_setPolarMotion(double xPolar_as,
                               double yPolar_as,
                               double t_cy,
                               Sky_PolarMot *polar)
/*! Update polar motion correction.
 \param[in] xPolar_as   Polar motion in x (arcseconds)
 \param[in] yPolar_as   Polar motion in y (arcseconds)
 \param[in] t_cy        Julian centuries since J2000.0, TT timescale
 \param[out] polar      Polar motion params and rotation matrix, calculated from
                            R1(-y) * R2(-x) * R3(s')

 \par Reference
    _Astronomical Almanac_ 2007, page B80

 \par When to call this function
    This routine needs to be called only when polar motion parameters
    change - which is no more frequently than once per day. This effect is so
    small that it can be ignored altogether with only a tiny loss of accuracy.

 \note Every time this routine is called, the routine
    sky_adjustSiteForPolarMotion() must be called for every site. (Typically
    this will be only one site.)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double      sPrime_as;          // Secular drift
    V3D_Matrix  yM, xM, sM, tempM;  // temporary rotation matrices

    REQUIRE_NOT_NULL(polar);

    polar->xPolar_as = xPolar_as;
    polar->yPolar_as = yPolar_as;

    if ((xPolar_as == 0.0) && (yPolar_as == 0.0)) {
        /* Yes I know one is not supposed to compare floats for equality, but
           these two values are entered directly by a user, not calculated. So
           being exactly zero is quite likely. */

        polar->correctionInUse = false;

    } else {
        polar->correctionInUse = true;

        /* Secular drift of TIO. This is so tiny it can be ignored altogether.
           But if we are to take it into account (as we do here) it only needs
           to be updated very infrequently. */
        sPrime_as = -0.000047 * t_cy;

        /* matrix is R1(-y) * R2(-x) * R3(s') */
        v3d_createRotationMatrix(&sM, Zaxis, arcsecToRad(sPrime_as));
        v3d_createRotationMatrix(&xM, Yaxis, arcsecToRad(-xPolar_as));
        v3d_createRotationMatrix(&yM, Xaxis, arcsecToRad(-yPolar_as));
        v3d_multMxM(&polar->corrM, &yM, v3d_multMxM(&tempM, &xM, &sM));
    }
}



/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules).
 *
 *------------------------------------------------------------------------------
 */
LOCAL double getDeltaUT1_s(double mjdUtc,
                           double usnoMjdBase,
                           double usnoCoeffC11,
                           double usnoCoeffC12)
/* Calculate deltaUT1 using the prediction formula from the US Naval Observatory
   deltaUT1 is the difference between UTC and the earth's actual rotational
   position, and it is subject to variations which are difficult to predict.
   Therefore this prediction is approximate.
 Returns         - deltaUT1 (= UT1 - UTC) (seconds)
 Inputs
    mjdUtc       - Modified Julian Date (= JD - 2 400 000.5), UTC timescale
    usnoMjdBase  - US Naval Observatory prediction formula values. These three
    usnoCoeffC11      values are updated at the USNO website from time to time
    usnoCoeffC12
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    const double A = 0.022;     // Markowitz time series constants (units = s)
    const double B = -0.012;
    const double C = -0.006;
    const double D = 0.007;
    const double bYear = TROP_CENT / 100.0;     // No of days in Besselian year
                                                // (= tropical year at B1900.0)
    double n;
    double mjDateBYR;
    double delBYR;
    double seasonalVariation;
    double deltaUT_s;

    /* Calculate the MJD of the beginning of the current Besselian year and
       the fractional part of the Besselian year to date */
    n = floor((mjdUtc - MJD_B1900) / bYear);
    mjDateBYR = MJD_B1900 + n * bYear;
    delBYR = (mjdUtc - mjDateBYR) / bYear;

    /* Calculate the seasonal variation of UT in seconds, using the standard
       Markowitz time series (this is UT2 - UT1) */
    seasonalVariation = A * sin(TWOPI * delBYR)
                        + B * cos(TWOPI * delBYR)
                        + C * sin(4.0 * PI * delBYR)
                        + D * cos(4.0 * PI * delBYR);

    /* Calculate the current value of Delta UT1 (= UT1 - UTC) in seconds from
       the seasonal variations and the USNO prediction formula.
       UT1 - UTC = (UT2 - UTC) - Seasonal_Variation */
    deltaUT_s = usnoCoeffC11 + usnoCoeffC12 * (mjdUtc - usnoMjdBase)
                 - seasonalVariation;

    return deltaUT_s;
}



LOCAL double calendarToJ2kd(int year, int month, int day)
/*  Return the number of days elapsed since Greenwich noon on 2000-01-01 of a
    calendar date, as at midnight at the beginning of the day, UTC.
 Returns - the days elapsed since Julian date 2 451 545.0 (2000 Jan 1.5).
           Because of the way this is defined, the result will always be an
           integer plus 0.5
 Inputs
    year, month, day - calendar date

 Reference
    The algorithm is inspired by a method of D A Hatcher, _Quarterly Journal of
    the Royal Astronomical Society_ 1984, Vol 25, pp 53-55.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
   /* Force 32-bit integer arithmetic, in case type int is only a 16-bit type */
    int32_t year32 = year;
    int32_t month32 = month - 3;
    int32_t day32 = day;
    int32_t j2kdi_d;        // Integer part of the J2KD date.


    if (month32 < 0) {
        month32 += 12;
        year32--;
    }
    j2kdi_d = ((1461 * (year32 + 4712)) >> 2) - 2451487;
    j2kdi_d += (306 * month32 + 5) / 10 + day32;
    if (j2kdi_d > -152386) {
        /* Supplied date must have been later than 4-Oct-1582. So use Gregorian
           calendar */
        j2kdi_d -= ((3 * (year32 / 100 - 11)) >> 2) + 7;
    }
    return (double)j2kdi_d + 0.5;
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/*! \page page-timescales   Timescales, and converting between them

    There are five timescales of interest:
    - TAI - International Atomic Time, as maintained by a network of atomic
            clocks
    - TT  - Terrestrial Time. Used for astronomical observations. It is designed
            to be free of the irregularities of the rotation of the earth. It is
            based on TAI, and for historical reasons has the following
            relationship
              TT = TAI + 32.184 seconds
    - UT1 - Universal Time. This is the timescale that essentially represents
            the rotation of the earth, and therefore is subject to the
            fluctuations that occur in that rotation. (Also called simply UT)
    - UTC - Coordinated Universal Time. This is the time we all use. It is
            based on TAI, but occasionally a leap second is inserted (or, in
            theory, removed) to keep UTC within 0.9 seconds of UT1.
    - local civil time - This is UTC + (local timezone offset)

    What this means for astronomical observations is that we need to obtain the
    time in the UTC timescale (either directly, or from the local civil time)
    and calculate TT and UT1 from it.
    TT is needed to calculate the celestial positions of objects in the sky, and
    UT1 is needed for converting that celestial coordinate into a direction as
    seen from a site on the surface of the earth, as it rotates.

    We need three pieces of information.
    1. ΔTAI = TT - TAI = 32.184 s (a fixed value, as described above)
    2. ΔAT = TAI - UTC.  This is the number of leap seconds that have been
            added in the past. This figure is available from the International
            Earth Rotation Service (IERS) at the Observatoire de Paris, via
            its regular Bulletin C, available at
            https://hpiers.obspm.fr/iers/bul/bulc/bulletinc.dat
            When accessed in August 2019, this value was 37 seconds. (Note that
            Bulletin C expresses the relation as UTC - TAI, so it said -37 s.)
    3. ΔUT = UT1 - UTC.   This is available from the IERS, at the US Naval
            Observatory website, in its Bulletin A.
            https://maia.usno.navy.mil/ser7/ser7.dat or from
            https://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html
            It also gives the leap seconds, this time in the form TAI - UTC.
            The value of deltaUT varies from day to day, but because its value
            is always within the range -0.9..+0.9, you have four choices as to
            what to do about it.
        - \b A. Ignore it altogether, and assume it is equal to zero. This will
                have no effect on the calculation of celestial coordinates, but
                could introduce an error of up to 0.00375° in the direction from
                the ground, if you are unlucky. For most purposes, this is good
                enough. Initialise your code with the function
                sky_initTimeSimple(), or call sky_initTime() with \a deltaUT_s
                set to 0.0
        - \b B. Check Bulletin A every few months, and call sky_initTime() with
                \a deltaUT_s set to the value of DUT1 given near the top of the
                file, described as "transmitted with time signals". This will
                give a better approximation.
        - \b C. Check Bulletin A every few months, and find the PREDICTIONS
                section, and look for a line that looks like the following:

                    UT1-UTC = -0.1733 - 0.00056 (MJD - 58725) - (UT2-UT1)
                              ^^^^^^^ ^^^^^^^^^        ^^^^^
                             CoeffC11  CoeffC12       MJDbase

                Take the three numbers and pass them to the function
                sky_initTimeDetailed()
                which will use this formula to predict the value. This method
                is good enough be used by large astronomical telescopes.
        - \b D. Check Bulletin A every week, and use the values in the table for
                each day. Almost nobody should need to do this.
        .
    There are two more such quantities, but we can calculate them from the above
    4. ΔT = TT - UT1. This is also obtainable from
            ΔT = ΔTAI + ΔAT - ΔUT.
    5. ΔTT = TT - UTC. This is also obtainable from
            ΔTT = ΔTAI + ΔAT.

    Once you have set the above with sky_initTime(), sky_initTimeSimple() or
    sky_initTimeDetailed(), then calling sky_updateTimes() will convert a
    time in the UTC timescale to various times in the other timescales, ready
    for use.

    Pursuit of really accurate UT1 (as outlined in 3C and 3D above) really only
    makes sense if your system clock is accurately synchronised to UTC.

    But wait! There's more! If you want really accurate coordinates, you have to
    take into account another variation in the earth's rotation. The polar axis
    itself moves around a tiny amount. This is called polar motion. If this
    matters to you, the values of polar motion can be obtained from the same
    Bulletin A described above. But for normal purposes, you can just enter zero
    for these parameters (or just not call sky_setPolarMotion() at all).
*/

/*
        ~~~~~~~~~~~~~~
Eh? 2007 edition of Astronomical almanac, page B8, suggests an error of 1 s in
deltaT (and therefore in deltaUT) here gives an error of 1.5e-6 arcseconds or
0.1e-6 seconds of time in the sidereal time calculations. This is clearly a
factor of 10^7 smaller than I calculated above. What have I done wrong?
Maybe it is referring to an error in TT when calculating GMST, rather than an
error in Du. If Du has the error, it gives an error in earth rotation angle.

Yes, on further examination, this is assuming you have the correct UT1, but
your TT is wrong. That is, deltaT is wrong, but deltaUT is correct, so it
must be deltaTT (or deltaAT) which is wrong by some number of whole seconds.
        ~~~~~~~~~~~~~~
*/

/*! \page page-time-variables Time variables
 *  Astronomical algorithms typically use a simple number as the time variable.
 *  This can take one of several forms:
 *      - Julian Date (JD). The number of days since Greenwich noon, 1 January
 *        4713 BC (on the Julian calendar, as if it had existed back then). So
 *        for example 2020-04-13 9:00 am is 2458952.875
 *      - Modified Julian Date (MJD). The Julian Date minus 2 400 000.5. This
 *        makes it the number of days since Greenwich midnight at the start of
 *        1858-11-17. So for example, 2020-04-13 9:00 am is 58952.375
 *      - J2KD (My abbreviation). Days since Greenwich noon, 2000-01-01.
 *        Equals JD - 2 451 545.0 (= MJD - 51544.5). So for example,
 *        2020-04-13 9:00 am is 7407.875
 *
 *  The three above could easily represent times in the TT or UT1 timescales,
 *  so the relevant timescale should be specified. (See \ref page-timescales)
 *      .
 *      - J2KC (My abbreviation). Julian centuries since Greenwich noon,
 *        2000-01-01. Equals J2KD ÷ 36525
 *      - J2KM (My abbreviation). Julian millennia since Greenwich noon,
 *        2000-01-01. Equals J2KC ÷ 10
 *
 *  Typically the two above will be used for times in the TT timescale.
 *      .
 *      - Julian epoch (e.g. J2000.0 or J2020.5). Calculated from
 *        J[2000.0 + J2KD / 365.25], TT timescale. This is used
 *        for star catalogue positions. For example, J2000.0 is JD 2 451 545.0.
 *        And 2-Jul-2020 3:00:00 (TT) is J2020.5
 *      - Besselian epoch (e.g. B1950.0). This is calculated from
 *        B[1900.0 + JD - (2 415 020.313 52) / 365.242 198 781].
 *        So B1950.0 is JD 2 433 282.423. Besselian epochs are often associated
 *        with FK4 catalogue positions. These are not supported by this
 *        software (yet).
 *
 *  Many astronomical algorithms use one of the forms based on time since
 *  J2000.0 (i.e. either J2KD, J2KC or J2KM) as their time variable.
 *  For this reason, the conversion routines from calendar date/time
 *  (sky_calTimeToJ2kd())and from the unix \c time_t and \c timespec formats
 *  (sky_unixTimeToJ2kd() and sky_unixTimespecToJ2kd()) return a
 *  variable in the J2KD form (rather than the JD or MJD forms).
 */

