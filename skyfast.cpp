/*==============================================================================
 * skyfast.c - set up and use interpolation for rapid calculation of a celestial
 *             object's apparent coordinates.
 *
 * Author:  David Hoadley
 *
 * Description: (see skyfast.h)
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

/* ANSI includes etc. */
#include <math.h>

/* Local and project includes */
#include "skyfast.h"

#include "sky0.h"
#include "astron.h"
#include "general.h"
#include "sun.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       // For use by REQUIRE() - assertions.

/*      --- Definitions that you may need to or wish to modify --- */
//#define BARE_METAL_THREADS
//#define POSIX_THREADS
//#define NO_THREADS
#if defined(BARE_METAL_THREADS)
#define startCriticalSection()      disableInterrupts()
#define endCriticalSection()        enableInterrupts()

#elif defined(POSIX_THREADS)
#include <pthread.h>

#define startCriticalSection()      pthread_mutex_lock(&mutex)
#define endCriticalSection()        pthread_mutex_unlock(&mutex)

#else /* Must be NO_THREADS */
#define startCriticalSection()      ((void)0)
#define endCriticalSection()        ((void)0)

#endif



/*
 * Prototypes for local functions (not called from other modules)
 */


/*
 * Global variables accessible by other modules
 */


/*
 * Local variables (not accessed by other modules)
 */
LOCAL Sky_TrueEquatorial  lfiA, lfiB, lfiC;
LOCAL Sky_TrueEquatorial  *last = &lfiA;    // items calculated for time in past
LOCAL Sky_TrueEquatorial  *next = &lfiB;    // items calculated for time ahead
LOCAL Sky_TrueEquatorial  *oneAfter = &lfiC;// ditto for time after next
LOCAL volatile bool       oneAfterIsValid = false;

LOCAL void (*callback)(double j2kTT_cy, Sky_TrueEquatorial *pos);
LOCAL double        recalcInterval_cy = 0.0; // time between full recalculations

#ifdef POSIX_THREADS
LOCAL pthread_mutex_t mutex;
#endif

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
GLOBAL void skyfast_init(double            tStartUtc_d,
                         int               fullRecalcInterval_mins,
                         const Sky_DeltaTs *deltas,
                         void (*getApparent)(double j2kTT_cy,
                                             Sky_TrueEquatorial *pos)
                         )
/*! Initialise those items that take a long time to calculate, but which do not
    need to be recalculated frequently. This routine calls a function that you
    supply to calculate the apparent coordinates of a celestial object, its
    distance, and the Equation of the Equinoxes, as derived from nutation
    calculations. This routine calls that function
        1.  for the time specified by \a tStartUtc_d,
        2.  for time \a tStartUtc_d + \a fullRecalcalcInterval_mins, and
        3.  for time \a tStartUtc_d + 2 x \a fullRecalcalcInterval_mins.
        .
    For example, to track the Sun, specify the function sun_nrelApparent() when
    calling this routine. To track the Moon, specify the function
    moon_nrelApparent().

    The routine skyfast_getApprox() can then be called
    (at a high frequency if required) to calculate the current position of the
    object, using these values to do it.
 \param[in]  tStartUtc_d  Time for first full calculation using function
                          \a getApparent(). UTC time in "J2KD" form - i.e days
                          since J2000.0 (= JD - 2 451 545.0)
 \param[in]  fullRecalcInterval_mins
                          Interval of time between full recalculation of the
                          object's position using the function supplied to
                          \a getApparent (minutes). This value must be greater
                          than zero. (Otherwise you will get a precondition
                          failure.)
 \param[in]  deltas       Delta T values, as set by the sky_initTime() (or
                          sky_initTimeSimple() or sky_initTimeDetailed())
                          routines
 \param      getApparent  Function to get the position of a celestial object in
                          apparent coordinates (i.e. referred to the true
                          equinox and equator at the time), in rectangular form.

 \par When to call this function
    At program initialisation time.
 \note
    Although the parameter name \a getApparent suggests that you must supply a
    function which calculates an object's position in apparent coordinates, it
    is also possible to supply a function which returns the position in
    Celestial Intermediate coordinates (i.e. referred to the true equator and
    the Celestial Intermediate Origin (CIO) at time \a t_cy) instead. If so,
    the function does not need to fill in the \a eqEq_rad field of struct
    Sky_TrueEquatorial.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_Times   atime;              // time, in various timescales
    double      calcTimeTT_cy;
#ifdef POSIX_THREADS
    int         ret;
#endif

    REQUIRE_NOT_NULL(deltas);
    REQUIRE_NOT_NULL(getApparent);
    REQUIRE(fullRecalcInterval_mins > 0);


#ifdef POSIX_THREADS
    ret = pthread_mutex_init(&mutex, NULL);
    ASSERT(ret == 0);   // There is no possible recovery from an error here.
#endif
    /* Save the function address for later call by skyfast_backgroundUpdate() */
    callback = getApparent;

    sky_updateTimes(tStartUtc_d, deltas, &atime);

    /* Save the recalculation rate, converted from minutes to centuries. */
    recalcInterval_cy = fullRecalcInterval_mins / (1440.0 * JUL_CENT);

    calcTimeTT_cy = atime.j2kTT_cy;
    getApparent(calcTimeTT_cy, last);

    /* Now do the same for the next time (e.g. next hour) */
    calcTimeTT_cy += recalcInterval_cy;
    getApparent(calcTimeTT_cy, next);

    /* And again for the time after */
    calcTimeTT_cy += recalcInterval_cy;
    getApparent(calcTimeTT_cy, oneAfter);
    oneAfterIsValid = true;
}



GLOBAL void skyfast_backgroundUpdate(void)
/*! Recalculation of the low frequency quantities.
    Checks whether the calculations of the celestial object's apparent
    coordinates and the equation of the equinoxes have been performed for the
    time called "oneAfter" (i.e. the time after the time called "next"). If they
    have not, this function performs those calculations.

    This function calls the function that you previously supplied to function
    skyfast_init() (parameter \a getApparent).

 \par When to call this function
    This function is designed to be called in a low priority loop, or at
    background level, using any available leftover processor time to slowly
    precalculate values for the "time after next". It needs to have been called
    before the high frequency routine skyfast_getApprox() (which is
    interpolating between time "last" and time "next") arrives at time "next"
    and therefore needs to access the data for time "oneAfter".
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double  t_cy;

    REQUIRE(recalcInterval_cy > 0.0);  // skyfast_init() not called before now?

    if (!oneAfterIsValid) {
        t_cy = next->timestamp_cy + recalcInterval_cy;
        callback(t_cy, oneAfter);

        startCriticalSection();
        oneAfterIsValid = true;
        endCriticalSection();
    }
}



GLOBAL void skyfast_getApprox(double t_cy,
                              Sky_TrueEquatorial *approx)
/*! Get the best approximation to the celestial object's apparent coordinates
    and distance, and the equation of the equinoxes, based on an interpolation
    between two sets of such data that we have previously calculated.
 \param[in]  t_cy     Julian centuries since J2000.0, TT timescale. This must
                      specify a time no earlier than the time specified in
                      argument \a tStartUtc_d in the call to skyfast_init().
 \param[out] approx   position vector, distance, etc, obtained by interpolation

    Although the position is described as approximate, the position returned can
    be very accurate, as shown in \ref page-interpolation, depending on the
    interpolation interval that was specified to routine skyfast_init().

 \note
    This function will return an approximate position in apparent coordinates
    (i.e. referred to the true equator and equinox at time \a t_cy) if the
    function that you passed to the skyfast_init() function returned apparent
    coordinates. This function will return an approximate position in
    Celestial Intermediate coordinates (i.e. referred to the true equator and
    the Celestial Intermediate Origin (CIO)) at time \a t_cy) if the
    function that you passed to the skyfast_init() function returned CIRS
    coordinates.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    Sky_TrueEquatorial *temp;
    double a;
    double b;

    if (t_cy > next->timestamp_cy) {
        /* Time t_cy is no longer between last and next, so we need to make next
           and oneAfter become the new last and next respectively. But this
           requires that our low frequency/low priority routine has completed
           filling in all the data for oneAfter. */
        REQUIRE(oneAfterIsValid);

        startCriticalSection();
        temp = last;
        last = next;
        next = oneAfter;
        oneAfter = temp;
        oneAfterIsValid = false;
        endCriticalSection();
    }

    /* It is a programming error if time t_cy is not between last and next */
    REQUIRE(t_cy >= last->timestamp_cy);
    REQUIRE(next->timestamp_cy >= t_cy);

    if ((next->timestamp_cy - last->timestamp_cy) < SFA) {
        a = 0.0;
        b = 1.0;
    } else {
        a = (t_cy - last->timestamp_cy)
            / (next->timestamp_cy - last->timestamp_cy);
        b = 1.0 - a;
    }
    /* Do a simple linear interpolation between the two appCirsV position
     * vectors last and next. Unlike the two appCirsV vectors, the resulting
     * vector will not be exactly of unit magnitude. But if the two appCirsV
     * vectors are less than one degree apart, the resulting position error is
     * very small (< 0.3′). If the two appCirsV are a few arcminutes
     * apart, the magnitude error of the resulting vector is negligible.  */
    approx->appCirsV.a[0] =  a * next->appCirsV.a[0] + b * last->appCirsV.a[0];
    approx->appCirsV.a[1] =  a * next->appCirsV.a[1] + b * last->appCirsV.a[1];
    approx->appCirsV.a[2] =  a * next->appCirsV.a[2] + b * last->appCirsV.a[2];
    /* And a linear interpolation of the other two quantities also. */
    approx->distance_au = a * next->distance_au + b * last->distance_au;
    approx->eqEq_rad    = a * next->eqEq_rad + b * last->eqEq_rad;
}


/*
 *------------------------------------------------------------------------------
 *
 * Local functions (not called from other modules)
 *
 *------------------------------------------------------------------------------
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/*! \page page-skyfast-c Edits you may want to make to skyfast.c
 *
 *  There are four different ways that you might want to use the skyfast
 *  module. It can be configured to run in two threads or one. To
 *  control this, define one or none of the following three macros. (Don't
 *  define more than one.) They are BARE_METAL_THREADS, POSIX_THREADS or
 *  NO_THREADS.
 *
 *  The four different approaches are:
 *      1. The simplest approach is to define NO_THREADS (or not to define any
 *          macro at all), and not to call the skyfast_backgroundUpdate()
 *          routine at all. The limitation is that you will only be able to
 *          track your object for a period that is twice as long as the number
 *          that you pass to the \a fullRecalcInterval_mins parameter of routine
 *          skyfast_init(), in minutes. If you try to track longer than this,
 *          the program will abort.
 *      2. The next simplest is once again to define NO_THREADS (or not to
 *          define any macro at all), and to put calls to
 *          skyfast_backgroundUpdate() in the same loop as your calls to
 *          skyfast_getApprox(). You will benefit from the faster execution of
 *          skyfast_getApprox(), but every \a fullRecalcInterval_mins minutes,
 *          a full recalculation will take place. That is, the loop will take
 *          much longer to execute than it does all the rest of the time.
 *      3. If the occasional long loop execution time is not acceptable to you,
 *          (say, it might upset your tracking control loop), consider using two
 *          threads. In this approach, you set up a high-frequency,
 *          high-priority task to drive the control calculations. This task
 *          calls skyfast_getApprox(). A background task calls
 *          skyfast_backgroundUpdate() to perform the slow calculations at very
 *          low priority. To use this approach, define either BARE_METAL_THREADS
 *          or POSIX_THREADS.
 *
 *          BARE_METAL_THREADS is intended for small embedded processors without
 *          an operating system. You set your high-priority task to be triggered
 *          by a timer interrupt. To control access to shared data, this makes
 *          use of two other macros, which you will need to define for your
 *          processor. They are called \c disableInterrupts() and
 *          \c enableInterrupts(), and it should be obvious from those names
 *          what it is that you need to make them do.
 *      4. The fourth approach is basically the same as the third, but it
 *          enables implementation on processors running a POSIX-compliant
 *          operating system (such as Linux). You will need to define the
 *          POSIX_THREADS macro, and then create the posix
 *          threads yourself; this macro simply causes this module to use the
 *          pthreads Mutex mechanism to control access to shared data.
 */
