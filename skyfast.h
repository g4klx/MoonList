#ifndef SKYFAST_H
#define SKYFAST_H
/*============================================================================*/
/*! \file
 * \brief
 * skyfast.h - set up and use interpolation for rapid calculation of a celestial
 *             object's apparent coordinates.
 *
 * \author  David Hoadley
 *
 * \details
 *          Routines to set up interpolation for celestial tracking, get the
 *          interpolated position at a given time, and to update the endpoints
 *          used by the interpolation algorithm. The error introduced by using
 *          interpolation rather than fully calculating each and every position
 *          can be very small - see \ref page-interpolation (the end of this
 *          source file)
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
#include "sky.h"

/*
 * Global #defines and typedefs
 */


#ifdef __cplusplus
extern "C" {
#endif
/*
 * Global functions available to be called by other modules
 */
void skyfast_init(double            tStartUtc_d,
                  int               fullRecalcInterval_mins,
                  const Sky_DeltaTs *deltas,
                  void (*getApparent)(double j2kTT_cy, Sky_TrueEquatorial *pos)
                  );
void skyfast_backgroundUpdate(void);
void skyfast_getApprox(double t_cy, Sky_TrueEquatorial *approx);

/*
 * Global variables accessible by other modules
 */

#ifdef __cplusplus
}
#endif

/*! \page page-interpolation Interpolation and its errors
 *  The interpolation process obtains an estimate of the apparent position (or 
 *  alternatively, the celestial intermediate position) of a celestial object
 *  by interpolating between previously calculated position vectors. The amount
 *  of error introduced by this process depends of course on how much the
 *  position changed between the two previously calculated times. Here is a
 *  table of the errors when this process is used to obtain the position of the
 *  Sun, Moon, planets and a star. The errors are calculated from
 *      sqrt((azimuth_error * cos(elevation))^2 + elevation_error^2)
 *
 *  The interpolation interval (in hours) is the time between the two previously
 *  fully calculated vectors. This value (multiplied by 60) is supplied to
 *  parameter \a fullRecalcInterval_mins of function skyfast_init().
 *  
 *  ###Maximum absolute position error (arcseconds) for different interpolation intervals (hours)

    hours   |1      |2      |3      |6      |9      |12     |15     |18     |24    |
    :-------|------:|------:|------:|------:|------:|------:|------:|------:|-----:|
    Sun     |0.0005 |0.0020 |0.0046 |0.0184 |0.0413 |0.0736 |0.1151 |0.1658 |0.2957|
    -       |       |       |       |       |       |       |       |       |      |
    Moon    |0.3381 |1.3536 |3.0455 |12.07  |27.62  |49.55  |77.74  |112.1  |185.4 |
    -       |       |       |       |       |       |       |       |       |      |
    Mercury |0.1533 |0.6132 |1.3797 |5.5177 |12.42  |22.07  |34.44  |49.64  |88.17 |
    Venus   |0.0324 |0.1295 |0.2914 |1.1656 |3.2375 |4.6621 |7.2847 |10.49  |18.64 |
    Mars    |0.0115 |0.0459 |0.1032 |0.4127 |0.9287 |1.6510 |2.5796 |3.7146 |6.6035|
    Jupiter |0.0025 |0.0101 |0.0227 |0.0907 |0.2040 |0.3627 |0.5667 |0.8161 |1.4508|
    Saturn  |0.0013 |0.0052 |0.0117 |0.0468 |0.1054 |0.1874 |0.2927 |0.4216 |0.7494|
    Uranus  |0.0007 |0.0027 |0.0062 |0.0247 |0.0556 |0.0988 |0.1544 |0.2223 |0.3951|
    Neptune |0.0005 |0.0018 |0.0041 |0.0163 |0.0366 |0.0651 |0.1017 |0.1465 |0.2604|
    -       |       |       |       |       |       |       |       |       |      |
    Antares |<0.0001|<0.0001|<0.0001|0.0001 |0.0003 |0.0005 |0.0007 |0.0010 |0.0019|

 *  As can be seen, the errors are very small for the Sun even with 24 hours
 *  between full calculations (i.e. a value of 1440 minutes supplied to
 *  parameter \a fullRecalcInterval_mins of function skyfast_init()).
 *  Actually a full twenty-four hours of tracking of the Sun can be obtained by
 *  supplying 720 minutes to this parameter, without ever having to call the
 *  function skyfast_backgroundUpdate(). This is because skyfast_init()
 *  calculates three apparent positions, not just two.
 *
 *  Errors for the Moon are larger. The position algorithm itself is accurate
 *  to about 5 arcseconds, it seems, so this suggests that an interpolation
 *  interval of no more than about 2 hours should be used for the Moon.
 *
 *  The maximum errors tracking planets are larger than for the Sun. These
 *  maximum errors seem to occur around the time the planet enters or leaves
 *  apparent retrograde motion.
 *
 *  The errors when tracking a star are almost negligible.
 *
 *  In all cases, the error appears to approximately quadruple for every
 *  doubling of the interpolation interval.
 */

#endif /* SKYFAST_H */
