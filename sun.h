#ifndef SUN_H
#define SUN_H
/*============================================================================*/
/*!\file
 * \brief sun.h - routines to calculate the Sun's position
 *
 * \author  David Hoadley
 *
 * \details
 *      Alternative routines which calculate the sun's position with differing
 *      levels of accuracy. The alternatives are
 *      - a low precision routine, taken from the _Astronomical Almanac_
 *      - a much more detailed, high precision routine, based on the Solar
 *        Position Algorithm, as published by the National Renewable Energy
 *        Laboratory (referred to as the NREL SPA throughout this software)
 * 
 *      Rectangular coordinates are used wherever possible, to minimise the
 *      unnecessary recalculation of trigonometric functions.
 *
 * \par References:
 *  The _Astronomical Almanac_, 2007, page C24\n
 *  Reda, Ibrahim and Andreas, Afshin,
 *  _Solar Position Algorithm for Solar Radiation Applications._
 *  National Renewable Energy Laboratory publication no. NREL/TP-560-34302,
 *  June 2003, revised 2008.
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
#include "sky0.h"
#include "vectors3d.h"

/*
 * Global #defines and typedefs
 */


/*
 * Global functions available to be called by other modules
 */
#ifdef __cplusplus
extern "C" {
#endif

void sun_aaApparentApprox(double n,
                          V3D_Vector *appV,
                          double     *dist_au);

void sun_nrelApp2(double             t_cy,
                  const Sky0_Nut1980 *nut,
                  V3D_Vector *appV,
                  double     *dist_au);
void sun_nrelApparent(double j2kTT_cy, Sky_TrueEquatorial *pos);
void sun_nrelTopocentric(double             j2kUtc_d,
                         const Sky_DeltaTs  *deltas,
                         const Sky_SiteProp *site,
                         Sky_SiteHorizon *topo);

double sun_solarNoon(int                year,
                     int                month,
                     int                day,
                     const Sky_DeltaTs  *deltas,
                     const Sky_SiteProp *site,
                     Sky_SiteHorizon *topo);
double sun_riseSet(int                year,
                   int                month,
                   int                day,
                   bool               getSunrise,
                   const Sky_DeltaTs  *deltas,
                   const Sky_SiteProp *site,
                   Sky_SiteHorizon *topo);

/*
 * Global variables accessible by other modules
 */
/*      (none) */

#ifdef __cplusplus
}
#endif

#endif /*SUN_H*/
