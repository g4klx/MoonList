#ifndef PLANET_H
#define PLANET_H
/*============================================================================*/
/*! \file
 * \brief planet.h - Astronomical routines to get the positions of planets
 *
 * \author  David Hoadley, based on a routine from the IAU's SOFA collection
 *
 * \details
 *          Routines to obtain a specified planet's position and convert that
 *          position to apparent coordinates and/or topocentric coordinates.
 *
 *==============================================================================
 */
/*
 *
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
#include "sky.h"
#include "sky1.h"
#include "vectors3d.h"

/*
 * Global #defines and typedefs
 */


/*
 * Global functions available to be called by other modules
 */
void planet_setCurrent(int np);

void planet_getApp2(double             t_cy,
                    int                np,
                    const Sky1_Nut1980 *nut,
                    V3D_Vector *appV,
                    double     *dist_au);
void planet_getApparent(double j2kTT_cy, Sky_TrueEquatorial *pos);
void planet_getTopocentric(double             j2kUtc_d,
                           const Sky_DeltaTs  *deltas,
                           const Sky_SiteProp *site,
                           Sky_SiteHorizon *topo);

void planet_getGeocentric(double t_cy,
                          int np,
                          V3D_Vector *p2V,
                          double *dist_au);

int planet_getHeliocentric(double t_cy,
                           int np,
                           V3D_Vector *j2kV_au,
                           V3D_Vector *velV_aupd);

/*
 * Global variables accessible by other modules
 */

#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif

#endif /* PLANET_H */

