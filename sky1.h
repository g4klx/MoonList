#ifndef SKY1_H
#define SKY1_H
/*============================================================================*/
/*!\file
 * \brief sky1.h - astronomical coordinate conversion routines, IAU 1980
 *
 * \author  David Hoadley
 *
 * \details
 *          This is one of three alternative modules: sky0.h / sky0.c,
 *          sky1.h / sky1.c and sky2.h / sky2.c. They contain routines for
 *          transforming astronomical positions from frame to another:
 *          precession, nutation, sidereal time etc. and they reflect changes
 *          in the International Astronomical Union's precession and nutation
 *          theory.
 *          The differences are:
 *          - sky0.h / sky0.c: nutation, obliquity  and sidereal time routines
 *            from the NREL Solar Position Algorithm document. These are based
 *            on the IAU 1980 nutation theory.
 *          - sky1.h / sky1.c: precession, nutation and their associated
 *            rotation matrices, obliquity and sidereal time. These are the
 *            IAU 1980 precession and nutation theory.
 *          - sky2.h / sky2.c: precession, nutation and their associated
 *            rotation matrices, obliquity and sidereal time. These are the
 *            newer IAU 2000 precession and nutation theory.
 * 
 *          This module (sky1.h / sky1.c) contains precession and nutation
 *          routines. The precession routine uses the IAU 1976 algorithm.
 *          The nutation routine uses the IAU 1980 algorithm.
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
 *      Character set: UTF-8. (Non-ASCII characters appear in this file)
 *----------------------------------------------------------------------------*/
#include "vectors3d.h"

/*
 * Global #defines and typedefs
 */
/*!     Precession angles.*/
typedef struct {
    double zeta_rad;        //!< Precession angle (ζ) (radian)
    double zed_rad;         //!< Precession angle (z) (radian)
    double theta_rad;       //!< Precession angle (θ) (radian)
} Sky1_Prec1976;

/*!      Nutation angles and obliquity */
typedef struct {
    double  dPsi_rad;       //!< Nutation in longitude (Δψ) (radian)
    double  dEps_rad;       //!< Nutation in obliquity (Δε) (radian)
    double  eps0_rad;       //!< Mean obliquity of ecliptic at date (ε0)(radian)
    double  eqEq_rad;       //!< Equation of the Equinoxes (radian)
} Sky1_Nut1980;


/*
 * Global functions available to be called by other modules
 */
#ifdef __cplusplus
extern "C" {
#endif

void sky1_frameBiasFK5(V3D_Matrix *biasM);

void sky1_precessionIAU1976(double t0, double t1, Sky1_Prec1976 *terms);
void sky1_createPrec1976Matrix(const Sky1_Prec1976 *terms, V3D_Matrix *precM);

void sky1_nutationIAU1980(double t_cy, int precision, Sky1_Nut1980 *nut);
void sky1_epsilon1980(double t_cy, Sky1_Nut1980 *nut);
void sky1_createNut1980Matrix(const Sky1_Nut1980 *nut, V3D_Matrix *nutM);

void sky1_createNPmatrix(double t0_cy,
                         double t1_cy,
                         int    precision,
                         V3D_Matrix* npM);

double sky1_gmSiderealTimeIAU1982(double du);
void sky1_appToTirs(const V3D_Vector *appV,
                    double           j2kUT1_d,
                    double           eqEq_rad,
                    V3D_Vector *terInterV);

/*
 * Global variables accessible by other modules
 */

#ifdef __cplusplus
}
#endif

#endif /* SKY1_H */
