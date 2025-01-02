#ifndef STAR_H
#define STAR_H
/*============================================================================*/
/*!\file
 * \brief star.h - Astronomical conversion routines for stars and other objects
 *                 beyond the Solar System
 *
 * \author  David Hoadley
 *
 * \details
 *          Routines to convert a celestial object's catalogue position (as
 *          seen from the Solar System Barycentre) to a geocentric apparent
 *          position (as seen from the centre of the Earth). This means
 *          - find the position and velocity of the Earth with respect to the
 *            Solar System Barycentre.
 *          - apply proper motion from the catalogue epoch to the desired epoch
 *          - offset the position by the position of the Earth (apply Annual
 *            Parallax)
 *          - use the velocity of the Earth to calculate the Annual Aberration
 *            correction
 * 
 *          The geocentric position calculated is then converted from the
 *          coordinate system of the catalogue position to apparent coordinates.
 *          (Exception: if the catalogue position is already in Celestial
 *          Intermediate coordinates, then the resulting vector will still be
 *          in Celestial Intermediate coordinates.)
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

#include "sky1.h"
#include "general.h"
#include "sky.h"
#include "vectors3d.h"

/*
 * Global #defines and typedefs
 */
/*!     Coordinate system equator and origin */
typedef enum {
    APPARENT,       //!< True equator and equinox of date
    INTERMEDIATE,   //!< True equator and CIO of date
    FK4,            //!< Mean equator and equinox of Besselian epoch
    FK5,            //!< Mean equator and equinox of Julian epoch
    ICRS            //!< International Celestial Reference System
} Star_CoordSys;


/*! Catalogue position of a celestial object. */
typedef struct {
    char   objectName[80];//!< Name of the celestial object (optional)
    double ra_rad;        //!< α - Right Ascension (radian)
    double dec_rad;       //!< δ - Declination (radian)
    Star_CoordSys cSys;   //!< equator and origin used for α & δ
    double eqnxT_cy;      /*!< time of equinox for FK4 or FK5 positions (Julian
                           *   centuries since J2000, TT timescale) */
    double epochT_cy;     /*!< time zero for proper motion (Julian centuries 
                           *   since J2000, TT timescale) */
    bool   epochSpecified;//!< equinox and epoch were separately specified
    double muRA_radpcy;   //!< μα - Proper motion in RA (radian/Julian century)
    double muDec_radpcy;  //!< μδ - Proper motion in Dec (radian/Julian century)
    double parallax_rad;  //!< π - annual parallax (radian)
    double radVel_aupcy;  //!< ν - radial velocity (AU/Julian century)
    bool   muRaInclCosDec;//!< μα term already includes the cosδ factor
} Star_CatalogPosn;

/*! Errors detected when decoding a text coordinate string or when filling a
    block of type Star_CatalogPosn with numerical data */
typedef enum {
    STAR_NORMAL,        /*!< Normal successful completion */
    STAR_IVEPOCH,       /*!< Invalid epoch or equinox specification */
    STAR_ERRINOBJ,      /*!< Error extracting object name from the string.
                         *   (comma missing at end?, or name started with a 
                         *    double quote char but there was no closing
                         *    double-quote?) */
    STAR_NORA,          /*!< Specification of Right Ascension is missing */
    STAR_ERRINRA,       /*!< Right Ascension can't be decoded */
    STAR_RARANERR,      /*!< Right Ascension out of range */
    STAR_NODEC,         /*!< Specification of Declination is missing */
    STAR_ERRINDEC,      /*!< Declination can't be decoded */
    STAR_DECRANERR,     /*!< Declination out of range */
    STAR_ERRINMURA,     /*!< Proper motion in RA can't be decoded */
    STAR_ERRINMUDEC,    /*!< Proper motion in Declination can't be decoded */
    STAR_NOMUDEC,       /*!< Specification of proper motion in Dec is missing
                         *    when mu_ra was provided. Check commas in string */
    STAR_NOMURA,        /*!< Specification of proper motion in RA is missing
                         *   when mu_dec was provided. Check commas in string */
    STAR_ERRINPARAL,    /*!< Parallax can't be decoded */
    STAR_PINEEDSPM,     /*!< Proper motion must be specified if specifying
                         *   parallax */
    STAR_ERRINVR,       /*!< Radial Velocity can't be decoded */
    STAR_VRNEEDSPM      /*!< Proper motion and parallax must be specified if
                         *   specifying Radial Velocity */
} Star_CoordErrors;

/*
 * Global functions available to be called by other modules
 */
void star_setCurrentObject(const Star_CatalogPosn *c);

/*      Functions to get the current position of the star */
void star_catalogToApp(const Star_CatalogPosn *c,
                       double                 j2kTT_cy,
                       const Sky1_Nut1980    *nut,
                       V3D_Vector *appV,
                       double     *dist_au);
void star_getApparent(double j2kTT_cy, Sky_TrueEquatorial *pos);
void star_getTopocentric(double             j2kUtc_d,
                         const Sky_DeltaTs  *deltas,
                         const Sky_SiteProp *site,
                         Sky_SiteHorizon *topo);

/*      Alternative functions to set up a coordinate block */
int star_setCatalogPosn(const char    objectName[],
                        double        ra_h,
                        double        dec_deg,
                        Star_CoordSys coordSys,
                        double        equinoxYear,
                        double        epochYear,
                        double        muRaSplat_maspa,
                        double        muDec_maspa,
                        double        annParallax_as,
                        double        radVel_kmps,
                        Star_CatalogPosn *coordBlock);
int star_setCatalogOldStyle(const char    objectName[],
                            double        ra_h,
                            double        dec_deg,
                            Star_CoordSys coordSys,
                            double        equinoxYear,
                            double        epochYear,
                            double        muRa_spcy,
                            double        muDec_aspcy,
                            double        annParallax_as,
                            double        radVel_kmps,
                            Star_CatalogPosn *coordBlock);
int star_parseCoordString(const char inStr[],
                          Star_CatalogPosn *coordBlock);

/*      Functions called by the above functions */
void star_catalogToVectors(const Star_CatalogPosn *c,
                            V3D_Vector *pV,
                            V3D_Vector *vV_radpcy);
void star_earth(double           t1_cy,
                const V3D_Matrix *npM,
                V3D_Vector *ebV_au,
                V3D_Vector *ebdotV_aupd,
                V3D_Vector *sbV_au);
void star_annAberr(const V3D_Vector *p1V,
                   const V3D_Vector *earthVelV_aupd,
                   V3D_Vector *p2V);

/*      Writing out equinox and epoch in standard form */
char *star_equinoxToStr(const Star_CatalogPosn *coordBlock,
                        char   equinoxStr[],
                        size_t eqnxStrSize);
char *star_epochToStr(const Star_CatalogPosn *coordBlock,
                      char   epochStr[],
                      size_t epochStrSize);

#if 0
void star_lightDeflection(const V3D_Vector *pV,
                          const V3D_Vector *earthBV_au,
                          const V3D_Vector *earthVelV_aupd,
                          const V3D_Vector *sunBV_au,
                          V3D_Vector *p2V);
#endif

/*
 * Global variables accessible by other modules
 */

#endif /*STAR_H*/
