#ifndef ASTRON_H
#define ASTRON_H
/*============================================================================*/
/*!\file
 * \brief astron.h - assorted definitions useful for astronomy
 *
 * \author  David Hoadley
 *
 * \details
 *          Some astronomical constants, and angle conversions that are used in
 *          astronomical software.
 * 
 *==============================================================================
 */
#include "general.h"        /* For PI, degToRad(), radToDeg() etc., C standard*/
#include "vectors3d.h"


/*
 * Global #defines and typedefs
 */
/* Useful astronomical constants */
#define MJD_B1900   15019.81352     /*!< MJD of Besselian std epoch B1900.0 */
#define MJD_B1950   33281.92346     /*!< MJD of Besselian std epoch B1950.0 */
#define MJD_J2000   51544.5         /*!< MJD of Fundamental Epoch J2000.0 */
#define TROP_CENT   36524.2198781   /*!< Length of Tropical Century in days */
#define JUL_CENT    36525.0         /*!< Length of Julian Century in days */

#define ARCSEC2RAD  (PI / 648000.0)     /*!< arcseconds to radians */
#define RAD2ARCSEC  (648000.0 / PI)     /*!< radians to arcseconds */
#define SEC2RAD     (PI / 43200.0)      /*!< seconds(time) to radians */
#define RAD2SEC     (43200.0 / PI)      /*!< radians to seconds(time) */
#define HRS2RAD     (PI / 12.0)         /*!< hours to radians */
#define RAD2HRS     (12.0 / PI)         /*!< radians to hours */

/*
 * Functions
 */
/*      Convert angle from one unit to another */
#ifdef PREDEF_STANDARD_C_1999
/*          Compiler supports inline functions */
/*! Returns \a angle_as converted from arcseconds to radians */
static inline double arcsecToRad(double angle_as) {
                                                return angle_as * ARCSEC2RAD; }
/*! Returns \a angle_rad converted from radians to arcseconds */
static inline double radToArcsec(double angle_rad) {
                                                return angle_rad * RAD2ARCSEC; }

/*! Returns \a angle_s converted from seconds to radians */
static inline double secToRad(double angle_s)   { return angle_s * SEC2RAD; }
/*! Returns \a angle_rad converted from radians to seconds */
static inline double radToSec(double angle_rad) { return angle_rad * RAD2SEC; }

/*! Returns \a angle_h converted from hours to radians */
static inline double hrsToRad(double angle_h)   { return angle_h * HRS2RAD; }
/*! Returns \a angle_rad converted from radians to hours */
static inline double radToHrs(double angle_rad) { return angle_rad * RAD2HRS; }

#else
/*          C89/C90 compiler only - no inline functions. Need macros instead */
#define arcsecToRad(angle_arcsec__) ((angle_arcsec__) * ARCSEC2RAD)
#define radToArcsec(angle_rad__)    ((angle_rad__) * RAD2ARCSEC)

#define secToRad(angle_s__)         ((angle_s__) * SEC2RAD)
#define radToSec(angle_rad__)       ((angle_rad__) * RAD2SEC)

#define hrsToRad(angle_h__)         ((angle_h__) * HRS2RAD)
#define radToHrs(angle_rad__)       ((angle_rad__) * RAD2HRS)

#endif


#endif /*ASTRON_H*/

