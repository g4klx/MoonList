#ifndef SKYIO_H
#define SKYIO_H
/*============================================================================*/
/*!\file
 * \brief skyio.h - output and formatting routines and a read routine
 *
 * \author  David Hoadley
 *
 * \details
 *          Routines for reading and writing out assorted astronomical things,
 *          such as:
 *          - formatting a radian angle, or an angle expressed in hours, into
 *            a string in Hours Minutes Seconds form, something only an
 *            astronomer would think of doing.
 *          - formatting a radian angle, or an angle expressed in degrees into
 *            a string in Degrees Minutes Seconds form, and doing it correctly.
 *            Quite a few people need this one.
 *          - reading an angle from a string
 *          - write out a Julian Date (in J2KD form - see 
 *            \ref page-time-variables) as a calendar date and time.
 * 
 *          These routines were mainly developed for debugging.
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

#include "general.h"
#include "astron.h"

/*
 * Global #defines and typedefs
 */
#define NO_ANGLE                (-1)
#define INVALID_ANGLE           (-2)


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Global functions available to be called by other modules
 */
/*      Convert from an angle (or a time) into a sexagesimal text string */
char *skyio_degToDmsStr(char destStr[],
                        size_t   destStrSize,
                        double   angle_deg,
                        unsigned decimals);

char *skyio_hrsToHmsStr(char destStr[],
                        size_t   destStrSize,
                        double   angle_h,
                        unsigned decimals);

/*      Convert from a sexagesimal text string to an angle */
double skyio_sxStrToAng(const char angleStr[],
                        const char **endPtr,
                        int *error);



/*! Routine to take an angle in radian and return a string in degrees, 
    arcminutes and arcseconds form - [±]DDD°MM′SS.sss″ - 
    correctly rounding according to the number of decimal places to be shown.
    The angle is assumed to be within the range (-2*Pi, 2*Pi). Angles which
    round to ±360° will be written out as 0°.
 \returns                Pointer to \a destStr
 \param[out] destStr     Destination character string
 \param[in]  destStrSize Size of destination string (max available length + 1)
 \param[in]  angle_rad   The angle to be written out (radian).
                           Valid range:(-2*Pi, 2*Pi); larger numbers may well be
                           written OK, but there is a risk of overflowing an
                           intermediate variable. No error will be detected,
                           just a wrong answer will be written.
 \param[in]  decimals    Number of digits after the decimal point to display.
                           Valid range: [0,9] (or [0,3] if long int is only a
                           32-bit number); numbers outside this range will be
                           clamped to this range.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined (PREDEF_STANDARD_C_1999)
/*          Compiler supports inline functions */
static inline char *skyio_radToDmsStr(char destStr[],
                                      size_t   destStrSize,
                                      double   angle_rad,
                                      unsigned decimals)
{
   return skyio_degToDmsStr(destStr, destStrSize, radToDeg(angle_rad),decimals);
}
#else
 /*          C89/C90 compiler - no inline functions. Need macros instead */
 #define skyio_radToDmsStr(destStr__, destStrSize__, angle_rad__, decimals__) \
    skyio_degToDmsStr(destStr__, destStrSize__,                               \
                      radToDeg(angle_rad__), decimals__)
#endif


/*! Routine to take an angle in radian and return a string in hours, minutes and
    seconds form - "±HH:MM:SS.sss" - 
    correctly rounding according to the number of decimal places to be shown.
    The angle is assumed to be within the range (-2*Pi, 2*Pi). Angles which
    round to ±24:00:00 will be written out as 0:00:00.
 \returns                Pointer to \a destStr
 \param[out] destStr     Destination character string
 \param[in]  destStrSize Size of destination string (max available length + 1)
 \param[in]  angle_rad   The angle to be written out (radian).
                         Valid range:(-2*Pi, 2*Pi); larger numbers may well be
                         written OK, but there is a risk of overflowing an
                         intermediate variable. No error will be detected,
                         just a wrong answer will be written.
 \param[in]  decimals    Number of digits after the decimal point to display.
                         Valid range: [0,9] (or [0,3] if long int is only a
                         32-bit number); numbers outside this range will be
                         clamped to this range.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined (PREDEF_STANDARD_C_1999)
/*          Compiler supports inline functions */
static inline char *skyio_radToHmsStr(char destStr[],
                                      size_t   destStrSize,
                                      double   angle_rad,
                                      unsigned decimals)
{
   return skyio_hrsToHmsStr(destStr, destStrSize, radToHrs(angle_rad),decimals);
}
#else
 /*          C89/C90 compiler - no inline functions. Need macros instead */
 #define skyio_radToHmsStr(destStr__, destStrSize__, angle_rad__, decimals__) \
    skyio_hrsToHmsStr(destStr__, destStrSize__,                               \
                      radToHrs(angle_rad__), decimals__)
#endif


void skyio_printJ2kd(double j2kd);

/*
 * Global variables accessible by other modules
 */


#ifdef __cplusplus
}
#endif

#endif /* SKYIO_H */

