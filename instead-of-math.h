#ifndef INSTEAD_OF_MATH_H
#define INSTEAD_OF_MATH_H
/*============================================================================*/
/*!\file
 * \brief instead-of-math.h - header to be included instead of math.h
 *
 * \author  David Hoadley <vcrumble@westnet.com.au>
 *
 * \details
 *          This header needs to be included instead of the standard C library
 *          header math.h, in order to provide some routines that are missing
 *          or possibly missing from math.h. This header will include math.h
 *          for you.
 * 
 *          The reason it must be included instead of math.h is a little
 *          obscure. One of the routines we want is the sincos() routine. Two
 *          versions of math.h actually do have it: The GNU C library does, and
 *          the Apple Clang library also has it, but under a different name.
 *
 *          The problem is, to gain access to the function, we need to include
 *          the appropriate incantation BEFORE including math.h if we are using
 *          the GNU C library, but we need to add definitions AFTER including
 *          math.h if we are using the Apple Clang library.
 *
 *          So use this header instead of math.h to look after this problem.
 *
 *          This header also provides the following missing function:
 *              - normalize()   
 *                      A proper modulo function. Unlike fmod(), result is
 *                      always positive. Used to to normalize a cyclic variable
 *                      to a range, e.g. an angle in radian to the range
 *                      [0, 2Pi), or an angle in degrees to the range [0, 360.0)
 *
 *
 * \copyright
 * \parblock
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
 * \endparblock
 *
 *==============================================================================
 */
#include "general.h"        /* Is this a C99 or later compiler? */

/*
 * Global #defines and typedefs
 */


/*
 * Global functions available to be called by other modules
 */
#ifdef __cplusplus
extern "C" {
#endif

/* -------------- The possible missing sincos() function -------------- */
/*! Calculates sine and cosine of an angle. Where the math library supports it,
    this is more efficient than calling \c sin() and \c cos() separately.
 \param[in]  angle_rad  Angle whose sine and cosine is wanted (radian)
 \param[out] sinA       sine of angle \a angle_rad
 \param[out] cosA       cosine of angle \a angle_rad
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined(__APPLE__)
    /* This is the CLANG compiler on an Apple computer (probably a mac). There
       is a sincos() function available, but it is called __sincos() instead.
       Re-direct all calls to sincos() to use __sincos() instead. */
#  include <math.h>             /* Must be done before defining fn below */
#  if defined (PREDEF_STANDARD_C_1999)
    /*          Compiler supports inline functions */
    static inline void sincos(double angle_rad, double *sinA, double *cosA)
    {
        __sincos(angle_rad, sinA, cosA);
    }
#  else
    /*          C89/C90 compiler - no inline functions. Need macros instead */
    #define sincos(angle_rad_, sinA_, cosA_) __sincos(angle_rad_, sinA_, cosA_)
#  endif

#elif defined(__GNUC__)
    /* This is the GNU C compiler. There is a sincos() function available in
       the GNU C library, but it is only accessible if the following macro is
       defined BEFORE including <math.h>. So do it here. */
#   define _GNU_SOURCE
#   include <math.h>

#else
    /* No sincos() function available. Oh well. Get sine and cosine separately.
       This is less efficient than a true sincos(), but at least it will work.*/
#  include <math.h>
#  if defined (PREDEF_STANDARD_C_1999)
    /*          Compiler supports inline functions */
    static inline void sincos(double angle_rad, double *sinA, double *cosA)
    {
        *sinA = sin(angle_rad);
        *cosA = cos(angle_rad);
    }
#  else
    /*          C89/C90 compiler - no inline functions. Need macros instead */
    #define sincos(angle_rad_, sinA_, cosA_) \
            { *sinA_ = sin(angle_rad_);      \
              *cosA_ = cos(angle_rad_); }
#  endif
#endif


/* -------------- The definitely missing normalize() function -------------- */

/*! Normalizes a cyclic double precision floating point variable \a x
    to the interval [0, range), assuming \a range is > 0. (If \a range is
    negative, this will normalize \a x to the interval (range, 0], but this is
    not the main use case for this function.)
 \returns The value of \a x, within the range [0, range)
 \param[in] x       The variable (e.g. an angle in degrees or radian)
 \param[in] range   The range to normalize \a x to (e.g. 2Pi or 360.0). It must
                    be non-zero.

    This function returns the same results as \c fmod() for positive \a x and
    \a range. Where it differs is in its handling of negative values of \a x.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined (PREDEF_STANDARD_C_1999)
/*          Compiler supports inline functions */
static inline double normalize(double x, double range)
{
    return x - floor(x / range) * range;
}
#else
 /*          C89/C90 compiler - no inline functions. Need macros instead */
 #define normalize(x__, range__)  ((x__) - floor((x__) / (range__)) * (range__))
#endif


/*
 * Global variables accessible by other modules
 */


#ifdef __cplusplus
}
#endif

#endif /* INSTEAD_OF_MATH_H */

