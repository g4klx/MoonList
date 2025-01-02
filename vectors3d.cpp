/*==============================================================================
 * vectors3d.c - Three dimensional geometry, vectors and matrices
 *
 * Author:  David Hoadley
 *
 * Description: (see vectors3d.h)
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
 *      Character set: UTF-8. (Non-ASCII characters appear in this file)
 *----------------------------------------------------------------------------*/

/* ANSI includes etc. */
#include "instead-of-math.h"                /* for sincos() */
#include <math.h>
#include <stdlib.h>                     /* for NULL */

/* Local and project includes */
#include "general.h"
#include "vectors3d.h"

/*
 * Local #defines and typedefs
 */
DEFINE_THIS_FILE;                       /* For use by REQUIRE() - assertions */

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
GLOBAL V3D_Vector *v3d_addV(V3D_Vector *destV,
                            const V3D_Vector *srcV1,
                            const V3D_Vector *srcV2)
/*! Add two vectors:
        [destV] = [srcV1] + [srcV2]
 \returns            pointer to \a destV
 \param[out] destV   vector which will contain the result
 \param[in]  srcV1   first vector
 \param[in]  srcV2   second vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(srcV1);
    REQUIRE_NOT_NULL(srcV2);
    REQUIRE_NOT_NULL(destV);

    destV->a[0] = srcV1->a[0] + srcV2->a[0];
    destV->a[1] = srcV1->a[1] + srcV2->a[1];
    destV->a[2] = srcV1->a[2] + srcV2->a[2];
    return destV;
}



GLOBAL V3D_Vector *v3d_addToV(V3D_Vector *modV1, const V3D_Vector *srcV2)
/*! Modify first vector by adding a second one:
        [modV1] += [srcV2]
 \returns               pointer to \a modV1
 \param[in,out] modV1   first vector, to which the second will be added
 \param[in]     srcV2   second vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(modV1);
    REQUIRE_NOT_NULL(srcV2);

    modV1->a[0] += srcV2->a[0];
    modV1->a[1] += srcV2->a[1];
    modV1->a[2] += srcV2->a[2];
    return modV1;
}



GLOBAL V3D_Vector *v3d_subtractV(V3D_Vector *destV,
                                 const V3D_Vector *srcV1,
                                 const V3D_Vector *srcV2)
/*! Vector subtraction:
        [destV] = [srcV1] - [srcV2]
 \returns            pointer to \a destV
 \param[out] destV   vector which will contain the result
 \param[in]  srcV1   first vector
 \param[in]  srcV2   second vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(srcV1);
    REQUIRE_NOT_NULL(srcV2);
    REQUIRE_NOT_NULL(destV);

    destV->a[0] = srcV1->a[0] - srcV2->a[0];
    destV->a[1] = srcV1->a[1] - srcV2->a[1];
    destV->a[2] = srcV1->a[2] - srcV2->a[2];
    return destV;
}



GLOBAL V3D_Vector *v3d_subFromV(V3D_Vector *modV1, const V3D_Vector *srcV2)
/*! Modify first vector by subtracting a second one:
        [modV1] -= [srcV2]
 \returns               pointer to \a modV1
 \param[in,out] modV1   first vector, from which the second will be subtracted
 \param[in]     srcV2   second vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(modV1);
    REQUIRE_NOT_NULL(srcV2);

    modV1->a[0] -= srcV2->a[0];
    modV1->a[1] -= srcV2->a[1];
    modV1->a[2] -= srcV2->a[2];
    return modV1;
}



GLOBAL double v3d_dotProductV(const V3D_Vector *srcV1, const V3D_Vector *srcV2)
/*! Return the dot product of the two vectors:
        res = [srcV1] · [srcV2]
    or  res = ||srcV1|| * ||srcV2|| * cos(θ)
 \returns            scalar result of the dot product
 \param[in]  srcV1   first vector
 \param[in]  srcV2   second vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(srcV1);
    REQUIRE_NOT_NULL(srcV2);

    return srcV1->a[0] * srcV2->a[0]
           + srcV1->a[1] * srcV2->a[1]
           + srcV1->a[2] * srcV2->a[2];
}



GLOBAL V3D_Vector *v3d_crossProductV(V3D_Vector *destV,
                                     const V3D_Vector *srcV1,
                                     const V3D_Vector *srcV2)
/*! Return the cross product of the two vectors:
        [destV] = [srcV1] x [srcV2]
 \returns            pointer to vector result of the cross product (\a destV)
 \param[out] destV   vector which will contain the result
 \param[in]  srcV1   first vector
 \param[in]  srcV2   second vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    V3D_Vector temp;            /* Allow caller to pass same vector for destV
                                 * and srcV1 or srcV2 */
    REQUIRE_NOT_NULL(srcV1);
    REQUIRE_NOT_NULL(srcV2);
    REQUIRE_NOT_NULL(destV);

    temp.a[0] = (srcV1->a[1] * srcV2->a[2]) - (srcV1->a[2] * srcV2->a[1]);
    temp.a[1] = (srcV1->a[2] * srcV2->a[0]) - (srcV1->a[0] * srcV2->a[2]);
    temp.a[2] = (srcV1->a[0] * srcV2->a[1]) - (srcV1->a[1] * srcV2->a[0]);
    *destV = temp;
    return destV;
}



GLOBAL V3D_Vector *v3d_addToUV(V3D_Vector *modV1, const V3D_Vector *srcV2)
/*! Modify a unit length rectangular position vector \a modV1 by adding a
    correction vector \a srcV2 to it and re-normalizing to unit length.
 \returns       Pointer to \a modV1
 \param[in,out] modV1   First vector, to which \a srcV2 will be added, and
                        then re-normalized to unit length.
 \param[in]     srcV2   vector to be added to \a modV1

 \par When to call this function
    If your vector \a modV1 is a unit vector, and you want it to be a unit
    vector after the addition of \a srcV2, then you can use this function. But
    if \a srcV2 is very small, you will be better off calling v3d_addToUVfast()
    instead.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double  mag;

    REQUIRE_NOT_NULL(modV1);
    REQUIRE_NOT_NULL(srcV2);

    v3d_addToV(modV1, srcV2);
    mag = v3d_magV(modV1);
    modV1->a[0] /= mag;
    modV1->a[1] /= mag;
    modV1->a[2] /= mag;
    return modV1;
}



GLOBAL V3D_Vector *v3d_addToUVfast(V3D_Vector *modV1, const V3D_Vector *srcV2)
/*! Modify a unit length rectangular position vector \a modV1 by adding a small
    correction vector \a srcV2 to it and re-normalizing to unit length. The
    scaling factor for the length is approximated by the expression
        1 - [modV1] · [srcV2],    where the dot indicates the Dot Product.
    The correction vector \a srcV2 must be very small.
 \returns       Pointer to \a modV1
 \param[in,out] modV1   First vector, to which \a srcV2 will be added, and
                        then re-normalized to unit length.
                        Valid range: |modV1| == 1
 \param[in]     srcV2   vector to be added to \a modV1.
                        Valid range: |srcV2| << 1

 \par When to call this function
    This routine performs the same function as v3d_addToUV(), but without
    calling \c sqrt(), or performing divisions. This can make it quicker than
    v3d_addToUV(), depending on what kind of floating point processor you have,
    if any. It is only valid if \a srcV2 is very small.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double  scale;

    REQUIRE_NOT_NULL(modV1);
    REQUIRE_NOT_NULL(srcV2);

    scale = 1.0 - v3d_dotProductV(modV1, srcV2);
    v3d_addToV(modV1, srcV2);
    modV1->a[0] *= scale;
    modV1->a[1] *= scale;
    modV1->a[2] *= scale;
    return modV1;
}



GLOBAL double v3d_magVSq(const V3D_Vector *srcV)
/*! Return the square of the magnitude of the specified vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(srcV);

    return srcV->a[0] * srcV->a[0]
           + srcV->a[1] * srcV->a[1]
           + srcV->a[2] * srcV->a[2];
}



GLOBAL double v3d_magV(const V3D_Vector *srcV)
/*! Return the magnitude of the specified vector
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    return sqrt(v3d_magVSq(srcV));
}



GLOBAL V3D_Vector *v3d_polarToRect(V3D_Vector *destV,
                                   double alpha_rad,
                                   double delta_rad)
/*! Converts polar (curvilinear) coordinates to equivalent rectangular
    (Cartesian) coordinates. The vector returned is of unit length and in a
    right or left-handed system according to the convention for measuring the
    angle alpha_rad: for example left-handed for Hour Angle/Declination or
    Azimuth/Elevation, and right-handed for Right Ascension/Declination.
 \returns               Pointer to \a destV, the resultant vector
 \param[out] destV      Resultant 3D unit vector in rectangular coordinates
 \param[in]  alpha_rad  Secondary angle coordinate, e.g. RA or Azimuth (radians)
 \param[in]  delta_rad  Primary angle coordinate, e.g. declination or elevation
                           (radians)
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double sinA, cosA;
    double sinD, cosD;

    REQUIRE_NOT_NULL(destV);

    sincos(alpha_rad, &sinA, &cosA);
    sincos(delta_rad, &sinD, &cosD);

    destV->a[0] = cosD * cosA;
    destV->a[1] = cosD * sinA;
    destV->a[2] = sinD;
    return destV;
}



GLOBAL void v3d_rectToPolar(double *alpha_rad,
                            double *delta_rad,
                            const V3D_Vector *srcV)
/*! Converts rectangular (Cartesian) coordinates to the equivalent polar
    (curvilinear) ones.
 \param[out] alpha_rad  Secondary angle coordinate (radians), range [-Pi, +Pi]
 \param[out] delta_rad  Primary angle coordinate (radians), range [-Pi/2, +Pi/2]
 \param[in]  srcV       3D unit vector in rectangular coordinates (direction
                           cosines)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    REQUIRE_NOT_NULL(alpha_rad);
    REQUIRE_NOT_NULL(delta_rad);
    REQUIRE_NOT_NULL(srcV);

    *alpha_rad = atan2(srcV->a[1], srcV->a[0]);
    *delta_rad = atan2(srcV->a[2],
                       sqrt(srcV->a[0] * srcV->a[0] + srcV->a[1] * srcV->a[1]));
}



#if 0
void v3d_polar2Rect(double alpha_rad, double delta_rad, V3D_Vector *R)
/*  Converts polar (curvilinear) coordinates to equivalent rectangular
    (cartesian) coordinates. The vector returned is of unit length and in a
    right or left-handed system according to the convention for measuring the
    angle alpha_rad: for example left-handed for HA/Dec or Az/El and
    right-handed for RA/Dec.
 Inputs
    alpha_rad - Secondary angle coordinate in radians
    delta_rad - Primary angle coordinate in radians
 Outputs
    R         - Resultant 3D unit vector in rectangular coordinates
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    R->a[0] = cos(delta_rad) * cos(alpha_rad);
    R->a[1] = cos(delta_rad) * sin(alpha_rad);
    R->a[2] = sin(delta_rad);
}



void v3d_rect2Polar(const V3D_Vector *R,
                    bool  alphaNonNegative,
                    double *alpha_rad,
                    double *delta_rad)
/*  Converts rectangular (cartesian) coordinates to the equivalent polar
    (curvilinear) ones.
 Inputs
    R           - 3D unit vector in rectangular coordinates (direction cosines)
    alphaNonNegative - if true, return alpha_rad in range [0, 2Pi)
 Outputs
    alpha_rad   - Secondary angle result in radians, range (-Pi, +Pi]
                    (or in range [0, 2Pi) if alphaNonNegative == TRUE)
    delta_rad   - Primary angle result in radians, range [-Pi/2, +Pi/2]
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    *alpha_rad = atan2(R->a[1], R->a[0]);
    if (alphaNonNegative && (*alpha_rad < 0.0)) { *alpha_rad += TWOPI; }
    *delta_rad = atan2(R->a[2], sqrt(R->a[0] * R->a[0] + R->a[1] * R->a[1]));
}
#endif



GLOBAL V3D_Matrix *v3d_createRotationMatrix(V3D_Matrix *destM,
                                            V3D_AxisNames axis,
                                            double        angle_rad)
/*! Creates a matrix to rotate a coordinate system about an axis. This matrix
    can then be used to convert a rectangular position vector from one
    coordinate system to another
 \returns               pointer to \a destM, the resultant matrix
 \param[out] destM      3x3 rotation matrix
 \param[in]  axis       axis of rotation.
                          Xaxis (=0), Yaxis (=1), Zaxis (=2) No other values
                          allowed.
 \param[in]  angle_rad  angle to rotate around the axis of rotation (radian).
                          When looking down the axis of rotation toward
                          the origin, positive angles of rotation are:
                          -  clockwise for left-handed systems (like
                             Azimuth/Elevation)
                          -  anticlockwise for right-handed systems (like
                             Right Ascension/Declination)

 \verbatim
                                                ┌  1     0     0  ┐
   createRotationMatrix(M, Xaxis, θ) is R1(θ) = │  0    cosθ  sinθ│
                                                └  0   -sinθ  cosθ┘

                                                ┌ cosθ   0   -sinθ┐
   createRotationMatrix(M, Yaxis, θ) is R2(θ) = │  0     1     0  │
                                                └ sinθ   0    cosθ┘

                                                ┌ cosθ  sinθ   0  ┐
   createRotationMatrix(M, Zaxis, θ) is R3(θ) = │-sinθ  cosθ   0  │
                                                └  0     0     1  ┘
 \endverbatim
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    double          cosA;
    double          sinA;
    unsigned int    i, j, k;

    REQUIRE((axis == Xaxis) || (axis == Yaxis) || (axis == Zaxis));
    REQUIRE_NOT_NULL(destM);

    i = axis;
    j = ((axis + 1) % 3);
    k = ((axis + 2) % 3);

    destM->a[i][i] = 1.0;
    destM->a[i][j] = 0.0;
    destM->a[i][k] = 0.0;
    destM->a[j][i] = 0.0;
    destM->a[k][i] = 0.0;

    sincos(angle_rad, &sinA, &cosA);

    destM->a[j][j] = cosA;
    destM->a[j][k] = sinA;
    destM->a[k][j] = -sinA;
    destM->a[k][k] = cosA;
    return destM;
}



GLOBAL V3D_Vector *v3d_multMxV(V3D_Vector *destV,
                               const V3D_Matrix *srcM,
                               const V3D_Vector *srcV)
/*! Multiply 3x3 matrix by 3x1 vector to give a new 3x1 vector, as per equation
        [destV] = [srcM] * [srcV]
 \returns            pointer to \a destV
 \param[out] destV   vector which will contain the result
 \param[in]  srcM    the matrix
 \param[in]  srcV    vector which will be multiplied by the matrix

 \note
    The same vector may be passed to \a destV and \a srcV, in which case the
    contents of \a srcV will be replaced.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int     i;
    /* Allow the caller to pass same vector for destV and srcV. Also allows the
     * compiler to optimise without us having to put "restrict" on the
     * arguments. */
    V3D_Vector temp;

    REQUIRE_NOT_NULL(srcV);
    REQUIRE_NOT_NULL(srcM);
    REQUIRE_NOT_NULL(destV);

    for (i = 0; i < 3; i++) {
        temp.a[i] =  srcM->a[i][0] * srcV->a[0]
                   + srcM->a[i][1] * srcV->a[1]
                   + srcM->a[i][2] * srcV->a[2];
    }
    *destV = temp;
    return destV;
}



GLOBAL V3D_Vector *v3d_multMtransxV(V3D_Vector *destV,
                                    const V3D_Matrix *srcM,
                                    const V3D_Vector *srcV)
/*! Multiply 3x1 vector by the transpose of the 3x3 matrix to give a new 3x1
   vector, as per equation
        [destV] = TRANSPOSE([srcM]) * [srcV]
 \returns            pointer to \a destV
 \param[out] destV   vector which will contain the result
 \param[in]  srcM    the matrix
 \param[in]  srcV    vector which will be multiplied by the matrix

 Note that for an orthogonal matrix which represents a rotation about
    axis j by angle θ, i.e. Matrix = Rj(θ), the following identity holds:
        TRANSPOSE(Matrix) = INVERSE(Matrix) = Rj(-θ).
    This is why this is a useful routine

 \note
    The same vector may be passed to \a destV and \a srcV, in which case the
    contents of \a srcV will be replaced.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int     i;
    /* Allow the caller to pass same vector for destV and srcV. Also allows the
     * compiler to optimise without us having to put "restrict" on the
     * arguments. */
    V3D_Vector temp;

    REQUIRE_NOT_NULL(srcV);
    REQUIRE_NOT_NULL(srcM);
    REQUIRE_NOT_NULL(destV);

    for (i = 0; i < 3; i++) {
        temp.a[i] =  srcM->a[0][i] * srcV->a[0]
                   + srcM->a[1][i] * srcV->a[1]
                   + srcM->a[2][i] * srcV->a[2];
    }
    *destV = temp;
    return destV;
}



GLOBAL V3D_Matrix *v3d_multMxM(V3D_Matrix *destM,
                               const V3D_Matrix *srcM1,
                               const V3D_Matrix *srcM2)
/*! Multiply 3x3 matrix by 3x3 matrix to give a new 3x3 matrix, as per
        [destM] = [srcM1]*[srcM2]
 \returns            pointer to \a destM
 \param[out] destM   matrix which will contain the result
 \param[in]  srcM1   first matrix
 \param[in]  srcM2   second matrix

 \note
    The same matrix may be passed to \a destM and either of the source matrices,
    in which case the contents of that source matrix will be replaced.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
{
    int         i, j;
    /* Allow the caller to pass same matrix for destM and srcM1 or srcM2. Also
     * allows the compiler to optimise without us having to put "restrict" on
     * the arguments. */
    V3D_Matrix  temp;

    REQUIRE_NOT_NULL(srcM1);
    REQUIRE_NOT_NULL(srcM2);
    REQUIRE_NOT_NULL(destM);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            temp.a[i][j] =  srcM1->a[i][0] * srcM2->a[0][j]
                          + srcM1->a[i][1] * srcM2->a[1][j]
                          + srcM1->a[i][2] * srcM2->a[2][j];
        }
    }
    *destM = temp;
    return destM;
}

/*! \page page-why-struct-array Why define vector and matrix arrays as structs?
    This has been done to assist compiler error checking, and reduce the chances
    of bugs that can potentially be rather difficult to diagnose.

    Take for example the following code. It has a vector of rectangular
    coordinates \c vec[], which contains the values one gets from converting
    polar angles of 45° and 45° to rectangular form. This program simply
    converts this rectangular vector back to polar form and writes out the two
    angles, after converting from radian back to degrees.

        #include <math.h>
        #include <stdio.h>

        void rectToPolar(double *alpha, double *delta, double rect[])
        {
            *alpha = atan2(rect[1], rect[0]);
            *delta = atan2(rect[2], sqrt(rect[0] * rect[0] + rect[1] * rect[1]));
        }

        int main(void)
        {
            double azimuth = 0.5;       // Dummy initial value, will be replaced
            double elevation = 0.5;     // Ditto
            double vec[3] = { 0.5, 0.5, 0.70710678119 };

            rectToPolar(vec, &azimuth, &elevation);
            printf("Azimuth = %f, Elevation = %f\n",
                   azimuth * 180.0 / 3.141592653589793,
                   elevation * 180.0 / 3.141592653589793);
        }

    There is a bug in the code, which in this short example is rather easy to
    see, but in a large program might not be. The expected output from this
    code is

        Azimuth = 45.000000, Elevation = 45.000000

    But what do you get if you compile and run it? On my running this program,
    I got

        Azimuth = 48.002776, Elevation = 28.647890

    You may get something else. But in any case, the results may look
    sufficiently plausible that you might not immediately notice that you have
    got a garbage answer. Depending on your program, this might end up being
    quite time-consuming to debug.

    Now instead, compile the following code, which uses our V3D_Vector struct
    instead of using a simple array:

        #include <math.h>
        #include <stdio.h>

        #include "vectors3d.h"

        int main(void)
        {
            double azimuth = 0.5;       // Dummy initial value, will be replaced
            double elevation = 0.5;     // Ditto
            V3D_Vector vec = {{ 0.5, 0.5, 0.70710678119 }};

            v3d_rectToPolar(&vec, &azimuth, &elevation);
            printf("Azimuth = %f, Elevation = %f\n",
                   azimuth * 180.0 / 3.141592653589793,
                   elevation * 180.0 / 3.141592653589793);
        }

    This time, you get two warning messages from the compiler. With the clang
    compiler, the messages are as follows:

        vec2.c:12:21: warning: incompatible pointer types passing 'V3D_Vector *'
              to parameter of type 'double *' [-Wincompatible-pointer-types]
            v3d_rectToPolar(&vec, &azimuth, &elevation);
                            ^~~~
        ./vectors3d.h:74:30: note: passing argument to parameter 'alpha_rad' here
        void v3d_rectToPolar(double *alpha_rad,
                                     ^
        vec2.c:12:37: warning: incompatible pointer types passing 'double *' to
              parameter of type 'const V3D_Vector *' [-Wincompatible-pointer-types]
            v3d_rectToPolar(&vec, &azimuth, &elevation);
                                            ^~~~~~~~~~
        ./vectors3d.h:76:40: note: passing argument to parameter 'srcV' here
                             const V3D_Vector *srcV);

    (Actually, there is a case here for elevating this warning to an error that
    stops compilation, using -Werror=incompatible-pointer-types)

    I believe that this justifies using this slightly more awkward way
    of implementing the vectors and arrays needed for 3-D geometry.

    It also allows one to use a simple assignment statement to copy a vector or
    matrix, as opposed to needing to write a function to do it. (Or in the case
    of the IAU SOFA routines for matrix copying - two functions!)
 */
