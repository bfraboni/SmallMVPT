/*
 * Copyright (C) 2012, Tomas Davidovic (http://www.davidovic.cz)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * (The above is MIT License: http://en.wikipedia.org/wiki/MIT_License)
 */

#ifndef __UTILS_HXX__
#define __UTILS_HXX__

#include <vector>
#include <cmath>
#include "math.hxx"

#define EPS_COSINE 1e-6f
#define EPS_RAY    1e-3f

// sRGB luminance
float Luminance(const Vec3f& aRGB)
{
    return 0.212671f * aRGB.x +
        0.715160f * aRGB.y +
        0.072169f * aRGB.z;
}

float FresnelDielectric(
    float aCosInc,
    float mIOR)
{
    if(mIOR < 0)
        return 1.f;

    float etaIncOverEtaTrans;

    if(aCosInc < 0.f)
    {
        aCosInc = -aCosInc;
        etaIncOverEtaTrans = mIOR;
    }
    else
    {
        etaIncOverEtaTrans = 1.f / mIOR;
    }

    const float sinTrans2 = Sqr(etaIncOverEtaTrans) * (1.f - Sqr(aCosInc));
    const float cosTrans = std::sqrt(std::max(0.f, 1.f - sinTrans2));

    const float term1 = etaIncOverEtaTrans * cosTrans;
    const float rParallel =
        (aCosInc - term1) / (aCosInc + term1);

    const float term2 = etaIncOverEtaTrans * aCosInc;
    const float rPerpendicular =
        (term2 - cosTrans) / (term2 + cosTrans);

    return 0.5f * (Sqr(rParallel) + Sqr(rPerpendicular));
}

// reflect vector through (0,0,1)
Vec3f ReflectLocal(const Vec3f& aVector)
{
    return Vec3f(-aVector.x, -aVector.y, aVector.z);
}

Vec3f Reflect(const Vec3f& aVector, const Vec3f& aNormal = Vec3f(0,0,1))
{
    return Normalize(-aVector + 2.f * Dot(aVector, aNormal) * aNormal);
}

//////////////////////////////////////////////////////////////////////////
// Cosine lobe hemisphere sampling

Vec3f SamplePowerCosHemisphereW(
    const Vec2f  &aSamples,
    const float  aPower,
    float        *oPdfW)
{
    const float term1 = 2.f * PI_F * aSamples.x;
    const float term2 = std::pow(aSamples.y, 1.f / (aPower + 1.f));
    const float term3 = std::sqrt(1.f - term2 * term2);

    if(oPdfW)
    {
        *oPdfW = (aPower + 1.f) * std::pow(term2, aPower) * (0.5f * INV_PI_F);
    }

    return Vec3f(
        std::cos(term1) * term3,
        std::sin(term1) * term3,
        term2);
}

float PowerCosHemispherePdfW(
    const Vec3f  &aNormal,
    const Vec3f  &aDirection,
    const float  aPower)
{
    const float cosTheta = std::max(0.f, Dot(aNormal, aDirection));

    return (aPower + 1.f) * std::pow(cosTheta, aPower) * (INV_PI_F * 0.5f);
}


//////////////////////////////////////////////////////////////////////////
// Disc sampling

Vec2f SampleConcentricDisc(
    const Vec2f &aSamples)
{
    float phi, r;

    float a = 2*aSamples.x - 1;   /* (a,b) is now on [-1,1]^2 */
    float b = 2*aSamples.y - 1;

    if(a > -b)      /* region 1 or 2 */
    {
        if(a > b)   /* region 1, also |a| > |b| */
        {
            r = a;
            phi = (PI_F/4.f) * (b/a);
        }
        else        /* region 2, also |b| > |a| */
        {
            r = b;
            phi = (PI_F/4.f) * (2.f - (a/b));
        }
    }
    else            /* region 3 or 4 */
    {
        if(a < b)   /* region 3, also |a| >= |b|, a != 0 */
        {
            r = -a;
            phi = (PI_F/4.f) * (4.f + (b/a));
        }
        else        /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
        {
            r = -b;

            if (b != 0)
                phi = (PI_F/4.f) * (6.f - (a/b));
            else
                phi = 0;
        }
    }

    Vec2f res;
    res.x = r * std::cos(phi);
    res.y = r * std::sin(phi);
    return res;
}

float ConcentricDiscPdfA()
{
    return INV_PI_F;
}


//////////////////////////////////////////////////////////////////////////
/// Sample direction in the upper hemisphere with cosine-proportional pdf
/** The returned PDF is with respect to solid angle measure */
Vec3f SampleCosHemisphereW(
    const Vec2f  &aSamples,
    float        *oPdfW)
{
    const float term1 = 2.f * PI_F * aSamples.x;
    const float term2 = std::sqrt(1.f - aSamples.y);

    const Vec3f ret(
        std::cos(term1) * term2,
        std::sin(term1) * term2,
        std::sqrt(aSamples.y));

    if(oPdfW)
    {
        *oPdfW = ret.z * INV_PI_F;
    }

    return ret;
}

float CosHemispherePdfW(
    const Vec3f  &aNormal,
    const Vec3f  &aDirection)
{
    return std::max(0.f, Dot(aNormal, aDirection)) * INV_PI_F;
}

// Sample Triangle
// returns barycentric coordinates
Vec2f SampleUniformTriangle(const Vec2f &aSamples)
{
    const float term = std::sqrt(aSamples.x);

    return Vec2f(1.f - term, aSamples.y * term);
}

//////////////////////////////////////////////////////////////////////////
// Sphere sampling

Vec3f SampleUniformSphereW(
    const Vec2f  &aSamples,
    float        *oPdfSA)
{
    const float term1 = 2.f * PI_F * aSamples.x;
    const float term2 = 2.f * std::sqrt(aSamples.y - aSamples.y * aSamples.y);

    const Vec3f ret(
        std::cos(term1) * term2,
        std::sin(term1) * term2,
        1.f - 2.f * aSamples.y);

    if(oPdfSA)
    {
        //*oPdfSA = 1.f / (4.f * PI_F);
        *oPdfSA = INV_PI_F * 0.25f;
    }

    return ret;
}

float UniformSpherePdfW()
{
    //return (1.f / (4.f * PI_F));
    return INV_PI_F * 0.25f;
}


//////////////////////////////////////////////////////////////////////////
// Utilities for converting PDF between Area (A) and Solid angle (W)
// WtoA = PdfW * cosine / distance_squared
// AtoW = PdfA * distance_squared / cosine

float PdfWtoA(
    const float aPdfW,
    const float aDist,
    const float aCosThere)
{
    return aPdfW * std::abs(aCosThere) / Sqr(aDist);
}

float PdfAtoW(
    const float aPdfA,
    const float aDist,
    const float aCosThere)
{
    return aPdfA * Sqr(aDist) / std::abs(aCosThere);
}

#endif //__UTILS_HXX__
