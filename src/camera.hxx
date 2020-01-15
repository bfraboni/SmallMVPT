/*
 * Copyright (C) 2019, Basile Fraboni
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

#ifndef __CAMERA_HXX__
#define __CAMERA_HXX__

#include <vector>
#include <cmath>
#include "math.hxx"
#include "ray.hxx"

class Camera
{
public:

    void Setup(
        const Vec3f &aPosition,
        const Vec3f &aForward,
        const Vec3f &aUp,
        const Vec2f &aResolution,
        float       aHorizontalFOV)
    {
        const Vec3f forward = Normalize(aForward);
        const Vec3f up      = Normalize(Cross(aUp, -forward));
        const Vec3f left    = Cross(-forward, up);

        mPosition   = aPosition;
        mForward    = forward;
        mResolution = aResolution;

        const Vec3f pos(
            Dot(up, aPosition),
            Dot(left, aPosition),
            Dot(-forward, aPosition));

        Mat4f worldToCamera = Mat4f::Identity();
        worldToCamera.SetRow(0, up,       -pos.x);
        worldToCamera.SetRow(1, left,     -pos.y);
        worldToCamera.SetRow(2, -forward, -pos.z);

        const Mat4f perspective = Mat4f::Perspective(aHorizontalFOV, 0.1f, 10000.f);
        const Mat4f worldToNScreen = perspective * worldToCamera;
        const Mat4f nscreenToWorld = Invert(worldToNScreen);

        mWorldToRaster  =
            Mat4f::Scale(Vec3f(aResolution.x * 0.5f, aResolution.y * 0.5f, 0)) *
            Mat4f::Translate(Vec3f(1.f, 1.f, 0)) * worldToNScreen;

        mRasterToWorld  = nscreenToWorld *
            Mat4f::Translate(Vec3f(-1.f, -1.f, 0)) *
            Mat4f::Scale(Vec3f(2.f / aResolution.x, 2.f / aResolution.y, 0));

        const float tanHalfAngle = std::tan(aHorizontalFOV * PI_F / 360.f);
        mImagePlaneDist = aResolution.x / (2.f * tanHalfAngle);

        // we need to find raster bounds at distance mImagePlaneDist
        // compute raster bounds
        Vec3f p0 = mRasterToWorld.TransformPoint(Vec3f(0, 0, 1));
        Vec3f px = mRasterToWorld.TransformPoint(Vec3f(aResolution.x, 0, 1));
        Vec3f py = mRasterToWorld.TransformPoint(Vec3f(0, aResolution.y, 1));

        // compute directions between camera origin and raster bounds
        Vec3f v0 = Normalize(mPosition - p0);
        Vec3f vx = Normalize(mPosition - px);
        Vec3f vy = Normalize(mPosition - py);
        
        // compute raster bounds in camera space at distance mImagePlaneDist
        p0 = Vec3f(v0 * mImagePlaneDist / Dot(v0, mForward));
        px = Vec3f(vx * mImagePlaneDist / Dot(vx, mForward));
        py = Vec3f(vy * mImagePlaneDist / Dot(vy, mForward));

        // compute film area
        mImagePlaneArea = Cross(Vec3f(px - p0), Vec3f(py - p0)).Length();
    }

    int RasterToIndex(const Vec2f &aPixelCoords) const
    {
        return int(std::floor(aPixelCoords.x) + std::floor(aPixelCoords.y) * mResolution.x);
    }

    Vec2f IndexToRaster(const int &aPixelIndex) const
    {
        const float y = std::floor(aPixelIndex / mResolution.x);
        const float x = float(aPixelIndex) - y * mResolution.x;
        return Vec2f(x, y);
    }

    Vec3f RasterToWorld(const Vec2f &aRasterXY) const
    {
        return mRasterToWorld.TransformPoint(Vec3f(aRasterXY.x, aRasterXY.y, 0));
    }

    Vec2f WorldToRaster(const Vec3f &aWorldPos) const
    {
        Vec3f temp = mWorldToRaster.TransformPoint(aWorldPos);
        return Vec2f(temp.x, temp.y);
    }

    // returns false when raster position is outside screen space
    bool CheckRaster(const Vec2f &aRasterPos) const
    {
        return aRasterPos.x >= 0 && aRasterPos.y >= 0 &&
            aRasterPos.x < mResolution.x && aRasterPos.y < mResolution.y;
    }

    Ray GenerateRay(const Vec2f &aRasterXY) const
    {
        const Vec3f worldRaster = RasterToWorld(aRasterXY);

        Ray res;
        res.org  = mPosition;
        res.dir  = Normalize(worldRaster - mPosition);
        res.tmin = 0;
        return res;
    }

    ///
    /// @ brief Returns the pdf of a ray
    ///
    float Pdf(const Ray& aRay) const
    {
        float cosTheta = Dot( aRay.dir, mForward );
        float lenSqr = mImagePlaneDist * mImagePlaneDist;
        return lenSqr / (cosTheta * cosTheta * cosTheta * mImagePlaneArea);
    }

    ///
    /// @ brief Returns the derivatives of the film point w.r.t the lens point
    ///
    float DfilmDlens(const Ray& aRay) const
    {
        float cosTheta = Dot( aRay.dir, mForward );
        float lenSqr = mImagePlaneDist * mImagePlaneDist;
        return lenSqr / (cosTheta * cosTheta * cosTheta);
    }

    ///
    /// @ brief Returns a ray that connects the camera origin to a world position
    ///
    Ray ConnectionRay(const Vec3f &worldPosition) const
    {
        Ray res;
        res.org  = mPosition;
        res.dir  = Normalize(worldPosition - mPosition);
        res.tmin = (worldPosition - mPosition).Length();
        return res;
    }

public:

    Vec3f mPosition;
    Vec3f mForward;
    Vec2f mResolution;
    Mat4f mRasterToWorld;
    Mat4f mWorldToRaster;
    float mImagePlaneDist;
    float mImagePlaneArea;
};

#endif //__CAMERA_HXX__