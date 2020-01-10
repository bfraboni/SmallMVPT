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

#ifndef __SCENE_HXX__
#define __SCENE_HXX__

#include <vector>
#include <map>
#include <cmath>
#include "math.hxx"
#include "geometry.hxx"
#include "camera.hxx"
#include "materials.hxx"
#include "lights.hxx"

class Scene
{
public:
    Scene() :
        mGeometry(NULL),
        mBackground(NULL)
    {}

    ~Scene()
    {
        delete mGeometry;

        for(size_t i=0; i<mLights.size(); i++)
            delete mLights[i];
    }

    bool Intersect(
        const Ray &aRay,
        Isect     &oResult) const
    {
        bool hit = mGeometry->Intersect(aRay, oResult);

        if(hit)
        {
            oResult.lightID = -1;
            std::map<int, int>::const_iterator it =
                mMaterial2Light.find(oResult.matID);

            if(it != mMaterial2Light.end())
                oResult.lightID = it->second;
        }

        return hit;
    }

    bool Occluded(
        const Vec3f &aPoint,
        const Vec3f &aDir,
        float aTMax) const
    {
        Ray ray;
        ray.org  = aPoint + aDir * EPS_RAY;
        ray.dir  = aDir;
        ray.tmin = 0;
        Isect isect;
        isect.dist = aTMax - 2*EPS_RAY;

        return mGeometry->IntersectP(ray, isect);
    }

    const Material& GetMaterial(const int aMaterialIdx) const
    {
        return mMaterials[aMaterialIdx];
    }

    int GetMaterialCount() const
    {
        return (int)mMaterials.size();
    }


    const AbstractLight* GetLightPtr(int aLightIdx) const
    {
        aLightIdx = std::min<int>(aLightIdx, mLights.size()-1);
        return mLights[aLightIdx];
    }

    int GetLightCount() const
    {
        return (int)mLights.size();
    }

    const BackgroundLight* GetBackground() const
    {
        return mBackground;
    }

    //////////////////////////////////////////////////////////////////////////
    // Loads a Cornell Box scene
    enum BoxMask
    {
        kLightCeiling      = 1,
        kLightSun          = 2,
        kLightPoint        = 4,
        kLightBackground   = 8,
        kLargeMirrorSphere = 16,
        kLargeGlassSphere  = 32,
        kSmallMirrorSphere = 64,
        kSmallGlassSphere  = 128,
        kGlossyFloor       = 256,
        kLargeGlossySphere = 512,
        kBothSmallSpheres  = (kSmallMirrorSphere | kSmallGlassSphere),
        kBothLargeSpheres  = (kLargeMirrorSphere | kLargeGlassSphere),
        kDefault           = (kLightCeiling | kBothSmallSpheres),
    };

    void LoadCornellBox(
        const Vec2i &aResolution,
        uint aBoxMask = kDefault, 
        float aPhongExp = 1)
    {
        mSceneName = GetSceneName(aBoxMask, &mSceneAcronym);

        if((aBoxMask & kBothLargeSpheres) == kBothLargeSpheres)
        {
            printf("Cannot have both large balls, using mirror\n\n");
            aBoxMask &= ~kLargeGlassSphere;
        }

        bool light_ceiling    = (aBoxMask & kLightCeiling)    != 0;
        bool light_sun        = (aBoxMask & kLightSun)        != 0;
        bool light_point      = (aBoxMask & kLightPoint)      != 0;
        bool light_background = (aBoxMask & kLightBackground) != 0;

        bool light_box = true;

        bool full_glossy = (aBoxMask & kGlossyFloor) && (aBoxMask & kLargeGlossySphere);

        // because it looks really weird with it
        if(light_point)
            light_box = false;

        // Camera
        mCamera.Setup(
            Vec3f(-0.0439815f, -4.12529f, 0.222539f),
            Vec3f(0.00688625f, 0.998505f, -0.0542161f),
            Vec3f(3.73896e-4f, 0.0542148f, 0.998529f),
            Vec2f(float(aResolution.x), float(aResolution.y)), 45);

        // Material metallic specular tint
        // https://seblagarde.wordpress.com/2011/08/17/feeding-a-physical-based-lighting-mode/
        const Vec3f Silver    = Vec3f(0.971519,0.959915,0.915324);
        const Vec3f Aluminium = Vec3f(0.913183,0.921494,0.924524);
        const Vec3f Gold      = Vec3f(1.000000,0.765557,0.336057);
        const Vec3f Copper    = Vec3f(0.955008,0.637427,0.538163);
        const Vec3f Chromium  = Vec3f(0.549585,0.556114,0.554256);
        const Vec3f Nickel    = Vec3f(0.659777,0.608679,0.525649);
        const Vec3f Titanium  = Vec3f(0.541931,0.496791,0.449419);
        const Vec3f Cobalt    = Vec3f(0.662124,0.654864,0.633732);
        const Vec3f Platinum  = Vec3f(0.672411,0.637331,0.585456);

        // Materials
        Material mat;
        // 0) light1, will only emit
        mMaterials.push_back(mat);
        // 1) light2, will only emit
        mMaterials.push_back(mat);

        // 2) glossy white floor
        mat.Reset();
        mat.mDiffuseReflectance = Vec3f(0.1f);
        mat.mPhongReflectance   = Vec3f(0.7f);
        mat.mPhongExponent      = 90.f;
        if( full_glossy )
        {
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        // 3) diffuse green left wall
        mat.Reset();
        mat.mDiffuseReflectance = Vec3f(0.156863f, 0.803922f, 0.172549f);
        if( full_glossy )
        {
            mat.mPhongReflectance   = Vec3f(0.4f);
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        // 4) diffuse red right wall
        mat.Reset();
        mat.mDiffuseReflectance = Vec3f(0.803922f, 0.152941f, 0.152941f);
        if( full_glossy )
        {
            mat.mPhongReflectance   = Vec3f(0.4f);
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        // 5) diffuse white back wall
        mat.Reset();
        mat.mDiffuseReflectance = Vec3f(0.803922f, 0.803922f, 0.803922f);
        if( full_glossy )
        {
            mat.mPhongReflectance   = Vec3f(0.4f);
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        // 6) mirror ball
        mat.Reset();
        if( !full_glossy )
        {
            mat.mMirrorReflectance = Vec3f(1.f);
        }
        else
        {
            mat.mDiffuseReflectance = Copper / 15.f;
            mat.mPhongReflectance   = Copper;
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        // 7) glass ball
        mat.Reset();
        if( !full_glossy )
        {
            mat.mMirrorReflectance  = Vec3f(1.f);
            mat.mIOR                = 1.6f;
        }
        else
        {
            mat.mDiffuseReflectance = Copper / 15.f;
            mat.mPhongReflectance   = Copper;
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        // 8) diffuse blue wall (back wall for glossy floor)
        mat.Reset();
        mat.mDiffuseReflectance = Vec3f(0.156863f, 0.172549f, 0.803922f);
        if( full_glossy )
        {
            mat.mPhongReflectance   = Vec3f(0.4f);
            mat.mPhongExponent      = aPhongExp;
        }
        mMaterials.push_back(mat);

        delete mGeometry;

        //////////////////////////////////////////////////////////////////////////
        // Cornell box
        Vec3f cb[8] = {
            Vec3f(-1.27029f,  1.30455f, -1.28002f),
            Vec3f( 1.28975f,  1.30455f, -1.28002f),
            Vec3f( 1.28975f,  1.30455f,  1.28002f),
            Vec3f(-1.27029f,  1.30455f,  1.28002f),
            Vec3f(-1.27029f, -1.25549f, -1.28002f),
            Vec3f( 1.28975f, -1.25549f, -1.28002f),
            Vec3f( 1.28975f, -1.25549f,  1.28002f),
            Vec3f(-1.27029f, -1.25549f,  1.28002f)
        };

        GeometryList *geometryList = new GeometryList;
        mGeometry = geometryList;

        if((aBoxMask & kGlossyFloor) != 0)
        {
            // Floor
            geometryList->mGeometry.push_back(new Triangle(cb[0], cb[4], cb[5], 2));
            geometryList->mGeometry.push_back(new Triangle(cb[5], cb[1], cb[0], 2));
            // Back wall
            geometryList->mGeometry.push_back(new Triangle(cb[0], cb[1], cb[2], 8));
            geometryList->mGeometry.push_back(new Triangle(cb[2], cb[3], cb[0], 8));
        }
        else
        {
            // Floor
            geometryList->mGeometry.push_back(new Triangle(cb[0], cb[4], cb[5], 5));
            geometryList->mGeometry.push_back(new Triangle(cb[5], cb[1], cb[0], 5));
            // Back wall
            geometryList->mGeometry.push_back(new Triangle(cb[0], cb[1], cb[2], 5));
            geometryList->mGeometry.push_back(new Triangle(cb[2], cb[3], cb[0], 5));
        }


        // Ceiling
        if(light_ceiling && !light_box)
        {
            geometryList->mGeometry.push_back(new Triangle(cb[2], cb[6], cb[7], 0));
            geometryList->mGeometry.push_back(new Triangle(cb[7], cb[3], cb[2], 1));
        }
        else
        {
            geometryList->mGeometry.push_back(new Triangle(cb[2], cb[6], cb[7], 5));
            geometryList->mGeometry.push_back(new Triangle(cb[7], cb[3], cb[2], 5));
        }

        // Left wall
        geometryList->mGeometry.push_back(new Triangle(cb[3], cb[7], cb[4], 3));
        geometryList->mGeometry.push_back(new Triangle(cb[4], cb[0], cb[3], 3));

        // Right wall
        geometryList->mGeometry.push_back(new Triangle(cb[1], cb[5], cb[6], 4));
        geometryList->mGeometry.push_back(new Triangle(cb[6], cb[2], cb[1], 4));

        // Ball - central
        float largeRadius = 0.8f;
        Vec3f center = (cb[0] + cb[1] + cb[4] + cb[5]) * (1.f / 4.f) + Vec3f(0, 0, largeRadius);

        if((aBoxMask & kLargeMirrorSphere) != 0)
            geometryList->mGeometry.push_back(new Sphere(center, largeRadius, 6));

        if((aBoxMask & kLargeGlassSphere) != 0)
            geometryList->mGeometry.push_back(new Sphere(center, largeRadius, 7));

        if((aBoxMask & kLargeGlossySphere) != 0)
            geometryList->mGeometry.push_back(new Sphere(center, largeRadius, 7));

        // Balls - left and right
        float smallRadius = 0.5f;
        Vec3f leftWallCenter  = (cb[0] + cb[4]) * (1.f / 2.f) + Vec3f(0, 0, smallRadius);
        Vec3f rightWallCenter = (cb[1] + cb[5]) * (1.f / 2.f) + Vec3f(0, 0, smallRadius);
        float xlen = rightWallCenter.x - leftWallCenter.x;
        Vec3f leftBallCenter  = leftWallCenter  + Vec3f(2.f * xlen / 7.f, 0, 0);
        Vec3f rightBallCenter = rightWallCenter - Vec3f(2.f * xlen / 7.f, 0, 0);

        if((aBoxMask & kSmallMirrorSphere) != 0)
            geometryList->mGeometry.push_back(new Sphere(leftBallCenter,  smallRadius, 6));

        if((aBoxMask & kSmallGlassSphere) != 0)
            geometryList->mGeometry.push_back(new Sphere(rightBallCenter, smallRadius, 7));

        //////////////////////////////////////////////////////////////////////////
        // Light box at the ceiling
        Vec3f lb[8] = {
            Vec3f(-0.25f,  0.25f, 1.26002f),
            Vec3f( 0.25f,  0.25f, 1.26002f),
            Vec3f( 0.25f,  0.25f, 1.28002f),
            Vec3f(-0.25f,  0.25f, 1.28002f),
            Vec3f(-0.25f, -0.25f, 1.26002f),
            Vec3f( 0.25f, -0.25f, 1.26002f),
            Vec3f( 0.25f, -0.25f, 1.28002f),
            Vec3f(-0.25f, -0.25f, 1.28002f)
        };

        if(light_box)
        {
            // Back wall
            geometryList->mGeometry.push_back(new Triangle(lb[0], lb[2], lb[1], 5));
            geometryList->mGeometry.push_back(new Triangle(lb[2], lb[0], lb[3], 5));
            // Left wall
            geometryList->mGeometry.push_back(new Triangle(lb[3], lb[4], lb[7], 5));
            geometryList->mGeometry.push_back(new Triangle(lb[4], lb[3], lb[0], 5));
            // Right wall
            geometryList->mGeometry.push_back(new Triangle(lb[1], lb[6], lb[5], 5));
            geometryList->mGeometry.push_back(new Triangle(lb[6], lb[1], lb[2], 5));
            // Front wall
            geometryList->mGeometry.push_back(new Triangle(lb[4], lb[5], lb[6], 5));
            geometryList->mGeometry.push_back(new Triangle(lb[6], lb[7], lb[4], 5));

            if(light_ceiling)
            {
                // Floor
                geometryList->mGeometry.push_back(new Triangle(lb[0], lb[5], lb[4], 0));
                geometryList->mGeometry.push_back(new Triangle(lb[5], lb[0], lb[1], 1));
            }
            else
            {
                // Floor
                geometryList->mGeometry.push_back(new Triangle(lb[0], lb[5], lb[4], 5));
                geometryList->mGeometry.push_back(new Triangle(lb[5], lb[0], lb[1], 5));
            }
        }

        //////////////////////////////////////////////////////////////////////////
        // Lights
        if(light_ceiling && !light_box)
        {
            // Without light box, whole ceiling is light
            mLights.resize(2);
            AreaLight *l = new AreaLight(cb[2], cb[6], cb[7]);
            l->mIntensity = Vec3f(0.95492965f);
            mLights[0] = l;
            mMaterial2Light.insert(std::make_pair(0, 0));

            l = new AreaLight(cb[7], cb[3], cb[2]);
            l->mIntensity = Vec3f(0.95492965f);
            mLights[1] = l;
            mMaterial2Light.insert(std::make_pair(1, 1));
        }
        else if(light_ceiling && light_box)
        {
            // With light box
            mLights.resize(2);
            AreaLight *l = new AreaLight(lb[0], lb[5], lb[4]);
            //l->mIntensity = Vec3f(0.95492965f);
            l->mIntensity = Vec3f(25.03329895614464f);
            mLights[0] = l;
            mMaterial2Light.insert(std::make_pair(0, 0));

            l = new AreaLight(lb[5], lb[0], lb[1]);
            //l->mIntensity = Vec3f(0.95492965f);
            l->mIntensity = Vec3f(25.03329895614464f);
            mLights[1] = l;
            mMaterial2Light.insert(std::make_pair(1, 1));
        }

        if(light_sun)
        {
            DirectionalLight *l = new DirectionalLight(Vec3f(-1.f, 1.5f, -1.f));
            l->mIntensity = Vec3f(0.5f, 0.2f, 0.f) * 20.f;
            mLights.push_back(l);
        }

        if(light_point)
        {
            PointLight *l = new PointLight(Vec3f(0.0, -0.5, 1.0));
            l->mIntensity = Vec3f(70.f * (INV_PI_F * 0.25f));
            mLights.push_back(l);
        }

        if(light_background)
        {
            BackgroundLight *l = new BackgroundLight;
            l->mScale = 1.f;
            mLights.push_back(l);
            mBackground = l;
        }
    }

    // add support for multi view path tracing
    void LoadMultiView( const Vec2i &aResolution )
    {
        // printf("loading cameras...\n");
        
        // initial camera parameters (middle camera)
        Vec3f position(-0.0439815f, -4.12529f, 0.222539f);
        Vec3f forward(0.00688625f, 0.998505f, -0.0542161f);
        Vec3f up(3.73896e-4f, 0.0542148f, 0.998529f);

        // find point of focus, just a hack...
        Vec3f focus = 4 * forward + position;

        // number of cameras
        const int nb = 3; 
        const float d = 0.25; 

        // setup translating cameras
        for(int i = -(nb-1)/2; i <= (nb-1)/2; i++)
        {
            Vec3f p = position + Vec3f(i * d, 0, 0);
            Camera c;
            c.Setup(p, focus - p, up, Vec2f(float(aResolution.x), float(aResolution.y)), 45);
            mCameras.push_back(c);
        }

        mCameraNumber = 3;
    }

    void BuildSceneSphere()
    {
        Vec3f bboxMin( 1e36f);
        Vec3f bboxMax(-1e36f);
        mGeometry->GrowBBox(bboxMin, bboxMax);

        const float radius2 = (bboxMax - bboxMin).LenSqr();

        mSceneSphere.mSceneCenter = (bboxMax + bboxMin) * 0.5f;
        mSceneSphere.mSceneRadius = std::sqrt(radius2) * 0.5f;
        mSceneSphere.mInvSceneRadiusSqr = 1.f / Sqr(mSceneSphere.mSceneRadius);
    }

    static std::string GetSceneName(
        uint        aBoxMask,
        std::string *oAcronym = NULL)
    {
        std::string name;
        std::string acronym;

        // Floor type
        if((aBoxMask & kGlossyFloor) == kGlossyFloor)
        {
            name    += "glossy ";
            acronym += "g";
        }

        // Box content
        if((aBoxMask & kBothSmallSpheres) == kBothSmallSpheres)
        {
            name    += "small spheres";
            acronym += "bs";
        }
        else if((aBoxMask & kSmallMirrorSphere) == kSmallMirrorSphere)
        {
            name    += "small mirror sphere";
            acronym += "sm";
        }
        else if((aBoxMask & kSmallGlassSphere) == kSmallGlassSphere)
        {
            name    += "small glass sphere";
            acronym += "sg";
        }
        else if((aBoxMask & kLargeMirrorSphere) == kLargeMirrorSphere)
        {
            name    += "large mirror sphere";
            acronym += "lm";
        }
        else if((aBoxMask & kLargeGlassSphere) == kLargeGlassSphere)
        {
            name    += "large glass sphere";
            acronym += "lg";
        }
        else if((aBoxMask & kLargeGlossySphere) == kLargeGlossySphere)
        {
            name    += "large glossy sphere";
            acronym += "ls";
        }
        else
        {
            name    += "empty";
            acronym += "e";
        }

        acronym += "_";

        // Lighting
        if((aBoxMask & kLightCeiling) == kLightCeiling)
        {
            name    += " + ceiling (area)";
            acronym += "c";
        }
        else if((aBoxMask & kLightSun) == kLightSun)
        {
            name    += " + sun (directional)";
            acronym += "s";
        }
        else if((aBoxMask & kLightPoint) == kLightPoint)
        {
            name    += " + point";
            acronym += "p";
        }
        else if((aBoxMask & kLightBackground) == kLightBackground)
        {
            name    += " + background (env. lighting)";
            acronym += "b";
        }

        if(oAcronym) *oAcronym = acronym;
        return name;
    }

public:

    AbstractGeometry      *mGeometry;
    Camera                mCamera;
    std::vector<Material> mMaterials;
    std::vector<AbstractLight*>   mLights;
    std::map<int, int>    mMaterial2Light;
    SceneSphere           mSceneSphere;
    BackgroundLight*      mBackground;

    std::string           mSceneName;
    std::string           mSceneAcronym;
    
    // add support for multi view
    std::vector<Camera>   mCameras;
    int                   mCameraNumber;
};

#endif //__SCENE_HXX__
