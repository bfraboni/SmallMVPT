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

#ifndef __MULTIVIEW_HXX__
#define __MULTIVIEW_HXX__

#include <vector>
#include <cmath>
#include "renderer.hxx"
#include "bsdf.hxx"
#include "rng.hxx"

class MultiViewPathTracer : public AbstractRenderer
{
public:

    MultiViewPathTracer(
        const Scene& aScene,
        int aSeed = 1234,
        bool aReuse = false
    ) :
        AbstractRenderer(aScene), mRng(aSeed), mFramebuffers(aScene.mCameraNumber), mReuse(aReuse)
    {
        for( int fbID = 0; fbID < (int)mFramebuffers.size(); ++fbID )
        {
            mFramebuffers[fbID].Setup(aScene.mCameras[fbID].mResolution);
        }
    }

    void GetFramebuffer(Framebuffer& oFramebuffer, int aFbID = -1) override
    {
        oFramebuffer = mFramebuffers[aFbID];

        // if(mIterations > 0)
        //     oFramebuffer.Scale(1.f / mIterations);
    }

    /*
        We want to determine how far a BRDF distribution is from 
        another one, the same material properties for 2 different 
        outgoing directions. We need on-the-fly distance 
        evaluation between our BDRF distributions. 
        And the distance value have to be normalized in [0,1] 
        to serve as a probability. Most statistical distances 
        rely on integration or are not normalized, thus do 
        not suit here.

        The Total Variation distance is convenient since we
        only have to evaluate the maximum difference between 
        the distributions, and is easily normalized.

        But the maximum of a function is difficult to find
        if we can not analytically derive a close formulae for
        the roots of the derivative. Which our Mathematica skills
        could not resolve yet.

        Hence we assume that the maximum difference between 
        distributions occurs at the mirror direction of the 
        first distribution.

        This assumption is true when glossy materials tends
        towards mirrors, and false when materials tends towards 
        diffuse. However the maximum difference between distributions
        becomes smaller when materials tends towards diffuse.
        For that reason this approximation is not far from the 
        real distance.
    */  
    float TotalVariationDistance(
        const Isect& aIsect, 
        const Scene& aScene,
        const Frame& frame, 
        const Ray& aRay1, 
        const Ray& aRay2,
        const BSDF<false>& bsdf)
    {
        // the BSDF can not be a Dirac Delta, we ensured this before
        // if the BSDF is diffuse, namely view-independent, the TV distance is 0
        if( bsdf.PhongProb() == 0 ) 
            return 0;
        
        // construct mirror directions
        const Vec3f worldRefl1 = Reflect(-aRay1.dir, aIsect.normal);
        const Vec3f worldRefl2 = Reflect(-aRay2.dir, aIsect.normal);
        // compute pdf w.r.t ray1 direction
        float pdf1 = bsdf.Pdf(aScene, worldRefl1);
        if( pdf1 == 0 ) 
            return 1;
        
        // compute pdf w.r.t ray2 direction 
        // expected max of the distribution for glossy materials
        float pdf2 = bsdf.Pdf(aScene, worldRefl2);
        if( pdf2 == 0 ) 
            return 1;
        
        // compute ratio
        float q = pdf2 > pdf1 ? pdf1 / pdf2 : pdf2 / pdf1;
            
        // return Total Variation distance approximation
        return std::abs(1.f - q);
    }

    // Compute the Jacobian of the transformation from one camera to another
    // for Pinhole camera model without lens.
    // |T'| = |d(pl)/d(pf)| * |d(ps)/d(pl)| * |d(pl')/d(pf')|^{-1} * |d(ps)/d(pl')|^{-1}
    // For thin lens cameras model please refer to the paper formulas.
    float Jacobian(
        const Vec3f& aHitPoint, 
        const Vec3f& aHitNormal, 
        const Ray& aRay1, 
        const int aCam1, 
        const Ray& aRay2,
        const int aCam2)
    {
        float lenSqr1 = (aHitPoint - aRay1.org).LenSqr();
        float cosTheta1 = std::max(0.f, Dot(aHitNormal, -aRay1.dir));
        float dpdl1 = cosTheta1 / lenSqr1;

        const Camera& cam1 = mScene.mCameras[aCam1];
        float cosThetaFilm1 = Dot( aRay1.dir, cam1.mForward );
        float lenSqrFilm1 = cam1.mImagePlaneDist * cam1.mImagePlaneDist / (cosThetaFilm1 * cosThetaFilm1);
        float dldf1 = lenSqrFilm1 / cosThetaFilm1;
        
        float lenSqr2 = (aHitPoint - aRay2.org).LenSqr();
        float cosTheta2 = std::max(0.f, Dot(aHitNormal, -aRay2.dir));
        float dpdl2 = cosTheta2 / lenSqr2;

        const Camera& cam2 = mScene.mCameras[aCam2];
        float cosThetaFilm2 = Dot( aRay2.dir, cam2.mForward );
        float lenSqrFilm2 = cam2.mImagePlaneDist * cam2.mImagePlaneDist / (cosThetaFilm2 * cosThetaFilm2);
        float dldf2 = lenSqrFilm2 / cosThetaFilm2;
    
        return dpdl1 * dldf1 / (dpdl2 * dldf2);        
    }

    // Utility structure to store information of a sample / reuse ray
    struct Data
    {
        Ray ray;
        Vec3f color = Vec3f(0.f);
        Vec3f bsdfDirect = Vec3f(0.f);
        Vec3f bsdfIndirect = Vec3f(0.f);
        Vec2f raster = Vec2f(-1, -1);
        float   weight = 1.f,
                pdf = 1.f,
                pdfTransformed = 1.f,
                pdfIndirect = 0,
                J = 1;
        int camID = -1;
        BSDF<false> bsdf;
    };

    // Performs camera selection for path reusing in multi view path tracing.
    void selection(std::vector<Data>& oData, const Isect& aIsect, const Vec3f& aPoint, const Vec3f& aNormal)
    {
        Data& initialData = oData.back(); 
        Frame frame( aNormal );

        // camera subset selection
        for(int otherCamID = 0; otherCamID < (int)mScene.mCameras.size(); otherCamID++)
        {
            if( otherCamID == initialData.camID ) 
                continue;
            const Camera& otherCamera = mScene.mCameras[otherCamID];
            Data data;

            // compute projection ray
            data.ray = otherCamera.ConnectionRay( aPoint );

            // check camera orientation
            if( Dot( data.ray.dir, aNormal ) >= 0 ) 
                continue;

            // check projection on image plane
            data.raster = otherCamera.WorldToRaster( aPoint ); 
            if( !otherCamera.CheckRaster(data.raster) ) 
                continue;

            // compute the Jacobian determinant of the transformation between cameras
            float jacobian = Jacobian(aPoint, aNormal, initialData.ray, initialData.camID, data.ray, otherCamID);
            if( jacobian <= 0 ) 
                continue;

            // compute the selection probability w.r.t Jacobian
            float jacobianProb = std::min(1.f, jacobian);

            // compute the selection probability w.r.t material
            data.bsdf = BSDF<false>(data.ray, aIsect, mScene);
            float tv = TotalVariationDistance(aIsect, mScene, frame, initialData.ray, data.ray, data.bsdf);
            if( tv >= 1 ) 
                continue;
            const Material &mat = mScene.GetMaterial( aIsect.matID );
            float materialProb = std::pow(1.f - tv, mat.mPhongExponent);

            // Russian Roulette
            float alpha = jacobianProb * materialProb;
            if( mRng.GetFloat() > alpha ) 
                continue;

            // check visibility
            if( mScene.Occluded(data.ray.org, data.ray.dir, data.ray.tmin) ) 
                continue;
            
            data.pdf = otherCamera.Pdf( data.ray );
            data.pdfTransformed = initialData.pdf * jacobian * alpha;
            data.camID = otherCamID;

            oData.push_back( data );
        } 
        // end camera selection

        // MIS weight computation for all visible cameras
        for(Data& datak : oData)
        {
            float pdfSum = 0;
            for(Data& dataj : oData)
            {
                // we already computed and stored the transformed pdf from the initial ray
                if( dataj.camID == initialData.camID )
                {
                    pdfSum += datak.pdfTransformed;
                    continue;
                }
                
                // compute the Jacobian determinant of the transformation between cameras
                float jacobian = Jacobian(aPoint, aNormal, dataj.ray, dataj.camID, datak.ray, datak.camID );

                // compute the selection probability w.r.t Jacobian
                float jacobianProb = std::min(1.f, jacobian);

                // compute the selection probability w.r.t material distance
                float tv = TotalVariationDistance(aIsect, mScene, frame, dataj.ray, datak.ray, datak.bsdf);
                const Material &mat = mScene.GetMaterial(aIsect.matID);
                float materialProb = std::pow(1.f - tv, mat.mPhongExponent);

                // selection
                float alpha = jacobianProb * materialProb;
                
                // compute pdf of reuse from j to k
                float pdf_j_to_k = 
                    // ray pdf
                    dataj.pdf 
                    // Jacobian of transformation
                    * jacobian
                    // selection probability
                    * alpha;

                // add to sum
                pdfSum += pdf_j_to_k;
            }

            // mis weight
            datak.weight = datak.pdfTransformed / pdfSum;    
        } 
        // end mis weights computation     
    }

    // Compute and add direct lightning for next event estimation
    void direct(std::vector<Data>& oData, const Vec3f& aPoint, const Vec3f& aNormal)
    {
        static const int   lightCount    = mScene.GetLightCount();
        static const float lightPickProb = 1.f / lightCount;

        // random light pick
        int lightID = int(mRng.GetFloat() * lightCount);
        const AbstractLight *light = mScene.GetLightPtr(lightID);

        // compute incoming radiance
        Vec3f directionToLight;
        float distance, directPdfW;
        Vec3f radiance = light->Illuminate(mScene.mSceneSphere, aPoint,
            mRng.GetVec2f(), directionToLight, distance, directPdfW);

        // add direct contribution
        bool visibility = !mScene.Occluded(aPoint, directionToLight, distance);
        if(!radiance.IsZero() && visibility)
        {
            // compute the mixture of selected observers pdf 
            float mixturePdfW = 0, count = 0;
            for( Data& data : oData )
            {
                float bsdfPdfW, cosThetaOut;
                data.bsdfDirect = data.bsdf.Evaluate(mScene, directionToLight, cosThetaOut, &bsdfPdfW);
                if( bsdfPdfW > 0 )
                {
                    mixturePdfW += bsdfPdfW;
                    count += 1;
                }
            }
            if( count > 0 ) mixturePdfW /= count;

            // compute MIS weight for direct illumination
            float weight = !light->IsDelta() ? Mis2(directPdfW * lightPickProb, mixturePdfW) : 1;

            // compute common quantity
            Vec3f common = (weight * Dot(aNormal, directionToLight)) / (lightPickProb * directPdfW) * radiance;
            
            // accumulate direct contributions for each visible camera
            for( Data& data : oData )
                data.color += data.J * data.bsdfDirect * common;
        }
    }

    void ambient(std::vector<Data>& oData, const Vec3f& aPoint, const Vec3f& aNormal)
    {
        // ambient occlusion
        static const int nb = 8;
        static const float invNb = 1.f / nb;
        static const float invPi = 1.f / M_PI;

        Frame frame( aNormal );

        for( int occID = 0; occID < nb; ++occID)
        {
            float pdf = 0;
            Vec3f localDir = SampleCosHemisphereW( mRng.GetVec2f(), &pdf );
            Vec3f worldDir = frame.ToWorld(localDir);

            if( !mScene.Occluded( aPoint, worldDir, mScene.mSceneSphere.mSceneRadius * 0.2 ))
            {   
                float common = localDir.z * invNb * invPi;
                for( Data& data : oData )
                {
                    data.color += common;
                }
            }
        }
    }

    virtual void RunIteration(int aIteration)
    {
        // We sample lights uniformly
        const int   lightCount    = mScene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        // Sample all cameras
        std::vector<Data> visibleData;
        visibleData.reserve( mScene.mCameras.size() );
        for(int camID = 0; camID < (int)mScene.mCameras.size(); camID++)
        {
            const Camera& camera = mScene.mCameras[camID];

            const int resX = int(camera.mResolution.x);
            const int resY = int(camera.mResolution.y);

            for(int pixID = 0; pixID < resX * resY; pixID++)
            {
                const int x = pixID % resX;
                const int y = pixID / resX;

                const Vec2f sample = Vec2f(float(x), float(y)) + mRng.GetVec2f();

                Isect isect;
                isect.dist = 1e36f;

                Data primary;
                primary.raster = sample;
                primary.color = Vec3f(0.f);
                primary.ray = camera.GenerateRay(sample);
                primary.pdf = camera.Pdf( primary.ray );
                primary.pdfTransformed = primary.pdf;
                primary.camID = camID;

                visibleData.push_back(primary);

                Ray& ray = visibleData.back().ray;

                if( mScene.Intersect(ray, isect) )
                {
                    Vec3f hitPoint = ray.org + ray.dir * isect.dist;
                    Vec3f hitNormal = isect.normal;

                    if( Dot(ray.dir, hitNormal) >= 0 ) 
                        continue;

                    isect.dist += EPS_RAY;
                    visibleData.back().bsdf = BSDF<false>(ray, isect, mScene);
                    const BSDF<false>& bsdf = visibleData.back().bsdf; 
                    if(!bsdf.IsValid())
                        continue;

                    // Directly hit some light, lights do not reflect
                    if( isect.lightID >= 0 )
                    {
                        const AbstractLight *light = mScene.GetLightPtr(isect.lightID);
                        Vec3f radiance = light->Radiance();
                        // do not reuse directly hit light can lead to artefacts at light edges
                        // and does not reduce variance since radiance is constant   
                        visibleData.back().color += radiance;
                    }
                    else
                    {
                        // Compute multi view projection and associated MIS weights
                        if((!bsdf.IsDelta()) && mReuse)
                        {
                            selection(visibleData, isect, hitPoint, hitNormal);
                        }

                        // Next event estimation
                        if(!bsdf.IsDelta())
                        {
                            direct(visibleData, hitPoint, hitNormal);
                        }

                        // Continue random walk
                        if( 1 )
                        {
                            int pickID = int(mRng.GetFloat() * visibleData.size());
                            
                            Vec3f rndTriplet = mRng.GetVec3f();
                            float pdf, cosThetaOut;
                            uint  sampledEvent;
                            Vec3f direction;
                            Vec3f sample = visibleData[pickID].bsdf.Sample(mScene, rndTriplet, direction, pdf, cosThetaOut, &sampledEvent);

                            bool specular = (sampledEvent & BSDF<true>::kSpecular) != 0;

                            float pdfW = 0;
                            if( specular )
                            {
                                visibleData.back().bsdfIndirect = sample;
                                pdfW = pdf;
                            }
                            else
                            {
                                // Compute the mixture pdf foreach selected observers
                                float count = 0;
                                for( Data& data : visibleData )
                                {   
                                    // Frame frame( hitNormal );

                                    // float mixture = 0, alphas = 0;
                                    // // compute the mixture
                                    // for( Data& data2 : visibleData )
                                    // {
                                    //     // compute bsdf value and pdf;
                                    //     float bsdfPdfW, cosThetaOut;
                                    //     data2.bsdfIndirect = data2.bsdf.Evaluate(mScene, direction, cosThetaOut, &bsdfPdfW);

                                    //     // compute alphas
                                    //     float tv = TotalVariationDistance(isect, mScene, frame, data.ray, data2.ray, data2.bsdf);
                                    //     const Material &mat = mScene.GetMaterial(isect.matID);
                                    //     float pm = std::pow(1.f - tv, mat.mPhongExponent);

                                    //     // compute mixture
                                    //     mixture += pm * bsdfPdfW;
                                    //     alphas += pm;
                                    // }

                                    // data.pdfIndirect = alphas > 0 ? mixture / alphas : 0;

                                    if( 1 )
                                    {
                                        float bsdfPdfW, cosThetaOut;
                                        data.bsdfIndirect = data.bsdf.Evaluate(mScene, direction, cosThetaOut, &bsdfPdfW);
                                        if( bsdfPdfW > 0 )
                                        {
                                            pdfW += bsdfPdfW;
                                            count += 1;
                                        }
                                    }
                                }
                                if( count > 0 ) pdfW /= count;
                            }

                            if( pdfW > 0 )
                            {
                                // We offset ray origin instead of setting tmin due to numeric
                                // issues in ray-sphere intersection. The isect.dist has to be
                                // extended by this EPS_RAY after hitpoint is determined
                                Ray nextRay(hitPoint + EPS_RAY * direction, direction, 0);
                                isect.dist = 1e36f;

                                // compute incoming radiance
                                Vec3f contrib = indirect(nextRay, isect, Vec3f(1), 2, pdfW, specular);
                                
                                // accumulate indirect illumination
                                for( Data& data : visibleData )
                                {
                                    // if( data.pdfIndirect <= 0 && !specular ) continue;
                                    // float common = specular ? cosThetaOut / pdfW : cosThetaOut / data.pdfIndirect;
                                    float common = cosThetaOut / pdfW;
                                    data.color += data.J * common * contrib * data.bsdfIndirect;
                                }
                            }
                        }
                    }
                }

                // accumulate contributions
                for(const Data& data : visibleData)
                {
                    mFramebuffers[data.camID].AddColor(data.raster, data.color, data.weight);
                }

                visibleData.clear();
            }
        }
        mIterations++;
    }

private:

    // Original code from SmallVCM "pathtracer.hxx"
    Vec3f indirect(Ray ray, Isect isect, Vec3f pathWeight, int pathLength, float lastPdfW, bool lastSpecular)
    {
        static const int   lightCount    = mScene.GetLightCount();
        static const float lightPickProb = 1.f / lightCount;

        Vec3f color(0.f);

        for(;; ++pathLength)
        {
            if(!mScene.Intersect(ray, isect))
            {
                if(pathLength < mMinPathLength)
                    break;

                const BackgroundLight* background = mScene.GetBackground();
                if(!background)
                    break;
                // For background we cheat with the A/W suffixes,
                // and GetRadiance actually returns W instead of A
                float directPdfW;
                Vec3f contrib = background->GetRadiance(mScene.mSceneSphere,
                    ray.dir, Vec3f(0), &directPdfW);
                if(contrib.IsZero())
                    break;

                float misWeight = 1.f;
                if(pathLength > 1 && !lastSpecular)
                {
                    misWeight = Mis2(lastPdfW, directPdfW * lightPickProb);
                }

                color += pathWeight * misWeight * contrib;
                break;
            }

            Vec3f hitPoint = ray.org + ray.dir * isect.dist;
            Vec3f hitNormal = isect.normal;

            isect.dist += EPS_RAY;

            BSDF<false> bsdf(ray, isect, mScene);
            if(!bsdf.IsValid())
                break;

            // directly hit some light, lights do not reflect
            if(isect.lightID >= 0)
            {
                if(pathLength < mMinPathLength)
                    break;

                const AbstractLight *light = mScene.GetLightPtr(isect.lightID);
                float directPdfA;
                Vec3f contrib = light->GetRadiance(mScene.mSceneSphere,
                    ray.dir, hitPoint, &directPdfA);
                if(contrib.IsZero())
                    break;

                float misWeight = 1.f;
                if(pathLength > 1 && !lastSpecular)
                {
                    const float directPdfW = PdfAtoW(directPdfA, isect.dist,
                        bsdf.CosThetaFix());
                    misWeight = Mis2(lastPdfW, directPdfW * lightPickProb);
                }

                color += pathWeight * misWeight * contrib;
                break;  
            }

            if(pathLength >= mMaxPathLength)
                break;

            if(bsdf.ContinuationProb() == 0)
                break;

            // next event estimation
            if(!bsdf.IsDelta() && pathLength + 1 >= mMinPathLength)
            {
                int lightID = int(mRng.GetFloat() * lightCount);
                const AbstractLight *light = mScene.GetLightPtr(lightID);

                Vec3f directionToLight;
                float distance, directPdfW;
                Vec3f radiance = light->Illuminate(mScene.mSceneSphere, hitPoint,
                    mRng.GetVec2f(), directionToLight, distance, directPdfW);

                if(!radiance.IsZero())
                {
                    float bsdfPdfW, cosThetaOut;
                    const Vec3f factor = bsdf.Evaluate(mScene,
                        directionToLight, cosThetaOut, &bsdfPdfW);

                    if(!factor.IsZero())
                    {
                        float weight = 1.f;
                        if(!light->IsDelta())
                        {
                            const float contProb = bsdf.ContinuationProb();
                            bsdfPdfW *= contProb;
                            weight = Mis2(directPdfW * lightPickProb, bsdfPdfW);
                        }

                        Vec3f contrib = (weight * cosThetaOut / (lightPickProb * directPdfW)) *
                            (radiance * factor);

                        if(!mScene.Occluded(hitPoint, directionToLight, distance))
                        {
                            color += pathWeight * contrib;
                        }
                    }
                }
            }

            // continue random walk
            {
                Vec3f rndTriplet = mRng.GetVec3f();
                float pdf, cosThetaOut;
                uint  sampledEvent;

                Vec3f factor = bsdf.Sample(mScene, rndTriplet, ray.dir,
                    pdf, cosThetaOut, &sampledEvent);

                if(factor.IsZero())
                    break;

                // Russian roulette
                const float contProb = bsdf.ContinuationProb();

                lastSpecular = (sampledEvent & BSDF<true>::kSpecular) != 0;
                lastPdfW     = pdf * contProb;

                if(contProb < 1.f)
                {
                    if(mRng.GetFloat() > contProb)
                    {
                        break;
                    }
                    pdf *= contProb;
                }

                pathWeight *= factor * (cosThetaOut / pdf);
                
                // We offset ray origin instead of setting tmin due to numeric
                // issues in ray-sphere intersection. The isect.dist has to be
                // extended by this EPS_RAY after hitpoint is determined
                ray.org    = hitPoint + EPS_RAY * ray.dir;
                ray.tmin   = 0.f;
                isect.dist = 1e36f;
            }
        }
        
        return color;
    }

    // Mis power (1 for balance heuristic)
    float Mis(float aPdf) const
    {
        return aPdf;
    }

    // Mis weight for 2 pdfs
    float Mis2(
        float aSamplePdf,
        float aOtherPdf) const
    {
        return Mis(aSamplePdf) / (Mis(aSamplePdf) + Mis(aOtherPdf));
    }

    // Add multi-view support
    std::vector<Framebuffer> mFramebuffers;
private:

    Rng mRng;
    bool mReuse = true;
};

#endif //__MULTIVIEW_HXX__
