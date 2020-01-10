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

#include <vector>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <sstream>
#include <chrono>
#include <iomanip>
#include "math.hxx"
#include "ray.hxx"
#include "geometry.hxx"
#include "camera.hxx"
#include "framebuffer.hxx"
#include "scene.hxx"
#include "bsdf.hxx"
#include "html_writer.hxx"
#include "config.hxx"

#ifndef NO_OMP
#include <omp.h>
#endif
#include <string>
#include <set>
#include <sstream>

//////////////////////////////////////////////////////////////////////////
// The main rendering function, renders what is in aConfig

float multirender(
    Config &aConfig,
    int *oUsedIterations = NULL)
{
    // Set number of used threads
#ifndef NO_OMP
    omp_set_num_threads(aConfig.mNumThreads);
#endif

    // Create 1 renderer per thread
    typedef AbstractRenderer* AbstractRendererPtr;
    AbstractRendererPtr *renderers;
    renderers = new AbstractRendererPtr[aConfig.mNumThreads];

    for(int i=0; i<aConfig.mNumThreads; i++)
    {
        renderers[i] = CreateRenderer(aConfig, aConfig.mBaseSeed + i);
        renderers[i]->mMaxPathLength = aConfig.mMaxPathLength;
        renderers[i]->mMinPathLength = aConfig.mMinPathLength;
    }

    clock_t startT = clock();
    auto start = std::chrono::system_clock::now();
    int iter = 0;

    // Rendering loop, when we have any time limit, use time-based loop,
    // otherwise go with required iterations
    if(aConfig.mMaxTime > 0)
    {
        // Time based loop
#pragma omp parallel
        // while(clock() < startT + aConfig.mMaxTime*CLOCKS_PER_SEC)
        while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < aConfig.mMaxTime)
        {
#ifndef NO_OMP
            int threadId = omp_get_thread_num();
#else
            int threadId = 0;
#endif
            renderers[threadId]->RunIteration(iter);

#pragma omp atomic
            iter++; // counts number of iterations
        }
    }
    else
    {
        // Iterations based loop
#pragma omp parallel for
        for(iter=0; iter < aConfig.mIterations; iter++)
        {
#ifndef NO_OMP
            int threadId = omp_get_thread_num();
#else
            int threadId = 0;
#endif
            renderers[threadId]->RunIteration(iter);
        }
    }

    clock_t endT = clock();
    auto end = std::chrono::system_clock::now();

    if(oUsedIterations)
        *oUsedIterations = iter+1;

    // Accumulate from all renderers into a common framebuffer
    int usedRenderers = 0;

    // With very low number of iterations and high number of threads
    // not all created renderers had to have been used.
    // Those must not participate in accumulation.
    for(int trheadID=0; trheadID<aConfig.mNumThreads; trheadID++)
    {
        if(!renderers[trheadID]->WasUsed())
            continue;

        // printf("renderer %d\n", trheadID);

        if(usedRenderers == 0)
        {
            for( int fbID = 0; fbID < (int)aConfig.mFramebuffers.size(); ++fbID )
            {
                renderers[trheadID]->GetFramebuffer(aConfig.mFramebuffers[fbID], fbID);
            }
        }
        else
        { 
            for( int fbID = 0; fbID < (int)aConfig.mFramebuffers.size(); ++fbID )
            {
                Framebuffer tmp;
                renderers[trheadID]->GetFramebuffer(tmp, fbID);
                aConfig.mFramebuffers[fbID].Add(tmp);
            }
        }

        usedRenderers++;
    }

    // Scale framebuffer by the number of used renderers
    for( int fbID = 0; fbID < (int)aConfig.mFramebuffers.size(); ++fbID )
    {
        aConfig.mFramebuffers[fbID].Scale(1.f / usedRenderers);
        // aConfig.mFramebuffers[fbID].Normalize();
    }

    // Clean up renderers
    for(int i=0; i<aConfig.mNumThreads; i++)
        delete renderers[i];

    delete [] renderers;

    int elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    // return float(endT - startT) / CLOCKS_PER_SEC;
    return elapsed;
}

//////////////////////////////////////////////////////////////////////////
// Main

int main(int argc, const char *argv[])
{
    // Setups config based on command line
    Config config;
    ParseCommandline(argc, argv, config);

    // If number of threads is invalid, set 1 thread per processor
    if(config.mNumThreads <= 0)
#ifndef NO_OMP
        config.mNumThreads = std::max(1, omp_get_num_procs());
#else
        config.mNumThreads = 1;
#endif

    // When some error has been encountered, exits
    if(config.mScene == NULL)
        return 1;

    // Sets up framebuffer and number of threads
    Framebuffer fbuffer;
    config.mFramebuffer = &fbuffer;

    // Prints what we are doing
    // printf("Scene:   %s\n", config.mScene->mSceneName.c_str());
    // if(config.mMaxTime > 0)
    //     printf("Target:  %g seconds render time\n", config.mMaxTime);
    // else
    //     printf("Target:  %d iteration(s)\n", config.mIterations);

    // // Renders the image
    // printf("Running: %s ", config.GetName(config.mAlgorithm));
    // if( config.mReuse )
    //     printf("+ path reusing");
    // printf("...\n");
    // fflush(stdout);
    float time = multirender( config );
    // printf("Time:    %.0f s\n", time);
    
    printf("%.0f\n", time);

    // Saves the image
    std::string extension = config.mOutputName.substr(config.mOutputName.length() - 3, 3);
    for( int fbID = 0; fbID < (int)config.mFramebuffers.size(); ++fbID )
    {
        std::string filename = config.mOutputName;
        std::stringstream ss;

        if(config.mSceneID == 4)
            ss << "_p" << std::setfill('0') << std::setw(4) << config.mPhongExp;
        
        if(config.mMaxTime > 0)
            ss << "_t" << std::setfill('0') << std::setw(4) << config.mMaxTime;
        else
            ss << "_i" << std::setfill('0') << std::setw(4) << config.mIterations;

        ss << "_f" << fbID;
        filename.insert(filename.length() - 4, ss.str());
        // printf("File:    %s\n", filename.c_str());

        if(extension == "bmp")
            config.mFramebuffers[fbID].SaveBMP(filename.c_str(), 2.2f /*gamma*/);
        else if(extension == "hdr")
            config.mFramebuffers[fbID].SaveHDR(filename.c_str());
        else
            printf("Used unknown extension %s\n", extension.c_str());
    }

    // Scene cleanup
    delete config.mScene;

    return 0;
}
