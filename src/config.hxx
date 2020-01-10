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
#ifndef __CONFIG_HXX__
#define __CONFIG_HXX__

#include <vector>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include "math.hxx"
#include "ray.hxx"
#include "geometry.hxx"
#include "camera.hxx"
#include "framebuffer.hxx"
#include "scene.hxx"
#include "multiview.hxx"
#include "bsdf.hxx"
#include "html_writer.hxx"

#ifndef NO_OMP
#include <omp.h>
#endif
#include <string>
#include <set>
#include <sstream>

// Renderer configuration, holds algorithm, scene, and all other settings
struct Config
{
    enum Algorithm
    {
        kMultiViewPathTracing,
        kAlgorithmMax
    };

    static const char* GetName(Algorithm aAlgorithm)
    {
        static const char* algorithmNames[1] = {"multi-view path tracing"};
        if(aAlgorithm < 0 || aAlgorithm > 1)
            return "unknown algorithm";
        return algorithmNames[aAlgorithm];
    }

    static const char* GetAcronym(Algorithm aAlgorithm)
    {
        static const char* algorithmNames[1] = { "mvpt" };
        if(aAlgorithm < 0 || aAlgorithm > 1)
            return "unknown";
        return algorithmNames[aAlgorithm];
    }

    const Scene *mScene;
    Algorithm   mAlgorithm;
    int         mIterations;
    float       mMaxTime;
    float       mRadiusFactor;
    float       mRadiusAlpha;
    Framebuffer *mFramebuffer;
    std::vector<Framebuffer> mFramebuffers;
    int         mNumThreads;
    int         mBaseSeed;
    uint        mMaxPathLength;
    uint        mMinPathLength;
    std::string mOutputName;
    Vec2i       mResolution;
    bool        mReuse;
    int         mCameraNumber;
    int         mSceneID = 0;
    float       mPhongExp = -1;
};

// Utility function, essentially a renderer factory
AbstractRenderer* CreateRenderer(
    const Config& aConfig,
    const int     aSeed)
{
    const Scene& scene = *aConfig.mScene;

    switch(aConfig.mAlgorithm)
    {
    // Add multi-view support
    case Config::kMultiViewPathTracing:
        return new MultiViewPathTracer(scene, aSeed, aConfig.mReuse);
    default:
        printf("Unknown algorithm!!\n");
        exit(2);
    }
}

// Scene configurations
uint g_SceneConfigs[] = {
    Scene::kGlossyFloor | Scene::kBothSmallSpheres  | Scene::kLightSun,
    Scene::kGlossyFloor | Scene::kLargeMirrorSphere | Scene::kLightCeiling,
    Scene::kGlossyFloor | Scene::kBothSmallSpheres  | Scene::kLightPoint,
    Scene::kGlossyFloor | Scene::kBothSmallSpheres  | Scene::kLightBackground,
    Scene::kGlossyFloor | Scene::kLargeGlossySphere | Scene::kLightCeiling
};

std::string DefaultFilename(
    const uint              aSceneConfig,
    const Scene             &aScene,
    const Config::Algorithm aAlgorithm,
    bool aReuse)
{
    std::string filename;
    // if scene has glossy floor, name will be prefixed by g
    if((aSceneConfig & Scene::kGlossyFloor) != 0)
        filename = "g";

    // We use scene acronym
    filename += aScene.mSceneAcronym;

    // We add acronym of the used algorithm
    filename += "_";
    filename += Config::GetAcronym(aAlgorithm);

    // We add acronym of the used algorithm
    if( aAlgorithm == Config::kMultiViewPathTracing && aReuse )
        filename += "_r";

    // And it will be written as bmp
    filename += ".bmp";

    return filename;
}

// Utility function, gives length of array
template <typename T, size_t N>
inline int SizeOfArray( const T(&)[ N ] )
{
    return int(N);
}

void PrintHelp(const char *argv[])
{
    printf("\n");
    printf("Usage: %s [ -s <scene_id> | -a <algorithm> |\n", argv[0]);
    printf("           -t <time> | -i <iteration> | -o <output_name> | --report ]\n\n");
    printf("    -s  Selects the scene (default 0):\n");

    for(int i = 0; i < SizeOfArray(g_SceneConfigs); i++)
        printf("          %d    %s\n", i, Scene::GetSceneName(g_SceneConfigs[i]).c_str());

    printf("    -a  Selects the rendering algorithm (default vcm):\n");

    for(int i = 0; i < (int)Config::kAlgorithmMax; i++)
        printf("          %-3s  %s\n",
            Config::GetAcronym(Config::Algorithm(i)),
            Config::GetName(Config::Algorithm(i)));

    printf("    -t  Number of seconds to run the algorithm\n");
    printf("    -i  Number of iterations to run the algorithm (default 1)\n");
    printf("    -o  User specified output name, with extension .bmp or .hdr (default .bmp)\n");
    printf("    -r  Activate path reusing for MVPT\n");
    printf("    --report\n");
    printf("        Renders all scenes using all algorithms and generates an index.html file\n");
    printf("        that displays all images. Obeys the -t and -i options, ignores the rest.\n");
    printf("        Recommended usage: --report -i 1   (fastest preview)\n");
    printf("        Recommended usage: --report -t 10  (takes 5.5 mins)\n");
    printf("        Recommended usage: --report -t 60  (takes 30 mins)\n");
    printf("\n    Note: Time (-t) takes precedence over iterations (-i) if both are defined\n");
}

// Parses command line, setting up config
void ParseCommandline(int argc, const char *argv[], Config &oConfig)
{
    // Parameters marked with [cmd] can be change from command line
    oConfig.mScene         = NULL;                          // [cmd] When NULL, renderer will not run
    oConfig.mAlgorithm     = Config::kMultiViewPathTracing; // [cmd]
    oConfig.mIterations    = 1;                             // [cmd]
    oConfig.mMaxTime       = -1.f;                          // [cmd]
    oConfig.mOutputName    = "";                            // [cmd]
    oConfig.mNumThreads    = 0;
    oConfig.mBaseSeed      = 1234;
    oConfig.mMaxPathLength = 30;
    oConfig.mMinPathLength = 0;
    oConfig.mResolution    = Vec2i(1024, 1024);
    oConfig.mRadiusFactor  = 0.003f;
    oConfig.mRadiusAlpha   = 0.75f;
    oConfig.mReuse         = false;
    oConfig.mCameraNumber  = 9;

    // Load arguments
    for(int i=1; i<argc; i++)
    {
        std::string arg(argv[i]);

        // print help string (at any position)
        if(arg == "-h" || arg == "--help" || arg == "/?")
        {
            PrintHelp(argv);
            return;
        }

        if(arg[0] != '-') // all our commands start with -
        {
            continue;
        }
        else if(arg == "-s") // scene to load
        {
            if(++i == argc)
            {
                printf("Missing <SceneID> argument, please see help (-h)\n");
                return;
            }

            std::istringstream iss(argv[i]);
            iss >> oConfig.mSceneID;

            if(iss.fail() || oConfig.mSceneID < 0 || oConfig.mSceneID >= SizeOfArray(g_SceneConfigs))
            {
                printf("Invalid <SceneID> argument, please see help (-h)\n");
                return;
            }

            if(oConfig.mSceneID == 4)
            {
                if(++i == argc)
                {
                    printf("Missing <SceneID> argument, please see help (-h)\n");
                    return;
                }

                iss = std::istringstream(argv[i]);
                iss >> oConfig.mPhongExp;
            }
        }
        else if(arg == "-i") // number of iterations to run
        {
            if(++i == argc)
            {
                printf("Missing <iteration> argument, please see help (-h)\n");
                return;
            }

            std::istringstream iss(argv[i]);
            iss >> oConfig.mIterations;

            if(iss.fail() || oConfig.mIterations < 1)
            {
                printf("Invalid <iteration> argument, please see help (-h)\n");
                return;
            }
        }
        else if(arg == "-t") // number of seconds to run
        {
            if(++i == argc)
            {
                printf("Missing <time> argument, please see help (-h)\n");
                return;
            }

            std::istringstream iss(argv[i]);
            iss >> oConfig.mMaxTime;

            if(iss.fail() || oConfig.mMaxTime < 0)
            {
                printf("Invalid <time> argument, please see help (-h)\n");
                return;
            }

            oConfig.mIterations = -1; // time has precedence
        }
        else if(arg == "-o") // output filename
        {
            if(++i == argc)
            {
                printf("Missing <output_name> argument, please see help (-h)\n");
                return;
            }

            oConfig.mOutputName = argv[i];

            if(oConfig.mOutputName.length() == 0)
            {
                printf("Invalid <output_name> argument, please see help (-h)\n");
                return;
            }
        }
        else if(arg == "-r") // output filename
        {
            oConfig.mReuse = true;
        }
    }

    // Load scene
    Scene *scene = new Scene;
    scene->LoadCornellBox(oConfig.mResolution, g_SceneConfigs[oConfig.mSceneID], oConfig.mPhongExp);

    // Load multiple cameras and multiple framebuffers for multi-view
    if(oConfig.mAlgorithm == Config::kMultiViewPathTracing)
    {
        scene->LoadMultiView(oConfig.mResolution, oConfig.mCameraNumber);
        oConfig.mFramebuffers.resize(oConfig.mCameraNumber);
    }

    scene->BuildSceneSphere();

    oConfig.mScene = scene;

    // If no output name is chosen, create a default one
    if(oConfig.mOutputName.length() == 0)
    {
        oConfig.mOutputName = DefaultFilename(g_SceneConfigs[oConfig.mSceneID],
            *oConfig.mScene, oConfig.mAlgorithm, oConfig.mReuse);
    }

    // Check if output name has valid extension (.bmp or .hdr) and if not add .bmp
    std::string extension = "";

    if(oConfig.mOutputName.length() > 4) // must be at least 1 character before .bmp
        extension = oConfig.mOutputName.substr(
            oConfig.mOutputName.length() - 4, 4);

    if(extension != ".bmp" && extension != ".hdr")
        oConfig.mOutputName += ".bmp";
}

#endif  //__CONFIG_HXX__
