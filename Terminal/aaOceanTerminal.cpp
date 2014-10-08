// aaOceanTerminal.cpp : standalone aaOcean console/terminal application
// Outputs RGBA, with Vector Displacement in RGB, and foam in Alpha
// Recommend OpenEXR 2.xx to avoid namespace runtime conflicts
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

// usage:
// ./aaOcean --resolution=256 --seed=1 --oceanscale=100 --oceandepth=10000 
// --surfacetension=0.0 --velocity=15 --smooth=0 --winddir=45 --windalign=1 
// --reflectedwaves=0.3 --speed=1.0 --waveheight=1.0 --wavechop=1.0 -
// -startframe=101 --endframe=101 --outputfolder=/tmp --postfix=Your_Postfix

// not all of the above arguments need to be specified. See default values in 
// input.h in input constructor

#include <stdio.h>
#include <iostream>
#include <iomanip> 
#include <omp.h>
#include <string.h>
#include <math.h> 

#include "functionLib.h"
#include "input.h"
#include "timer/Timer.cpp"
#include "aaOceanClass.cpp"
#include "openEXROutput.h"

int main(int argc, char* argv[])
{
    Timer t;
    char msg[256];
    input oceanInput;

    if(!processInput(argc, argv, oceanInput))
        return 1;

    int dimension = oceanInput.resolution;
    LOG(logINFO) << "Starting ocean evaluation for " << dimension << "x" << dimension;

    // modify dimension to account of aaOcean's
    // arbitrary resolution scaling of 4
    dimension = (int)(log((float)dimension)/ log(2.0f)) - 4;
    aaOcean *pOcean = new aaOcean;
    
    float timestep = 1.0f/oceanInput.fps;
    int currentFrame = oceanInput.startFrame;

    while(currentFrame <= oceanInput.endFrame)
    {
        t.start();

        int absoluteFrame = currentFrame - oceanInput.startFrame;
        float time = float(absoluteFrame) * timestep;

        sprintf(msg,"Evaluating for time %0.3f seconds and frame %d", time, currentFrame);
        LOG(logINFO) << msg;

        pOcean->input(
            dimension,                  // resolution 
            oceanInput.seed,            // seed
            oceanInput.oceanScale,      // ocean scale
            oceanInput.oceanDepth,      // ocean depth
            oceanInput.surfaceTension,  // surface tension
            oceanInput.velocity,        // velocity
            oceanInput.smooth,          // cutoff/smooth
            oceanInput.windDir,         // wind dir
            oceanInput.windAlign,       // wind align
            oceanInput.reflectedWaves,  // damp
            oceanInput.waveSpeed,		// wave speed
            oceanInput.waveHeight,		// wave height
            oceanInput.waveChop,		// chop amount
            time,                       // time in seconds
            100000.f,                   // repeat/loop time
            TRUE,                       // calculate foam
            FALSE);                     // calculate normals
        
        t.stop();
        
        LOG(logDEBUG) << "Logging Ocean Core messages\n" << pOcean->m_state;
        sprintf(msg,"Time to create ocean grid: %0.2f seconds", t.getElapsedTimeInSec());
        LOG(logINFO) << msg;
        
        char outputFileName[512];
        oceanDataToEXR(pOcean, 
                       &oceanInput.outputFolder[0], 
                       &oceanInput.postfix[0], 
                       currentFrame, 
                       &outputFileName[0]);

        sprintf(msg,"Written OpenEXR object-space vector-displacement map");
        LOG(logDEBUG) << msg;
        sprintf(msg,"OpenEXR image location: %s", &outputFileName[0]);
        LOG(logINFO) << msg;
        sprintf(msg,"OpenEXR RGB contains position, Alpha contains raw foam/spray emission data");
        LOG(logDEBUG) << msg;

        currentFrame++;
    }

    delete pOcean;
    return 0;
}

void sressTest(aaOcean *pOcean)
{
    Timer t;
     // stress-test
    int steps = INT_MAX / 10;
    char msg[512];
    sprintf(msg,"Evaluating Ocean for sample UV's, total samples %d", steps);
    LOG(logINFO) << msg;
    
    t.start();
    #pragma omp parallel for 
    for(int i = 0; i < steps; ++i)
    {
        // generate random UVs between -5.0 to 5.0
        float u = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * 10.f) - 5.f; 
        float v = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) * 10.f) - 5.f; 

        float y = pOcean->getOceanData(u, v, aaOcean::eHEIGHTFIELD);
        float x = pOcean->getOceanData(u, v, aaOcean::eCHOPX);
        float z = pOcean->getOceanData(u, v, aaOcean::eCHOPZ);
    }
    t.stop();
    sprintf(msg,"Time to randomly sample ocean grid: %0.2f seconds", t.getElapsedTimeInSec());
    LOG(logINFO) << msg;
}

