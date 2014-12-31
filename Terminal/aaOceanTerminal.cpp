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
#include "aaOceanClass.cpp"
#include "openEXROutput.h"

#include "timer/Timer.cpp"

int main(int argc, char* argv[])
{
    char msg[512];
    input oceanInput;

    Timer t;
    t.start();

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
            oceanInput.waveSpeed,       // wave speed
            oceanInput.waveHeight,      // wave height
            oceanInput.waveChop,        // chop amount
            time,                       // time in seconds
            100000.f,                   // repeat/loop time
            TRUE);                       // calculate foam
        
        
        LOG(logDEBUG) << "Logging Ocean Core messages\n" << pOcean->m_state;
        LOG(logINFO) << msg;
        
        /*char outputFileName[512];
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
        LOG(logDEBUG) << msg;*/

        t.stop();
        sprintf(msg,"Elapsed time: %f secs", t.getElapsedTimeInSec());
        LOG(logDEBUG) << msg;


        currentFrame++;
    }

    delete pOcean;
    return 0;
}
