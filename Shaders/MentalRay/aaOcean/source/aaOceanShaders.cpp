// aaOcean Mental Ray Shaders
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

// aaOcean Mental Ray Shader, main compile file


#include <omp.h>
#include <shader.h>
#include <geoshader.h>

#include "timer/Timer.cpp"
#include "fftw3.h"
#include "aaOceanClass.cpp"
#include "shaders/oceanStore.h"
#include "shaders/shader_funcs.h"
#include "shaders/aaOceanData.cpp"

// The shaders below are now redundant
// #include "shaders/aaOceanImgVecDisplace.cpp"
// #include "shaders/aaOceanDisplace.cpp"
// #include "shaders/aaOceanFoam.cpp"



