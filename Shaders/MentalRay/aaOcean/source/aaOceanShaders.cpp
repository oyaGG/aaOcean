//////////////////////////////////
//aaOcean Mental Ray Shaders, main compile file
/////////////////////////////////

#include <omp.h>
#include <shader.h>
#include <geoshader.h>

#include "timer\Timer.cpp"
#include "fftw3.h"
#include "aaOceanClass.cpp"
#include "shaders/oceanStore.h"
#include "shaders/shader_funcs.h"
#include "shaders/aaOceanData.cpp"
#include "shaders/aaOceanDisplace.cpp"
#include "shaders/aaOceanFoam.cpp"
#include "shaders/aaOceanImgVecDisplace.cpp"
//#include "shaders/aaOceanNormals.cpp"


