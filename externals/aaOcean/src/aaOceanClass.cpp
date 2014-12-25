// aaOcean v2.5
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef AAOCEANCLASS_CPP
#define AAOCEANCLASS_CPP

// Test for gcc
// gcc 4.2.0+ required for OpenMP
// 4.2.0 is ABI-compatible with 4.1.x.
// recommended compiler: gcc-4.2.4
#ifdef __GNUC__
    #if __GNUC_MINOR__ < 2 
        #warning [aaOcean] GNU C compiler 4.2.x, which is ABI-compatible with 4.1.0, is required
    #endif
#endif

#include <cmath>
#include <omp.h>
#include <climits>
#include <float.h>
#include <string.h>
#include "dSFMT/dSFMT.h"
#include "dSFMT/dSFMT.c"
#include "constants.h"
#include "functionLib.h"

#include "aaOceanClass.h"

aaOcean::aaOcean() :
    // input variables
    m_resolution(-1),
    m_seed(1),
    m_windAlign(0),         
    m_velocity(15.0f),
    m_windDir(45.0f),
    m_cutoff(-1.0f),
    m_damp(0.985f),
    m_oceanScale(100.0f),
    m_oceanDepth(10000.f),
    m_surfaceTension(0.0f),
    m_chopAmount(1.0f),
    m_waveHeight(1.0f),
    m_waveSpeed(1.0f),
    m_time(1.0f),
    m_loopTime(10000.0f),
    m_foamBoundmin(FLT_MAX),
    m_foamBoundmax(-FLT_MAX),
    m_randWeight(0.f),

    // working arrays
    m_xCoord(0),
    m_zCoord(0),
    m_hokReal(0),
    m_hokImag(0),
    m_hktReal(0),
    m_hktImag(0),
    m_kX(0),
    m_kZ(0),
    m_omega(0),
    m_rand1(0),
    m_rand2(0),

    // output arrays
    m_out_fft_htField(0),
    m_out_fft_chopX(0),
    m_out_fft_chopZ(0),
    m_out_fft_jxxX(0),
    m_out_fft_jxxZ(0),
    m_out_fft_jzzX(0),
    m_out_fft_jzzZ(0),
    m_out_fft_jxz(0),

    // bools to check ocean state
    m_isAllocated(0),
    m_isFoamAllocated(0),
    m_doHoK(0),
    m_doSetup(0),
    m_doChop(0),
    m_doFoam(0),

    // memory tracking
    m_memory(0),
    
    // fftw arrays
    m_fft_htField(0),
    m_fft_chopX(0),
    m_fft_chopZ(0),
    m_fft_jxx(0),
    m_fft_jzz(0),
    m_fft_jxz(0)
{
    strcpy (m_state, "");
}

aaOcean::aaOcean(const aaOcean &cpy)
{
    // empty copy constructor
    input(  cpy.m_resolution,
            cpy.m_seed,
            cpy.m_oceanScale,
            cpy.m_oceanDepth,
            cpy.m_surfaceTension,
            cpy.m_velocity,
            cpy.m_cutoff,
            cpy.m_windDir,
            cpy.m_windAlign,
            cpy.m_damp,
            cpy.m_waveSpeed,
            cpy.m_waveHeight,
            cpy.m_chopAmount,
            cpy.m_time,
            cpy.m_loopTime,
            cpy.m_doFoam
            );
}

aaOcean::~aaOcean()
{
    clearArrays();
}

bool aaOcean::isChoppy()
{
    return m_doChop;
}

char* aaOcean::getState()
{
    return &m_state[0];
}

int aaOcean::getResolution()
{
    return m_resolution;
}

void aaOcean::input(int resolution, unsigned int seed, float oceanScale, float oceanDepth, float surfaceTension, 
                    float velocity, float cutoff, float windDir, int windAlign, float damp, float waveSpeed, 
                    float waveHeight, float chopAmount, float time, float loopTime, bool doFoam, float randWeight)
{
    // forcing to be power of two, setting minimum resolution of 2^4
    resolution  = (int)pow(2.0f, (4 + abs(resolution))); 
    if(m_resolution != resolution)
    {
        m_resolution = resolution;
        allocateBaseArrays();
    }

    // scaled for better UI control
    m_waveHeight    = waveHeight * 0.01f;
    m_time          = time;
    m_waveSpeed     = waveSpeed;
    m_doFoam        = doFoam;

    if(chopAmount > 1.0e-6f)
    {
        // scaled for better UI control
        m_chopAmount = chopAmount * 0.01f;
        m_doChop = 1;
    }
    else
    {
        m_doChop = 0;
        m_chopAmount = 0.0f;
    }

    // clamping to minimum value
    oceanScale  = maximum<float>(oceanScale, 1.0e-6f);  
    velocity    = maximum<float>(velocity,   1.0e-6f); 
    oceanDepth  = maximum<float>(oceanDepth, 1.0e-6f); 
    // scaling for better UI control
    cutoff      = fabs(cutoff * 0.01f);
    // to radians
    windDir     = windDir * aa_PIBY180; 
    // forcing to even numbers
    windAlign   = maximum<int>(((windAlign + 1) * 2), 2); 
    // clamping to a maximum value of 1
    damp        = minimum<float>(damp, 1.0f); 

    if( m_oceanScale        != oceanScale       ||
        m_oceanDepth        != oceanDepth       ||
        m_surfaceTension    != surfaceTension   || 
        m_windDir           != windDir          ||
        m_cutoff            != cutoff           ||
        m_velocity          != velocity         ||
        m_windAlign         != windAlign        ||
        m_damp              != damp             ||
        m_loopTime          != loopTime)
    {
        m_oceanScale        = oceanScale;
        m_oceanDepth        = oceanDepth;
        m_surfaceTension    = surfaceTension;
        m_windDir           = windDir;
        m_cutoff            = cutoff;
        m_velocity          = velocity;
        m_windAlign         = windAlign;
        m_damp              = damp;
        m_loopTime          = loopTime;

        // re-evaluate HoK if any of these vars change
        m_doHoK = 1;
    }

    if(m_seed != seed || m_randWeight != randWeight)
    {
        m_seed = seed;
        m_randWeight = randWeight;
        m_doSetup = m_doHoK = 1;
    }

    if(!m_doHoK || !m_doSetup)
        sprintf(m_state,"[aaOcean Core] Ocean base state unchanged. Re-evaluating ocean with cached data");
    
    // we have our inputs. start preparing ocean arrays
    prepareOcean();
}

void aaOcean::prepareOcean()
{
    if(m_doSetup)
        setupGrid();
    
    if(m_doHoK)
        evaluateHokData();
    
    evaluateHieghtField();

    if(m_doChop)
        evaluateChopField();

    if(m_doFoam)
    {
        if(!m_isFoamAllocated)
            allocateFoamArrays();
        evaluateJacobians();
    }
    sprintf(m_state,"%s\n[aaOcean Core] Working memory allocated: %.2f megabytes", m_state, float(m_memory)/1048576.f);
}

void aaOcean::allocateBaseArrays()
{
    if(m_isAllocated) 
    {
        sprintf(m_state,"[aaOcean Core] Reallocating memory for ocean data structures for resolution %dx%d", m_resolution, m_resolution);
        clearArrays();
    }
    else
        sprintf(m_state,"[aaOcean Core] Allocating memory for ocean data structures for resolution %dx%d", m_resolution, m_resolution);

    int size = m_resolution * m_resolution;
    int dims[2] = {m_resolution, m_resolution};

    m_memory = size * sizeof(int) * 2 + size * sizeof(float) * 9 + size * sizeof(kiss_fft_cpx) * 3;

    m_fft_htField = (kiss_fft_cpx *) KISS_FFT_MALLOC(size * sizeof(kiss_fft_cpx));
    m_fft_chopX   = (kiss_fft_cpx *) KISS_FFT_MALLOC(size * sizeof(kiss_fft_cpx));
    m_fft_chopZ   = (kiss_fft_cpx *) KISS_FFT_MALLOC(size * sizeof(kiss_fft_cpx));

    m_xCoord    = (int*)   malloc(size * sizeof(int)); 
    m_zCoord    = (int*)   malloc(size * sizeof(int)); 
    
    m_hokReal   = (float*) malloc(size * sizeof(float)); 
    m_hokImag   = (float*) malloc(size * sizeof(float)); 
    m_hktReal   = (float*) malloc(size * sizeof(float)); 
    m_hktImag   = (float*) malloc(size * sizeof(float)); 
    m_kX        = (float*) malloc(size * sizeof(float)); 
    m_kZ        = (float*) malloc(size * sizeof(float)); 
    m_omega     = (float*) malloc(size * sizeof(float)); 
    m_rand1     = (float*) malloc(size * sizeof(float)); 
    m_rand2     = (float*) malloc(size * sizeof(float)); 

    m_out_fft_htField  = (float*) malloc(size * sizeof(float)); 
    m_out_fft_chopX    = (float*) malloc(size * sizeof(float)); 
    m_out_fft_chopZ    = (float*) malloc(size * sizeof(float)); 

    m_planHeightField   = kiss_fftnd_alloc(dims, 2, 1, 0, 0);
    m_planChopX         = kiss_fftnd_alloc(dims, 2, 1, 0, 0);
    m_planChopZ         = kiss_fftnd_alloc(dims, 2, 1, 0, 0);

    m_arrayPointer[eHEIGHTFIELD] = m_out_fft_htField;
    m_arrayPointer[eCHOPX]       = m_out_fft_chopX;
    m_arrayPointer[eCHOPZ]       = m_out_fft_chopZ;

    m_doHoK         = 1;
    m_doSetup       = 1;
    m_isAllocated   = 1;
}

void aaOcean::allocateFoamArrays()
{
    int size = m_resolution * m_resolution;
    int dims[2] = {m_resolution, m_resolution};
    m_memory += size * sizeof(kiss_fft_cpx) * 3;

    m_fft_jxx   = (kiss_fft_cpx *) KISS_FFT_MALLOC(size * sizeof(kiss_fft_cpx));
    m_fft_jzz   = (kiss_fft_cpx *) KISS_FFT_MALLOC(size * sizeof(kiss_fft_cpx));
    m_fft_jxz   = (kiss_fft_cpx *) KISS_FFT_MALLOC(size * sizeof(kiss_fft_cpx));

    m_out_fft_jxxX = (float*) malloc(size * sizeof(float)); 
    m_out_fft_jxxZ = (float*) malloc(size * sizeof(float)); 
    m_out_fft_jzzX = (float*) malloc(size * sizeof(float)); 
    m_out_fft_jzzZ = (float*) malloc(size * sizeof(float)); 
    m_out_fft_jxz  = (float*) malloc(size * sizeof(float)); 

    m_planJxx = kiss_fftnd_alloc(dims, 2, 1, 0, 0);
    m_planJxz = kiss_fftnd_alloc(dims, 2, 1, 0, 0);
    m_planJzz = kiss_fftnd_alloc(dims, 2, 1, 0, 0);

    m_arrayPointer[eFOAM]        = m_out_fft_jxz;
    m_arrayPointer[eEIGENPLUSX]  = m_out_fft_jxxX;
    m_arrayPointer[eEIGENPLUSZ]  = m_out_fft_jxxZ;
    m_arrayPointer[eEIGENMINUSX] = m_out_fft_jzzX;
    m_arrayPointer[eEIGENMINUSZ] = m_out_fft_jzzZ;

    m_isFoamAllocated = 1;
}


void aaOcean::clearResidualArrays()
{
    bool isResidualAllocated = 1;

    if(m_rand2)
    {
        free(m_rand2); 
        m_rand2 = 0;
    }
    else
        isResidualAllocated = 0;

    if(m_rand1)
    {
        free(m_rand1); 
        m_rand1 = 0;
    }
    if(m_omega)
    {
        free(m_omega); 
        m_omega = 0;
    }
    if(m_kZ)
    {
        free(m_kZ); 
        m_kZ = 0;
    }
    if(m_kX)
    {
        free(m_kX); 
        m_kX = 0;
    }
    if(m_hktImag)
    {
        free(m_hktImag); 
        m_hktImag = 0;
    }
    if(m_hktReal)
    {
        free(m_hktReal); 
        m_hktReal = 0;
    }
    if(m_hokImag)
    {
        free(m_hokImag); 
        m_hokImag = 0;
    }
    if(m_hokReal)
    {
        free(m_hokReal); 
        m_hokReal = 0;
    }
    if(m_zCoord)
    {
        free(m_zCoord); 
        m_zCoord = 0;
    }
    if(m_xCoord)
    {
        free(m_xCoord); 
        m_xCoord = 0;
    }

    if(isResidualAllocated)
    {
        int size = m_resolution * m_resolution;
        float cleared_memory = float(m_memory);
        m_memory = m_memory - (size * sizeof(float) * 9 + size * sizeof(int) * 2);
        cleared_memory = (cleared_memory - float(m_memory)) / 1048576.f;
        float working_memory = (float)m_memory/1048576.f;
        sprintf(m_state,"%s\n[aaOcean Core] Clearing %.2f megabytes of working memory. Current usage %.2f megabytes", \
            m_state, cleared_memory, working_memory);
    }
}

void aaOcean::clearArrays()
{
    if(m_isAllocated)
    {
        if(m_isFoamAllocated)
        {
            if(m_fft_jxx)
            {
                KISS_FFT_FREE(m_fft_jxx); 
                free(m_out_fft_jxxX);
                free(m_out_fft_jxxZ);
                free(m_planJxx);
                m_out_fft_jxxZ = m_out_fft_jxxX = 0;
                m_fft_jxx = 0;
            }
            if(m_fft_jzz)
            {
                KISS_FFT_FREE(m_fft_jzz);  
                free(m_out_fft_jzzX);
                free(m_out_fft_jzzZ);
                free(m_planJzz);
                m_out_fft_jzzX = m_out_fft_jzzZ = 0;
                m_fft_jzz = 0;
            }
            if(m_fft_jxz)
            {
                KISS_FFT_FREE(m_fft_jxz); 
                free(m_out_fft_jxz);
                free(m_planJxz);
                m_out_fft_jxz = 0;
                m_fft_jxz = 0;
            }
            m_isFoamAllocated = 0;
        }
        if(m_fft_chopZ)
        {
            KISS_FFT_FREE(m_fft_chopZ); 
            free(m_out_fft_chopZ);
            free(m_planChopZ);
            m_out_fft_chopZ = 0;
            m_fft_chopZ = 0;
        }
        if(m_fft_chopX)
        {
            KISS_FFT_FREE(m_fft_chopX); 
            free(m_out_fft_chopX);
            free(m_planChopX);
            m_out_fft_chopX = 0;
            m_fft_chopX = 0;
        }
        if(m_fft_htField)
        {
            KISS_FFT_FREE(m_fft_htField); 
            free(m_out_fft_htField);
            free(m_planHeightField);
            m_out_fft_htField = 0;
            m_fft_htField = 0;
        }
        m_isAllocated = 0;
    }
    clearResidualArrays();
}

unsigned int aaOcean::generateUID(const float xCoord, const float zCoord) const
{
    // a very simple hash function. should probably do a better one at some point
    register float angle;
    register float length;
    register float coordSq;
    register float id_out;
    unsigned int returnVal = 1;

    if (zCoord != 0.0f && xCoord != 0.0f)
    {
        coordSq = zCoord * zCoord + xCoord * xCoord;
        length = sqrt(coordSq);
        angle = xCoord / length; 
        angle = acos(angle);
        angle = RadsToDegs(angle);
        
        if (angle > 180.0f)
            angle = 360.0f - angle;
        
        id_out = coordSq + (length * angle) + 0.5f;
        
        if(angle == 0.0f)
            returnVal = (unsigned int)coordSq;
        else if (zCoord <= 0.0f)
            returnVal = (unsigned int)floor(id_out);
        else
            returnVal = INT_MAX - (unsigned int)floor(id_out) ;
    }
    return returnVal;
}

void aaOcean::setupGrid()
{
    const int n = m_resolution;
    const int half_n = (-n / 2) - ((n-1) / 2);
    unsigned int index, uID; 

    #pragma omp parallel for private(index, uID)
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            index = (i * n) + j;

            m_xCoord[index] = half_n + i * 2 ;
            m_zCoord[index] = half_n + j * 2 ;
            float x = m_xCoord[index];
            float z = m_zCoord[index];
            uID = (unsigned int)generateUID(x, z);

            // generate random numbers
            dsfmt_t dsfmt1;
            dsfmt_t dsfmt2;
            dsfmt_init_gen_rand(&dsfmt1, uID + m_seed);
            dsfmt_init_gen_rand(&dsfmt2, uID + m_seed + 1);
            
            float g1 = (float)gaussian(dsfmt1);
            float g2 = (float)gaussian(dsfmt2);

            float u1 = (float)uniform(dsfmt1);
            float u2 = (float)uniform(dsfmt2);

            m_rand1[index] = lerp(m_randWeight, g1, u1);
            m_rand2[index] = lerp(m_randWeight, g2, u2);
        }
    }
    m_doSetup = 0;
}

 void aaOcean::evaluateHokData()
{
    float k_sq, k_mag, k_dot_w, philips;

    const int      n        = m_resolution * m_resolution;
    const float    k_mult   = aa_TWOPI / m_oceanScale;
    const float    L        = m_velocity * m_velocity / aa_GRAVITY;
    const float    L_sq     = L * L;
    const float    windx    = cos(m_windDir);
    const float    windz    = sin(m_windDir);
    const float    omega0   = aa_TWOPI / m_loopTime;
    
    bool bDamp  = 0;
    if (m_damp > 0.0f)
        bDamp = 1;

    #pragma omp parallel for private( k_sq, k_mag, k_dot_w, philips)  
    for(int index = 0; index < n; ++index)
    {
        // build Kx and Kz
        m_kX[index] =  (float)m_xCoord[index] * k_mult; 
        m_kZ[index] =  (float)m_zCoord[index] * k_mult;

        // philips spectrum vars
        k_sq        = (m_kX[index] * m_kX[index]) + (m_kZ[index] * m_kZ[index]);
        k_mag       = 1.0f / sqrt( k_sq );
        k_dot_w     = (m_kX[index] * k_mag) * windx + (m_kZ[index] * k_mag) * windz;

        // calculate philips spectrum
        philips     = sqrt((( exp(-1.0f / ( L_sq * k_sq)) * pow(k_dot_w, m_windAlign)) / 
                      (k_sq * k_sq)) * exp(-k_sq * m_cutoff));

        // reduce reflected waves
        if(bDamp)
        {
            if(k_dot_w < 0.0f)
                philips *= (1.0f - m_damp);
        }       

        // build dispersion relationship with oceanDepth relationship
        m_omega[index]   = (aa_GRAVITY / k_mag) * tanh( sqrt( k_sq ) * m_oceanDepth);

        // modifying dispersion for capillary waves
        m_omega[index] = m_omega[index] * (1.0f + k_sq * m_surfaceTension * m_surfaceTension);

        m_omega[index] = sqrt(m_omega[index]);

        // add time looping support with OmegaNought
        m_omega[index]   = (int(m_omega[index] / omega0)) * omega0;

        m_hokReal[index] = (aa_INV_SQRTTWO) * (m_rand1[index]) * philips;
        m_hokImag[index] = (aa_INV_SQRTTWO) * (m_rand2[index]) * philips;
    }

    sprintf(m_state,"\n[aaOcean Core] Finished initializing all ocean data");
    m_doHoK = 0;
}

 void aaOcean::evaluateHieghtField()
{
    int  i,j,index, index_rev;
    register float hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt;

    const float wt  = m_waveSpeed * m_time;
    const int n = m_resolution;
    const int nn = n * n;
    register const int n_sq = n * n - 1;
    const float signs[2] = { 1.0f, -1.0f };

    #pragma omp parallel for private(index, index_rev, hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt)  
    for(index = 0; index < nn; ++index)
    {
        index_rev   = n_sq - index; //tail end 
        hokReal     = m_hokReal[index];
        hokImag     = m_hokImag[index];
        hokRealOpp  = m_hokReal[index_rev];
        hokImagOpp  = m_hokImag[index_rev];

        coswt = cos( m_omega[index] * wt);
        sinwt = sin( m_omega[index] * wt);

        m_hktReal[index]  = ( hokReal    *  coswt )  +  ( hokImag    *  sinwt )  + 
                            ( hokRealOpp *  coswt )  -  ( hokImagOpp *  sinwt )  ;  //complex conjugage
                
        m_hktImag[index]  = (-hokReal    *  sinwt )  +  ( hokImag    *  coswt )  +
                            ( hokRealOpp *  sinwt )  +  ( hokImagOpp *  coswt )  ;  //complex conjugage
        
        m_fft_htField[index].r = m_hktReal[index];
        m_fft_htField[index].i = m_hktImag[index];
    }

    kiss_fftnd(m_planHeightField, m_fft_htField, m_fft_htField);

    for(i = 0; i < n; ++i)
    {
        #pragma omp parallel for private(j)
        for(j = 0; j < n; ++j)
            m_out_fft_htField[(i*n) + j] = m_fft_htField[(i*n) + j].r * signs[(i + j) & 1]  * m_waveHeight;
    }
}

 void aaOcean::evaluateChopField()
{
    int  i, j, index;
    register float  kX, kZ, kMag;
    int n = m_resolution * m_resolution;
    const float signs[2] = { 1.0f, -1.0f };
    float multiplier;

    #pragma omp parallel for private( index,  kX,  kZ, kMag)  
    for(index = 0; index < n; ++index)
    {           
        kMag = sqrt(m_kX[index] * m_kX[index] + m_kZ[index] * m_kZ[index]);
        kX = m_kX[index] / kMag;
        kZ = m_kZ[index] / kMag;
        
        m_fft_chopX[index].r =  m_hktImag[index] * kX;
        m_fft_chopX[index].i = -m_hktReal[index] * kX;

        m_fft_chopZ[index].r =  m_hktImag[index] * kZ;
        m_fft_chopZ[index].i = -m_hktReal[index] * kZ;
    }

    kiss_fftnd(m_planChopX, m_fft_chopX, m_fft_chopX);
    kiss_fftnd(m_planChopZ, m_fft_chopZ, m_fft_chopZ);

    n = m_resolution;
    for(i = 0; i < n; ++i)
    {
        #pragma omp parallel for private(multiplier, j, index)  
        for(j = 0; j < n; ++j)
        {
            index = (i*n) + j;
            multiplier = m_chopAmount * signs[(i + j) & 1] * -1.0f;
            m_out_fft_chopX[index] = m_fft_chopX[index].r * multiplier;
            m_out_fft_chopZ[index] = m_fft_chopZ[index].r * multiplier;
        }
    }
}

void aaOcean::evaluateJacobians()
{
    int  i, j, index;
    register float kX, kZ, kMag, kXZ, multiplier;
    const float signs[2] = { 1.0f, -1.0f };
    int n = m_resolution * m_resolution;

    #pragma omp parallel for private( index, kX, kZ, kXZ, kMag) 
    for(index = 0; index < n; ++index)
    {           
        kMag = 1.0f / sqrt(m_kX[index] * m_kX[index] + m_kZ[index] * m_kZ[index]);
        kX   = (m_kX[index] * m_kX[index]) * kMag;
        kZ   = (m_kZ[index] * m_kZ[index]) * kMag;
        kXZ  = (m_kX[index] * m_kZ[index]) * kMag;

        m_fft_jxx[index].r =  m_hktReal[index] * kX;
        m_fft_jxx[index].i =  m_hktImag[index] * kX;

        m_fft_jzz[index].r =  m_hktReal[index] * kZ;
        m_fft_jzz[index].i =  m_hktImag[index] * kZ;

        m_fft_jxz[index].r =  m_hktReal[index] * kXZ;
        m_fft_jxz[index].i =  m_hktImag[index] * kXZ;
    }

    kiss_fftnd(m_planJxx, m_fft_jxx, m_fft_jxx);
    kiss_fftnd(m_planJzz, m_fft_jzz, m_fft_jzz);
    kiss_fftnd(m_planJxz, m_fft_jxz, m_fft_jxz);

    n = m_resolution;
    for(i = 0; i < n; ++i)
    {
        #pragma omp parallel for private(multiplier, j, index)  
        for(j = 0; j < n; ++j)
        {
            index = (i*n) + j;
            multiplier = -m_chopAmount * signs[(i + j) & 1];
            m_fft_jxx[index].r = m_fft_jxx[index].r * multiplier + 1.0f;
            m_fft_jzz[index].r = m_fft_jzz[index].r * multiplier + 1.0f;
            m_fft_jxz[index].r = m_fft_jxz[index].r * multiplier;
        }
    }

    register float jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz, temp;
    #pragma omp parallel for private(index, jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz, temp)  
    for(index = 0; index < n*n; ++index)
    {
        Jxx = m_fft_jxx[index].r;
        Jzz = m_fft_jzz[index].r;
        Jxz = m_fft_jxz[index].r;

        temp = (0.5f * sqrt( ((Jxx - Jzz) * (Jxx - Jzz)) + 4.0f * (Jxz*Jxz)));
        jPlus   = (0.5f * (Jxx + Jzz))  +  temp;
        jMinus  = (0.5f * (Jxx + Jzz))  -  temp;
        qPlus   = (jPlus  - Jxx) / Jxz;
        qMinus  = (jMinus - Jxx) / Jxz;

        temp = sqrt( 1.0f + qPlus * qPlus);
        m_out_fft_jxxX[index] = 1.0f  / temp;
        m_out_fft_jxxZ[index] = qPlus / temp;

        temp = sqrt( 1.0f + qMinus * qMinus);
        m_out_fft_jzzX[index] = 1.0f   / temp;
        m_out_fft_jzzZ[index] = qMinus / temp;

        //store foam back in this array for convenience
        m_out_fft_jxz[index] = (Jxx * Jzz) - (Jxz * Jxz); 
    }
}

void aaOcean::getFoamBounds(float& outBoundsMin, float& outBoundsMax) const
{
    outBoundsMax = -FLT_MAX;
    outBoundsMin =  FLT_MAX;

    int index, n;
    n = m_resolution * m_resolution;
    for(index = 0; index < n; index++)
    {
        if(outBoundsMax < m_out_fft_jxz[index])
            outBoundsMax = m_out_fft_jxz[index];

        if(outBoundsMin > m_out_fft_jxz[index]) 
            outBoundsMin = m_out_fft_jxz[index];
    }
}

void aaOcean::getOceanArray(float *&outArray, aaOcean::arrayType type) const
{
    // get the pointer to the aaOcean array that we want to pull data from
    float *arrayPointer = m_arrayPointer[type];

    const int arraySize = m_resolution * m_resolution;

    #pragma omp parallel for 
    for(int i = 0; i < arraySize; ++i)
            outArray[i] = arrayPointer[i];
}

float aaOcean::getOceanData(float uCoord, float vCoord, aaOcean::arrayType type) const
{
    float u, v, du, dv = 0;
    int xMinus1, yMinus1, x, y, xPlus1, yPlus1, xPlus2, yPlus2;

    // get the pointer to the aaOcean array that we want to pull data from
    float *arrayPointer = m_arrayPointer[type];

    // maya and softimage V axis runs along negative z axis
    // aaOcean uses convention of V axis along positive z axis
    vCoord = -vCoord;

    // begin UV coordinate wrapping to [0-1] interval
    uCoord = fmod(uCoord, 1.0f);
    vCoord = fmod(vCoord, 1.0f);
    if(uCoord < 0.0f)
        uCoord = 1.0f + uCoord;
    if(vCoord < 0.0f)
        vCoord =  1.0f + vCoord;

    // use UV coordinates to work out ocean array indeces
    u = uCoord * float(m_resolution);
    v = vCoord * float(m_resolution);
    x = (int)floor(u);
    y = (int)floor(v);
    du = u - x; 
    dv = v - y; 

    // prepare catmul-rom end points for interpolation
    // wrap any indices that are outside the array boundaries
    xMinus1 = wrap(x - 1) * m_resolution;
    xPlus1  = wrap(x + 1) * m_resolution;
    xPlus2  = wrap(x + 2) * m_resolution;
    x       = wrap(x)   * m_resolution;
    y       = wrap(y);
    yMinus1 = wrap(y - 1);
    yPlus1  = wrap(y + 1);
    yPlus2  = wrap(y + 2);
    
    // prepare for catmul-rom interpolation
    const float a1 = catmullRom(du, 
                            arrayPointer[xMinus1    + yMinus1],
                            arrayPointer[x          + yMinus1],
                            arrayPointer[xPlus1     + yMinus1],
                            arrayPointer[xPlus2     + yMinus1]);

    const float b1 = catmullRom(du, 
                            arrayPointer[xMinus1    +   y],
                            arrayPointer[x          +   y],
                            arrayPointer[xPlus1     +   y],
                            arrayPointer[xPlus2     +   y]);

    const float c1 = catmullRom(du, 
                            arrayPointer[xMinus1    + yPlus1],
                            arrayPointer[x          + yPlus1],
                            arrayPointer[xPlus1     + yPlus1],
                            arrayPointer[xPlus2     + yPlus1]);

    const float d1 = catmullRom(du, 
                            arrayPointer[xMinus1    + yPlus2],
                            arrayPointer[x          + yPlus2],
                            arrayPointer[xPlus1     + yPlus2],
                            arrayPointer[xPlus2     + yPlus2]);

    return catmullRom(dv, a1, b1, c1, d1);
}

inline float aaOcean::catmullRom(const float t, const float a, const float b, const float c, const float d) const
{
    return  0.5f * ( ( 2.0f * b ) + ( -a + c ) * t + 
            ( 2.0f * a - 5.0f * b + 4.0f * c - d ) * t * t + 
            ( -a + 3.0f * b - 3.0f * c + d )* t * t * t );
}

inline int aaOcean::wrap(int x) const
{
    // fast modulo for power of 2 numbers
    x = x & (m_resolution - 1);
    return x;
}

#endif  /* AAOCEANCLASS_CPP */
