#ifndef AAOCEAN_H
#define AAOCEAN_H

#include <float.h>
#include <cmath>
#include <omp.h>
#include <climits>
#include "fftw3.h"
#include "constants.h"
#include "functionLib.h"
#include "alignedMalloc.h"
#include "vectorSSE.h"
#include "agnerFog\sfmt.cpp" 
#include "agnerFog\stocc.h"
#include "agnerFog\stoc1.cpp"
#include "agnerFog\userintf.cpp"

class aaOcean
{ 
public:
	int		m_pointCount;
	int		m_resolution;
	ULONG	m_seed;
	int		m_windAlign;
	float	m_velocity;
	float	m_windDir;
	float	m_cutoff;
	float	m_damp;
	float	m_oceanScale;
	float	m_chopAmount;
	float	m_waveHeight;
	float	m_waveSpeed;
	float	m_time;
	float	m_3DGridULength;
	float	m_3DGridVLength;
	float	m_fmin, m_fmax; //for holding min/max foam
	char	m_state[256];

	//ocean array pointers
	int		*m_xCoord;
	int		*m_zCoord;
	float	*m_hokReal;
	float	*m_hokImag;
	float	*m_hktReal;
	float	*m_hktImag;
	float	*m_kX;
	float	*m_kZ;
	float	*m_omega;
	float	*m_rand1;
	float	*m_rand2;
	float	*m_eigenPlusX;
	float	*m_eigenPlusZ;
	float	*m_eigenMinusX;
	float	*m_eigenMinusZ;

	float m_randomType;

	//bool types for various checks during run-time
	bool	m_renderReady;
	bool	m_isAllocated;
	bool	m_isValid;
	bool	m_isNormalsAllocated;
	bool	m_isFoamAllocated;
	bool	m_isSplashAllocated;
	bool	m_redoHoK;
	bool	m_isShader;

	fftwf_complex *m_fft_htField;
	fftwf_complex *m_fft_chopX;
	fftwf_complex *m_fft_chopZ;
	fftwf_complex *m_fft_jxx;
	fftwf_complex *m_fft_jzz;
	fftwf_complex *m_fft_jxz;
	fftwf_complex *m_fft_normX;
	fftwf_complex *m_fft_normY;
	fftwf_complex *m_fft_normZ;

	fftwf_plan m_planHeightField;
	fftwf_plan m_planChopX;
	fftwf_plan m_planChopZ;
	fftwf_plan m_planJxx;
	fftwf_plan m_planJxz;
	fftwf_plan m_planJzz;

	aaOcean();
	aaOcean(const aaOcean &cpy);
	~aaOcean();

	void input(	int 	resolution,
				ULONG 	seed, 
				float 	oceanScale, 
				float 	velocity, 
				float 	cutoff, 
				float 	windDir, 
				int 	windAlign, 
				float 	damp,  
				float 	waveSpeed, 
				float 	waveHeight,
				float 	chopAmount,
				float 	time,
				float	randomType);

	void init();
	void allocateBaseArrays();
	void allocateFoamArrays();
	void allocateSplashArrays();
	void allocateNormalsArrays();
	void clearArrays();
	void clearResidualArrays();
	bool reInit(int data_size);
	ULONG get_uID(float, float);
	void setup_grid();
	void evaluateHokData();
	void evaluateHieghtField();
	void evaluateChopField();
	void evaluateJacobians();
	void evaluateNormalsFinDiff();
	void prepareOcean(bool doHeightField, bool doChopField, bool doJacobians, bool doNormals);
	void makeTileable(fftwf_complex *&fft_array);
};

#endif  /* AAOCEAN_H */

