// aaOcean v2.5 header file
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef AAOCEANCLASS_H
#define AAOCEANCLASS_H

#include "fftw3.h"

class aaOcean
{ 
public:
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
				bool	doFoam,
				bool    powTwoConversion);

	float getOceanData(float uCoord, float vCoord, int TYPE, bool rotateUV);
	void clearResidualArrays();

	bool isValid();
	bool isChoppy();
	int getResolution();
	char* getState();
	void getFoamBounds(float inBoundsMin, float inBoundsMax, float& outBoundsMin, float& outBoundsMax);

private:
	int		m_resolution;
	int		m_seed;
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
	float	m_foamBoundmin; //for holding min/max foam
	float   m_foamBoundmax; //for holding min/max foam
	
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

	//bool types for various checks during run-time
	bool	m_isAllocated;
	bool	m_isValid;
	bool	m_isFoamAllocated;
	bool	m_doHoK;
	bool	m_doSetup;
	bool	m_doChop;
	bool	m_doFoam;

	fftwf_complex *m_fft_htField;
	fftwf_complex *m_fft_chopX;
	fftwf_complex *m_fft_chopZ;
	fftwf_complex *m_fft_jxx;
	fftwf_complex *m_fft_jzz;
	fftwf_complex *m_fft_jxz;

	fftwf_plan m_planHeightField;
	fftwf_plan m_planChopX;
	fftwf_plan m_planChopZ;
	fftwf_plan m_planJxx;
	fftwf_plan m_planJxz;
	fftwf_plan m_planJzz;

	char m_state[256]; // array for holding the current state of aaOcean object
	
	// memory management functions
	void allocateBaseArrays();
	void allocateFoamArrays();
	void clearArrays();
	
	// initialization functions
	bool reInit(int data_size);
	ULONG generateUID(float, float);
	void prepareOcean();
	void setupGrid();

	// tessendorf ocean functions
	void evaluateHokData();
	void evaluateHieghtField();
	void evaluateChopField();
	void evaluateJacobians();
	
	// interpolation functions
	inline float catmullRom(float t, float a, float b, float c, float d);
	inline int wrap(int x, int n);

	void makeTileable(fftwf_complex *&fft_array);
};

#endif  /* AAOCEANCLASS_H */

