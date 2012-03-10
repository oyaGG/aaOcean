#ifndef AAOCEANCLASS_H
#define AAOCEANCLASS_H

//////////////////////////////////
//aaOcean custom includes and defines
//////////////////////////////////
#include "fftw3.h"
#include <float.h>
#include "3rdParty/RW_MersenneTwister.h"
#include "3rdParty/Timer.cpp"
#include "3rdParty/vectorClass.h"
#include <omp.h>
#include <stdio.h>
/////////////////////////////////

inline double isEven(int);
inline long int_sqrt(long);
inline int wrap(int ,int );
inline bool isEqual(double x, double y, const double epsilon);

class aaOceanClass
{ 
public:
	//ocean variables
	int	   resolution;
	int	   normalsResolution;
	int	   windAlign;
	double velocity;
	double windDir;
	double cutoff;
	double damp;
	double oceanScale;
	double chopAmount;
	double waveHeight;
	double waveSpeed;

	//foam array bounds
	float foamMin;
	float foamMax;

	//Softimage time information
	double time;
	float  simStep;
	int	   frame;
	
	//misc variables
	int	   pointCount;
	int	   seed;
	int	   shaderResolution;

	//ocean array pointers
	int *xCoord;
	int *zCoord;
	double *hokReal;
	double *hokImag;
	double *hktReal;
	double *hktImag;
	double *kX;
	double *kZ;
	double *omega;
	double *rand1;
	double *rand2;
	double *eigenPlusX;
	double *eigenPlusZ;
	double *eigenMinusX;
	double *eigenMinusZ;

	//bool types for various checks during run-time
	bool renderReady;
	bool isAllocated;
	bool isValid;
	bool isNormalsAllocated;
	bool isFoamAllocated;
	bool isSplashAllocated;
	bool redoHoK;
	bool isChoppy;	
	bool isArnoldInit;

	AtCritSec   critical_section;
	char name[50];

	fftw_complex *fft_htField;
	fftw_complex *fft_chopX;
	fftw_complex *fft_chopZ;

	fftw_complex *fft_jxx;
	fftw_complex *fft_jzz;
	fftw_complex *fft_jxz;

	fftw_complex *fft_normX;
	fftw_complex *fft_normY;
	fftw_complex *fft_normZ;

	fftw_plan planHeightField;
	fftw_plan planChopX;
	fftw_plan planChopZ;

	fftw_plan planJxx;
	fftw_plan planJxz;
	fftw_plan planJzz;

	aaOceanClass();

	~aaOceanClass()
	{
		clearArrays();
	}

	inline void init();
	inline void allocateBaseArrays();
	inline void allocateFoamArrays();
	inline void allocateSplashArrays();
	inline void allocateNormalsArrays();
	inline void clearArrays();

	inline void copyInput(aaOceanClass *&input);
	inline bool reInit(int data_size);

	inline ULONG get_uID(double, double);
	inline void setup_grid();
	inline void evaluateHokData();
	inline void evaluateHieghtField();
	inline void evaluateChopField();
	inline void evaluateJacobians();
	inline void evaluateNormalsFinDiff();
	inline void prepareOcean(bool doHoK, bool doHeightField, bool doJacobians, bool doNormals);
	inline void makeTileable(fftw_complex *&fft_array);

};

inline aaOceanClass::aaOceanClass()
{
	normalsResolution = frame = shaderResolution = pointCount = resolution = windAlign = seed = -1;
	velocity = windDir = cutoff = damp = oceanScale = waveHeight = time = -1.0;
	isAllocated = isValid = isNormalsAllocated = isFoamAllocated = isChoppy = isArnoldInit = 0;
	isSplashAllocated = redoHoK = 0;
	sprintf(name,"");
	foamMin = foamMax = 0.0f;

	xCoord = zCoord = 0;
	hokReal=0;
	hokImag=0;
	hktReal=0;
	hktImag=0;
	kX=0;
	kZ=0;
	omega=0;
	rand1=0;
	rand2=0;
	eigenPlusX=0;
	eigenPlusZ=0;
	eigenMinusX=0;
	eigenMinusZ=0;
}
inline bool aaOceanClass::reInit(int data_size)
{
	if(((data_size & (data_size - 1)) != 0) || data_size == 0) //	not power of 2
	{		
		AiMsgWarning("[aaOcean %s] invalid point resolution of %d. Please select a power-of-2 subdivision value",name,data_size );
		isValid = 0;
	}
	else
	{
		if(resolution != data_size || !isAllocated)
		{
			resolution = data_size;
			allocateBaseArrays();				
			redoHoK  = 1;
			setup_grid();
			char msg[200];
			sprintf(msg,"[aaOcean %s] Ocean data (%dx%d) allocated ",name,resolution,resolution);
			AiMsgWarning(msg);
		}
		else
			resolution =  data_size;
		isValid = 1;
	}
	return isValid;
}


inline void aaOceanClass::prepareOcean(bool doHoK, bool doHeightField, bool doJacobians, bool doNormals)
{
	if(doHoK || redoHoK )
	{
		evaluateHokData();
		redoHoK = 0;
	}

	if(doHeightField)
		evaluateHieghtField();

	if(isChoppy)
		evaluateChopField();

	if(doJacobians && isAllocated)
	{
		if(!isFoamAllocated)
			allocateFoamArrays();
		if(!isSplashAllocated)
			allocateSplashArrays();
		evaluateJacobians();
	}
	if(doNormals)
	{
		if(!isNormalsAllocated)
			allocateNormalsArrays();
		evaluateNormalsFinDiff();
	}
}

inline void aaOceanClass::allocateBaseArrays()
{
	if(isAllocated) 
		clearArrays();

	kX		= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	kZ		= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	omega	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	hokReal	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	hokImag	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	hktReal	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	hktImag	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	rand1	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	rand2	= (double*) _aligned_malloc(resolution * resolution * sizeof(double), 16);
	xCoord	= (int*)	_aligned_malloc(resolution * resolution * sizeof(int)	, 16);
	zCoord	= (int*)	_aligned_malloc(resolution * resolution * sizeof(int)	, 16);

	if(resolution >= 1024)
	{
		if(!fftw_init_threads())
			AiMsgWarning("[aaOcean %s] FAILED to allocate FFTW threads!",name);
		else
		{
			int num_procs = omp_get_num_procs();
			fftw_plan_with_nthreads(num_procs);
			AiMsgWarning("[aaOcean %s] enabling FFTW threading with %d threads",name, num_procs);
		}
	}

	fft_htField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((resolution+1)*(resolution+1)));
	fft_chopX	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((resolution+1)*(resolution+1)));
	fft_chopZ	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((resolution+1)*(resolution+1)));

	planHeightField = fftw_plan_dft_2d(resolution,resolution,fft_htField,fft_htField,1,FFTW_MEASURE);
	planChopX		= fftw_plan_dft_2d(resolution,resolution,fft_chopX	,fft_chopX	,1,FFTW_MEASURE);
	planChopZ		= fftw_plan_dft_2d(resolution,resolution,fft_chopZ	,fft_chopZ	,1,FFTW_MEASURE);
	isAllocated = 1;
	
	
}

inline void aaOceanClass::allocateNormalsArrays()
{
	int arSize = normalsResolution + 1;
	fft_normX	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((arSize)*(arSize)));
	fft_normY	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((arSize)*(arSize)));
	fft_normZ	= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((arSize)*(arSize)));
	isNormalsAllocated = 1;
}

inline void aaOceanClass::allocateFoamArrays()
{
	if(resolution >= 1024)
		fftw_plan_with_nthreads(omp_get_num_procs());

	fft_jxx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((resolution+1)*(resolution+1)));
	fft_jxz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((resolution+1)*(resolution+1)));
	fft_jzz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ((resolution+1)*(resolution+1)));
	planJxx = fftw_plan_dft_2d(resolution,resolution,fft_jxx,fft_jxx,1,FFTW_MEASURE);
	planJxz = fftw_plan_dft_2d(resolution,resolution,fft_jxz,fft_jxz,1,FFTW_MEASURE);
	planJzz = fftw_plan_dft_2d(resolution,resolution,fft_jzz,fft_jzz,1,FFTW_MEASURE);
	isFoamAllocated = 1;
}
inline void aaOceanClass::allocateSplashArrays()
{
	eigenPlusX	= (double*) _aligned_malloc((resolution+1)*(resolution+1) * sizeof(double), 16);
	eigenPlusZ	= (double*) _aligned_malloc((resolution+1)*(resolution+1) * sizeof(double), 16);
	eigenMinusX	= (double*) _aligned_malloc((resolution+1)*(resolution+1) * sizeof(double), 16);
	eigenMinusZ	= (double*) _aligned_malloc((resolution+1)*(resolution+1) * sizeof(double), 16);
	isSplashAllocated = 1;
}

inline void aaOceanClass::clearArrays()
{
	if(isAllocated)
	{
		_aligned_free(kX);
		_aligned_free(kZ);
		_aligned_free(omega);
		_aligned_free(hokReal);
		_aligned_free(hokImag);
		_aligned_free(hktReal);
		_aligned_free(hktImag);
		_aligned_free(rand1);
		_aligned_free(rand2);
		_aligned_free(xCoord);
		_aligned_free(zCoord);
		fftw_free(fft_htField);
		fftw_destroy_plan(planHeightField);
		fftw_free(fft_chopX);
		fftw_destroy_plan(planChopX);
		fftw_free(fft_chopZ);
		fftw_destroy_plan(planChopZ);
		isAllocated = 0;
	}
	if(isFoamAllocated)
	{
		fftw_free(fft_jxx);
		fftw_destroy_plan(planJxx);
		fftw_free(fft_jxz);
		fftw_destroy_plan(planJxz);
		fftw_free(fft_jzz);
		fftw_destroy_plan(planJzz);
		isFoamAllocated = 0;
	}
	if(isSplashAllocated)
	{
		_aligned_free(eigenPlusX);
		_aligned_free(eigenPlusZ);
		_aligned_free(eigenMinusX);
		_aligned_free(eigenMinusZ);
		isSplashAllocated = 0;
	}
	if(isNormalsAllocated)
	{
		fftw_free(fft_normX);
		fftw_free(fft_normY);
		fftw_free(fft_normZ);
		isNormalsAllocated = 0;
	}
	fftw_cleanup_threads();
	fftw_cleanup();
	
}

inline ULONG aaOceanClass::get_uID(double xCoord, double zCoord)
{
	double angle=0;
	double length=0;
	double id_out;
	ULONG returnVal=0;

	if (zCoord ==0 && xCoord == 0)
		return 1;
	else
	{
		angle = xCoord/sqrt(zCoord*zCoord + xCoord*xCoord); 
		angle = acos(angle);
		angle = ( angle * (180/3.141592653589793238462643383279502884) ); 
		
		if (angle > 180)
			angle = 360 - angle;
		
		length = sqrt(zCoord*zCoord + xCoord*xCoord); 
		id_out = (length * length)  +  (length * angle);
		
		if( angle == 0)
			returnVal = (ULONG)length * (ULONG)length;
		else if (zCoord <= 0)
			returnVal = (ULONG)floor(id_out + 0.5f);
		else
			returnVal = ULONG_MAX - (ULONG)floor(id_out + 0.5f) ;
	}
	return returnVal;
}

inline void aaOceanClass::setup_grid()
{
	if(!isAllocated)
		return;
	const int n = resolution;
	int index; ULONG uID;

	#pragma omp parallel for private(index, uID)
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			index = (i*n) + j;

			xCoord[index] = ((-n/2)+i*2-(n-1)/2);
			zCoord[index] = ((-n/2)+j*2-(n-1)/2);

			uID = get_uID(xCoord[index], zCoord[index]);

			MTRand randNum(uID + seed); 
			rand1[index] = randNum.randNorm(); 
			rand2[index] = randNum.randNorm();
		}
	}
}
inline void aaOceanClass::evaluateHokData()
{
	double k_sq, k_mag,  k_dot_w, low_freq, exp_term, philips;
	const int	n			 = resolution * resolution;
	const double	k_mult	 = (2 * 3.14159265358979323846264338327950288) / oceanScale;
	const double	L		 = velocity;
	const double	L_sq	 = L*L;
	const double	windx	 = cos(windDir);
	const double	windz	 = sin(windDir);
	const double    invSqRt2 = 0.70710678118654752440084436210485;
	const double	gravity	 = 9.81;

	bool bDamp	= 0;
	if (damp > 0.0f)
		bDamp = 1;

	#pragma omp parallel for private( k_sq, k_mag, k_dot_w, low_freq, exp_term, philips) 
	for(int index = 0; index < n; index++)
	{
		kX[index] =  xCoord[index] * k_mult; 
		kZ[index] =  zCoord[index] * k_mult;
		//philips spectrum vars
		k_sq		= (kX[index] * kX[index]) + (kZ[index] * kZ[index]);
		k_mag		= sqrt( k_sq );
		k_dot_w		= (kX[index]/k_mag) * windx + (kZ[index]/k_mag) * windz;
		low_freq	= -k_sq * cutoff;
		exp_term	= -1.0 / ( L_sq * k_sq);
		philips		= sqrt((( exp(exp_term) * pow(k_dot_w, windAlign)) / (k_sq * k_sq)) * exp(low_freq));

		if(bDamp)
		{
			if(k_dot_w < 0.0f)
				philips *= (1.0f - damp);
		}		
		omega[index]   = sqrt(gravity * k_mag );
		hokReal[index] = (invSqRt2) * (rand1[index]) * philips;
		hokImag[index] = (invSqRt2) * (rand2[index]) * philips;
	}
}

inline void aaOceanClass::evaluateHieghtField()
{
	int  i,j,index, index_rev;
	double  hokRealOpp, hokImagOpp, sinwt, coswt;
	const int n = resolution;
	const int n_sq = n*n - 1;

	#pragma omp parallel for  private( index, index_rev ) 
	for(index = 0; index < n*n; index++)
	{
		index_rev = n_sq - index;//tail end  
		hokRealOpp	=  hokReal[index_rev];
		hokImagOpp	=  hokImag[index_rev];

		coswt = cos( omega[index] * time * waveSpeed);
		sinwt = sin( omega[index] * time * waveSpeed);

		hktReal[index]  =	( hokReal[index] *	coswt )  + 	( hokImag[index] *	sinwt )  + 
							( hokRealOpp	 *	coswt )  -  ( hokImagOpp	 *	sinwt )  ;  //complex conjugage
				
		hktImag[index]  =	(-hokReal[index] *	sinwt )  + 	( hokImag[index] *	coswt )  +
							( hokRealOpp	 *	sinwt )  +  ( hokImagOpp	 *	coswt )  ;  //complex conjugage
		
		fft_htField[index][0] = hktReal[index];
		fft_htField[index][1] = hktImag[index];
	}

	fftw_execute(planHeightField);

	#pragma omp parallel for private(i,j)
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			fft_htField[(i*n) + j][0] *= isEven(i+j)  * waveHeight;
		}
	}
}

inline void aaOceanClass::evaluateChopField()
{
	int  i,j,index;
	double _kX,_kZ, kMag;
	int n = resolution * resolution;

	#pragma omp parallel for private( index, _kX, _kZ, kMag) 
	for(index = 0; index < n; index++)
	{			
		kMag = sqrt(kX[index]*kX[index] + kZ[index]*kZ[index]);
		_kX = kX[index]/kMag;
		_kZ = kZ[index]/kMag;
		
		fft_chopX[index][0] =  hktImag[index] * _kX;
		fft_chopX[index][1] = -hktReal[index] * _kX;

		fft_chopZ[index][0] =  hktImag[index] * _kZ;
		fft_chopZ[index][1] = -hktReal[index] * _kZ;
	}

	#pragma omp sections
	{
		#pragma omp section
		{
			fftw_execute(planChopX);
		}
		#pragma omp section
		{
			fftw_execute(planChopZ);
		}
	}

	n = resolution;
	#pragma omp parallel for private(i,j, index)
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			index = (i*n) + j;
			fft_chopX[index][0] *= chopAmount * isEven(i+j) ;
			fft_chopZ[index][0] *= chopAmount * isEven(i+j) ;
		}
	}
}

inline void aaOceanClass::evaluateNormalsFinDiff()
{
	int i,j, nextI, prevI, nextJ, prevJ;
	double X,Z, nextX, prevX, nextZ, prevZ;
	mVector current,prev_iPos,next_iPos,prev_jPos,next_jPos,a,b,c,d,v1,v2;

	const int n = normalsResolution;	
	const int halfRes = n/2;	
	//const double mult1 = float(xsiGridULength) / float(n);
	//const double mult2 = float(xsiGridVLength) / float(n);
	const double mult1 = oceanScale / double(n);
	const double mult2 = oceanScale / double(n);

	#pragma omp parallel for private( i,j, nextI, prevI, nextJ, prevJ, nextX, prevX, nextZ, prevZ, current, prev_iPos, next_iPos, prev_jPos, next_jPos,a,b,c,d,v1,v2) 
	for(i = 0; i< n; i++)
	{					
		for(j = 0; j< n; j++)
		{
			nextI = i+1;
			prevI = i-1;
			nextJ = j+1;
			prevJ = j-1;

			X = (-halfRes + i) * mult1 ;
			Z =	(-halfRes + j) * mult2 ;
			nextX = (-halfRes + nextI) * mult1 ;
			nextZ =	(-halfRes + nextJ) * mult2 ;
			prevX = (-halfRes + prevI) * mult1 ;
			prevZ =	(-halfRes + prevJ) * mult2 ;

			nextI = wrap(i+1,n-1);
			prevI = wrap(i-1,n-1);
			nextJ = wrap(j+1,n-1);
			prevJ = wrap(j-1,n-1);
			
			current.i = (X - fft_chopX	[i*n+j][0]);
			current.j = (fft_htField	[i*n+j][0]);
			current.k = (Z - fft_chopZ	[i*n+j][0]);
			
			prev_iPos.i = (X - fft_chopX[		prevI*n+j][0]);
			prev_iPos.j = (fft_htField[			i*n+ prevJ][0]);
			prev_iPos.k = (prevZ - fft_chopZ[	prevI*n+j][0]);					
			
			next_iPos.i = (X - fft_chopX[		nextI*n+j][0]);
			next_iPos.j = (fft_htField[			i*n+nextJ][0]);
			next_iPos.k = (nextZ - fft_chopZ[	nextI*n+j][0]);
			
			prev_jPos.i = (prevX - fft_chopX[	i*n+ prevJ][0]);
			prev_jPos.j = (fft_htField[			prevI*n+ j][0]);
			prev_jPos.k = (Z - fft_chopZ[		i*n+ prevJ][0]);

			next_jPos.i = (nextX - fft_chopX[	i*n+nextJ][0]);
			next_jPos.j = (fft_htField[			nextI*n+j][0]);
			next_jPos.k = (Z - fft_chopZ[		i*n+nextJ][0]);
			
			a = prev_iPos - current;
			b = next_iPos - current;

			c = prev_jPos - current;
			d = next_jPos - current ;

			v1 = a/c;			
			v2 = b/d;
			v1 = ((v1+v2)/2);
			v1 = v1.unit();

			fft_normX[i*n+j][0] = v1.i; 
			fft_normY[i*n+j][0] = v1.j; 
			fft_normZ[i*n+j][0] = v1.k; 
		}
	}
}

inline void aaOceanClass::evaluateJacobians()
{
	int  i,j,index;
	double _kX,_kZ, kMag, kXZ;
	int n = resolution*resolution;

	#pragma omp parallel for private( index, _kX, _kZ, kXZ, kMag) 
	for(index = 0; index < n; index++)
	{			
		kMag = sqrt(kX[index]*kX[index] + kZ[index]*kZ[index]);
		_kX  = (kX[index] * kX[index]) / kMag;
		_kZ  = (kZ[index] * kZ[index]) / kMag;
		kXZ  = (kX[index] * kZ[index]) / kMag;

		fft_jxx[index][0] =  hktReal[index] * _kX;
		fft_jxx[index][1] =  hktImag[index] * _kX;

		fft_jzz[index][0] =  hktReal[index] * _kZ;
		fft_jzz[index][1] =  hktImag[index] * _kZ;

		fft_jxz[index][0] =  hktReal[index] * kXZ;
		fft_jxz[index][1] =  hktImag[index] * kXZ;
	}

	#pragma omp sections
	{
		#pragma omp section
		{
			fftw_execute(planJxx);
		}
		#pragma omp section
		{
			fftw_execute(planJzz);
		}
		#pragma omp section
		{
			fftw_execute(planJxz);
		}
	}

	n = resolution;
	#pragma omp parallel for private(i,j, index)
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			index = (i*n) + j;
			fft_jxx[index][0] = (fft_jxx[index][0] * -chopAmount * isEven(i+j)) + 1;
			fft_jzz[index][0] = (fft_jzz[index][0] * -chopAmount * isEven(i+j)) + 1;
			fft_jxz[index][0] =  fft_jxz[index][0] * -chopAmount * isEven(i+j);
		}
	}

	double jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz;
	#pragma omp parallel for private(index, jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz)
	for(index = 0; index < n*n; index++)
	{
		Jxx = fft_jxx[index][0];
		Jzz = fft_jzz[index][0];
		Jxz = fft_jxz[index][0];

		jPlus	= (0.5 * (Jxx + Jzz))  +  (0.5 * sqrt( ((Jxx - Jzz)*(Jxx - Jzz)) + 4 * (Jxz*Jxz) ));
		jMinus	= (0.5 * (Jxx + Jzz))  -  (0.5 * sqrt( ((Jxx - Jzz)*(Jxx - Jzz)) + 4 * (Jxz*Jxz) ));
		qPlus	= (jPlus  - Jxx) / Jxz;
		qMinus	= (jMinus - Jxx) / Jxz;

		eigenPlusX[index] = 1.0	    /  sqrt( 1.0 + qPlus * qPlus);
		eigenPlusZ[index] = qPlus   /  sqrt( 1.0 + qPlus * qPlus);

		eigenMinusX[index] = 1.0	/  sqrt( 1.0 + qMinus * qMinus);
		eigenMinusZ[index] = qMinus /  sqrt( 1.0 + qMinus * qMinus);

		//store foam back in this array for convenience
		//fft_jxz[index][0] =   (Jxx * Jzz) - (Jxz * Jxz); //original jacobian.
		fft_jxz[index][0] =   (Jxz * Jxz) - (Jxx * Jzz) + 1; // invert hack
	}
}

#endif  /* AAOCEANCLASS_H */

