#ifndef AAOCEAN_CPP
#define AAOCEAN_CPP

#include "aaOceanClass.h"

aaOcean::aaOcean() :	
	m_pointCount(0),			m_resolution(0),			m_windAlign(0),			m_seed(1),			
	m_3DGridULength(-1.0f),		m_3DGridVLength(-1.0f),		m_velocity(-1.0f),		m_windDir(-1.0f),
	m_cutoff(-1.0f),			m_damp(-1.0f),				m_oceanScale(-1.0f),	m_waveHeight(-1.0f), 
	m_time(-1.0),				m_isAllocated(FALSE),		m_isValid(FALSE),		m_isNormalsAllocated(FALSE),
	m_isFoamAllocated(FALSE),	m_isSplashAllocated(FALSE),	m_redoHoK(FALSE),		m_isDisplacementDirty(FALSE),
	m_isNormalsDirty(FALSE),	m_isFoamDirty(FALSE),		m_isShader(FALSE),		m_fmin(FLT_MAX), 
	m_fmax(FLT_MAX),			m_zCoord(0),				m_xCoord(0),			m_hokReal(0),
	m_hokImag(0),				m_hktReal(0),				m_hktImag(0),			m_kX(0),
	m_kZ(0),					m_omega(0),					m_rand1(0),				m_rand2(0),
	m_eigenPlusX(0),			m_eigenPlusZ(0),			m_eigenMinusX(0),		m_eigenMinusZ(0),
	m_fft_htField(0),			m_fft_chopX(0),				m_fft_chopZ(0),			m_fft_jxx(0),
	m_fft_jzz(0),				m_fft_jxz(0),				m_fft_normX(0),			m_fft_normY(0),
	m_fft_normZ(0),				m_randomType(0.f)
{
	fftwf_init_threads();
}

aaOcean::aaOcean(const aaOcean &cpy)
{

}

aaOcean::~aaOcean()
{
	clearArrays();
	fftwf_cleanup_threads();
	fftwf_cleanup();
}

void aaOcean::input(int resolution, ULONG seed, float oceanScale, float velocity, float cutoff, float windDir, 
			int windAlign, float damp, float	waveSpeed, float waveHeight, float chopAmount, float time, float randomType)
{
	bool isDirty = FALSE; 
	resolution	= (int)intpow(2.0f,  (4 + resolution)); //pow(2.0f, (4 + resolution));
	oceanScale	= maximum<float>(oceanScale, 0.00001f);
	velocity	= maximum<float>(((velocity  * velocity) / aa_GRAVITY), 0.00001f);
	cutoff		= fabs(cutoff * 0.01f);
	windDir		= windDir * aa_PIBY180;
	windAlign	= maximum<int>(((windAlign + 1) * 2),2); 
	damp		= minimum<float>(damp, 1.0f);

	waveHeight *= 0.01f;
	chopAmount *= 0.01f;

	if(m_time != time || m_waveSpeed != waveSpeed || m_waveHeight != waveHeight || m_chopAmount != chopAmount )
		isDirty = TRUE;

	m_time			= time;
	m_waveSpeed		= waveSpeed;
	m_waveHeight	= waveHeight;
	m_chopAmount	= chopAmount;

	if( m_oceanScale	!= oceanScale	||
		m_windDir		!= windDir		||
		m_cutoff		!= cutoff		||
		m_velocity		!= velocity	||
		m_windAlign		!= windAlign	||
		m_damp			!= damp		)
	{
		m_oceanScale	= oceanScale;
		m_windDir		= windDir;
		m_cutoff		= cutoff;
		m_velocity		= velocity;
		m_windAlign		= windAlign;
		m_damp			= damp;
		m_redoHoK		= TRUE;
		isDirty			= TRUE;
	}

	if(m_seed != seed)
	{
		m_seed		= seed;
		m_redoHoK	= TRUE;
		if(m_resolution == resolution)
			setup_grid(); 
		isDirty = TRUE;
	}
	if(isDirty)
		m_isDisplacementDirty = m_isNormalsDirty = m_isFoamDirty = TRUE;
	else
		m_isDisplacementDirty = m_isNormalsDirty = m_isFoamDirty = FALSE;

	reInit(resolution);

}

bool aaOcean::reInit(int data_size)
{
	if(((data_size & (data_size - 1)) != 0) || data_size == 0) //	not power of 2
	{	
		sprintf_s(m_state,"[aaOcean] :  invalid point resolution of %d. Please select a power-of-2 subdivision value", data_size);
		m_isValid = FALSE;
	}
	else
	{
		if(m_resolution != data_size || !m_isAllocated )
		{
			if(m_isShader)
			{
				sprintf_s(m_state,"[aaOcean] : Building ocean shader with resolution %dx%d", data_size, data_size);
			}
			else
			{
				sprintf_s(m_state,"[aaOcean ICE] : Building ocean grid with resolution %dx%d", data_size, data_size);
			}

			m_resolution = data_size;
			allocateBaseArrays();				
			m_redoHoK  = TRUE;
			setup_grid();
		}
		m_isValid = TRUE;
	}
	return m_isValid;
}

void aaOcean::prepareOcean(bool doHeightField, bool doChopField, bool doJacobians, bool doNormals)
{
	if( m_redoHoK )
	{
		evaluateHokData();
		m_redoHoK = FALSE;
		m_isDisplacementDirty = m_isNormalsDirty = m_isFoamDirty = TRUE;
	}

	if(doHeightField)
	{
		evaluateHieghtField();
		m_isDisplacementDirty = m_isNormalsDirty = m_isFoamDirty = TRUE;
	}

	if(doChopField && m_chopAmount > 0.0f)
	{
		evaluateChopField();
		m_isDisplacementDirty = m_isNormalsDirty = m_isFoamDirty = TRUE;
	}

	if(doJacobians && m_isAllocated)
	{
		if(!m_isFoamAllocated)
			allocateFoamArrays();
		if(!m_isSplashAllocated)
			allocateSplashArrays();
		evaluateJacobians();
		m_isDisplacementDirty = m_isNormalsDirty = m_isFoamDirty = TRUE;
	}
	if(doNormals)
	{
		if(!m_isNormalsAllocated)
			allocateNormalsArrays();
		evaluateNormalsFinDiff();
	}
}

void aaOcean::allocateBaseArrays()
{
	if(m_isAllocated) 
		clearArrays();

	m_kX		= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_kZ		= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_omega		= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_hokReal	= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_hokImag	= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_hktReal	= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_hktImag	= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_rand1		= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_rand2		= (float*) aligned_malloc(m_resolution * m_resolution * sizeof(float), 16);
	m_xCoord	= (int*)   aligned_malloc(m_resolution * m_resolution * sizeof(int),   16);
	m_zCoord	= (int*)   aligned_malloc(m_resolution * m_resolution * sizeof(int),   16);

	if(m_resolution > 254)
	{
		int threads = omp_get_num_procs();
		fftwf_plan_with_nthreads(threads);
		char msg[100];
		sprintf_s(msg,"[aaOcean] : Launching threaded FFTW with %d threads", threads);
	}
	else
		fftwf_plan_with_nthreads(1);

	m_fft_htField	= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((m_resolution+1) * (m_resolution+1)));
	m_fft_chopX		= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((m_resolution+1) * (m_resolution+1)));
	m_fft_chopZ		= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((m_resolution+1) * (m_resolution+1)));

	m_planHeightField	= fftwf_plan_dft_2d(m_resolution,m_resolution,m_fft_htField,m_fft_htField,1, FFTW_ESTIMATE);
	m_planChopX			= fftwf_plan_dft_2d(m_resolution,m_resolution,m_fft_chopX ,m_fft_chopX	,1, FFTW_ESTIMATE);
	m_planChopZ			= fftwf_plan_dft_2d(m_resolution,m_resolution,m_fft_chopZ ,m_fft_chopZ	,1, FFTW_ESTIMATE);
	m_isAllocated		= TRUE;
}

void aaOcean::allocateNormalsArrays()
{
	int arSize = m_resolution + 1;
	m_fft_normX	= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((arSize)*(arSize)));
	m_fft_normY	= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((arSize)*(arSize)));
	m_fft_normZ	= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((arSize)*(arSize)));
	m_isNormalsAllocated = TRUE;
}

 void aaOcean::allocateFoamArrays()
{
	m_fft_jxx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((m_resolution+1) * (m_resolution+1)));
	m_fft_jxz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((m_resolution+1) * (m_resolution+1)));
	m_fft_jzz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((m_resolution+1) * (m_resolution+1)));
	m_planJxx = fftwf_plan_dft_2d(m_resolution, m_resolution, m_fft_jxx, m_fft_jxx, 1, FFTW_MEASURE);
	m_planJxz = fftwf_plan_dft_2d(m_resolution, m_resolution, m_fft_jxz, m_fft_jxz, 1, FFTW_MEASURE);
	m_planJzz = fftwf_plan_dft_2d(m_resolution, m_resolution, m_fft_jzz, m_fft_jzz, 1, FFTW_MEASURE);
	m_isFoamAllocated = TRUE;
}

void aaOcean::allocateSplashArrays()
{
	m_eigenPlusX	= (float*) aligned_malloc((m_resolution+1) * (m_resolution+1) * sizeof(float), 16);
	m_eigenPlusZ	= (float*) aligned_malloc((m_resolution+1) * (m_resolution+1) * sizeof(float), 16);
	m_eigenMinusX	= (float*) aligned_malloc((m_resolution+1) * (m_resolution+1) * sizeof(float), 16);
	m_eigenMinusZ	= (float*) aligned_malloc((m_resolution+1) * (m_resolution+1) * sizeof(float), 16);
	m_isSplashAllocated = TRUE;
}

void aaOcean::clearResidualArrays()
{
	if(m_kX)
	{
		aligned_free(m_kX);
		m_kX = FALSE;
	}
	if(m_kZ)
	{
		aligned_free(m_kZ);
		m_kZ = FALSE;
	}
	if(m_omega)
	{
		aligned_free(m_omega);
		m_omega = FALSE;
	}
	if(m_hokReal)
	{
		aligned_free(m_hokReal);
		m_hokReal = FALSE;
	}
	if(m_hokImag)
	{
		aligned_free(m_hokImag);
		m_hokImag = FALSE;
	}
	if(m_hktReal)
	{
		aligned_free(m_hktReal);
		m_hktReal = FALSE;
	}
	if(m_hktImag)
	{
		aligned_free(m_hktImag);
		m_hktImag = FALSE;
	}
	if(m_rand1)
	{
		aligned_free(m_rand1);
		m_rand1 = FALSE;
	}
	if(m_rand2)
	{
		aligned_free(m_rand2);
		m_rand2 = FALSE;
	}
	if(m_xCoord)
	{
		aligned_free(m_xCoord);
		m_xCoord = FALSE;
	}
	if(m_zCoord)
	{
		aligned_free(m_zCoord);
		m_zCoord = FALSE;
	}
	if(m_fft_jxx)
	{
		fftwf_free(m_fft_jxx);
		fftwf_destroy_plan(m_planJxx);
		m_fft_jxx=0;
	}
	if(m_fft_jzz)
	{
		fftwf_free(m_fft_jzz);
		fftwf_destroy_plan(m_planJzz);
		m_fft_jzz = FALSE;
	}
}

void aaOcean::clearArrays()
{
	if(m_isAllocated)
	{
		clearResidualArrays();

		if(m_fft_htField)
		{
			fftwf_free(m_fft_htField);
			fftwf_destroy_plan(m_planHeightField);
			m_fft_htField = FALSE;
		}
		if(m_fft_chopX)
		{
			fftwf_free(m_fft_chopX);
			fftwf_destroy_plan(m_planChopX);
			m_fft_chopX = FALSE;
		}
		if(m_fft_chopZ)
		{
			fftwf_free(m_fft_chopZ);
			fftwf_destroy_plan(m_planChopZ);
			m_fft_chopZ = FALSE;
		}
		m_isAllocated = FALSE;
	}
	if(m_isFoamAllocated)
	{
		if(m_fft_jxz)
		{
			fftwf_free(m_fft_jxz);
			fftwf_destroy_plan(m_planJxz);
			m_fft_jxz = FALSE;
		}
		m_isFoamAllocated = FALSE;
	}
	if(m_isSplashAllocated)
	{
		if(m_eigenPlusX)
		{
			aligned_free(m_eigenPlusX);
			m_eigenPlusX = FALSE;
		}
		if(m_eigenPlusZ)
		{
			aligned_free(m_eigenPlusZ);
			m_eigenPlusZ = FALSE;
		}
		if(m_eigenMinusX)
		{
			aligned_free(m_eigenMinusX);
			m_eigenMinusX = FALSE;
		}
		if(m_eigenMinusZ)
		{
			aligned_free(m_eigenMinusZ);
			m_eigenMinusZ = FALSE;
		}

		m_isSplashAllocated = FALSE;
	}
	if(m_isNormalsAllocated)
	{
		if(m_fft_normX)
		{
			fftwf_free(m_fft_normX);
			m_fft_normX = FALSE;
		}
		if(m_fft_normZ)
		{
			fftwf_free(m_fft_normY);
			m_fft_normY = FALSE;
		}
		if(m_fft_normZ)
		{
			fftwf_free(m_fft_normZ);
			m_fft_normZ = FALSE;
		}
		m_isNormalsAllocated = FALSE;
	}
}

ULONG aaOcean::get_uID(float xCoord, float zCoord)
{
	// a very simple hash function. should probably do a better one at some point
	register float angle = 0.0f;
	register float length = 0.0f;
	register float id_out;
	ULONG returnVal = 1;

	if (zCoord != 0.0f && xCoord != 0.0f)
	{
		angle = xCoord / sqrt(zCoord * zCoord + xCoord * xCoord); 
		angle = acos(angle);
		angle = RadsToDegs(angle);
		
		if (angle > 180.0f)
			angle = 360.0f - angle;
		
		length = sqrt(zCoord * zCoord + xCoord * xCoord); 
		id_out = (length * length)  +  (length * angle);
		
		if( angle == 0.0f)
			returnVal = (ULONG)(length * length);
		else if (zCoord <= 0.0f)
			returnVal = (ULONG)floor(id_out + 0.5f);
		else
			returnVal = INT_MAX - (ULONG)floor(id_out + 0.5) ;
	}
	return returnVal;
}

void aaOcean::setup_grid()
{
	if(!m_isAllocated)
		return;
	register const int n = m_resolution;
	register const int half_n = (-n / 2) - ((n-1) / 2);
	register ULONG index, uID; 

	Timer timer;
	timer.start();
	#pragma omp parallel for private(index, uID)
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			index = (i * n) + j;

			m_xCoord[index] = half_n + i * 2 ;
			m_zCoord[index] = half_n + j * 2 ;

			uID = (ULONG)get_uID((float)m_xCoord[index], (float)m_zCoord[index]);

			StochasticLib1 sto(uID + m_seed);
			m_rand1[index] = (float)sto.Normal(0.0f, 1.0f);
			m_rand2[index] = (float)sto.Normal(0.0f, 1.0f);
		}
	}
	timer.stop();
	Application().LogMessage("RanGen Time : " + CString(timer.getElapsedTimeInMilliSec()) + "ms");
}

 void aaOcean::evaluateHokData()
{
	register float k_sq, k_mag,  k_dot_w, philips, x, z;
	register const int	n			 = m_resolution * m_resolution;
	register const float	k_mult	 = (2.0f * aa_PI) / m_oceanScale;
	register const float	L		 = m_velocity;
	register const float	L_sq	 = L*L;
	register const float	windx	 = cos(m_windDir);
	register const float	windz	 = sin(m_windDir);
	
	bool bDamp	= FALSE;
	if (m_damp > 0.0f)
		bDamp = true;

	#pragma omp parallel for private( k_sq, k_mag, k_dot_w, philips, x, z)  
	for(int index = 0; index < n; ++index)
	{
		x = m_kX[index] =  m_xCoord[index] * k_mult; 
		z = m_kZ[index] =  m_zCoord[index] * k_mult;

		//philips spectrum vars
		k_sq		= (x * x) + (z * z);
		k_mag		= 1.0f / sqrt( k_sq );
		k_dot_w		= (x * k_mag) * windx + (z * k_mag) * windz;
		philips		= sqrt((( exp(-1.0f / ( L_sq * k_sq)) * pow(k_dot_w, m_windAlign)) / 
					  (k_sq * k_sq)) * exp(-k_sq * m_cutoff));

		if(bDamp)
		{
			if(k_dot_w < 0.0f)
				philips *= (1.0f - m_damp);
		}		
		m_omega[index]   = sqrt(aa_GRAVITY / k_mag );
		m_hokReal[index] = (aa_INV_SQRTTWO) * (m_rand1[index]) * philips;
		m_hokImag[index] = (aa_INV_SQRTTWO) * (m_rand2[index]) * philips;
	}
}

 void aaOcean::evaluateHieghtField()
{
	int  i,j,index, index_rev;
	register float  hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt;
	float wt = m_waveSpeed * m_time;
	const int n = m_resolution;
	const int nn = n * n;
	register const int n_sq = n * n - 1;

	#pragma omp parallel for private( index, index_rev, hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt )  
	for(index = 0; index < nn; ++index)
	{
		index_rev   = n_sq - index; //tail end 
		hokReal		= m_hokReal[index];
		hokImag		= m_hokImag[index];
		hokRealOpp	= m_hokReal[index_rev];
		hokImagOpp	= m_hokImag[index_rev];

		coswt = cos( m_omega[index] * wt);
		sinwt = sin( m_omega[index] * wt);

		m_hktReal[index]  =	( hokReal    *	coswt )  + 	( hokImag	 *	sinwt )  + 
							( hokRealOpp *	coswt )  -  ( hokImagOpp *	sinwt )  ;  //complex conjugage
				
		m_hktImag[index]  =	(-hokReal	 *	sinwt )  + 	( hokImag	 *	coswt )  +
							( hokRealOpp *	sinwt )  +  ( hokImagOpp *	coswt )  ;  //complex conjugage
		
		m_fft_htField[index][0] = m_hktReal[index];
		m_fft_htField[index][1] = m_hktImag[index];
	}

	fftwf_execute(m_planHeightField);

	#pragma omp parallel for private(i,j)
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			m_fft_htField[(i*n) + j][0] *= isEven(i+j)  * m_waveHeight;
		}
	}
}

 void aaOcean::evaluateChopField()
{
	int  i,j,index;
	register float _kX,_kZ, kMag;
	int n = m_resolution * m_resolution;
	#pragma omp parallel for private( index, _kX, _kZ, kMag)  
	for(index = 0; index < n; ++index)
	{			
		kMag = sqrt(m_kX[index] * m_kX[index] + m_kZ[index] * m_kZ[index]);
		_kX = m_kX[index] / kMag;
		_kZ = m_kZ[index] / kMag;
		
		m_fft_chopX[index][0] =  m_hktImag[index] * _kX;
		m_fft_chopX[index][1] = -m_hktReal[index] * _kX;

		m_fft_chopZ[index][0] =  m_hktImag[index] * _kZ;
		m_fft_chopZ[index][1] = -m_hktReal[index] * _kZ;
	}

	fftwf_execute(m_planChopX);
	fftwf_execute(m_planChopZ);

	n = m_resolution;
	#pragma omp parallel for private(i,j, index)  
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			index = (i*n) + j;
			m_fft_chopX[index][0] *= m_chopAmount * isEven(i+j) ;
			m_fft_chopZ[index][0] *= m_chopAmount * isEven(i+j) ;
		}
	}
}

void aaOcean::evaluateNormalsFinDiff()
{
	int i,j, nextI, prevI, nextJ, prevJ;
	float X,Z, nextX, prevX, nextZ, prevZ;
	vector3 current, prev_iPos, next_iPos, prev_jPos, next_jPos, a, b, c, d, v1, v2;

	const int n = m_resolution;	
	const int halfRes = n / 2;	
	register const float mult1 = float(m_oceanScale) / float(n);
	register const float mult2 = float(m_oceanScale) / float(n);

	// #pragma omp parallel for private( i,j, nextI, prevI, nextJ, prevJ, nextX, prevX, nextZ, prevZ, current, prev_iPos, next_iPos, prev_jPos, next_jPos,a,b,c,d,v1,v2)  
	for(i = 0; i< n; ++i)
	{					
		for(j = 0; j< n; ++j)
		{
			nextI = i+1;
			prevI = i-1;
			nextJ = j+1;
			prevJ = j-1;

			X = (-halfRes + i) * mult1;
			Z =	(-halfRes + j) * mult2;
			nextX = (-halfRes + nextI) * mult1;
			nextZ =	(-halfRes + nextJ) * mult2;
			prevX = (-halfRes + prevI) * mult1;
			prevZ =	(-halfRes + prevJ) * mult2;

			nextI = wrap(i+1, n-1);
			prevI = wrap(i-1, n-1);
			nextJ = wrap(j+1, n-1);
			prevJ = wrap(j-1, n-1);
			
			current.x = (X - m_fft_chopX		[i * n + j][0]);
			current.y = (m_fft_htField			[i * n + j][0]);
			current.z = (Z - m_fft_chopZ		[i * n + j][0]);
			
			prev_iPos.x = (X - m_fft_chopX		[prevI * n + j][0]);
			prev_iPos.y = (m_fft_htField		[i * n + prevJ][0]);
			prev_iPos.z = (prevZ - m_fft_chopZ	[prevI * n + j][0]);					
			
			next_iPos.x = (X - m_fft_chopX		[nextI * n + j][0]);
			next_iPos.y = (m_fft_htField		[i * n + nextJ][0]);
			next_iPos.z = (nextZ - m_fft_chopZ	[nextI * n + j][0]);
			
			prev_jPos.x = (prevX - m_fft_chopX	[i * n + prevJ][0]);
			prev_jPos.y = (m_fft_htField		[prevI * n + j][0]);
			prev_jPos.z = (Z - m_fft_chopZ		[i * n + prevJ][0]);

			next_jPos.x = (nextX - m_fft_chopX	[i * n + nextJ][0]);
			next_jPos.y = (m_fft_htField		[nextI * n + j][0]);
			next_jPos.z = (Z - m_fft_chopZ		[i * n + nextJ][0]);
			
			a = prev_iPos - current;
			b = next_iPos - current;

			c = prev_jPos - current;
			d = next_jPos - current ;

			v1 = a / c;			
			v2 = b / d;
			v1 = ((v1 + v2) / 2.0f);
			v1 = v1.normalize();

			m_fft_normX[i * n + j][0] = v1.x; 
			m_fft_normY[i * n + j][0] = v1.y; 
			m_fft_normZ[i * n + j][0] = v1.z; 
		}
	}
}

void aaOcean::evaluateJacobians()
{
	int  i, j, index;
	register float _kX, _kZ, kMag, kXZ;
	int n = m_resolution * m_resolution;

	#pragma omp parallel for private( index, _kX, _kZ, kXZ, kMag) 
	for(index = 0; index < n; ++index)
	{			
		kMag = sqrt(m_kX[index] * m_kX[index] + m_kZ[index] * m_kZ[index]);
		_kX  = (m_kX[index] * m_kX[index]) / kMag;
		_kZ  = (m_kZ[index] * m_kZ[index]) / kMag;
		kXZ  = (m_kX[index] * m_kZ[index]) / kMag;

		m_fft_jxx[index][0] =  m_hktReal[index] * _kX;
		m_fft_jxx[index][1] =  m_hktImag[index] * _kX;

		m_fft_jzz[index][0] =  m_hktReal[index] * _kZ;
		m_fft_jzz[index][1] =  m_hktImag[index] * _kZ;

		m_fft_jxz[index][0] =  m_hktReal[index] * kXZ;
		m_fft_jxz[index][1] =  m_hktImag[index] * kXZ;
	}

	fftwf_execute(m_planJxx);
	fftwf_execute(m_planJzz);
	fftwf_execute(m_planJxz);

	n = m_resolution;
	#pragma omp parallel for private(i, j, index)  
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			index = (i*n) + j;
			m_fft_jxx[index][0] = (m_fft_jxx[index][0] * -m_chopAmount * isEven(i+j)) + 1.0f;
			m_fft_jzz[index][0] = (m_fft_jzz[index][0] * -m_chopAmount * isEven(i+j)) + 1.0f;
			m_fft_jxz[index][0] =  m_fft_jxz[index][0] * -m_chopAmount * isEven(i+j);
		}
	}

	register float jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz;
	#pragma omp parallel for private(index, jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz)  
	for(index = 0; index < n*n; ++index)
	{
		Jxx = m_fft_jxx[index][0];
		Jzz = m_fft_jzz[index][0];
		Jxz = m_fft_jxz[index][0];

		jPlus	= (0.5f * (Jxx + Jzz))  +  (0.5f * sqrt( ((Jxx - Jzz)*(Jxx - Jzz)) + 4.0f * (Jxz*Jxz) ));
		jMinus	= (0.5f * (Jxx + Jzz))  -  (0.5f * sqrt( ((Jxx - Jzz)*(Jxx - Jzz)) + 4.0f * (Jxz*Jxz) ));
		qPlus	= (jPlus  - Jxx) / Jxz;
		qMinus	= (jMinus - Jxx) / Jxz;

		m_eigenPlusX[index] = 1.0f	    /  sqrt( 1.0f + qPlus * qPlus);
		m_eigenPlusZ[index] = qPlus		/  sqrt( 1.0f + qPlus * qPlus);

		m_eigenMinusX[index] = 1.0f		/  sqrt( 1.0f + qMinus * qMinus);
		m_eigenMinusZ[index] = qMinus	/  sqrt( 1.0f + qMinus * qMinus);

		//store foam back in this array for convenience
		//fft_jxz[index][0] =   (Jxx * Jzz) - (Jxz * Jxz); //original jacobian.
		m_fft_jxz[index][0] =   (Jxz * Jxz) - (Jxx * Jzz) + 1.0f; // invert hack
	}
}
#endif  /* AAOCEANCLASS_CPP */