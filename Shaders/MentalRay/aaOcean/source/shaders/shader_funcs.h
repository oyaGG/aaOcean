#ifndef SHADER_FUNCS_H
#define SHADER_FUNCS_H

#define maxim(a,b) (((a) > (b)) ? (a) : (b))

void shaderCleanup(aaOcean *&ocean)
{
	if(ocean->m_kX)
	{
		aligned_free(ocean->m_kX);
		ocean->m_kX=0;
	}
	if(ocean->m_kZ)
	{
		aligned_free(ocean->m_kZ);
		ocean->m_kZ=0;
	}
	if(ocean->m_omega)
	{
		aligned_free(ocean->m_omega);
		ocean->m_omega=0;
	}
	if(ocean->m_hokReal)
	{
		aligned_free(ocean->m_hokReal);
		ocean->m_hokReal=0;
	}
	if(ocean->m_hokImag)
	{
		aligned_free(ocean->m_hokImag);
		ocean->m_hokImag=0;
	}
	if(ocean->m_hktReal)
	{
		aligned_free(ocean->m_hktReal);
		ocean->m_hktReal=0;
	}
	if(ocean->m_hktImag)
	{
		aligned_free(ocean->m_hktImag);
		ocean->m_hktImag=0;
	}
	if(ocean->m_rand1)
	{
		aligned_free(ocean->m_rand1);
		ocean->m_rand1=0;
	}
	if(ocean->m_rand2)
	{
		aligned_free(ocean->m_rand2);
		ocean->m_rand2=0;
	}
	if(ocean->m_xCoord)
	{
		aligned_free(ocean->m_xCoord);
		ocean->m_xCoord=0;
	}
	if(ocean->m_zCoord)
	{
		aligned_free(ocean->m_zCoord);
		ocean->m_zCoord=0;
	}
	if(ocean->m_fft_jxx)
	{
		fftwf_free(ocean->m_fft_jxx);
		fftwf_destroy_plan(ocean->m_planJxx);
		ocean->m_fft_jxx=0;
	}
	if(ocean->m_fft_jzz)
	{
		fftwf_free(ocean->m_fft_jzz);
		fftwf_destroy_plan(ocean->m_planJzz);
		ocean->m_fft_jzz=0;
	}
}
/*
int getShaderInstanceID(miState *&state)
{
	int id=0;
	miTag shaderInst;
	const char  *shaderName;
	if ( mi_query(miQ_FUNC_TAG, state, 0, &shaderInst) )
		shaderName = mi_api_tag_lookup(shaderInst);
	XSI::CString name(shaderName);
	CString num = name.GetSubString(name.Length()-1);
	if(num != L"m")
		id = atoi(num.GetAsciiString());
	return id;
}*/

void getArrayBounds(fftwf_complex *&fftw_array, int idx, int n, float &min, float &max)
{
	max = -FLT_MAX;
	min =  FLT_MAX;

	int i, j, n1, index;
	n1 = n+1;
	//if fft_array has been copied and tiled, set idx = 1, else 0
	//#pragma omp parallel for private(index,i,j) shared(min, max)
	for(i = 0; i< n1; i++)
	{
		for(j = 0; j< n1; j++)
		{
			index = i*(n1)+j;
			//#pragma omp critical  // not cool
			{
				if(max < fftw_array[index][idx]) //precision -- epsilon check
					max = fftw_array[index][idx];

				if(min > fftw_array[index][idx]) //precision -- epsilon check
					min = fftw_array[index][idx];
			}
		}
	}
}

bool shader_init(oceanStore **&os, miState *&state)
{
	if(mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os ))
		*os = (oceanStore *)mi_mem_allocate( sizeof( oceanStore ) );
	else
		return FALSE;	
	
	return TRUE;
}

inline float catrom(float t, float a, float b, float c, float d) 
{
   return  0.5f * ( ( 2.0f * b ) + ( -a + c ) * t + ( 2.0f * a - 5.0f * b + 4 * c - d ) * t*t + ( -a + 3.0f * b - 3.0f * c + d )* t*t*t );
}

float catromPrep(aaOcean *&ocean, fftwf_complex *&fftw_array, miState *&state, miVector *&coord)
{
	float u,v,du,dv=0;
	int xMinus1, yMinus1, x, y, xPlus1, yPlus1, xPlus2, yPlus2;
	int n = ocean->m_resolution;	

	miVector point;	
	mi_point_to_object(state, &point, coord); // not needed

	if(point.y > 1.0f)
		point.y = point.y - floor(point.y);
	if(point.x > 1.0f)
		point.x = point.x - floor(point.x);
	if(point.y < 0.0f)
		point.y = point.y - floor(point.y);
	if(point.x < 0.0f)
		point.x = point.x - floor(point.x);

	u = point.y * float(n);
	v = point.x * float(n);
	x = (int)floor(u);
	y = (int)floor(v);

	xMinus1 = wrap((x-1),n);
	yMinus1 = wrap((y-1),n);	
	x = wrap(x,n);
	y = wrap(y,n);		
	xPlus1 = wrap((x+1),n);	
	yPlus1 = wrap((y+1),n);	
	xPlus2 = wrap((x+2),n);	
	yPlus2 = wrap((y+2),n);	

	du = u - x; 
	dv = v - y;	

	float a1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) + yMinus1][1],
						fftw_array[x*(n+1)		 + yMinus1][1],
						fftw_array[xPlus1*(n+1)	 + yMinus1][1],
						fftw_array[xPlus2*(n+1)	 + yMinus1][1] 
						);

	float b1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) +	y][1],
						fftw_array[x*(n+1)		 +	y][1],
						fftw_array[xPlus1*(n+1)	 +	y][1],
						fftw_array[xPlus2*(n+1)	 +	y][1]
						);

	float c1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) + yPlus1][1],
						fftw_array[x*(n+1)		 + yPlus1][1],
						fftw_array[xPlus1*(n+1)	 + yPlus1][1],
						fftw_array[xPlus2*(n+1)	 + yPlus1][1] 
						);

	float d1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) + yPlus2][1],
						fftw_array[x*(n+1)		 + yPlus2][1],
						fftw_array[xPlus1*(n+1)	 + yPlus2][1],
						fftw_array[xPlus2*(n+1)	 + yPlus2][1]
						);

	return catrom(dv,a1,b1,c1,d1);
	
}

void rotateArray(fftwf_complex *&fft_array, aaOcean *&ocean)
{
	int n = ocean->m_resolution + 1;
	int halfn = n/2;
	float tmp;
	for (int i = 0; i < halfn; ++i)
	{
        for (int j = i; j < n-i-1; ++j)
		{
            tmp								 = fft_array[i*n+j][1];
            fft_array[i*n+j][1]				 = fft_array[j*n +(n-i-1)][1];
            fft_array[j*n +(n-i-1)][1]		 = fft_array[(n-i-1)*n +(n-j-1)][1];
            fft_array[(n-i-1)*n +(n-j-1)][1] = fft_array[(n-j-1)*n + i][1];
            fft_array[(n-j-1)*n + i][1]		 = tmp;
        }
	}
}

void copy_and_tile(fftwf_complex *&fft_array, aaOcean *&ocean)
{
	int n = ocean->m_resolution;
	int n1 = n + 1;
	int index, i, j;
	
	#pragma omp parallel for private(index, i, j)
	for(i = 0; i< n1; i++)
	{					
		for(j = 0; j< n1; j++)
		{
			index = i*n1 + j;
			if( i<n && j<n) // regular array copy
			{
				fft_array[index][1] = fft_array[i*n+j][0];
			}
			else
			{
				if(i==n)  //copy left-most col to right-most col
				{
					fft_array[index][1] = fft_array[j][0];
				}
				if(j==n) // copy top row to bottom row
				{
					fft_array[index][1] = fft_array[i*n][0];
				}
				if(i==n && j==n) // copy top-left corner to bottom-left
				{
					fft_array[index][1] = fft_array[0][0];
				}
			}
			
		}
	}
}

#endif //SHADER_FUNCS
