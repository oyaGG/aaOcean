#ifndef SHADER_FUNCS_H
#define SHADER_FUNCS_H

bool fetchInput(AtParamValue *&params, aaOcean *&ocean)
{
	float tempFLT;
	int   tempINT;
	bool  isDirty = FALSE;
	ocean->m_redoHoK = FALSE;

	tempFLT = maximum<float>(params[p_oceanScale].FLT, 0.00001f);
	if(ocean->m_oceanScale	!= tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("m_oceanScale = %f", params[p_oceanScale].FLT);
		#endif
		ocean->m_oceanScale = tempFLT;
		ocean->m_redoHoK = TRUE;
	}

	tempINT =  params[p_seed].INT;
	if(ocean->m_seed != tempINT)
	{
		#ifdef DEBUG 
		AiMsgWarning("seed = %d", tempINT); 
		#endif
		ocean->m_seed	= tempINT;
		ocean->m_redoHoK = TRUE;
		ocean->setup_grid();
	}

	tempFLT = maximum<float>(((params[p_velocity].FLT * params[p_velocity].FLT) / (aa_GRAVITY)), 0.00001f);
	if(ocean->m_velocity !=  tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("velocity = %f", params[p_velocity].FLT);
		#endif
		ocean->m_velocity = tempFLT;
		ocean->m_redoHoK = TRUE;
	}

	tempFLT = fabs(params[p_cutoff].FLT * 0.01f);
	if(ocean->m_cutoff != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("m_cutoff = %f", params[p_cutoff].FLT);
		#endif
		ocean->m_cutoff = tempFLT;
		ocean->m_redoHoK = TRUE;
	}

	tempFLT = DegsToRads(params[p_windDir].FLT);
	if(ocean->m_windDir != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("m_windDir = %f degrees, %f radians", params[p_windDir].FLT, tempFLT); 
		#endif
		ocean->m_windDir = tempFLT;
		ocean->m_redoHoK = TRUE;
	}
	
	tempINT = maximum<int>(((params[p_windAlign].INT + 1) * 2), 2); 
	if(ocean->m_windAlign != tempINT)
	{
		#ifdef DEBUG 
		AiMsgWarning("windAlign = %d", tempINT);
		#endif
		ocean->m_windAlign = tempINT;
		ocean->m_redoHoK = TRUE;
	}
	
	tempFLT = minimum<float>(params[p_damp].FLT,1);
	if(ocean->m_damp != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("damp = %f", tempFLT); 
		#endif
		ocean->m_damp = tempFLT;
		ocean->m_redoHoK = TRUE;
	}

	isDirty = ocean->m_redoHoK;
	
	tempFLT = params[p_chopAmount].FLT * 0.01f;
	if(ocean->m_chopAmount != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("chopAmount = %f", tempFLT); 
		#endif
		ocean->m_chopAmount = tempFLT;
		isDirty = TRUE;
	}

	tempFLT = params[p_waveHeight].FLT * 0.01f;
	if(ocean->m_waveHeight != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("m_waveHeight = %f", tempFLT); 
		#endif
		ocean->m_waveHeight = tempFLT;
		isDirty = TRUE;
	}

	tempFLT = params[p_waveSpeed].FLT;
	if(ocean->m_waveSpeed != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("m_waveSpeed = %f", tempFLT); 
		#endif
		ocean->m_waveSpeed = tempFLT;
		isDirty = TRUE;
	}

	tempFLT = params[p_time].FLT;
	if(ocean->m_time != tempFLT)
	{
		#ifdef DEBUG 
		AiMsgWarning("m_time = %f", tempFLT); 
		#endif
		ocean->m_time = tempFLT;
		isDirty = TRUE;
	}

	return isDirty;
}

void getArrayBounds(fftwf_complex *&fftw_array, int idx, int n, float &min, float &max)
{
	max = -FLT_MAX;
	min =  FLT_MAX;

	int i, j, n1, index;
	n1 = n+1;
	//if fft_array has been copied and tiled, set idx = 1, else 0
	//#pragma omp parallel for private(index,i,j)
	for(i = 0; i< n1; i++)
	{
		for(j = 0; j< n1; j++)
		{
			index = i*(n1)+j;
			if(max < fftw_array[index][idx]) //precision -- epsilon check
				max = fftw_array[index][idx];

			if(min > fftw_array[index][idx]) //precision -- epsilon check
				min = fftw_array[index][idx];
		}
	}
}

inline float catrom(float &t, float &a, float &b, float &c, float &d)
{
   return  0.5f * ( ( 2.0f * b ) + ( -a + c ) * t + ( 2.0f * a - 5.0f * b + 4 * c - d ) * t*t + ( -a + 3.0f * b - 3.0f * c + d )* t*t*t );
}


float catromPrep(aaOcean *&ocean, fftwf_complex *&fftw_array, AtPoint point)
{
	// prepares end points 
	float u, v, du, dv = 0.f;
	int xMinus1, yMinus1, x, y, xPlus1, yPlus1, xPlus2, yPlus2;
	int n = ocean->m_resolution;

	if(point.y > 1.00f)
		point.y = point.y - floor(point.y);
	if(point.x > 1.00f)
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


void copy_and_tile(fftwf_complex *&fft_array, aaOcean *&ocean)
{
	// This function is admittedly crap
	// It makes the ocean grid tileable (when it should already be tileable!!)
	// Needlessly inefficient. I should fix this.
	// Also, optimze this from O(N^2) to O(2N) later

	int n = ocean->m_resolution;
	int n1 = n + 1;
	int index, i, j;
	
	#pragma omp parallel for private(index,i,j)
	for(i = 0; i< n1; i++)
	{					
		for(j = 0; j< n1; j++)
		{
			index = i*n1 + j;
			if( i<n && j<n) // regular array copy -- what a waste of resources
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


void rotateArray(fftwf_complex *&fft_array, aaOcean *&ocean)
{
	int n = ocean->m_resolution + 1;
	float tmp;

	for (int i=0; i<n/2; i++)
	{
        for (int j=i; j<n-i-1; j++)
		{
            tmp								 = fft_array[i*n+j][1];
            fft_array[i*n+j][1]				 = fft_array[j*n +(n-i-1)][1];
            fft_array[j*n +(n-i-1)][1]		 = fft_array[(n-i-1)*n +(n-j-1)][1];
            fft_array[(n-i-1)*n +(n-j-1)][1] = fft_array[(n-j-1)*n + i][1];
            fft_array[(n-j-1)*n + i][1]		 = tmp;
        }
	}
}
#endif //SHADER_FUNCS
