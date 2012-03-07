#ifndef SHADER_FUNCS_H
#define SHADER_FUNCS_H

#if (_MSC_VER > 1000)
#pragma once
#endif

void arnoldDebugPrint(AtParamValue *&params)
{
	AiMsgWarning("[aaOcean] Render Resolution %d", params[p_renderResolution].INT);	
	AiMsgWarning("[aaOcean] Ocean Size %f", params[p_oceanSize].FLT);
	AiMsgWarning("[aaOcean] Seed %d", params[p_seed].INT);	
	AiMsgWarning("[aaOcean] Wave Height %f", params[p_waveHeight].FLT);
	AiMsgWarning("[aaOcean] Wave Size %f", params[p_waveSize].FLT);
	AiMsgWarning("[aaOcean] Wave Speed %f", params[p_waveSpeed].FLT);
	AiMsgWarning("[aaOcean] Wave Chop %f", params[p_waveChop].FLT);	
	AiMsgWarning("[aaOcean] Wave Smooth %f", params[p_waveSmooth].FLT);
	AiMsgWarning("[aaOcean] Wave Direction %f", params[p_waveDirection].FLT);
	AiMsgWarning("[aaOcean] Wave Reflection %f", params[p_waveReflection].FLT);	
	AiMsgWarning("[aaOcean] Wind Align %d", params[p_waveAlign].INT);		
	AiMsgWarning("[aaOcean] Time %f", params[p_time].FLT);
}


bool fetchInput(AtParamValue *&params, aaOceanClass *&pOcean)
{	
	float tempFLT;
	int tempINT;
	bool retVal = false;

	int renderRes = params[p_renderResolution].INT  + 4;
	renderRes = (int)pow(2.0f,renderRes);
	pOcean->redoHoK = false;

	tempFLT = ((params[p_waveDirection].FLT)/180) * 3.141592653589793238462643383;
	if((float)pOcean->windDir != (float)tempFLT)
	{
		pOcean->windDir = tempFLT;
		pOcean->redoHoK = true;
	}
	tempFLT = (params[p_waveSmooth].FLT / 100);
	if((float)pOcean->cutoff != (float)tempFLT)
	{
		pOcean->cutoff = tempFLT;
		pOcean->redoHoK = true;
	}
	tempFLT = maximum<float>(params[p_oceanSize].FLT,0.00001);
	if((float)pOcean->oceanScale	!= (float)tempFLT)
	{
		pOcean->oceanScale = tempFLT;
		pOcean->redoHoK = true;
	}
	tempFLT = maximum<float>(((params[p_waveSize].FLT  * params[p_waveSize].FLT) / (9.81)),0.00001);
	if((float)pOcean->velocity !=  (float)tempFLT)
	{
		pOcean->velocity = tempFLT;
		pOcean->redoHoK = true;
	}
	tempINT = maximum<int>	 (((params[p_waveAlign].INT + 1) * 2),2);
	if(pOcean->windAlign != tempINT)
	{
		pOcean->windAlign = tempINT;
		pOcean->redoHoK = true;
	}
	tempFLT = params[p_waveReflection].FLT;
	if((float)pOcean->damp != (float)tempFLT)
	{
		pOcean->damp = tempFLT;
		pOcean->redoHoK = true;
	}

	tempINT =  params[p_seed].INT;
	if(pOcean->seed	!= tempINT)
	{
		pOcean->seed	= tempINT;
		pOcean->redoHoK = true;
		pOcean->setup_grid(); 
	}

	retVal = pOcean->redoHoK;

	tempFLT = params[p_waveChop].FLT  * .01; //divided by scale for better ocean control;
	if((float)pOcean->chopAmount != (float)tempFLT)
	{
		pOcean->chopAmount = tempFLT;
		retVal = true;
	}

	tempFLT = params[p_waveHeight].FLT* .01;		//divided by scale for better ocean control;
	if((float)pOcean->waveHeight != (float)tempFLT)
	{
		pOcean->waveHeight = tempFLT;
		retVal = true;
	}

	tempFLT = params[p_waveSpeed].FLT;
	if((float)pOcean->waveSpeed != (float)tempFLT)
	{
		pOcean->waveSpeed = tempFLT;
		retVal = true;
	}

	tempFLT = params[p_time].FLT;
	if((float)pOcean->time != (float)tempFLT)
	{
		pOcean->time = tempFLT;
		retVal = true;
	}

	if(pOcean->resolution != renderRes)
		retVal = true;

	return retVal;
}


bool fetchInputNodeFinish(AtParamValue *&params, aaOceanClass *&pOcean)
{	
	float tempFLT;
	int tempINT;
	bool retVal = 0;

	tempFLT = ((params[p_waveDirection].FLT)/180) * 3.141592653589793238462643383;
	if((float)pOcean->windDir != (float)tempFLT)
	{
		pOcean->windDir = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = (params[p_waveSmooth].FLT / 100);
	if((float)pOcean->cutoff != (float)tempFLT)
	{
		pOcean->cutoff = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = maximum<float>(params[p_oceanSize].FLT,0.00001);
	if((float)pOcean->oceanScale	!= (float)tempFLT)
	{
		pOcean->oceanScale = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = maximum<float>(((params[p_waveSize].FLT  * params[p_waveSize].FLT) / (9.81)),0.00001);
	if((float)pOcean->velocity !=  (float)tempFLT)
	{
		pOcean->velocity = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempINT = maximum<int>	 (((params[p_waveAlign].INT + 1) * 2),2);
	if(pOcean->windAlign != tempINT)
	{
		pOcean->windAlign = tempINT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = params[p_waveReflection].FLT;
	if((float)pOcean->damp != (float)tempFLT)
	{
		pOcean->damp = tempFLT;
		retVal = pOcean->redoHoK = true;
	}

	tempINT =  params[p_seed].INT;
	if(pOcean->seed	!= tempINT)
	{
		pOcean->seed	= tempINT;
		retVal = pOcean->redoHoK = true;
		pOcean->setup_grid(); 
	}
	
	tempFLT = params[p_waveChop].FLT  * .01;
	if((float)pOcean->chopAmount != tempFLT)
	{
		pOcean->chopAmount	= tempFLT;		
		retVal = 1;
	}
	
	tempFLT = params[p_waveSpeed].FLT;
	if(	(float)pOcean->waveSpeed != tempFLT)
	{
		pOcean->waveSpeed = tempFLT;
		retVal = 1;
	}
	
	pOcean->waveHeight	= params[p_waveHeight].FLT* .01;		//divided by scale for better ocean control;
	pOcean->time		= params[p_time].FLT;

	return retVal;
}

bool fetchInputEval(AtShaderGlobals *&sg, AtNode *&node, aaOceanClass *&pOcean)
{	
	float tempFLT;
	int tempINT;
	bool retVal = false;

	tempFLT = ((AiShaderEvalParamFlt(p_waveDirection)/180) * 3.141592653589793238462643383);
	if((float)pOcean->windDir != (float)tempFLT)
	{
		pOcean->windDir = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = (AiShaderEvalParamFlt(p_waveSmooth) / 100);
	if((float)pOcean->cutoff != (float)tempFLT)
	{
		pOcean->cutoff = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = maximum<float>(AiShaderEvalParamFlt(p_oceanSize),0.00001);
	if((float)pOcean->oceanScale	!= (float)tempFLT)
	{
		pOcean->oceanScale = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = maximum<float>(((AiShaderEvalParamFlt(p_waveSize)  * AiShaderEvalParamFlt(p_waveSize)) / (9.81)),0.00001);
	if((float)pOcean->velocity !=  (float)tempFLT)
	{
		pOcean->velocity = tempFLT;
		retVal = pOcean->redoHoK = true;
	}
	tempINT = maximum<int>	 (((AiShaderEvalParamInt(p_waveAlign) + 1) * 2),2);
	if(pOcean->windAlign != tempINT)
	{
		pOcean->windAlign = tempINT;
		retVal = pOcean->redoHoK = true;
	}
	tempFLT = AiShaderEvalParamFlt(p_waveReflection);
	if((float)pOcean->damp != (float)tempFLT)
	{
		pOcean->damp = tempFLT;
		retVal = pOcean->redoHoK = true;
	}

	tempINT =  AiShaderEvalParamInt(p_seed);
	if(pOcean->seed	!= tempINT)
	{
		pOcean->seed	= tempINT;
		pOcean->redoHoK = 1;
		pOcean->setup_grid();
		retVal = true;
	}
	tempFLT	= AiShaderEvalParamFlt(p_waveChop)   * .01;	
	if(pOcean->chopAmount != tempFLT)
	{
		pOcean->chopAmount	= tempFLT;
		retVal = true;
	}
	
	tempFLT	=AiShaderEvalParamFlt(p_waveSpeed);
	if(pOcean->waveSpeed != tempFLT)
	{
		pOcean->waveSpeed	=  tempFLT;
		retVal = true;
	}

	tempFLT	=AiShaderEvalParamFlt(p_waveSpeed);
	if(pOcean->waveSpeed != tempFLT)
	{
		pOcean->waveSpeed	=  tempFLT;
		retVal = true;
	}
	return retVal;
}

void getArrayBounds(fftw_complex *&fftw_array, int idx, int n, float &min, float &max)
{
	max = FLT_MIN;
	min = FLT_MAX;

	int i,j,index;
	//if fft_array has been copied and tiled, set idx = 1, else 0
	#pragma omp parallel for private(index,i,j)
	for(i = 0; i< n+1; i++)
	{
		for(j = 0; j< n+1; j++)
		{
			index = i*(n+1)+j;
			if(i==n)  //copy left-most col to right-most col
			{
				if(max < fftw_array[index][idx]) //precision -- epsilon check
					max = fftw_array[index][idx];

				if(min > fftw_array[index][idx]) //precision -- epsilon check
					min = fftw_array[index][idx];
			}
		}
	}
}

inline int wrap(int x,int n)
{
	if(x > n)
		x =x%(n+1);
	else if(x < 0)
		x = (n+1) + x%(n+1);
	
	return x;
}

inline double catrom(double &t, double &a, double &b, double &c, double &d) 
{
   return  0.5f * ( ( 2.0f * b ) + ( -a + c ) * t + ( 2.0f * a - 5.0f * b + 4 * c - d ) * t*t + ( -a + 3.0f * b - 3.0f * c + d )* t*t*t );
}


double catromPrep(aaOceanClass *&ocean, fftw_complex *&fftw_array,  AtPoint point)
{
	double u,v,du,dv=0;
	int xMinus1, yMinus1, x, y, xPlus1, yPlus1, xPlus2, yPlus2;
	int n = ocean->resolution;	

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
	x = floor(u);
	y = floor(v);

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

	double a1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) + yMinus1][1],
						fftw_array[x*(n+1)		 + yMinus1][1],
						fftw_array[xPlus1*(n+1)	 + yMinus1][1],
						fftw_array[xPlus2*(n+1)	 + yMinus1][1] 
						);

	double b1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) +	y][1],
						fftw_array[x*(n+1)		 +	y][1],
						fftw_array[xPlus1*(n+1)	 +	y][1],
						fftw_array[xPlus2*(n+1)	 +	y][1]
						);

	double c1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) + yPlus1][1],
						fftw_array[x*(n+1)		 + yPlus1][1],
						fftw_array[xPlus1*(n+1)	 + yPlus1][1],
						fftw_array[xPlus2*(n+1)	 + yPlus1][1] 
						);

	double d1 = catrom(	du, 
						fftw_array[xMinus1*(n+1) + yPlus2][1],
						fftw_array[x*(n+1)		 + yPlus2][1],
						fftw_array[xPlus1*(n+1)	 + yPlus2][1],
						fftw_array[xPlus2*(n+1)	 + yPlus2][1]
						);

	return catrom(dv,a1,b1,c1,d1);
	
}
#endif //SHADER_FUNCS