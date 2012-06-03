// aaOcean Mental Ray Shader functions
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef SHADER_FUNCS_H
#define SHADER_FUNCS_H

char* miaux_tag_to_string(miTag tag, char *default_value) 
{ 
    char *result = default_value; 
    if (tag != 0) { 
	result = (char*)mi_db_access(tag); 
	mi_db_unpin(tag); 
    } 
    return result; 
}

void genFullFilePath(char* dest, const char* outputFolder, const char* postfix, const int frame)
{
	#ifdef _MSC_VER
	char* slash = "\\";
	#else
	char* slash = "/";
	#endif

	char cFrame[10];
	char prefixFrame[10];
	if ( frame < 10)
		strcpy(prefixFrame, "000");
	else if ( frame < 100) 
		strcpy(prefixFrame, "00");
	else if ( frame < 1000) 
		strcpy(prefixFrame, "0");
	sprintf(cFrame, "%s%d", prefixFrame, frame);
					
	strcpy(dest, outputFolder);
	strcat(dest, slash);

	if(strcmp(postfix,"") == 0)
		strcat(dest, "aaOceanData");
	else
		strcat(dest, "aaOceanData_");

	strcat(dest, postfix);
	strcat(dest, ".");
	strcat(dest, cFrame);
	strcat(dest, ".exr");
}

void writeExr(const char* outputFileName, aaOcean *&ocean)
{
	int res = ocean->m_resolution;
	Imf::Array2D<Imf::Rgba> pixels(res, res);

	#pragma omp parallel for
	for(int i = 0; i < res; i++)
	{
		for(int j = 0; j < res; j++)
		{
			pixels[j][i].r = ocean->m_fft_chopX[i*res+j][1];
			pixels[j][i].g = ocean->m_fft_htField[i*res+j][1];
			pixels[j][i].b = ocean->m_fft_chopZ[i*res+j][1];
			pixels[j][i].a = ocean->m_fft_jxz[i*res+j][1];
		}
	}
	Imf::RgbaOutputFile outA(&outputFileName[0], Imf::Header(ocean->m_resolution, ocean->m_resolution), Imf::WRITE_RGBA);
	outA.setFrameBuffer (&pixels[0][0], 1, ocean->m_resolution);
	outA.writePixels (ocean->m_resolution);
}

void shaderCleanup(aaOcean *&ocean)
{
	if(ocean->m_kX)
	{
		free(ocean->m_kX);
		ocean->m_kX=0;
	}
	if(ocean->m_kZ)
	{
		free(ocean->m_kZ);
		ocean->m_kZ=0;
	}
	if(ocean->m_omega)
	{
		free(ocean->m_omega);
		ocean->m_omega=0;
	}
	if(ocean->m_hokReal)
	{
		free(ocean->m_hokReal);
		ocean->m_hokReal=0;
	}
	if(ocean->m_hokImag)
	{
		free(ocean->m_hokImag);
		ocean->m_hokImag=0;
	}
	if(ocean->m_hktReal)
	{
		free(ocean->m_hktReal);
		ocean->m_hktReal=0;
	}
	if(ocean->m_hktImag)
	{
		free(ocean->m_hktImag);
		ocean->m_hktImag=0;
	}
	if(ocean->m_rand1)
	{
		free(ocean->m_rand1);
		ocean->m_rand1=0;
	}
	if(ocean->m_rand2)
	{
		free(ocean->m_rand2);
		ocean->m_rand2=0;
	}
	if(ocean->m_xCoord)
	{
		free(ocean->m_xCoord);
		ocean->m_xCoord=0;
	}
	if(ocean->m_zCoord)
	{
		free(ocean->m_zCoord);
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

#endif //SHADER_FUNCS
