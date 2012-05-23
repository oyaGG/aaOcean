#ifndef SHADER_FUNCS_H
#define SHADER_FUNCS_H

inline void print(MString str)
{
	MGlobal::displayInfo( str);
}

void copy_and_tile(fftwf_complex *&fft_array, aaOcean *&ocean)
{
	int n = ocean->m_resolution;
	int index, i, j;
	
	#pragma omp parallel for private(index,i,j)
	for(i = 0; i< n+1; i++)
	{					
		for(j = 0; j< n+1; j++)
		{
			index = i*(n+1)+j;
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
			if( i<n && j<n) // regular array copy
			{
				fft_array[index][1] = fft_array[i*n+j][0];
			}
			
		}
	}
}

#endif /*SHADER_FUNCS*/