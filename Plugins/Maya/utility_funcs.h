#ifndef SHADER_FUNCS_H
#define SHADER_FUNCS_H

inline void print(MString str)
{
	MGlobal::displayInfo( str);
}
void rotateArray(fftwf_complex *&fft_array, aaOcean *&ocean)
{
	int n = ocean->m_resolution + 1;
	double tmp;
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

#endif /*SHADER_FUNCS*/