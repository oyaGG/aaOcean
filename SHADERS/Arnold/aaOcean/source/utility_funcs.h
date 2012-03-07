#ifndef UTILITY_FUNCS_H
#define UTILITY_FUNCS_H

#if (_MSC_VER > 1000)
#pragma once
#endif

template <class T> const T& maximum ( const T& a, const T& b ) 
{
  return (b<a)?a:b;     
}
template <class T> const T& minimum ( const T& a, const T& b ) 
{
  return (b>a)?a:b;     
}

inline float clamp(float x, float a, float b)
{
    return x < a ? a : (x > b ? b : x);
}

inline bool isdEqual(double x, double y, const double epsilon)
{
   return std::abs(x - y) <= epsilon * std::abs(x);  // use this to compare float values
}

inline bool isfEqual(float x, float y, const double epsilon)
{
   return std::abs(x - y) <= epsilon * std::abs(x);  // use this to compare float values
}

inline double isEven(int x)
{
    if(!(x % 2))
		return 1.0;
	else
		return -1.0;
}

inline long int_sqrt(long r) // paul bourke
{
   long t,b,c=0;

   for (b=0x10000000; b!=0; b>>=2) 
   {
      t = c + b;
      c >>= 1;
      if (t <= r) 
	  {
         r -= t;
         c += b;
      }
   }
   return(c);
}

inline int round(double x)
{
   return  (int)(x > 0.0 ? x + 0.5 : x - 0.5);
}

inline bool isInt(float a)
{
	if(static_cast<int>(a) == a)
		return true;
	else
		return false;
}

inline double rescale(const double& value, const double& oldMin, const double& oldMax, const double& newMin, const double& newMax)
{
  const double oldDistance = oldMax - oldMin;
  const double newDistance = newMax - newMin;

  const double distance = (value - oldMin) / oldDistance;
  const double newValue = newMin + (distance * newDistance);

  return newValue;
}

void rotateArray(fftw_complex *&fft_array, aaOceanClass *&ocean)
{
	int n = ocean->resolution + 1;
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

void copy_and_tile(fftw_complex *&fft_array, aaOceanClass *&ocean)
{
	int n = ocean->resolution;
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

#endif /*UTILITY_FUNCS*/