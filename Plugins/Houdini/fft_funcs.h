//////FFT START

#include <stdlib.h>
#include <math.h>


typedef struct {
    double real, imag;
} COMPLEX;



bool PowerOfTwo(int n)
{
   return (n) && !(n & (n - 1)); //this checks if the integer n is a power of two or not
}

void create_2Dcomplex_array(COMPLEX **&comp, int size)
{
	comp=new COMPLEX*[size];

	for (int row = 0; row < size; row++)
		comp[row] = new COMPLEX[size];        
}

void delete_2Dcomplex_array(COMPLEX **&comp, int size)
{
	for (int row = 0; row < size; row++)
			delete [] comp[row];

	delete [] comp;
}


//FFT CODE STARTS HERE. BASED ON CODE FROM PAUL BOURKE AND FFTGURU.COM

void swap(long s1, long s2, double*& r, double*& i) 
{
	double t; // function swap is
	t = r[s1] ; r[s1] = r[s2] ; r[s2] = t ; // needed for bit reverse
	t = i[s1] ; i[s1] = i[s2] ; i[s2] = t ; // at the end of the FFT
}
void fft_1D(long N, double *&r, double *&i) 
{
	long lo = 0, hi, k, m, delta; 
	double t, x, y, twcos, twsin, pi = 3.1415926535897932 ;
	
	for (delta = N/2; delta > 0; delta = delta/2) 
	{
		x = pi/delta;
		for (k = 0; k < N/(delta*2); k++) 
		{
			hi = lo + delta;
			for (m = 0; m < delta; m++) 
			{
				t = r[lo] - r[hi] ; r[lo] = r[lo] + r[hi] ; r[hi] = t ;
				t = i[lo] - i[hi] ; i[lo] = i[lo] + i[hi] ; i[hi] = t ;
				if ( m > 0 && delta > 1 ) 
				{ // do twiddle multiply, but not for twiddles of 0
					y = m*x;
					twcos = cos (y) ; twsin = -sin (y) ;
					t = (r[hi] - i[hi])*twsin ;
					r[hi] = t + r[hi]*(twcos - twsin) ;
					i[hi] = t + i[hi]*(twcos + twsin) ;
				} 
				lo++; 
				hi++;
			} // end for over m
			
			lo = hi%N;
			
		} // end for over k
	} // end for over delta
	
	//******** bit reverse for radix 2 *****
	long J = 0, K, L, N2 = N/2;
	for (L = 0; L < (N-2); L++) 
	{
		if (L < J) 
			swap( L, J, r, i) ;
		K = N2;
		while (K <= J) {
			J = J - K; K = K/2; }
		J = J + K;
	} // end for ***** bit reverse done ******
} 

void fft_2D(COMPLEX **c, int nx, int ny, int direction)
{
	int row, col;
	int N = nx;
	double *real,*imag;

	//---------------------------
	// inverse fft rows

	real = (double *)malloc(nx * sizeof(double));
	imag = (double *)malloc(nx * sizeof(double));

	 
	for (row = 0; row < N; row++)
	{ 

		for (col = 0; col < N; col++)
		{
			real[col] = c[row][col].real; 
			imag[col] = c[row][col].imag;
		}
		if (direction==1) //forward FFT
			fft_1D(N,real,imag);

		if (direction == -1) //inverse FFT
			fft_1D(N,imag,real);

		for (col = 0; col < N; col++)
		{ 
			//put data back into arrays
			if (direction == -1) //inverse FFT
			{
				c[row][col].real = real[col]; 
				c[row][col].imag = imag[col];
			}
			if (direction == 1) //forward FFT
			{
				c[row][col].real = real[col]/N; 
				c[row][col].imag = imag[col]/N;
			}
		}
	}


	//---------------------------
	// inverse fft columns
	for (col = 0; col< N; col++)
	{ 
		//get 8 points (a column) from image arrays
		for (row = 0; row< N; row++)
		{
			real[row] = c[row][col].real; 
			imag[row] = c[row][col].imag;
		}

		if (direction == 1) //forward FFT
			fft_1D(N,real,imag);

		if (direction == -1) //inverse FFT
			fft_1D(N,imag,real);

		for (row = 0; row< N; row++)
		{ //put data back into image arrays
			
			if (direction == -1) //inverse FFT
			{
				c[row][col].real = real[row]; 
				c[row][col].imag = imag[row];
			}

		}
	}

	free(real);
	free(imag);

}