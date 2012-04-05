#ifndef FUNCTION_LIB_H
#define FUNCTION_LIB_H

#include "constants.h"

inline float DegsToRads(float degrees)	// Degrees to radians conversion...
{ 
	return(degrees * aa_PIBY180);  
}
inline float RadsToDegs(float rads)		// Radians to degrees conversion...
{ 
	return(rads * aa_180BYPI); 
}

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

inline bool isfEqual(float x, float y, const float epsilon)
{
   return fabs(x - y) <= epsilon * fabs(x);  // use this to compare float values
}

inline float isEven(int x)
{
	// used in FFT
	// (-1)^(x+y), shift back from centre representation for zero-frequencies to left
	// x, y being coordinates of 2d array holding fft data
    if(!(x % 2))
		return 1.0f;
	else
		return -1.0f;
}

inline long int_sqrt(long r) // paul bourke
{
   long t,b,c = 0;

   for (b = 0x10000000; b != 0; b >>= 2) 
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

inline int round(float x)
{
   return  (int)(x > 0.0f ? x + 0.5f : x - 0.5f);
}

inline bool isInt(float a)
{
	if(static_cast<int>(a) == a)
		return true;
	else
		return false;
}

inline float rescale(const float& value, const float& oldMin, const float& oldMax, const float& newMin, const float& newMax)
{
  const float oldDistance = oldMax - oldMin;
  const float newDistance = newMax - newMin;

  const float distance = (value - oldMin) / oldDistance;
  const float newValue = newMin + (distance * newDistance);

  return newValue;
}

inline int wrap(int x, int n)
{
	if(x > n)
		x = x % (n+1);
	else if(x < 0)
		x = (n+1) + x % (n+1);
	
	return x;
}

float intpow( float base, int exponent )
{
	float out = 1.f, curpwr = base;
	for( ; exponent > 0; exponent = exponent >> 1)
	{
		if ((exponent & 1) > 0)
			base *= curpwr;
		curpwr *= curpwr;
	}
	return out;
}


float fastPow(float a, float b) 
{
  union 
  {
    float d;
    int x[2];
  } 
  u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

inline float lerp(float t, float a, float b)
{ 
	// t = 0 retunes a
	// t = 1 returns b
	return a*(1.0f - t) + b * t;
}

//bool isAligned(void* data, int alignment = 16)
//{
//	// check that the alignment is a power of two
//	assert((alignment & (alignment-1)) == 0); 
//	return ((uintptr_t)data & (alignment-1)) == 0;
//}

//void* aligned_malloc(int size, int alignment)
//{
//	const int pointerSize = sizeof(void*);
//	const int requestedSize = size + alignment - 1 + pointerSize;
//	void* raw = std::malloc(requestedSize);
//	void* start = (char*)raw + pointerSize;
//	void* aligned = (void*)(((unsigned int)((char*)start+alignment-1)) & ~(alignment-1));
//	*(void**)((char*)aligned-pointerSize) = raw;
//	return aligned;
//}
//
//void aligned_free(void* aligned)
//{
//	void* raw = *(void**)((char*)aligned-sizeof(void*));
//	std::free(raw);
//}

#if defined(_MSC_VER)
	#define DATA_ALIGN(declaration, alignment) __declspec(align(alignment)) declaration
#elif defined(GCC)
	#define DATA_ALIGN(declaration, alignment) declaration __attribute__ ((aligned (alignment)))
#else
	#define DATA_ALIGN(declaration, alignment)
#endif

#endif /*FUNCTION_LIB_H*/
