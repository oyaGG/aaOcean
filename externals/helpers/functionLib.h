#ifndef FUNCTION_LIB_H
#define FUNCTION_LIB_H

#include <cmath>
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

inline int round(float x)
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
#endif /*FUNCTION_LIB_H*/
