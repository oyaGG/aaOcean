#ifndef aaOceanDisplaceShader_H
#define aaOceanDisplaceShader_H 

#include <shader.h>
#include <geoshader.h>

typedef struct
{
	miVector	uv_coords;			// uv_coords
	miScalar	fade;

	miInteger	resolution;
	miScalar	oceanScale;
	miInteger	seed;
	miScalar	waveHeight;
	miScalar	velocity;
	miScalar	waveSpeed;
	miScalar	chopAmount;
	miScalar	cutoff;
	miScalar	windDir;
	miScalar	damp;
	miInteger	windAlign;
	miScalar	time;
	miScalar	layerOcean;
} aaOceanDisplaceShader_t;

#endif // aaOceanDisplaceShader_H


