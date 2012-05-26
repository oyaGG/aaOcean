// aaOcean Mental Ray Main Shader
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

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


