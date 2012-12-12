// aaOcean Mental Ray Main Shader
// Outputs RGBA, with Vector Displacement in RGB, and foam in Alpha
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef aaOceanDataShader_H
#define aaOceanDataShader_H 

#include <shader.h>

typedef struct
{
	miVector	uv_coords;			// uv_coords
	miBoolean	use_uv_input;
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

	miScalar	gamma;
	miScalar	brightness;
	miBoolean	rawOutput;
	miScalar	fmin;
	miScalar	fmax;

	miScalar	layerOcean;
	miMatrix	transform;

	miBoolean	writeFile;
	miTag		outputFolder;
	miTag		postfix;
	miInteger	currentFrame;

} aaOceanDataShader_t;

#endif // aaOceanDataShader_H


