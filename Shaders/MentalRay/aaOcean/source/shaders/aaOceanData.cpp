// aaOcean Mental Ray Main Shader
// Outputs RGBA, with Vector Displacement in RGB, and foam in Alpha
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef aaOceanDataShader_CPP
#define aaOceanDataShader_CPP

#include <string>
#include <shader.h>
#include "aaOceanClass.cpp"
#include "oceanStore.h"
#include "aaOceanData.h"

extern "C" DLLEXPORT 
miBoolean aaOceanDataShader(miColor *result, miState *state, aaOceanDataShader_t *params)
{
	oceanStore** os;
	mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
	aaOcean *ocean = (*os)->ocean;
	
	miScalar	layerOcean	= *mi_eval_scalar(&params->layerOcean);
	miVector	*coord		=  mi_eval_vector(&params->uv_coords);
	miScalar	fade		= *mi_eval_scalar(&params->fade);
	
	if(*mi_eval_boolean(&params->use_uv_input) == FALSE)
	{
		coord->x = state->tex_list[0].x;
		coord->y = state->tex_list[0].y;
	}

	result->g = ocean->getOceanData(coord->x, coord->y, HEIGHTFIELD) * (1.0f - fade);
	
	if(ocean->isChoppy())
	{
		result->r = ocean->getOceanData(coord->x, coord->y, CHOPX)  *  (1.0f - fade);
		result->b = ocean->getOceanData(coord->x, coord->y, CHOPZ)  *  (1.0f - fade);

		miScalar	gamma		= *mi_eval_scalar(&params->gamma);
		miScalar	brightness  = *mi_eval_scalar(&params->brightness);
		miBoolean	rawOutput	= *mi_eval_boolean(&params->rawOutput);
		miScalar	fmin		= *mi_eval_scalar(&params->fmin);
		miScalar	fmax		= *mi_eval_scalar(&params->fmax);

		miScalar	foam		= ocean->getOceanData(coord->x, coord->y, FOAM);

		if(!rawOutput)
		{
			foam  = rescale(foam, fmin, fmax, 0.0f, 1.0f);
			foam  = maximum<float>(foam, 0.0f);
			foam  = pow(foam, gamma);
			foam *= brightness;
			foam *= (1.0f - fade);
		}
		result->a = foam;
	}
	else
		result->a = 0.0f;

	return( miTRUE );
}

extern "C" DLLEXPORT 
void aaOceanDataShader_init(miState *state, aaOceanDataShader_t *params, miBoolean *inst_init_req)
{
	if( !params )
	{
		*inst_init_req = miTRUE;
		fftwf_init_threads();
	}
	else
	{
		// evaluated any previously connected ocean shaders
		miScalar layerOcean	= *mi_eval_scalar(&params->layerOcean);

		// allocate memory for our ocean object
		oceanStore** os = NULL;
		mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
		*os = (oceanStore *)mi_mem_allocate( sizeof(oceanStore));	
		(*os)->ocean = new aaOcean;
		aaOcean *ocean = (*os)->ocean;
	
		// retrieve user input
		ocean->input(*mi_eval_integer(&params->resolution), 
		*mi_eval_integer(&params->seed),
		*mi_eval_scalar(&params->oceanScale), 
		*mi_eval_scalar(&params->velocity), 
		*mi_eval_scalar(&params->cutoff), 
		*mi_eval_scalar(&params->windDir), 
		*mi_eval_integer(&params->windAlign), 
		*mi_eval_scalar(&params->damp), 
		*mi_eval_scalar(&params->waveSpeed), 
		*mi_eval_scalar(&params->waveHeight),
		*mi_eval_scalar(&params->chopAmount), 
		*mi_eval_scalar(&params->time),
		miTRUE);

		// get the tag of the currently running shader so that we can call its name
		miTag shaderInst; 
		mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);

		if(!ocean->isValid())
		{
			mi_error("%s. Shader ID: %d", ocean->getState(), shaderInst);	
			return;
		}
		mi_info("%s. Shader ID: %d", ocean->getState(), shaderInst);	

		if(ocean->isChoppy())
		{
			float outMin, outMax;

			miScalar	fmin		= *mi_eval_scalar(&params->fmin);
			miScalar	fmax		= *mi_eval_scalar(&params->fmax);
			miBoolean	rawOutput	= *mi_eval_boolean(&params->rawOutput);

			ocean->getFoamBounds(*mi_eval_scalar(&params->fmin), *mi_eval_scalar(&params->fmax), outMin, outMax);

			float epsilon = 1e-3f;
			if(!rawOutput && ( !isfEqual(fmin, outMin, epsilon) || !isfEqual(fmax, outMax, epsilon) ))
				mi_warning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
							outMin, outMax);
		}

		//clear arrays that are not required
		ocean->clearResidualArrays();

		mi_info("[aaOcean Shader] Data shader initiated. Shader ID: %d, location: %p", shaderInst, ocean);	
	}
}


extern "C" DLLEXPORT 
void aaOceanDataShader_exit(miState	*state,	aaOceanDataShader_t *params)
{
	if( params )
	{
		oceanStore** os = NULL;
		if(mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os))
		{
			aaOcean *ocean = (*os)->ocean;

			if(ocean)
				delete ocean;
			mi_mem_release( (void *)(*os) );

			miTag shaderInst; // tag of the currently running shader
			mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);
			mi_info("[aaOcean Shader] Data Shader ID %d, terminated at location: %p", shaderInst, ocean);
		}
	}
	else
	{
		fftwf_cleanup_threads();
		fftwf_cleanup();
	}

}

extern "C" DLLEXPORT int aaOceanDataShader_version( )
{
	return( 1 );
}

#endif /* aaOceanDataShader_CPP */