// aaOcean Mental Ray Main Shader
// Outputs RGBA, with Vector Displacement in RGB, and foam in Alpha
// Author: Amaan Akram 
// www.amaanakram.com
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html

#ifndef aaOceanDataShader_CPP
#define aaOceanDataShader_CPP

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

	// rotating UV space by 90 degrees to align with ICE representation
	coord->z = coord->x;
	coord->x = coord->y * -1.0f;
	coord->y = coord->z;

	result->g = catromPrep(ocean, ocean->m_fft_htField,  state,  coord) * (1.0f - fade);
	
	if(ocean->m_chopAmount > 0.0f)
	{
		result->r = -catromPrep( ocean, ocean->m_fft_chopX,  state,  coord)  *  (1.0f - fade);
		result->b = -catromPrep( ocean, ocean->m_fft_chopZ,  state,  coord)  *  (1.0f - fade);

		miScalar	gamma		= *mi_eval_scalar(&params->gamma);
		miScalar	brightness  = *mi_eval_scalar(&params->brightness);
		miBoolean	rawOutput	= *mi_eval_boolean(&params->rawOutput);
		miScalar	fmin		= *mi_eval_scalar(&params->fmin);
		miScalar	fmax		= *mi_eval_scalar(&params->fmax);

		miScalar	foam		= catromPrep(ocean, ocean->m_fft_jxz,state,coord);

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
		Timer timer;
		timer.start();

		miScalar layerOcean	= *mi_eval_scalar(&params->layerOcean);

		oceanStore** os = NULL;
		mi_query( miQ_FUNC_USERPTR, state, 0, (void *)&os);
		*os = (oceanStore *)mi_mem_allocate( sizeof(oceanStore));	
		(*os)->ocean = new aaOcean;
		aaOcean *ocean = (*os)->ocean;
	
		int resolution			= (int)pow(2.0f,(4 + (*mi_eval_integer(&params->resolution))));
		ocean->m_oceanScale		= maximum<float>(*mi_eval_scalar(&params->oceanScale),0.00001f);
		ocean->m_seed			= *mi_eval_integer(&params->seed);
		ocean->m_waveHeight		= *mi_eval_scalar(&params->waveHeight) * .01f;
		ocean->m_velocity		= maximum<float>((((*mi_eval_scalar(&params->velocity)) * (*mi_eval_scalar(&params->velocity))) / (aa_GRAVITY)),0.00001f);
		ocean->m_waveSpeed		= *mi_eval_scalar(&params->waveSpeed);
		ocean->m_chopAmount		= *mi_eval_scalar(&params->chopAmount) * 0.01f;
		ocean->m_cutoff			= fabs(*mi_eval_scalar(&params->cutoff) * 0.01f);
		ocean->m_windDir		= DegsToRads(*mi_eval_scalar(&params->windDir));
		ocean->m_windAlign		= maximum<int>(((*mi_eval_integer(&params->windAlign) + 1) * 2), 2); 
		ocean->m_damp			= minimum<float>(*mi_eval_scalar(&params->damp),1);
		ocean->m_time			= *mi_eval_scalar(&params->time);

		miTag shaderInst; // tag of the currently running shader
		mi_query(miQ_FUNC_TAG, state, 0, &shaderInst);

		if(!ocean->reInit(resolution))
		{
			mi_error("%s. Shader ID: %d", ocean->m_state, shaderInst);	
			return;
		}
		mi_info("%s. Shader ID: %d", ocean->m_state, shaderInst);	

		ocean->prepareOcean(TRUE, TRUE, TRUE, FALSE, TRUE);

		bool doChop = FALSE;
		if(ocean->m_chopAmount > 0.0f)
			doChop = TRUE;
		
		if(doChop)
		{
			getArrayBounds(ocean->m_fft_jxz, 1, ocean->m_resolution, ocean->m_fmin, ocean->m_fmax);

			miScalar	fmin		= *mi_eval_scalar(&params->fmin);
			miScalar	fmax		= *mi_eval_scalar(&params->fmax);
			miBoolean	rawOutput	= *mi_eval_boolean(&params->rawOutput);

			float epsilon = 1e-3f;
			if(( !isfEqual(fmin, ocean->m_fmin, epsilon) || !isfEqual(fmax, ocean->m_fmax, epsilon) ) && !rawOutput)
				mi_warning("[aaOcean Shader] Foam Min/Max mismatch. Please set the Foam Min/Max values in foam shader to Min: %f, Max: %f", 
							ocean->m_fmin, ocean->m_fmax);
		}

		//clear arrays that are not required
		ocean->clearResidualArrays();
		timer.stop();

		mi_info("[aaOcean Shader] Data shader initiated in %f seconds. Shader ID: %d, location: %p", 
				timer.getElapsedTimeInSec(), shaderInst, ocean);	
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
		// commented out because of possible conflict between openmp and fftw
		//fftwf_cleanup_threads();
		//fftwf_cleanup();
	}

}

extern "C" DLLEXPORT int aaOceanDataShader_version( )
{
	return( 1 );
}

#endif /* aaOceanDataShader_CPP */